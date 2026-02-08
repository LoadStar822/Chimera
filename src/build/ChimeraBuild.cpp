/*
 * Chimera build entry point (split from monolithic ChimeraBuild.cpp)
 */
#include "ChimeraBuildCommon.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#if defined(__GLIBC__)
#include <malloc.h>
#endif
#include <queue>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>
#if defined(__linux__)
#include <sys/resource.h>
#include <unistd.h>
#endif

namespace ChimeraBuild {

namespace {

struct MemorySnapshot {
  uint64_t rss_kb{0};
  uint64_t hwm_kb{0};
};

inline MemorySnapshot read_memory_snapshot_kb() {
  MemorySnapshot out{};
#if defined(__linux__)
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line)) {
    if (line.rfind("VmRSS:", 0) == 0) {
      std::istringstream iss(line.substr(std::strlen("VmRSS:")));
      uint64_t value = 0;
      std::string unit;
      if (iss >> value >> unit) {
        out.rss_kb = value;
      }
    } else if (line.rfind("VmHWM:", 0) == 0) {
      std::istringstream iss(line.substr(std::strlen("VmHWM:")));
      uint64_t value = 0;
      std::string unit;
      if (iss >> value >> unit) {
        out.hwm_kb = value;
      }
    }
  }
  if (out.hwm_kb == 0) {
    struct rusage usage {};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
      out.hwm_kb = static_cast<uint64_t>(usage.ru_maxrss);
    }
  }
  if (out.rss_kb == 0) {
    out.rss_kb = out.hwm_kb;
  }
#endif
  return out;
}

inline void log_memory_checkpoint(const char *stage) {
  const auto mem = read_memory_snapshot_kb();
  if (mem.rss_kb == 0 && mem.hwm_kb == 0) {
    return;
  }
  auto kb_to_mib = [](uint64_t kb) {
    return static_cast<double>(kb) / 1024.0;
  };
  std::ostringstream oss;
  oss << "  [mem] " << stage << ": rss=" << std::fixed << std::setprecision(1)
      << kb_to_mib(mem.rss_kb) << " MiB, hwm=" << kb_to_mib(mem.hwm_kb)
      << " MiB";
  std::cout << oss.str() << std::endl;
}

} // namespace

static std::vector<chimera::imcf::Group> partitionHashCountNoSplit(
    const robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    int maxTaxidsPerGroup = 16) {
  if (hashCount.empty()) {
    return {};
  }
  if (maxTaxidsPerGroup <= 0) {
    throw std::invalid_argument(
        "IMCF partition: maxTaxidsPerGroup must be positive");
  }

  const size_t maxTaxids = static_cast<size_t>(maxTaxidsPerGroup);
  std::vector<std::pair<std::string, uint64_t>> items;
  items.reserve(hashCount.size());
  for (const auto &kv : hashCount) {
    items.emplace_back(kv.first, kv.second);
  }
  std::sort(items.begin(), items.end(),
            [](const auto &a, const auto &b) {
              if (a.second != b.second) {
                return a.second > b.second;
              }
              return a.first < b.first;
            });

  size_t groupCount = (items.size() + maxTaxids - 1) / maxTaxids;
  groupCount = std::max<size_t>(1, groupCount);
  std::vector<chimera::imcf::Group> groups(groupCount);

  struct QueueItem {
    uint64_t total{0};
    size_t size{0};
    size_t index{0};
  };
  struct QueueCmp {
    bool operator()(const QueueItem &a, const QueueItem &b) const {
      if (a.total != b.total) {
        return a.total > b.total;
      }
      if (a.size != b.size) {
        return a.size > b.size;
      }
      return a.index > b.index;
    }
  };

  std::priority_queue<QueueItem, std::vector<QueueItem>, QueueCmp> pq;
  for (size_t i = 0; i < groupCount; ++i) {
    pq.push({0, 0, i});
  }

  for (const auto &item : items) {
    if (pq.empty()) {
      groups.emplace_back();
      pq.push({0, 0, groups.size() - 1});
    }

    while (!pq.empty()) {
      QueueItem head = pq.top();
      pq.pop();

      auto &group = groups[head.index];
      if (head.total != group.totalHash || head.size != group.taxids.size()) {
        continue; // stale
      }
      if (group.taxids.size() >= maxTaxids) {
        continue;
      }

      group.taxids.push_back(item.first);
      group.assignedHashes.push_back(item.second);
      group.totalHash += item.second;
      pq.push({group.totalHash, group.taxids.size(), head.index});
      break;
    }
  }

  std::vector<chimera::imcf::Group> compact;
  compact.reserve(groups.size());
  for (auto &g : groups) {
    if (!g.taxids.empty()) {
      compact.push_back(std::move(g));
    }
  }
  return compact;
}

void run(BuildConfig config) {
  if (config.threads == 0) {
    unsigned int hardwareThreads = std::thread::hardware_concurrency();
    if (hardwareThreads == 0) {
      hardwareThreads = 1;
    }
    const auto maxThreads =
        static_cast<unsigned int>(std::numeric_limits<uint16_t>::max());
    if (hardwareThreads > maxThreads) {
      hardwareThreads = maxThreads;
    }
    config.threads = static_cast<uint16_t>(hardwareThreads);
  }
  const char *kTaxonomyMetaFilename = "taxonomy.meta";
  auto trim = [](std::string &value) {
    while (!value.empty() &&
           std::isspace(static_cast<unsigned char>(value.front()))) {
      value.erase(value.begin());
    }
    while (!value.empty() &&
           std::isspace(static_cast<unsigned char>(value.back()))) {
      value.pop_back();
    }
  };
  if ((config.taxonomy_kind == "auto" || config.taxonomy_kind.empty()) ||
      (config.taxonomy_version == "auto" ||
       config.taxonomy_version.empty())) {
    std::filesystem::path metaPath =
        std::filesystem::path(config.input_file).parent_path() /
        kTaxonomyMetaFilename;
    if (std::filesystem::exists(metaPath)) {
      std::ifstream metaStream(metaPath);
      std::string line;
      while (std::getline(metaStream, line)) {
        auto pos = line.find('=');
        if (pos == std::string::npos) {
          continue;
        }
        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);
        trim(key);
        trim(value);
        if (key == "taxonomy_kind") {
          if (config.taxonomy_kind == "auto" ||
              config.taxonomy_kind.empty()) {
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char ch) {
                             return static_cast<char>(std::tolower(ch));
                           });
            if (!value.empty()) {
              config.taxonomy_kind = value;
            }
          }
        } else if (key == "taxonomy_version") {
          if (config.taxonomy_version == "auto" ||
              config.taxonomy_version.empty()) {
            if (!value.empty()) {
              config.taxonomy_version = value;
            }
          }
        }
      }
    }
  }
  if (config.verbose) {
    std::cout << config << std::endl;
  }
  omp_set_num_threads(config.threads);
  auto build_start = std::chrono::high_resolution_clock::now();
  auto read_start = std::chrono::high_resolution_clock::now();
  std::cout << "Reading input files..." << std::endl;
  FileInfo fileInfo;
  if (config.smer_size == 0) {
    std::cerr << "Syncmer s-mer size must be greater than 0." << std::endl;
    return;
  }
  if (config.smer_size >= config.kmer_size) {
    std::cerr << "Syncmer s-mer size must be smaller than k-mer size."
              << std::endl;
    return;
  }
  const uint16_t syncmer_span =
      static_cast<uint16_t>(config.kmer_size - config.smer_size + 1);
  if (config.syncmer_position >= syncmer_span) {
    std::cerr << "Syncmer offset must satisfy 0 <= pos < k - s + 1."
              << std::endl;
    return;
  }
  robin_hood::unordered_flat_map<std::string, uint64_t> hashCount;
  robin_hood::unordered_flat_map<std::string, uint64_t> bpCount;
  robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
      inputFiles;
  parseInputFile(config.input_file, inputFiles, hashCount, fileInfo);
  auto read_end = std::chrono::high_resolution_clock::now();
  auto read_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                             read_end - read_start)
                             .count();
  if (config.verbose) {
    std::cout << "Read time: ";
    print_build_time(read_total_time);
    std::cout << std::endl;
    log_memory_checkpoint("after_read_inputs");
  }
  HashFrequencyContext hashFreqContext;
  build_hash_frequency_sketch(config, inputFiles, hashFreqContext);
  if (config.verbose) {
    log_memory_checkpoint("after_hash_frequency_sketch");
  }
  std::filesystem::path tmpRoot = "tmp";
  if (const char *envTmpRoot = std::getenv("CHIMERA_BUILD_TMP_ROOT")) {
    std::string value(envTmpRoot);
    if (!value.empty()) {
      tmpRoot = value;
    }
  }
#if defined(__linux__)
  const auto pid_part = static_cast<long long>(::getpid());
#else
  const auto pid_part = static_cast<long long>(std::hash<std::thread::id>{}(
      std::this_thread::get_id()));
#endif
  const auto tmpNonce = static_cast<unsigned long long>(
      std::chrono::steady_clock::now().time_since_epoch().count());
  std::filesystem::path dir =
      tmpRoot / ("build_" + std::to_string(pid_part) + "_" +
                 std::to_string(tmpNonce));
  set_tmp_work_dir(dir);
  createOrResetDirectory(dir.string(), config);
  if (config.verbose) {
    std::cout << "Build tmp directory: " << dir.string() << std::endl;
  }
  auto calculate_start = std::chrono::high_resolution_clock::now();
  std::cout << "Calculating feature hashes..." << std::endl;
  syncmer_count(config, inputFiles, hashCount, fileInfo,
                hashFreqContext.enabled() ? &hashFreqContext : nullptr,
                &bpCount);
  auto calculate_end = std::chrono::high_resolution_clock::now();
  auto calculate_total_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(calculate_end -
                                                            calculate_start)
          .count();
  if (config.verbose) {
    std::cout << "Calculate time: ";
    print_build_time(calculate_total_time);
    std::cout << "File information:" << std::endl;
    std::cout << "Number of files: " << fileInfo.fileNum << std::endl;
    std::cout << "Number of invalid files: " << fileInfo.invalidNum
              << std::endl;
    std::cout << "Number of skipped files: " << fileInfo.skippedNum
              << std::endl;
    std::cout << "Number of sequences: " << fileInfo.sequenceNum << std::endl;
    std::cout << "Number of skipped sequences: " << fileInfo.skippedSeqNum
              << std::endl;
    std::cout << "Total base pairs: " << fileInfo.bpLength << std::endl
              << std::endl;
    log_memory_checkpoint("after_feature_hash_calculation");
  }
  if (hashFreqContext.enabled()) {
    const uint64_t total_checked =
        hashFreqContext.passB_total_hashes.load(std::memory_order_relaxed);
    const uint64_t filtered = hashFreqContext.passB_filtered_hashes.load(
        std::memory_order_relaxed);
    const uint64_t kept =
        (total_checked >= filtered)
            ? (total_checked - filtered)
            : 0;

    auto fmt_ratio = [&](uint64_t num) -> std::string {
      if (total_checked == 0) {
        return "0.00";
      }
      std::ostringstream ss;
      double ratio = static_cast<double>(num) * 100.0 /
                     static_cast<double>(total_checked);
      ss << std::fixed << std::setprecision(2) << ratio;
      return ss.str();
    };

    std::cout << "Hash DF filter dropped " << filtered << " / " << total_checked
              << " (" << fmt_ratio(filtered) << "%)"
              << " with df >= " << hashFreqContext.stats.df_high_threshold
              << std::endl;
    std::cout << "Hash DF kept " << kept << " / " << total_checked << " ("
              << fmt_ratio(kept) << "%)" << std::endl;
  }
  // Input path lists are only needed for feature extraction.
  robin_hood::unordered_flat_map<std::string, std::vector<std::string>>()
      .swap(inputFiles);
#if defined(__GLIBC__)
  ::malloc_trim(0);
#endif
  if (config.verbose) {
    log_memory_checkpoint("after_release_input_paths");
  }
  chimera::feature::Method imcf_feature_method{};
  uint64_t imcf_feature_seed = 0;
  auto imcf_feature_params =
      make_feature_params(config, imcf_feature_method, imcf_feature_seed);
  const std::string feature_suffix =
      (imcf_feature_method == chimera::feature::Method::Strobemer) ? ".strb"
                                                                   : ".sync";
  const uint16_t effective_span = static_cast<uint16_t>(
      std::min<size_t>(std::numeric_limits<uint16_t>::max(),
                       chimera::feature::min_required_length(
                           imcf_feature_params)));
  constexpr uint16_t ref_read_len = 150;
  chimera::presence::CoverageMeta presence_meta;
  auto partition_start = std::chrono::high_resolution_clock::now();
  std::vector<chimera::imcf::Group> groups;
  if (config.max_hashes_per_taxid > 0) {
    std::cout << "Partitioning hashcout (no-split)..." << std::endl;
    groups = partitionHashCountNoSplit(hashCount, 16);
  } else {
    std::cout << "Partitioning hashcout..." << std::endl;
    struct PartitionEval {
      uint64_t area{std::numeric_limits<uint64_t>::max()};
      size_t hashSize{std::numeric_limits<size_t>::max()};
      uint64_t maxLoad{std::numeric_limits<uint64_t>::max()};
      size_t groupCount{std::numeric_limits<size_t>::max()};
      double imbalance{std::numeric_limits<double>::infinity()};
      int maxTaxids{16};
      long long elapsedMs{0};
    };
    auto next_pow2_local = [](size_t v) -> size_t {
      if (v <= 1) {
        return 1;
      }
      --v;
      v |= v >> 1;
      v |= v >> 2;
      v |= v >> 4;
      v |= v >> 8;
      v |= v >> 16;
#if SIZE_MAX > 0xFFFFFFFFu
      v |= v >> 32;
#endif
      return v + 1;
    };
    const uint64_t effectiveLoadDen = std::max<uint64_t>(
        1, static_cast<uint64_t>(config.load_factor * 4.0));
    auto score_groups = [&](const std::vector<chimera::imcf::Group> &cand,
                            int maxTaxids,
                            long long elapsedMs) -> PartitionEval {
      PartitionEval out;
      out.maxTaxids = maxTaxids;
      out.elapsedMs = elapsedMs;
      out.groupCount = cand.size();
      uint64_t totalHash = 0;
      uint64_t maxLoad = 0;
      for (const auto &g : cand) {
        totalHash += g.totalHash;
        if (g.totalHash > maxLoad) {
          maxLoad = g.totalHash;
        }
      }
      out.maxLoad = maxLoad;
      const double needBinSize =
          std::ceil(static_cast<double>(ceil_div_u64(maxLoad, effectiveLoadDen)) *
                    1.05);
      out.hashSize = next_pow2_local(
          static_cast<size_t>(std::max<double>(1.0, needBinSize)));
      out.area = static_cast<uint64_t>(out.groupCount) *
                 static_cast<uint64_t>(out.hashSize);
      if (out.groupCount > 0) {
        const double avgLoad =
            static_cast<double>(totalHash) / static_cast<double>(out.groupCount);
        out.imbalance =
            (avgLoad > 0.0) ? (static_cast<double>(maxLoad) / avgLoad) : 0.0;
      } else {
        out.imbalance = 0.0;
      }
      return out;
    };
    auto better = [](const PartitionEval &lhs, const PartitionEval &rhs) -> bool {
      if (lhs.area != rhs.area) {
        return lhs.area < rhs.area;
      }
      if (lhs.hashSize != rhs.hashSize) {
        return lhs.hashSize < rhs.hashSize;
      }
      if (lhs.maxLoad != rhs.maxLoad) {
        return lhs.maxLoad < rhs.maxLoad;
      }
      if (lhs.groupCount != rhs.groupCount) {
        return lhs.groupCount < rhs.groupCount;
      }
      if (lhs.imbalance != rhs.imbalance) {
        return lhs.imbalance < rhs.imbalance;
      }
      return lhs.elapsedMs < rhs.elapsedMs;
    };

    std::vector<int> candidateMaxTaxids = {12, 11, 10, 9};
    if (const char *env = std::getenv("CHIMERA_PARTITION_MAXTAXIDS")) {
      std::vector<int> parsed;
      std::stringstream ss(env);
      std::string token;
      while (std::getline(ss, token, ',')) {
        token.erase(std::remove_if(token.begin(), token.end(),
                                   [](unsigned char c) { return std::isspace(c); }),
                    token.end());
        if (token.empty()) {
          continue;
        }
        try {
          int value = std::stoi(token);
          if (value > 0) {
            parsed.push_back(value);
          }
        } catch (...) {
          // ignore invalid token
        }
      }
      if (!parsed.empty()) {
        std::sort(parsed.begin(), parsed.end(), std::greater<int>());
        parsed.erase(std::unique(parsed.begin(), parsed.end()), parsed.end());
        candidateMaxTaxids.swap(parsed);
      }
    }
    const size_t candidateCount = candidateMaxTaxids.size();
    std::vector<PartitionEval> evals(candidateCount);
    std::vector<uint8_t> evalReady(candidateCount, 0u);
    PartitionEval bestEval{};
    size_t bestIdx = std::numeric_limits<size_t>::max();
    bool hasBest = false;
    auto better_or_earlier = [&](const PartitionEval &lhs, size_t lhsIdx,
                                 const PartitionEval &rhs,
                                 size_t rhsIdx) -> bool {
      if (better(lhs, rhs)) {
        return true;
      }
      if (better(rhs, lhs)) {
        return false;
      }
      return lhsIdx < rhsIdx;
    };

#pragma omp parallel for schedule(dynamic, 1)
    for (ptrdiff_t idx = 0; idx < static_cast<ptrdiff_t>(candidateCount); ++idx) {
      const size_t candIdx = static_cast<size_t>(idx);
      const int maxTaxids = candidateMaxTaxids[candIdx];
      auto candStart = std::chrono::high_resolution_clock::now();
      auto cand = chimera::imcf::partitionHashCount(hashCount, maxTaxids,
                                                    config.load_factor);
      auto candEnd = std::chrono::high_resolution_clock::now();
      const auto candMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                              candEnd - candStart)
                              .count();
      PartitionEval eval = score_groups(cand, maxTaxids, candMs);
      evals[candIdx] = eval;
      evalReady[candIdx] = 1u;
#pragma omp critical(chimera_partition_pick_best)
      {
        if (!hasBest || better_or_earlier(eval, candIdx, bestEval, bestIdx)) {
          bestEval = eval;
          groups = std::move(cand);
          bestIdx = candIdx;
          hasBest = true;
        }
      }
    }

    for (size_t i = 0; i < candidateCount; ++i) {
      if (evalReady[i] == 0u) {
        continue;
      }
      const auto &eval = evals[i];
      if (config.verbose) {
        std::cout << "  partition candidate maxTaxids="
                  << candidateMaxTaxids[i]
                  << ": groups=" << eval.groupCount
                  << ", hashSize=" << eval.hashSize
                  << ", maxLoad=" << eval.maxLoad
                  << ", area=" << eval.area
                  << ", imbalance=" << eval.imbalance
                  << ", time=" << eval.elapsedMs << " ms" << std::endl;
      }
    }
    if (config.verbose && hasBest) {
      std::cout << "  selected partition maxTaxids=" << bestEval.maxTaxids
                << " (groups=" << bestEval.groupCount
                << ", hashSize=" << bestEval.hashSize
                << ", area=" << bestEval.area
                << ", time=" << bestEval.elapsedMs << " ms)" << std::endl;
    }
  }
  auto partition_end = std::chrono::high_resolution_clock::now();
  auto partition_total_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(partition_end -
                                                            partition_start)
          .count();
  if (config.verbose) {
    std::cout << "Partition time: ";
    print_build_time(partition_total_time);
    std::cout << std::endl;
    log_memory_checkpoint("after_partition");
  }
  auto imcf_build_start = std::chrono::high_resolution_clock::now();
  std::cout << "Building IMCF..." << std::endl;
  IMCFConfig imcfConfig;
  imcfConfig.loadFactor = config.load_factor;
  imcfConfig.kmerSize = config.kmer_size;
  imcfConfig.smerSize = config.smer_size;
  imcfConfig.syncmerPosition = config.syncmer_position;
  imcfConfig.seed64 = imcf_feature_seed;
  imcfConfig.fpSalt = IMCFConfig::DefaultFingerprintSalt;
  imcfConfig.hashVersion = IMCFConfig::CurrentHashVersion;
  imcfConfig.presenceUniqueDeg = config.presence_unique_deg;
  if (imcf_feature_method == chimera::feature::Method::Strobemer) {
    imcfConfig.featureMethod = 1;
    imcfConfig.strobeOrder = config.strobemer_order;
    imcfConfig.strobeWmin = config.strobemer_w_min;
    imcfConfig.strobeWmax = config.strobemer_w_max;
    imcfConfig.strobeK = config.strobemer_k;
  } else {
    imcfConfig.featureMethod = 0;
    imcfConfig.strobeOrder = 0;
    imcfConfig.strobeWmin = 0;
    imcfConfig.strobeWmax = 0;
    imcfConfig.strobeK = 0;
  }
  if (config.taxonomy_kind == "auto" || config.taxonomy_kind.empty()) {
    imcfConfig.taxonomyKind = "ncbi";
  } else {
    imcfConfig.taxonomyKind = config.taxonomy_kind;
  }
  if (config.taxonomy_version == "auto" || config.taxonomy_version.empty()) {
    imcfConfig.taxonomyVersion = "ncbi-taxdump";
  } else {
    imcfConfig.taxonomyVersion = config.taxonomy_version;
  }
  chimera::imcf::InterleavedMergedCuckooFilter imcf(groups, imcfConfig);
  std::vector<std::vector<std::string>> indexToTaxid =
      buildIMCF(imcf, groups, hashCount, feature_suffix,
                hashFreqContext.enabled() ? &hashFreqContext : nullptr,
                &bpCount, effective_span, ref_read_len,
                config.presence_unique_deg, &presence_meta);
  if (presence_meta.entries.empty()) {
    std::cout << "  coverage meta: empty (fallback to runtime heuristics)"
              << std::endl;
  } else {
    std::cout << "  coverage meta: " << presence_meta.entries.size()
              << " taxa cached" << std::endl;
  }
  auto imcf_build_end = std::chrono::high_resolution_clock::now();
  auto imcf_build_total_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(imcf_build_end -
                                                            imcf_build_start)
          .count();
  if (config.verbose) {
    std::cout << "Build time: ";
    print_build_time(imcf_build_total_time);
    std::cout << std::endl;
    std::cout << imcf << std::endl;
    log_memory_checkpoint("after_imcf_build");
  }
  if (config.verbose) {
    std::cout << "Releasing pre-QIMCF buffers..." << std::endl;
  }
  robin_hood::unordered_flat_map<std::string, std::vector<std::string>>()
      .swap(inputFiles);
  robin_hood::unordered_flat_map<std::string, uint64_t>().swap(hashCount);
  robin_hood::unordered_flat_map<std::string, uint64_t>().swap(bpCount);
  std::vector<chimera::imcf::Group>().swap(groups);
  hashFreqContext.sketch.reset();
#if defined(__GLIBC__)
  // Return free heap pages before QIMCF materialization to lower RSS peak.
  ::malloc_trim(0);
#endif
  if (config.verbose) {
    log_memory_checkpoint("before_qimcf_build");
  }
  auto qidx_start = std::chrono::high_resolution_clock::now();
  std::cout << "Building QIMCF..." << std::endl;
  constexpr bool kQimcfVerify = false;
  constexpr bool kQimcfLowPeak = true;
  bool drop_classic_before_materialize = kQimcfLowPeak && !kQimcfVerify;
  std::cout << "  QIMCF options: low_peak=on, verify=off (fixed)" << std::endl;
  imcf.build_query_index(/*include_stash=*/true, kQimcfVerify, kQimcfLowPeak,
                         drop_classic_before_materialize);
  imcf.set_storage_mode(1);
  imcf.release_classic_storage();
  auto qidx_end = std::chrono::high_resolution_clock::now();
  auto qidx_total_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(qidx_end -
                                                            qidx_start)
          .count();
  if (config.verbose) {
    std::cout << "QIMCF build time: ";
    print_build_time(qidx_total_time);
    std::cout << std::endl;
    log_memory_checkpoint("after_qimcf_build");
  }
  auto save_start = std::chrono::high_resolution_clock::now();
  std::cout << "Saving IMCF..." << std::endl;
  saveIMCF(imcf, config.output_file, indexToTaxid, imcfConfig, &presence_meta);
  auto save_end = std::chrono::high_resolution_clock::now();
  auto save_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                             save_end - save_start)
                             .count();
  if (config.verbose) {
    std::cout << "Save time: ";
    print_build_time(save_total_time);
    std::cout << std::endl;
    log_memory_checkpoint("after_save_imcf");
  }
  auto build_end = std::chrono::high_resolution_clock::now();
  auto build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                              build_end - build_start)
                              .count();
  if (config.verbose) {
    std::cout << "Total build time: ";
    print_build_time(build_total_time);
    std::cout << "Remove temporary files..." << std::endl;
  }
  std::filesystem::remove_all(dir);
}

} // namespace ChimeraBuild
