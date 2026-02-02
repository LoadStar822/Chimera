/*
 * Chimera build entry point (split from monolithic ChimeraBuild.cpp)
 */
#include "ChimeraBuildCommon.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace ChimeraBuild {

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
  }
  HashFrequencyContext hashFreqContext;
  build_hash_frequency_sketch(config, inputFiles, hashFreqContext);
  std::filesystem::path dir = "tmp";
  createOrResetDirectory(dir, config);
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
    groups = chimera::imcf::partitionHashCount(hashCount, 16);
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
