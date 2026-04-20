/*
 * Chimera build entry point (split from monolithic ChimeraBuild.cpp)
 */
#include "ChimeraBuildCommon.hpp"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#if defined(__GLIBC__)
#include <malloc.h>
#endif
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>
#if defined(__linux__)
#include <unistd.h>
#endif

namespace ChimeraBuild {

namespace {

constexpr uint64_t kBlockClassicTableTargetBytes =
    72ull * 1024ull * 1024ull * 1024ull;
constexpr uint8_t kBlockQidxGroupBits = 8;
constexpr uint32_t kBlockQidxGroupSize = 1u << kBlockQidxGroupBits;
constexpr uint32_t kBlockQidxStride = kBlockQidxGroupSize + 1;

size_t compute_imcf_bin_size(const std::vector<chimera::imcf::Group> &groups,
                             double loadFactor) {
  uint64_t maxTotalHash = 0;
  for (const auto &group : groups) {
    maxTotalHash = std::max<uint64_t>(maxTotalHash, group.totalHash);
  }
  const uint64_t denom =
      std::max<uint64_t>(1, static_cast<uint64_t>(loadFactor * 4.0));
  size_t needBinSize = ceil_div_u64(maxTotalHash, denom);
  needBinSize = static_cast<size_t>(
      std::ceil(static_cast<double>(std::max<size_t>(1, needBinSize)) *
                1.05));
  return next_pow2(std::max<size_t>(1, needBinSize));
}

std::vector<std::vector<std::string>>
build_index_to_taxid(const std::vector<chimera::imcf::Group> &groups) {
  std::vector<std::vector<std::string>> indexToTaxid(groups.size());
  for (size_t idx = 0; idx < groups.size(); ++idx) {
    indexToTaxid[idx] = groups[idx].taxids;
  }
  return indexToTaxid;
}

std::vector<chimera::imcf::Group>
slice_groups(const std::vector<chimera::imcf::Group> &groups, size_t begin,
             size_t end) {
  return std::vector<chimera::imcf::Group>(groups.begin() + begin,
                                           groups.begin() + end);
}

robin_hood::unordered_flat_map<std::string, uint64_t>
make_shard_offset_base(const std::vector<chimera::imcf::Group> &groups,
                       size_t blockBegin) {
  robin_hood::unordered_flat_map<std::string, uint64_t> offsets;
  for (size_t groupIdx = 0; groupIdx < blockBegin; ++groupIdx) {
    const auto &group = groups[groupIdx];
    if (group.taxids.size() != group.assignedHashes.size()) {
      throw std::runtime_error(
          "IMCF block: group taxids/assignedHashes size mismatch");
    }
    for (size_t slot = 0; slot < group.taxids.size(); ++slot) {
      offsets[group.taxids[slot]] += group.assignedHashes[slot];
    }
  }
  return offsets;
}

void write_prefix_spool_and_bucket_base(
    const std::filesystem::path &prefixSpoolPath,
    const std::vector<uint32_t> &groupCounts, size_t hashSize,
    std::vector<uint64_t> &bucketBase, uint64_t &entriesCount,
    uint32_t &maxBucketTotal) {
  std::ofstream prefix(prefixSpoolPath, std::ios::binary);
  if (!prefix.is_open()) {
    throw std::runtime_error("QIMCF block: failed to open prefix spool: " +
                             prefixSpoolPath.string());
  }
  bucketBase.assign(hashSize + 1, 0);
  entriesCount = 0;
  maxBucketTotal = 0;
  std::vector<uint32_t> prefixRow(kBlockQidxStride, 0);
  for (size_t bucket = 0; bucket < hashSize; ++bucket) {
    bucketBase[bucket] = entriesCount;
    prefixRow[0] = 0;
    uint32_t running = 0;
    const uint64_t offset =
        static_cast<uint64_t>(bucket) * kBlockQidxGroupSize;
    for (uint32_t lo = 0; lo < kBlockQidxGroupSize; ++lo) {
      running += groupCounts[offset + lo];
      prefixRow[lo + 1] = running;
    }
    maxBucketTotal = std::max<uint32_t>(maxBucketTotal, running);
    entriesCount += running;
    prefix.write(reinterpret_cast<const char *>(prefixRow.data()),
                 static_cast<std::streamsize>(prefixRow.size() *
                                              sizeof(uint32_t)));
    if (!prefix.good()) {
      throw std::runtime_error("QIMCF block: failed to write prefix spool");
    }
  }
  bucketBase[hashSize] = entriesCount;
  prefix.close();
  if (!prefix.good()) {
    throw std::runtime_error("QIMCF block: failed to close prefix spool");
  }
}

int open_sized_entries_spool(const std::filesystem::path &entriesSpoolPath,
                             uint64_t entriesCount) {
  int fd = ::open(entriesSpoolPath.c_str(), O_CREAT | O_TRUNC | O_WRONLY | O_CLOEXEC,
                  0644);
  if (fd < 0) {
    throw std::runtime_error("QIMCF block: failed to open entries spool: " +
                             std::string(std::strerror(errno)));
  }
  const uint64_t bytes = entriesCount * sizeof(uint32_t);
  if (::ftruncate(fd, static_cast<off_t>(bytes)) != 0) {
    const int err = errno;
    ::close(fd);
    throw std::runtime_error("QIMCF block: failed to size entries spool: " +
                             std::string(std::strerror(err)));
  }
  return fd;
}

void close_entries_spool(int fd) {
  if (fd >= 0 && ::close(fd) != 0) {
    throw std::runtime_error("QIMCF block: failed to close entries spool: " +
                             std::string(std::strerror(errno)));
  }
}

std::vector<uint64_t> make_bucket_base_from_group_counts(
    const std::vector<uint32_t> &groupCounts, size_t hashSize,
    uint64_t &entriesCount) {
  const uint64_t expectedSize =
      static_cast<uint64_t>(hashSize) * kBlockQidxGroupSize;
  if (groupCounts.size() != expectedSize) {
    throw std::runtime_error("QIMCF block: invalid group count size");
  }
  std::vector<uint64_t> bucketBase(hashSize + 1, 0);
  entriesCount = 0;
  for (size_t bucket = 0; bucket < hashSize; ++bucket) {
    bucketBase[bucket] = entriesCount;
    const uint64_t offset =
        static_cast<uint64_t>(bucket) * kBlockQidxGroupSize;
    uint32_t running = 0;
    for (uint32_t lo = 0; lo < kBlockQidxGroupSize; ++lo) {
      running += groupCounts[offset + lo];
    }
    entriesCount += running;
  }
  bucketBase[hashSize] = entriesCount;
  return bucketBase;
}

void write_all_u32_fd(int fd, const uint32_t *values, uint64_t count) {
  const char *data = reinterpret_cast<const char *>(values);
  uint64_t bytesRemaining = count * sizeof(uint32_t);
  while (bytesRemaining > 0) {
    ssize_t rc = ::write(fd, data, static_cast<size_t>(bytesRemaining));
    if (rc < 0) {
      if (errno == EINTR) {
        continue;
      }
      throw std::runtime_error("QIMCF block: failed to write merged entries: " +
                               std::string(std::strerror(errno)));
    }
    if (rc == 0) {
      throw std::runtime_error("QIMCF block: merged entry write returned zero");
    }
    data += rc;
    bytesRemaining -= static_cast<uint64_t>(rc);
  }
}

void read_exact_u32(std::ifstream &input, const std::filesystem::path &path,
                    uint32_t *values, uint64_t count) {
  if (count == 0) {
    return;
  }
  const auto bytes = static_cast<std::streamsize>(count * sizeof(uint32_t));
  input.read(reinterpret_cast<char *>(values), bytes);
  if (input.gcount() != bytes) {
    throw std::runtime_error("QIMCF block: truncated block entries spool: " +
                             path.string());
  }
}

void merge_block_entry_spools(
    const std::vector<std::filesystem::path> &blockEntryPaths,
    const std::vector<std::vector<uint32_t>> &blockGroupCounts,
    size_t hashSize, int finalEntriesFd, uint64_t entriesCount) {
  const size_t blockCount = blockEntryPaths.size();
  if (blockGroupCounts.size() != blockCount) {
    throw std::runtime_error("QIMCF block: block entry merge input mismatch");
  }
  const uint64_t expectedGroupSize =
      static_cast<uint64_t>(hashSize) * kBlockQidxGroupSize;
  std::vector<std::ifstream> inputs;
  inputs.reserve(blockCount);
  for (size_t blockIdx = 0; blockIdx < blockCount; ++blockIdx) {
    if (blockGroupCounts[blockIdx].size() != expectedGroupSize) {
      throw std::runtime_error("QIMCF block: invalid block count vector");
    }
    inputs.emplace_back(blockEntryPaths[blockIdx], std::ios::binary);
    if (!inputs.back().is_open()) {
      throw std::runtime_error("QIMCF block: failed to open block entries: " +
                               blockEntryPaths[blockIdx].string());
    }
  }

  std::vector<std::vector<uint32_t>> bucketBuffers(blockCount);
  std::vector<size_t> bucketOffsets(blockCount, 0);
  std::vector<uint32_t> output;
  constexpr size_t kOutputFlushValues = 1u << 22;
  output.reserve(kOutputFlushValues);
  uint64_t writtenValues = 0;

  auto flush_output = [&]() {
    if (output.empty()) {
      return;
    }
    write_all_u32_fd(finalEntriesFd, output.data(), output.size());
    writtenValues += static_cast<uint64_t>(output.size());
    output.clear();
  };

  for (size_t bucket = 0; bucket < hashSize; ++bucket) {
    const uint64_t groupBase =
        static_cast<uint64_t>(bucket) * kBlockQidxGroupSize;
    for (size_t blockIdx = 0; blockIdx < blockCount; ++blockIdx) {
      uint64_t bucketTotal = 0;
      const auto &counts = blockGroupCounts[blockIdx];
      for (uint32_t lo = 0; lo < kBlockQidxGroupSize; ++lo) {
        bucketTotal += counts[groupBase + lo];
      }
      bucketBuffers[blockIdx].resize(bucketTotal);
      read_exact_u32(inputs[blockIdx], blockEntryPaths[blockIdx],
                     bucketBuffers[blockIdx].data(), bucketTotal);
      bucketOffsets[blockIdx] = 0;
    }

    for (uint32_t lo = 0; lo < kBlockQidxGroupSize; ++lo) {
      for (size_t blockIdx = 0; blockIdx < blockCount; ++blockIdx) {
        const uint32_t n = blockGroupCounts[blockIdx][groupBase + lo];
        if (n == 0) {
          continue;
        }
        const size_t offset = bucketOffsets[blockIdx];
        const auto &buffer = bucketBuffers[blockIdx];
        if (offset + n > buffer.size()) {
          throw std::runtime_error("QIMCF block: block entry merge overflow");
        }
        output.insert(output.end(), buffer.begin() + offset,
                      buffer.begin() + offset + n);
        bucketOffsets[blockIdx] += n;
        if (output.size() >= kOutputFlushValues) {
          flush_output();
        }
      }
    }
    for (size_t blockIdx = 0; blockIdx < blockCount; ++blockIdx) {
      if (bucketOffsets[blockIdx] != bucketBuffers[blockIdx].size()) {
        throw std::runtime_error("QIMCF block: block entry merge underflow");
      }
    }
  }
  flush_output();
  if (writtenValues != entriesCount) {
    throw std::runtime_error("QIMCF block: merged entries count mismatch");
  }
}

} // namespace

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
  robin_hood::unordered_flat_map<std::string, uint64_t> hashCount;
  robin_hood::unordered_flat_map<std::string, uint64_t> bpCount;
  FeatureBuildLayout featureLayout;
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
  feature_count(config, inputFiles, hashCount, fileInfo,
                hashFreqContext.enabled() ? &hashFreqContext : nullptr,
                &bpCount, &featureLayout);
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
  // Input path lists are only needed for feature extraction.
  robin_hood::unordered_flat_map<std::string, std::vector<std::string>>()
      .swap(inputFiles);
#if defined(__GLIBC__)
  ::malloc_trim(0);
#endif
  uint64_t imcf_feature_seed = 0;
  auto imcf_feature_params = make_feature_params(config, imcf_feature_seed);
  const uint16_t effective_span = static_cast<uint16_t>(
      std::min<size_t>(std::numeric_limits<uint16_t>::max(),
                       chimera::feature::min_required_length(
                           imcf_feature_params)));
  constexpr uint16_t ref_read_len = 150;
  chimera::presence::CoverageMeta presence_meta;
  auto partition_start = std::chrono::high_resolution_clock::now();
  constexpr int kDefaultMaxTaxidsPerGroup = 16;
  std::cout << "Partitioning hashcout..." << std::endl;
  std::vector<chimera::imcf::Group> groups =
      chimera::imcf::partitionHashCount(hashCount, kDefaultMaxTaxidsPerGroup,
                                        config.load_factor);
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
  std::cout << "Building IMCF in blocks..." << std::endl;
  IMCFConfig imcfConfig;
  imcfConfig.loadFactor = config.load_factor;
  imcfConfig.kmerSize = config.strobemer_k;
  imcfConfig.smerSize = 0;
  imcfConfig.syncmerPosition = 0;
  imcfConfig.seed64 = imcf_feature_seed;
  imcfConfig.fpSalt = IMCFConfig::DefaultFingerprintSalt;
  imcfConfig.hashVersion = IMCFConfig::CurrentHashVersion;
  imcfConfig.presenceUniqueDeg = config.presence_unique_deg;
  imcfConfig.featureMethod = 1;
  imcfConfig.strobeOrder = config.strobemer_order;
  imcfConfig.strobeWmin = config.strobemer_w_min;
  imcfConfig.strobeWmax = config.strobemer_w_max;
  imcfConfig.strobeK = config.strobemer_k;
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
  const size_t globalBinNum = groups.size();
  const size_t globalBinSize = compute_imcf_bin_size(groups, config.load_factor);
  imcfConfig.binNum = globalBinNum;
  imcfConfig.binSize = globalBinSize;
  const size_t blockBinCapacity = std::max<size_t>(
      1, static_cast<size_t>(kBlockClassicTableTargetBytes /
                             (static_cast<uint64_t>(globalBinSize) * 8ull)));
  const size_t blockCount =
      (globalBinNum + blockBinCapacity - 1) / blockBinCapacity;

  std::vector<std::vector<std::string>> indexToTaxid =
      build_index_to_taxid(groups);
  populateCoverageMeta(hashFreqContext.enabled() ? &hashFreqContext : nullptr,
                       &featureLayout, effective_span, ref_read_len,
                       config.presence_unique_deg, presence_meta);
  if (presence_meta.entries.empty()) {
    std::cout << "  coverage meta: empty (fallback to runtime heuristics)"
              << std::endl;
  } else {
    std::cout << "  coverage meta: " << presence_meta.entries.size()
              << " taxa cached" << std::endl;
  }

  const uint64_t groupCountSize =
      static_cast<uint64_t>(globalBinSize) * kBlockQidxGroupSize;
  std::vector<uint32_t> qidxGroupCounts(groupCountSize, 0u);
  std::vector<std::vector<uint32_t>> blockQidxCounts;
  blockQidxCounts.reserve(blockCount);
  std::vector<std::vector<uint64_t>> blockBucketBases(blockCount);
  std::vector<uint64_t> blockEntryCounts(blockCount, 0);
  for (size_t blockIdx = 0; blockIdx < blockCount; ++blockIdx) {
    const size_t blockBegin = blockIdx * blockBinCapacity;
    const size_t blockEnd = std::min(globalBinNum, blockBegin + blockBinCapacity);
    auto blockGroups = slice_groups(groups, blockBegin, blockEnd);
    auto shardOffsetBase = make_shard_offset_base(groups, blockBegin);
    chimera::imcf::InterleavedMergedCuckooFilter blockImcf;
    blockImcf.initialize_classic_storage(blockGroups.size(), globalBinSize);
    buildIMCF(blockImcf, blockGroups, hashCount, nullptr, &featureLayout,
              effective_span, ref_read_len, config.presence_unique_deg,
              nullptr, blockBegin, /*verifyShardTotals=*/false,
              &shardOffsetBase);
    std::vector<uint32_t> blockCounts(groupCountSize, 0u);
    blockImcf.accumulate_qidx_group_counts(blockCounts,
                                           /*include_stash=*/true);
    uint64_t blockEntryCount = 0;
    blockBucketBases[blockIdx] = make_bucket_base_from_group_counts(
        blockCounts, globalBinSize, blockEntryCount);
    blockEntryCounts[blockIdx] = blockEntryCount;
#pragma omp parallel for schedule(static)
    for (int64_t idx = 0; idx < static_cast<int64_t>(groupCountSize); ++idx) {
      qidxGroupCounts[static_cast<size_t>(idx)] +=
          blockCounts[static_cast<size_t>(idx)];
    }
    blockQidxCounts.push_back(std::move(blockCounts));
#if defined(__GLIBC__)
    ::malloc_trim(0);
#endif
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
    std::cout << "Interleaved Merged Cuckoo Filter" << std::endl;
    std::cout << "Bin number: " << globalBinNum << std::endl;
    std::cout << "Bin size: " << globalBinSize << std::endl;
  }
  if (config.verbose) {
    std::cout << "Building block QIMCF..." << std::endl;
  }
#if defined(__GLIBC__)
  ::malloc_trim(0);
#endif
  auto qidx_start = std::chrono::high_resolution_clock::now();
  const std::filesystem::path prefixSpoolPath = dir / "qidx_prefix.bin";
  const std::filesystem::path entriesSpoolPath = dir / "qidx_entries.bin";
  std::vector<uint64_t> bucketBase;
  uint64_t entriesCount = 0;
  uint32_t maxBucketTotal = 0;
  write_prefix_spool_and_bucket_base(prefixSpoolPath, qidxGroupCounts,
                                     globalBinSize, bucketBase, entriesCount,
                                     maxBucketTotal);
  const uint64_t prefixSpoolCount =
      static_cast<uint64_t>(globalBinSize) * kBlockQidxStride;
  std::vector<uint32_t>().swap(qidxGroupCounts);
  std::vector<std::filesystem::path> blockEntryPaths(blockCount);
  try {
    for (size_t blockIdx = 0; blockIdx < blockCount; ++blockIdx) {
      const size_t blockBegin = blockIdx * blockBinCapacity;
      const size_t blockEnd =
          std::min(globalBinNum, blockBegin + blockBinCapacity);
      auto blockGroups = slice_groups(groups, blockBegin, blockEnd);
      auto shardOffsetBase = make_shard_offset_base(groups, blockBegin);
      chimera::imcf::InterleavedMergedCuckooFilter blockImcf;
      blockImcf.initialize_classic_storage(blockGroups.size(), globalBinSize);
      buildIMCF(blockImcf, blockGroups, hashCount, nullptr, &featureLayout,
                effective_span, ref_read_len, config.presence_unique_deg,
                nullptr, blockBegin, /*verifyShardTotals=*/false,
                &shardOffsetBase);
      blockEntryPaths[blockIdx] =
          dir / ("qidx_block_entries." + std::to_string(blockIdx) + ".bin");
      int blockEntriesFd =
          open_sized_entries_spool(blockEntryPaths[blockIdx],
                                   blockEntryCounts[blockIdx]);
      uint64_t writtenBlockEntries = 0;
      try {
        writtenBlockEntries = blockImcf.write_qidx_block_entries_by_bucket(
            blockEntriesFd, blockBucketBases[blockIdx], blockBegin,
            /*include_stash=*/true);
        close_entries_spool(blockEntriesFd);
        blockEntriesFd = -1;
      } catch (...) {
        if (blockEntriesFd >= 0) {
          ::close(blockEntriesFd);
        }
        throw;
      }
      if (writtenBlockEntries != blockEntryCounts[blockIdx]) {
        throw std::runtime_error("QIMCF block: block entry count mismatch");
      }
#if defined(__GLIBC__)
      ::malloc_trim(0);
#endif
    }
    int entriesFd = open_sized_entries_spool(entriesSpoolPath, entriesCount);
    try {
      merge_block_entry_spools(blockEntryPaths, blockQidxCounts, globalBinSize,
                               entriesFd, entriesCount);
      close_entries_spool(entriesFd);
      entriesFd = -1;
    } catch (...) {
      if (entriesFd >= 0) {
        ::close(entriesFd);
      }
      throw;
    }
  } catch (...) {
    throw;
  }
  std::vector<std::vector<uint32_t>>().swap(blockQidxCounts);
  std::vector<std::vector<uint64_t>>().swap(blockBucketBases);
  const uint8_t prefixBits = bits_required_u64(maxBucketTotal);
  const uint8_t entryBits = static_cast<uint8_t>(
      bits_required_u64(globalBinNum > 0 ? (globalBinNum - 1) : 0) +
      (12 - kBlockQidxGroupBits) + 4);
  chimera::imcf::InterleavedMergedCuckooFilter imcf;
  imcf.initialize_qidx_spool_only(globalBinNum, globalBinSize,
                                  std::move(bucketBase), prefixBits, entryBits,
                                  prefixSpoolPath, prefixSpoolCount,
                                  entriesSpoolPath, entriesCount);
  robin_hood::unordered_flat_map<std::string, uint64_t>().swap(hashCount);
  robin_hood::unordered_flat_map<std::string, uint64_t>().swap(bpCount);
  std::vector<chimera::imcf::Group>().swap(groups);
  hashFreqContext.sketch.reset();
#if defined(__GLIBC__)
  ::malloc_trim(0);
#endif
  auto qidx_end = std::chrono::high_resolution_clock::now();
  auto qidx_total_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(qidx_end -
                                                            qidx_start)
          .count();
  if (config.verbose) {
    std::cout << "QIMCF build time: ";
    print_build_time(qidx_total_time);
    std::cout << std::endl;
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
