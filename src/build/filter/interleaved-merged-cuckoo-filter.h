/*
 * -----------------------------------------------------------------------------
 * Filename:      interleaved-merged-cuckoo-filter.h
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-11-12
 *
 * Last Modified: 2024-11-18
 *
 * Description:
 *  This is the header file of the Interleaved Merged Cuckoo Filter,
 *	which contains the basic operations of the Interleaved Merged Cuckoo
 *Filter
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <bitset>
#include <buildConfig.hpp>
#include <cassert>
#include <cereal/archives/binary.hpp>
#include <cereal/details/helpers.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <kvec.h>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <ranges>
#include <robin_hood.h>
#include <sdsl/int_vector.hpp>
#include <simde/x86/avx2.h>
#include <simde/x86/sse2.h>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#include <xxhash.h>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif


static inline size_t ceil_div_u64(uint64_t a, uint64_t b) {
  return (size_t)((a + b - 1) / b);
}

static inline size_t next_pow2(size_t v) {
  if (v <= 1)
    return 1;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
#if SIZE_MAX > 0xFFFFFFFFu
  v |= v >> 32;
#endif
  return v + 1;
}

static inline uint8_t bits_required_u64(uint64_t v) {
  if (v == 0) {
    return 1;
  }
  return static_cast<uint8_t>(64 - std::countl_zero(v));
}

namespace chimera::imcf {
static inline uint64_t mix64(uint64_t value) {
  return XXH3_64bits(&value, sizeof(value));
}

struct Group {
  std::vector<std::string> taxids;
  std::vector<uint64_t> assignedHashes;
  uint64_t totalHash{0};
};

struct HashChunk {
  std::string taxid;
  uint64_t hashCount;
};

/**
 * @brief Partitions hash counts into groups using a greedy algorithm for the
 * Interleaved Merged Cuckoo Filter (IMCF).
 *
 * This function partitions the provided hash counts into groups to balance the
 * total hash count in each group as evenly as possible, while adhering to a
 * maximum number of taxids per group. The goal is to optimize the distribution
 * for efficient processing within the IMCF.
 *
 * @param hashCount A map containing taxid strings as keys and their
 * corresponding hash counts as values.
 * @param maxTaxidsPerGroup The maximum number of taxids allowed in a single
 * group. Default is 16.
 *
 * @return A vector of Group objects, where each Group contains a collection of
 * taxids and the total hash count for that group.
 *
 * @details
 * The function first calculates the median of the hash counts and sets a
 * threshold as `median * 64`. This threshold helps determine if a taxid's hash
 * count should be split into smaller chunks to distribute the load more evenly.
 *
 * - Step 1: Calculate the median of hash counts and establish a threshold.
 * - Step 2: Split any hash count exceeding the threshold into smaller chunks
 * and add them to `hashChunks`.
 * - Step 3: Sort the `hashChunks` vector in descending order of hash counts.
 * - Step 4: Use a priority queue to greedily assign chunks to groups, ensuring
 * the group with the current lowest total hash count receives the next chunk.
 *
 * The function ensures that each group receives chunks in a balanced manner,
 * maintaining the load across all groups.
 *
 * @note
 * - The function assumes that `HashChunk` and `Group` types are defined.
 * `HashChunk` should have `taxid` and `hashCount` fields, while `Group` should
 * contain a `taxids` vector and `totalHash` field.
 * - The `minHeap` priority queue helps maintain the group with the smallest
 * total hash count for efficient chunk assignment.
 */
inline std::vector<Group> partitionHashCount(
    const robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    int maxTaxidsPerGroup = 16, double loadFactor = 0.85) {
  if (hashCount.empty()) {
    return {};
  }
  if (maxTaxidsPerGroup <= 0) {
    throw std::invalid_argument(
        "IMCF partition: maxTaxidsPerGroup must be positive");
  }
  // Keep behavior aligned with c060 baseline: this function optimizes
  // load-balance ratio only and does not use loadFactor in partitioning.
  (void)loadFactor;
  size_t maxTaxids = static_cast<size_t>(maxTaxidsPerGroup);

  // 1) 计算全局统计量，用更小的阈值做细粒度切分
  std::vector<uint64_t> counts;
  counts.reserve(hashCount.size());
  uint64_t totalHashes = 0;
  for (const auto &[taxid, count] : hashCount) {
    counts.push_back(count);
    totalHashes += count;
  }
  std::sort(counts.begin(), counts.end());
  const uint64_t median = counts[counts.size() / 2];
  const uint64_t mean =
      counts.empty() ? 0 : (totalHashes / static_cast<uint64_t>(counts.size()));

  const uint64_t baseMedian = median == 0 ? 1 : median;
  const uint64_t thresholdFromMedian =
      (baseMedian > std::numeric_limits<uint64_t>::max() / 4ull)
          ? std::numeric_limits<uint64_t>::max()
          : baseMedian * 4ull;
  const uint64_t thresholdFromMean = (mean == 0) ? thresholdFromMedian : mean;
  uint64_t threshold =
      std::max<uint64_t>(1, std::min<uint64_t>(thresholdFromMedian,
                                               thresholdFromMean));
  const double kTargetImbalance = 1.05; // 目标 M/A

  auto makeChunks = [&](uint64_t splitThreshold) {
    std::vector<HashChunk> chunks;
    chunks.reserve(hashCount.size() * 2);
    for (const auto &[taxid, count] : hashCount) {
      if (count == 0) {
        chunks.push_back({taxid, 0});
        continue;
      }
      if (count <= splitThreshold) {
        chunks.push_back({taxid, count});
        continue;
      }
      uint64_t numChunks = (count + splitThreshold - 1) / splitThreshold;
      numChunks = std::max<uint64_t>(1, numChunks);
      uint64_t chunkSize = (count + numChunks - 1) / numChunks;
      uint64_t assigned = 0;
      for (uint64_t i = 0; i < numChunks; ++i) {
        uint64_t remaining = count - assigned;
        uint64_t current = std::min<uint64_t>(chunkSize, remaining);
        if (i == numChunks - 1) {
          current = remaining;
        }
        if (current == 0) {
          break;
        }
        chunks.push_back({taxid, current});
        assigned += current;
      }
    }
    std::sort(chunks.begin(), chunks.end(),
              [](const HashChunk &a, const HashChunk &b) {
                return a.hashCount > b.hashCount;
              });
    return chunks;
  };

  auto packChunks = [&](size_t seedGroupCount,
                        const std::vector<HashChunk> &chunks) {
    std::vector<Group> groups(seedGroupCount);
    std::vector<robin_hood::unordered_flat_set<std::string>> used(seedGroupCount);

    auto cmp = [&](size_t a, size_t b) {
      return groups[a].totalHash > groups[b].totalHash;
    };
    std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)> minHeap(cmp);
    for (size_t i = 0; i < seedGroupCount; ++i) {
      minHeap.push(i);
    }

    for (const auto &chunk : chunks) {
      size_t target = std::numeric_limits<size_t>::max();
      std::vector<size_t> popped;
      while (!minHeap.empty()) {
        size_t cand = minHeap.top();
        minHeap.pop();
        if (groups[cand].taxids.size() >= maxTaxids ||
            used[cand].contains(chunk.taxid)) {
          popped.push_back(cand);
          continue;
        }
        target = cand;
        break;
      }

      if (target == std::numeric_limits<size_t>::max()) {
        groups.emplace_back();
        used.emplace_back();
        target = groups.size() - 1;
      }
      for (size_t idx : popped) {
        minHeap.push(idx);
      }

      groups[target].taxids.push_back(chunk.taxid);
      groups[target].assignedHashes.push_back(chunk.hashCount);
      groups[target].totalHash += chunk.hashCount;
      used[target].insert(chunk.taxid);
      minHeap.push(target);
    }

    std::vector<Group> compact;
    compact.reserve(groups.size());
    for (auto &g : groups) {
      if (!g.taxids.empty()) {
        compact.push_back(std::move(g));
      }
    }
    return compact;
  };

  size_t seedGroups =
      std::max<size_t>(1, (hashCount.size() + maxTaxids - 1) / maxTaxids);

  std::vector<Group> bestGroups;
  double bestRatio = std::numeric_limits<double>::infinity();
  uint64_t currentThreshold = threshold;

  for (int iter = 0; iter < 6; ++iter) {
    auto hashChunks = makeChunks(currentThreshold);
    auto groups = packChunks(seedGroups, hashChunks);
    if (groups.empty()) {
      return groups;
    }
    double maxLoad = 0.0;
    for (const auto &g : groups) {
      if (g.totalHash > maxLoad) {
        maxLoad = static_cast<double>(g.totalHash);
      }
    }
    double avgLoad =
        groups.empty()
            ? 0.0
            : static_cast<double>(totalHashes) /
                  static_cast<double>(groups.size());
    double ratio = (avgLoad == 0.0) ? 0.0 : (maxLoad / avgLoad);
    if (ratio < bestRatio) {
      bestRatio = ratio;
      bestGroups = groups;
    }
    if (ratio <= kTargetImbalance || currentThreshold <= 1) {
      break;
    }
    currentThreshold = std::max<uint64_t>(1, currentThreshold / 2);
  }

  if (bestGroups.empty()) {
    bestGroups = packChunks(seedGroups, makeChunks(threshold));
  }
  return bestGroups;
}

class InterleavedMergedCuckooFilter {
  typedef kvec_t(int) kvector;
  typedef kvec_t(bool) kvectorBool;
  sdsl::bit_vector data;
  size_t binNum;
  size_t binSize;
  size_t tagNum{4};
  int MaxCuckooCount{500};
  size_t hashSize{};
  struct StashEntry {
    uint64_t bucket{0};
    uint16_t fingerprint{0};
    uint16_t speciesMask{0};

    template <class Archive> void serialize(Archive &ar) {
      ar(bucket, fingerprint, speciesMask);
    }
  };
  struct QueryIndex {
    static constexpr uint32_t FP_BITS = 12;
    std::vector<uint64_t> bucketBase;  // size = hashSize + 1
    sdsl::int_vector<0> prefix;        // width = prefix_bits
    sdsl::int_vector<0> entries;       // width = entry_bits
    uint8_t g{8};                      // group bits for fp_lo
    uint32_t groupSize{0};             // 1 << g
    uint32_t stride{0};                // groupSize + 1
    uint8_t fpHiBits{0};               // 12 - g
    uint8_t prefix_bits{0};
    uint8_t entry_bits{0};
    std::filesystem::path prefixSpoolPath;
    std::filesystem::path entriesSpoolPath;
    uint64_t prefixSpoolCount{0};
    uint64_t entriesSpoolCount{0};
    bool spoolBacked{false};

    void refresh() {
      groupSize = 1u << g;
      stride = groupSize + 1;
      fpHiBits = static_cast<uint8_t>(FP_BITS - g);
    }

    void set_spool_backing(const std::filesystem::path &prefix_path,
                           uint64_t prefix_count,
                           const std::filesystem::path &entries_path,
                           uint64_t entries_count) {
      prefixSpoolPath = prefix_path;
      entriesSpoolPath = entries_path;
      prefixSpoolCount = prefix_count;
      entriesSpoolCount = entries_count;
      spoolBacked = true;
    }

    void clear_spool_backing(bool remove_files) {
      if (remove_files) {
        std::error_code ec;
        if (!prefixSpoolPath.empty()) {
          std::filesystem::remove(prefixSpoolPath, ec);
          ec.clear();
        }
        if (!entriesSpoolPath.empty()) {
          std::filesystem::remove(entriesSpoolPath, ec);
        }
      }
      prefixSpoolPath.clear();
      entriesSpoolPath.clear();
      prefixSpoolCount = 0;
      entriesSpoolCount = 0;
      spoolBacked = false;
    }

    bool has_spool_backing() const { return spoolBacked; }

    template <class Archive> void serialize(Archive &ar) {
      ar(bucketBase, prefix, entries, g, prefix_bits, entry_bits);
      if constexpr (Archive::is_loading::value) {
        refresh();
        clear_spool_backing(false);
      }
    }
  };

  struct StashFlat {
    uint32_t bucket{0};
    uint16_t fp{0};
    uint16_t mask{0};
    uint32_t bin{0};
  };

  std::vector<std::vector<StashEntry>> stash;
  std::unique_ptr<QueryIndex> qidx;
  uint8_t storageMode{0}; // 0=classic, 1=qidx-only, 2=both
  std::atomic<uint64_t> insertFailureTotal{0};
  std::atomic<uint64_t> insertFailureSaturated{0};

  struct BucketStats {
    size_t bucketCount{0};
    std::array<uint64_t, 5> occupancy{};
    double meanLoad{0.0};
    double stddevLoad{0.0};
    double percentFull{0.0};
    double maxBucketLoad{0.0};
    size_t lowBitBins{0};
    double lowBitMinRatio{0.0};
    double lowBitMaxRatio{0.0};
    uint64_t stashEntries{0};
    size_t stashEntryCount{0};
    size_t stashMaxPerEntry{0};
  };

  static inline uint16_t popcount16(uint16_t value) {
#if defined(__cpp_lib_bitops)
    return static_cast<uint16_t>(std::popcount(value));
#else
    return static_cast<uint16_t>(__builtin_popcount(static_cast<unsigned int>(value)));
#endif
  }

  inline bool insertIntoStash(size_t binIndex, size_t bucket, uint16_t fingerprint,
                              size_t speciesIndex) {
    if (binIndex >= stash.size()) {
      return false;
    }
    auto &entries = stash[binIndex];
    uint16_t speciesBit = static_cast<uint16_t>(1u << speciesIndex);
    for (auto &entry : entries) {
      if (entry.bucket == bucket && entry.fingerprint == fingerprint) {
        if ((entry.speciesMask & speciesBit) != 0) {
          return true;
        }
        if (entry.speciesMask == 0xFFFFu) {
          return false;
        }
        entry.speciesMask = static_cast<uint16_t>(entry.speciesMask | speciesBit);
        return true;
      }
    }
    entries.push_back(
        StashEntry{static_cast<uint64_t>(bucket), fingerprint, speciesBit});
    return true;
  }

  inline size_t bucketLinearIndex(size_t bucket, size_t bin) const {
    return bucket * binNum + bin;
  }

  inline size_t bucketBitOffset(size_t bucketLinear) const {
    return bucketLinear * tagNum * 16;
  }

  inline void loadBucketWords(size_t bucket, std::vector<uint64_t> &out) const {
    out.resize(binNum);
    if (binNum == 0) {
      return;
    }
    const uint64_t start = static_cast<uint64_t>(bucket) * static_cast<uint64_t>(binNum);
    std::memcpy(out.data(), data.data() + start, binNum * sizeof(uint64_t));
  }

  inline uint64_t readBucketRaw(size_t bucketLinear) const {
    // 桶布局为 4×16bit，总长度固定 64bit，且 bit_offset 总是 64 的倍数。
    // 直接按 word 访问可避免 sdsl::int_vector 的 get_int/set_int 抽象开销。
    return data.data()[bucketLinear];
  }

  inline void writeBucketRaw(size_t bucketLinear, uint64_t value) {
    data.data()[bucketLinear] = value;
  }

  inline uint64_t readBucket64(size_t bucket, size_t bin) const {
    size_t linear = bucketLinearIndex(bucket, bin);
    return readBucketRaw(linear);
  }

  inline void writeBucket64(size_t bucket, size_t bin, uint64_t value) {
    size_t linear = bucketLinearIndex(bucket, bin);
    writeBucketRaw(linear, value);
  }

  inline bool has_qidx() const {
    return static_cast<bool>(qidx);
  }

  template <class KeepFn>
  static inline void collect_subset_bins(uint32_t n, size_t reserve_hint,
                                         KeepFn &&keep,
                                         std::vector<uint32_t> &out) {
    out.clear();
    out.reserve(reserve_hint);
    for (uint32_t bin = 0; bin < n; ++bin) {
      if (keep(bin)) {
        out.push_back(bin);
      }
    }
  }

  inline void fetchActiveBins(size_t bucket, uint16_t fingerprint,
                              std::vector<uint32_t> &out) const {
    out.clear();
    if (binNum == 0 || hashSize == 0 || bucket >= hashSize) {
      return;
    }
    out.reserve(16);
    for (uint32_t bin = 0; bin < binNum; ++bin) {
      uint64_t bucketValue = readBucket64(bucket, bin);
      bool match = false;
      if (bucketValue != 0) {
        for (int lane = 0; lane < static_cast<int>(tagNum); ++lane) {
          uint16_t tag = static_cast<uint16_t>((bucketValue >> (lane * 16)) & 0xFFFFu);
          if (tag == 0u) {
            continue;
          }
          if ((tag & 0x0FFFu) == fingerprint) {
            match = true;
            break;
          }
        }
      }
      if (!match) {
        bool stashHit = false;
        forEachStashMatch(bin, bucket, fingerprint,
                          [&](uint16_t) { stashHit = true; });
        if (!stashHit) {
          continue;
        }
        match = true;
      }
      if (match) {
        out.push_back(bin);
      }
    }
  }

  inline uint64_t qidx_prefix(size_t bucket, uint32_t group) const {
    return qidx->prefix[(uint64_t)bucket * qidx->stride + group];
  }

  template <class Archive>
  static inline void save_spooled_int_vector(Archive &ar,
                                             const std::filesystem::path &path,
                                             uint64_t valueCount,
                                             uint8_t width,
                                             const char *label) {
    if (width == 0 || width > 64) {
      throw std::runtime_error(std::string("QIMCF stream-save invalid width for ") +
                               label);
    }
    const uint64_t expectedBytes = valueCount * sizeof(uint32_t);
    const uint64_t actualBytes = std::filesystem::file_size(path);
    if (actualBytes != expectedBytes) {
      throw std::runtime_error(std::string("QIMCF stream-save spool size mismatch for ") +
                               label + " (expected " +
                               std::to_string(expectedBytes) + ", got " +
                               std::to_string(actualBytes) + ")");
    }

    using dynamic_iv = sdsl::int_vector<0>;
    using int_width_type = typename dynamic_iv::int_width_type;
    const int_width_type widthTag = static_cast<int_width_type>(width);
    ar(cereal::make_size_tag(widthTag));
    float growthFactor = 1.5f;
    ar(growthFactor);
    const uint64_t bitSize = valueCount * static_cast<uint64_t>(width);
    ar(cereal::make_size_tag(bitSize));
    if (bitSize == 0) {
      return;
    }

    constexpr size_t kValueChunkU32 = 1u << 20;
    constexpr size_t kWordFlush = 1u << 18;
    std::ifstream spool(path, std::ios::binary);
    if (!spool.is_open()) {
      throw std::runtime_error(std::string("QIMCF stream-save failed to open spool for ") +
                               label);
    }

    std::vector<uint32_t> values(kValueChunkU32, 0u);
    std::vector<uint64_t> words;
    words.reserve(kWordFlush);
    uint64_t currentWord = 0;
    uint8_t bitsUsed = 0;
    uint64_t readValues = 0;
    uint64_t writtenWords = 0;
    const uint64_t wordCount = (bitSize + 63u) >> 6;

    auto flush_words = [&]() {
      if (words.empty()) {
        return;
      }
      ar(cereal::binary_data(words.data(), words.size() * sizeof(uint64_t)));
      writtenWords += static_cast<uint64_t>(words.size());
      words.clear();
    };

    while (readValues < valueCount) {
      const uint64_t need =
          std::min<uint64_t>(kValueChunkU32, valueCount - readValues);
      const auto needBytes =
          static_cast<std::streamsize>(need * sizeof(uint32_t));
      spool.read(reinterpret_cast<char *>(values.data()), needBytes);
      if (spool.gcount() != needBytes) {
        throw std::runtime_error(std::string("QIMCF stream-save truncated spool for ") +
                                 label);
      }
      for (uint64_t i = 0; i < need; ++i) {
        uint64_t value = values[i];
        uint8_t remaining = width;
        while (remaining > 0) {
          const uint8_t space = static_cast<uint8_t>(64 - bitsUsed);
          const uint8_t take = std::min<uint8_t>(space, remaining);
          const uint64_t mask =
              (take == 64) ? std::numeric_limits<uint64_t>::max()
                           : ((uint64_t{1} << take) - 1u);
          currentWord |= (value & mask) << bitsUsed;
          if (take == 64) {
            value = 0;
          } else {
            value >>= take;
          }
          bitsUsed = static_cast<uint8_t>(bitsUsed + take);
          remaining = static_cast<uint8_t>(remaining - take);
          if (bitsUsed == 64) {
            words.push_back(currentWord);
            currentWord = 0;
            bitsUsed = 0;
            if (words.size() >= kWordFlush) {
              flush_words();
            }
          }
        }
      }
      readValues += need;
    }

    if (bitsUsed != 0) {
      words.push_back(currentWord);
    }
    flush_words();
    if (writtenWords != wordCount) {
      throw std::runtime_error(std::string("QIMCF stream-save packed word mismatch for ") +
                               label + " (expected " +
                               std::to_string(wordCount) + ", got " +
                               std::to_string(writtenWords) + ")");
    }
  }

  template <class Archive>
  inline void save_qidx_only_from_spool(Archive &ar) {
    if (storageMode != 1 || !qidx || !qidx->has_spool_backing()) {
      throw std::runtime_error("QIMCF stream-save requested without spool-backed qidx");
    }
    ar(storageMode, binNum, binSize, tagNum, MaxCuckooCount, hashSize);
    ar(qidx->bucketBase);
    save_spooled_int_vector(ar, qidx->prefixSpoolPath, qidx->prefixSpoolCount,
                            qidx->prefix_bits, "prefix");
    save_spooled_int_vector(ar, qidx->entriesSpoolPath, qidx->entriesSpoolCount,
                            qidx->entry_bits, "entries");
    ar(qidx->g, qidx->prefix_bits, qidx->entry_bits);
  }

  struct QidxLookupState {
    size_t hash1{0};
    size_t hash2{0};
    uint32_t group{0};
    uint16_t hi{0};
    uint64_t hiMask{0};
  };

  inline QidxLookupState make_qidx_lookup_state(uint64_t value) const;
  inline void qidx_group_range(size_t bucket, uint32_t group, uint64_t &start,
                               uint64_t &end) const;

  static inline void encodeVarint(std::vector<uint8_t> &out, uint32_t value) {
    while (value >= 0x80) {
      out.push_back(static_cast<uint8_t>((value & 0x7F) | 0x80));
      value >>= 7;
    }
    out.push_back(static_cast<uint8_t>(value));
  }

  static inline bool decodeVarint(const std::vector<uint8_t> &dataBytes,
                                  size_t &offset, uint32_t &value) {
    uint32_t result = 0;
    int shift = 0;
    while (offset < dataBytes.size()) {
      uint8_t byte = dataBytes[offset++];
      result |= static_cast<uint32_t>(byte & 0x7F) << shift;
      if ((byte & 0x80) == 0) {
        value = result;
        return true;
      }
      shift += 7;
      if (shift > 28) {
        return false;
      }
    }
    return false;
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const InterleavedMergedCuckooFilter &filter) {
    os << "Interleaved Merged Cuckoo Filter" << std::endl;
    os << "Bin number: " << filter.binNum << std::endl;
    os << "Bin size: " << filter.binSize << std::endl;

    return os;
  }

public:
  InterleavedMergedCuckooFilter() = default;
  ~InterleavedMergedCuckooFilter() = default;
  InterleavedMergedCuckooFilter(std::vector<Group> &groups,
                                ChimeraBuild::IMCFConfig &config) {
    uint64_t maxTotalHash = 0;
    for (const auto &group : groups) {
      if (group.taxids.size() != group.assignedHashes.size()) {
        throw std::runtime_error(
            "IMCF constructor: group taxids/assignedHashes size mismatch");
      }
      uint64_t assignedSum =
          std::accumulate(group.assignedHashes.begin(),
                          group.assignedHashes.end(), uint64_t{0});
      if (assignedSum != group.totalHash) {
        throw std::runtime_error(
            "IMCF constructor: group totalHash inconsistent with shard sums");
      }
      if (group.totalHash > maxTotalHash) {
        maxTotalHash = group.totalHash;
      }
    }
    binNum = groups.size();

    // 关键：向上取整 + 2 的幂，便于用按位掩码实现可逆备桶
    size_t needBinSize = ceil_div_u64(
        maxTotalHash,
        static_cast<uint64_t>(config.loadFactor * static_cast<double>(tagNum)));
    constexpr double safetyCoeff = 1.05;
    needBinSize = static_cast<size_t>(
        std::ceil(static_cast<double>(needBinSize) * safetyCoeff));
    hashSize = next_pow2(std::max<size_t>(1, needBinSize));
    binSize = hashSize;

    // 保持旧布局：每桶 4×16bit，总位数 = binNum * binSize * tagNum * 16
    data = sdsl::bit_vector(
        (uint64_t)binNum * (uint64_t)binSize * (uint64_t)tagNum * 16ull, 0);

    config.binNum = binNum;
    config.binSize = binSize;
    stash.assign(binNum, {});
  }

  template <class Archive> void serialize(Archive &ar) {
    ar(storageMode, binNum, binSize, tagNum, MaxCuckooCount, hashSize);
    if (storageMode == 0) {
      ar(data, stash);
      if constexpr (Archive::is_loading::value) {
        if (stash.size() != binNum) {
          stash.resize(binNum);
        }
        qidx.reset();
      }
    } else if (storageMode == 1) {
      if constexpr (Archive::is_loading::value) {
        qidx = std::make_unique<QueryIndex>();
      }
      ar(*qidx);
      if constexpr (Archive::is_loading::value) {
        data = sdsl::bit_vector();
        stash.clear();
      }
    } else { // storageMode == 2 (classic + qidx)
      if constexpr (Archive::is_loading::value) {
        qidx = std::make_unique<QueryIndex>();
      }
      ar(data, stash, *qidx);
      if constexpr (Archive::is_loading::value) {
        if (stash.size() != binNum) {
          stash.resize(binNum);
        }
      }
    }
  }

  template <class Archive> inline void save_for_archive(Archive &ar) {
    if constexpr (Archive::is_saving::value) {
      if (storageMode == 1 && qidx && qidx->has_spool_backing()) {
        save_qidx_only_from_spool(ar);
        return;
      }
    }
    ar(*this);
  }

  inline bool has_spooled_qidx() const {
    return qidx && qidx->has_spool_backing();
  }

  inline void cleanup_qidx_spool_files() {
    if (qidx) {
      qidx->clear_spool_backing(true);
    }
  }

  /**
   * @brief 轻量主哈希索引。
   *
   * 使用若干位运算完成快速搅拌，并利用 `hashSize` 为 2
   * 的幂这一前提直接按位掩码。
   * 相比重型哈希大幅减少运算量，同时仍保持良好的分布特性。
   */
inline size_t hashIndex(uint64_t value) const {
  uint64_t hashed = mix64(value);
  return (size_t)hashed & (hashSize - 1);
}

  /**
   * @brief 可逆备桶哈希。
   *
   * 基于 fingerprint 乘黄金常数得到的扰动值与主桶索引 XOR，并依赖 2^k
   * 大小做掩码， 保证 `altHash(altHash(b, fp), fp) ==
   * b`，便于在踢出链中双向定位。
   */
inline size_t altHash(size_t b, uint16_t fingerprint) const {
  uint64_t fp64 = (uint64_t)fingerprint ^ 0x9E3779B97F4A7C15ull;
  size_t alt = (b ^ (size_t)mix64(fp64)) & (hashSize - 1);
  if (alt == b) {
    alt = (b + 1) & (hashSize - 1);
  }
  return alt;
}

  /**
   * @brief 组合 12 bit 指纹与物种索引。
   *
   * 步骤：
   * 1. 调用 `reduceTo12bit` 生成非零的 12 bit 指纹。
   * 2. 将物种索引左移 12 位并与指纹按位或，构造完整标签。
   *
   * @param value 用于生成指纹的值（通常是 minimizer 或 k-mer）。
   * @param index 需要编码进高 4 bit 的物种索引，取值范围 [0, 15]。
   *
   * @return 编码好物种索引与指纹的 16 bit 标签。
   */
	inline uint16_t reduceTo12bitAndAddIndex(uint64_t value, size_t index) const {
    assert(index < 16);
    uint16_t reduced_value = reduceTo12bit(value);
    uint16_t combined = (index << 12) | reduced_value;
    return combined;
  }

  /**
   * @brief 轻量 12 bit 指纹化。
   *
   * 复用主哈希的位混合策略，只截取低 12 bit，并保证不返回 0。
   */
	inline uint16_t reduceTo12bit(uint64_t value) const {
		uint64_t mixed = mix64((uint64_t)value ^ ChimeraBuild::IMCFConfig::DefaultFingerprintSalt);
		uint16_t v = (uint16_t)(mixed & 0x0FFFu);
		return v ? v : 1;
	}

	inline void route(uint64_t value, std::vector<uint32_t> &bins) const {
    if (has_qidx()) {
      route_qidx(value, bins);
      return;
    }
    route_scan(value, bins);
  }

  inline void route_scan(uint64_t value, std::vector<uint32_t> &bins) const {
    bins.clear();

    uint16_t fingerprint = reduceTo12bit(value);
    size_t hash1 = hashIndex(value);
    size_t hash2 = altHash(hash1, fingerprint);

    static thread_local std::vector<uint32_t> routeBuf1;
    static thread_local std::vector<uint32_t> routeBuf2;

    fetchActiveBins(hash1, fingerprint, routeBuf1);
    if (hash2 != hash1) {
      fetchActiveBins(hash2, fingerprint, routeBuf2);
    } else {
      routeBuf2.clear();
    }

    if (routeBuf1.empty() && routeBuf2.empty()) {
      return;
    }

    bins.reserve(routeBuf1.size() + routeBuf2.size());
    size_t i = 0;
    size_t j = 0;
    while (i < routeBuf1.size() || j < routeBuf2.size()) {
      uint32_t valueOut;
      if (j >= routeBuf2.size() ||
          (i < routeBuf1.size() && routeBuf1[i] < routeBuf2[j])) {
        valueOut = routeBuf1[i++];
      } else if (i >= routeBuf1.size() || routeBuf2[j] < routeBuf1[i]) {
        valueOut = routeBuf2[j++];
      } else {
        valueOut = routeBuf1[i];
        ++i;
        ++j;
      }
      if (bins.empty() || bins.back() != valueOut) {
        bins.push_back(valueOut);
      }
    }
  }

  inline void route_qidx(uint64_t value, std::vector<uint32_t> &bins) const;

  /**
   * @brief 插入指纹。
   *
   * 构造标签后先命中主桶，再跳到可逆备桶，两桶任意存在空位或已有指纹即视为成功；
   * 若均告满载则进入踢出流程。
   */
	inline bool insertTag(size_t binIndex, uint64_t value, size_t index) {
    if (index >= 16) {
      assert(false &&
             "IMCF group index overflow (taxids per group must be <=16)");
      return false;
    }
    uint16_t tag = reduceTo12bitAndAddIndex(value, index);
    uint16_t fp = (uint16_t)(tag & 0x0FFFu);

    size_t b1 = hashIndex(value);
    size_t b2 = altHash(b1, fp);

    auto try_insert_bucket = [&](size_t b) -> bool {
      uint64_t q = readBucket64(b, binIndex);
      for (int i = 0; i < 4; ++i) {
        uint16_t chunk = (uint16_t)((q >> (i * 16)) & 0xFFFFu);
        if (chunk == 0u) {
          q &= ~((uint64_t)0xFFFFu << (i * 16));
          q |= ((uint64_t)tag << (i * 16));
          writeBucket64(b, binIndex, q);
          return true;
        }
        if (chunk == tag)
          return true;
      }
      return false;
    };

    if (try_insert_bucket(b1))
      return true;
    if (try_insert_bucket(b2))
      return true;

    auto countMatching = [&](uint64_t bucketValue) -> int {
      int matches = 0;
      for (int i = 0; i < 4; ++i) {
        uint16_t chunk = (uint16_t)((bucketValue >> (i * 16)) & 0xFFFFu);
        if ((chunk & 0x0FFFu) == fp) {
          ++matches;
        }
      }
      return matches;
    };
    int sameFingerprint = countMatching(readBucket64(b1, binIndex)) +
                          countMatching(readBucket64(b2, binIndex));
    if (sameFingerprint >= 8) {
      bool stored = insertIntoStash(binIndex, b1, fp, index);
      if (!stored && b2 != b1) {
        stored = insertIntoStash(binIndex, b2, fp, index);
      }
      if (stored) {
        return true;
      }
      insertFailureTotal.fetch_add(1, std::memory_order_relaxed);
      insertFailureSaturated.fetch_add(1, std::memory_order_relaxed);
      return false;
    }

    bool kicked = kickOut(binIndex, value, tag);
    if (!kicked) {
      insertFailureTotal.fetch_add(1, std::memory_order_relaxed);
    }
    return kicked;
  }

  inline uint64_t getInsertFailureTotal() const {
    return insertFailureTotal.load(std::memory_order_relaxed);
  }

  inline uint64_t getInsertFailureSaturatedFingerprint() const {
    return insertFailureSaturated.load(std::memory_order_relaxed);
  }

  BucketStats computeBucketStats(size_t binIndex,
                                 size_t lowBitBins = 256) const {
    if (binIndex >= binNum) {
      throw std::out_of_range("IMCF bucket stats: binIndex out of range");
    }
    if (lowBitBins == 0 || (lowBitBins & (lowBitBins - 1)) != 0) {
      throw std::invalid_argument(
          "IMCF bucket stats: lowBitBins must be power of two");
    }
    BucketStats stats;
    stats.bucketCount = hashSize;
    stats.lowBitBins = lowBitBins;
    double sum = 0.0;
    double sumSq = 0.0;
    uint64_t maxLoad = 0;
    std::vector<uint64_t> lowBitLoad(lowBitBins, 0);
    std::vector<uint64_t> lowBitBucketCount(lowBitBins, 0);
    size_t mask = lowBitBins - 1;
    for (size_t bucket = 0; bucket < hashSize; ++bucket) {
      uint64_t raw = readBucket64(bucket, binIndex);
      uint64_t load = 0;
      for (int i = 0; i < static_cast<int>(tagNum); ++i) {
        uint16_t chunk = (uint16_t)((raw >> (i * 16)) & 0xFFFFu);
        if (chunk != 0u) {
          ++load;
        }
      }
      if (load <= 4) {
        stats.occupancy[load]++;
      }
      sum += static_cast<double>(load);
      sumSq += static_cast<double>(load) * static_cast<double>(load);
      if (load > maxLoad) {
        maxLoad = load;
      }
      size_t lowIdx = bucket & mask;
      lowBitLoad[lowIdx] += load;
      lowBitBucketCount[lowIdx]++;
    }
    if (stats.bucketCount > 0) {
      stats.meanLoad = sum / static_cast<double>(stats.bucketCount);
      double variance = sumSq / static_cast<double>(stats.bucketCount) -
                        stats.meanLoad * stats.meanLoad;
      if (variance < 0.0) {
        variance = 0.0;
      }
      stats.stddevLoad = std::sqrt(variance);
      stats.percentFull = static_cast<double>(stats.occupancy[4]) /
                          static_cast<double>(stats.bucketCount);
    }
    stats.maxBucketLoad = static_cast<double>(maxLoad);
    double minRatio = std::numeric_limits<double>::infinity();
    double maxRatio = 0.0;
    for (size_t idx = 0; idx < lowBitBins; ++idx) {
      if (lowBitBucketCount[idx] == 0) {
        continue;
      }
      double avg = static_cast<double>(lowBitLoad[idx]) /
                   static_cast<double>(lowBitBucketCount[idx]);
      double ratio = stats.meanLoad > 0.0 ? avg / stats.meanLoad : 0.0;
      if (ratio > maxRatio) {
        maxRatio = ratio;
      }
      if (ratio < minRatio) {
        minRatio = ratio;
      }
    }
    if (!std::isfinite(minRatio)) {
      minRatio = 0.0;
    }
    stats.lowBitMinRatio = minRatio;
    stats.lowBitMaxRatio = maxRatio;
    if (binIndex < stash.size()) {
      const auto &entries = stash[binIndex];
      if (!entries.empty()) {
        stats.stashEntryCount = entries.size();
        uint64_t totalSpecies = 0;
        size_t maxPerEntry = 0;
        for (const auto &entry : entries) {
          uint16_t count = popcount16(entry.speciesMask);
          totalSpecies += count;
          if (count > maxPerEntry) {
            maxPerEntry = count;
          }
        }
        stats.stashEntries = totalSpecies;
        stats.stashMaxPerEntry = maxPerEntry;
      }
    }
    return stats;
  }

  inline void release_classic_storage() {
    sdsl::bit_vector().swap(data);
    stash.clear();
    stash.shrink_to_fit();
  }

  inline void set_storage_mode(uint8_t mode) {
    storageMode = mode;
  }

  inline uint8_t get_storage_mode() const {
    return storageMode;
  }

  inline void build_query_index(bool include_stash = true, bool verify = true,
                                bool low_peak_mode = false,
                                bool drop_classic_before_materialize = false);

  /**
   * @brief 踢出流程，只在两桶之间往返。
   *
   * 随机逐出当前桶一项换入新项，再借助可逆哈希跳往另一桶继续尝试，最多迭代
   * `MaxCuckooCount` 次。
   */
	inline bool kickOut(size_t binIndex, uint64_t value, uint16_t tag) {
    uint16_t cur = tag;
    uint16_t fp = (uint16_t)(cur & 0x0FFFu);
    size_t b = hashIndex(value);

    static thread_local std::mt19937_64 gen{std::random_device{}()};
    std::uniform_int_distribution<int> dis(0, (int)tagNum - 1);

    for (int cnt = 0; cnt < MaxCuckooCount; ++cnt) {
      uint64_t q = readBucket64(b, binIndex);

      for (int i = 0; i < 4; ++i) {
        uint16_t chunk = (uint16_t)((q >> (i * 16)) & 0xFFFFu);
        if (chunk == 0u) {
          q &= ~((uint64_t)0xFFFFu << (i * 16));
          q |= ((uint64_t)cur << (i * 16));
          writeBucket64(b, binIndex, q);
          return true;
        }
        if (chunk == cur)
          return true;
      }

      int rp = dis(gen);
      uint16_t victim = (uint16_t)((q >> (rp * 16)) & 0xFFFFu);
      q &= ~((uint64_t)0xFFFFu << (rp * 16));
      q |= ((uint64_t)cur << (rp * 16));
      writeBucket64(b, binIndex, q);

      cur = victim;
      fp = (uint16_t)(cur & 0x0FFFu);
      b = altHash(b, fp);
    }
    return false;
  }

  template <typename Emit>
  inline void forEachStashMatch(size_t binIndex, size_t bucket,
                                uint16_t fingerprint, Emit &&emit) const {
    if (binIndex >= stash.size()) {
      return;
    }
    const auto &entries = stash[binIndex];
    if (entries.empty()) {
      return;
    }
    for (const auto &entry : entries) {
      if (entry.bucket != bucket || entry.fingerprint != fingerprint) {
        continue;
      }
      unsigned int mask = entry.speciesMask;
      while (mask) {
        auto species = static_cast<uint16_t>(std::countr_zero(mask));
        emit(species);
        mask &= (mask - 1);
      }
      return;
    }
  }

  /**
   * @brief Checks for the presence of a given minimizer hash in the Interleaved
   * Merged Cuckoo Filter and updates the result vector.
   *
   * This function uses the SIMDe library to accelerate the process of checking
   * whether a given minimizer hash exists in the Interleaved Merged Cuckoo
   * Filter (IMCF). The hash value is reduced to 12 bits, and two hash positions
   * are checked for matches. If a match is found, the function extracts the
   * position index from the upper 4 bits and updates the result vector.
   *
   * @param value The input minimizer hash value to be checked.
   * @param result A vector of std::bitset<16> where each bitset corresponds to
   * a bin index. Each bit represents whether a particular position index
   * (species index) has been found in the filter.
   *
   * @details
   * The function first reduces the input value to a 12-bit fingerprint using
   * the reduceTo12bit function. It then computes two hash positions (hash1 and
   * hash2) to identify potential bucket locations in the filter. The function
   * iterates over each bin index, checks the data in these bucket positions,
   * and loads them into SIMD vectors for comparison.
   *
   * - SIMDe Vectorization:
   *   The function uses SIMD vectors to perform parallel comparisons:
   *   - fingerprint_vec: Contains the 12-bit fingerprint replicated across all
   * vector positions.
   *   - fingerprint_mask: A mask used to extract the 12-bit fingerprint from
   * each 16-bit entry.
   *   - species_shift: The number of bits to shift to extract the species index
   * from the upper 4 bits of each entry.
   *
   * The function updates the result vector by setting bits corresponding to
   * matching species indices.
   *
   * @note Ensure the data object provides a method to retrieve 64-bit data
   * (e.g., data.get_int). The function requires the SIMDe library for SIMD
   * operations.
   */
	inline void bulkContain(uint64_t value, std::vector<std::bitset<16>> &result) {
    // Step 1: Reduce the minimizer hash to a 12-bit fingerprint
    uint16_t fingerprint = reduceTo12bit(value);
    // Step 2: Calculate primary and alternative hash positions
    size_t hash1 = hashIndex(value);
    size_t hash2 = altHash(hash1, fingerprint);

    // Step 3: Prepare SIMD vectors for comparison (128-bit lanes for 8 entries)
    simde__m128i fingerprint_vec = simde_mm_set1_epi16(fingerprint);
    simde__m128i fingerprint_mask = simde_mm_set1_epi16(0x0FFF);
    const int species_shift = 12;

    // Step 4: Iterate over each bin index to check for matches
    for (size_t binIndex = 0; binIndex < binNum; binIndex++) {
      uint64_t bucketData1 = readBucket64(hash1, binIndex);
      uint64_t bucketData2 = readBucket64(hash2, binIndex);

      // Load entries into an array for SIMD processing
      uint16_t entries[8];
      entries[0] = (bucketData1 >> 0) & 0xFFFF;
      entries[1] = (bucketData1 >> 16) & 0xFFFF;
      entries[2] = (bucketData1 >> 32) & 0xFFFF;
      entries[3] = (bucketData1 >> 48) & 0xFFFF;

      entries[4] = (bucketData2 >> 0) & 0xFFFF;
      entries[5] = (bucketData2 >> 16) & 0xFFFF;
      entries[6] = (bucketData2 >> 32) & 0xFFFF;
      entries[7] = (bucketData2 >> 48) & 0xFFFF;

      // Load data into SIMD vectors (8 x 16-bit = 128 bits)
      simde__m128i bucket_vec =
          simde_mm_loadu_si128(reinterpret_cast<const simde__m128i *>(entries));

      // Step 5: Extract fingerprints and species indices
      simde__m128i fingerprints =
          simde_mm_and_si128(bucket_vec, fingerprint_mask);
      simde__m128i species_indices =
          simde_mm_srli_epi16(bucket_vec, species_shift);

      // Compare extracted fingerprints with the input fingerprint
      simde__m128i cmp = simde_mm_cmpeq_epi16(fingerprints, fingerprint_vec);

      // Step 6: Update the result vector with matching species indices
      // Store comparison results and species indices to arrays for variable
      // indexing
      uint16_t cmp_vals[8];
      uint16_t species_vals[8];
      simde_mm_storeu_si128(reinterpret_cast<simde__m128i *>(cmp_vals), cmp);
      simde_mm_storeu_si128(reinterpret_cast<simde__m128i *>(species_vals),
                            species_indices);

      for (int i = 0; i < 8; ++i) {
        if (cmp_vals[i] == 0xFFFF) {
          uint16_t speciesIndex = species_vals[i];
          result[binIndex].set(speciesIndex);
        }
      }
      forEachStashMatch(binIndex, hash1, fingerprint,
                        [&](uint16_t speciesIndex) {
                          result[binIndex].set(speciesIndex);
                        });
      if (hash2 != hash1) {
        forEachStashMatch(binIndex, hash2, fingerprint,
                          [&](uint16_t speciesIndex) {
                            result[binIndex].set(speciesIndex);
                          });
      }
    }
  }

  /**
   * @brief Performs bulk counting of input minimizers and updates the result
   * based on their presence in the filter.
   *
   * This function iterates over a range of input minimizers and checks their
   * presence in the Interleaved Merged Cuckoo Filter (IMCF). For each
   * minimizer, it performs a bulk containment query and updates the `result`
   * vector with the counts based on the corresponding positions in the filter.
   *
   * @tparam value_range_t A range type representing the collection of
   * minimizers to be checked.
   * @param values The input range of minimizers to query in the IMCF.
   * @param result A 2D vector that stores the count results. Each sub-vector
   * corresponds to a specific position index, and each element in the
   * sub-vector represents the count for a particular tag at that index.
   *
   * @details
   * The function uses a temporary vector of `std::bitset<16>` to store the
   * intermediate results of the containment check for each minimizer. For each
   * input value:
   * - The `bitset` vector is reset to ensure no residual data from previous
   * iterations.
   * - The `bulkContain` function is called to determine the presence of the
   * minimizer across multiple positions.
   * - The results in `tmpResult` are analyzed, and the corresponding positions
   * in the `result` vector are updated:
   *   - If a bit is set in `tmpResult`, the respective count in `result` is
   * incremented.
   *   - If a sub-vector in `result` does not have enough space for a new tag,
   * it is resized and initialized to zero.
   *
   * @note
   * Ensure the `result` vector is pre-initialized with sub-vectors
   * corresponding to the positions to be checked. The `bulkContain` function
   * should be defined to handle the minimizer containment query.
   */
  template <std::ranges::range value_range_t>
  inline void bulkCount(value_range_t &&values,
                        std::vector<std::vector<size_t>> &result) {
    bulkCount_sparse_compat(std::forward<value_range_t>(values), result);
  }

  template <class EmitFn>
  inline void bulkContain_events(uint64_t value, EmitFn &&emit) {
    if (has_qidx()) {
      bulkContain_events_qidx(value, std::forward<EmitFn>(emit));
      return;
    }
    bulkContain_events_scan(value, std::forward<EmitFn>(emit));
  }

  template <class EmitFn>
	inline void bulkContain_events_scan(uint64_t value, EmitFn &&emit) {
    uint16_t fingerprint = reduceTo12bit(value);
    size_t hash1 = hashIndex(value);
    size_t hash2 = altHash(hash1, fingerprint);

    const uint64_t fpBroadcast =
        static_cast<uint64_t>(fingerprint) * 0x0001000100010001ULL;
    const uint64_t fpMask = 0x0FFF0FFF0FFF0FFFULL;
    const uint64_t idxMask = 0x000F000F000F000FULL;

    constexpr size_t BLOCK = 16;
    constexpr size_t SUB_BATCH = 4;

    auto process64 = [&](uint64_t q, uint32_t binIndex) {
      uint64_t diff = (q ^ fpBroadcast) & fpMask;
      uint64_t idxBits = (q >> 12) & idxMask;
      for (int lane = 0; lane < 4; ++lane) {
        if (((diff >> (lane * 16)) & 0x0FFFULL) == 0ULL) {
          uint16_t speciesIndex =
              static_cast<uint16_t>((idxBits >> (lane * 16)) & 0xF);
          emit(static_cast<uint32_t>(binIndex), speciesIndex);
        }
      }
    };

    for (size_t base = 0; base < binNum; base += BLOCK) {
      size_t end = std::min(binNum, base + BLOCK);
      for (size_t batchStart = base; batchStart < end;
           batchStart += SUB_BATCH) {
        size_t batchEnd = std::min(end, batchStart + SUB_BATCH);
        uint64_t q1[SUB_BATCH];
        uint64_t q2[SUB_BATCH];
        size_t bins[SUB_BATCH];
        size_t cnt = 0;
        for (size_t binIndex = batchStart; binIndex < batchEnd; ++binIndex) {
          q1[cnt] = readBucket64(hash1, binIndex);
          q2[cnt] = readBucket64(hash2, binIndex);
          bins[cnt] = binIndex;
          ++cnt;
        }

        for (size_t i = 0; i < cnt; ++i) {
          uint32_t bin = static_cast<uint32_t>(bins[i]);
          if (q1[i] != 0) {
            process64(q1[i], bin);
          }
          if (q2[i] != 0) {
            process64(q2[i], bin);
          }
          forEachStashMatch(bin, hash1, fingerprint,
                            [&](uint16_t speciesIndex) {
                              emit(bin, speciesIndex);
                            });
          if (hash2 != hash1) {
            forEachStashMatch(bin, hash2, fingerprint,
                              [&](uint16_t speciesIndex) {
                                emit(bin, speciesIndex);
                              });
          }
        }
      }
    }
  }

  template <class EmitFn>
  inline void bulkContain_events_qidx(uint64_t value, EmitFn &&emit);

  template <class EmitFn>
	inline void bulkContain_events_subset(uint64_t value,
                                        const std::vector<uint32_t> &binSubset,
                                        EmitFn &&emit) {
    if (has_qidx()) {
      bulkContain_events_subset_qidx(value, binSubset, std::forward<EmitFn>(emit));
      return;
    }
    bulkContain_events_subset_scan(value, binSubset, std::forward<EmitFn>(emit));
  }

  template <class EmitFn>
  inline void bulkContain_events_subset_mask(
      uint64_t value, const std::vector<uint8_t> &binMask, EmitFn &&emit) {
    if (has_qidx()) {
      bulkContain_events_subset_qidx_mask(value, binMask,
                                          std::forward<EmitFn>(emit));
      return;
    }
    static thread_local std::vector<uint32_t> subsetBins;
    collect_subset_bins(
        static_cast<uint32_t>(binMask.size()), binMask.size(),
        [&](uint32_t bin) { return binMask[bin] != 0u; }, subsetBins);
    bulkContain_events_subset_scan(value, subsetBins, std::forward<EmitFn>(emit));
  }

  template <class EmitFn>
  inline void bulkContain_events_subset_marked(
      uint64_t value, const std::vector<uint32_t> &binMarks,
      uint32_t activeMark, EmitFn &&emit) {
    if (activeMark == 0u || binMarks.empty()) {
      return;
    }
    if (has_qidx()) {
      bulkContain_events_subset_qidx_marked(
          value, binMarks, activeMark, std::forward<EmitFn>(emit));
      return;
    }
    static thread_local std::vector<uint32_t> subsetBins;
    const size_t reserve_hint = std::min<size_t>(binMarks.size(), 512);
    collect_subset_bins(static_cast<uint32_t>(binMarks.size()), reserve_hint,
                        [&](uint32_t bin) { return binMarks[bin] == activeMark; },
                        subsetBins);
    bulkContain_events_subset_scan(value, subsetBins,
                                   std::forward<EmitFn>(emit));
  }

  template <class EmitFn>
	inline void bulkContain_events_subset_scan(uint64_t value,
                                        const std::vector<uint32_t> &binSubset,
                                        EmitFn &&emit) {
    uint16_t fingerprint = reduceTo12bit(value);
    size_t hash1 = hashIndex(value);
    size_t hash2 = altHash(hash1, fingerprint);

    const uint64_t fpBroadcast =
        static_cast<uint64_t>(fingerprint) * 0x0001000100010001ULL;
    const uint64_t fpMask = 0x0FFF0FFF0FFF0FFFULL;
    const uint64_t idxMask = 0x000F000F000F000FULL;

    constexpr size_t SUB_BATCH = 4;

    auto process64 = [&](uint64_t q, uint32_t binIndex) {
      uint64_t diff = (q ^ fpBroadcast) & fpMask;
      uint64_t idxBits = (q >> 12) & idxMask;
      for (int lane = 0; lane < 4; ++lane) {
        if (((diff >> (lane * 16)) & 0x0FFFULL) == 0ULL) {
          uint16_t speciesIndex =
              static_cast<uint16_t>((idxBits >> (lane * 16)) & 0xF);
          emit(static_cast<uint32_t>(binIndex), speciesIndex);
        }
      }
    };

    for (size_t offset = 0; offset < binSubset.size(); offset += SUB_BATCH) {
      size_t batchEnd = std::min(binSubset.size(), offset + SUB_BATCH);
      uint64_t q1[SUB_BATCH];
      uint64_t q2[SUB_BATCH];
      uint32_t bins[SUB_BATCH];
      size_t cnt = 0;
      for (size_t i = offset; i < batchEnd; ++i) {
        uint32_t binIndex = binSubset[i];
        if (binIndex >= binNum)
          continue;
        bins[cnt] = binIndex;
        q1[cnt] = readBucket64(hash1, binIndex);
        q2[cnt] = readBucket64(hash2, binIndex);
        ++cnt;
      }
      for (size_t i = 0; i < cnt; ++i) {
        uint32_t bin = bins[i];
        if (q1[i] != 0)
          process64(q1[i], bin);
        forEachStashMatch(bin, hash1, fingerprint,
                          [&](uint16_t speciesIndex) {
                            emit(bin, speciesIndex);
                          });
        if (q2[i] != 0)
          process64(q2[i], bin);
        if (hash2 != hash1) {
          forEachStashMatch(bin, hash2, fingerprint,
                            [&](uint16_t speciesIndex) {
                              emit(bin, speciesIndex);
                            });
        }
      }
    }
  }

  template <class EmitFn>
  inline void bulkContain_events_subset_qidx(
      uint64_t value, const std::vector<uint32_t> &binSubset, EmitFn &&emit);

  template <class EmitFn>
  inline void bulkContain_events_subset_qidx_mask(
      uint64_t value, const std::vector<uint8_t> &binMask, EmitFn &&emit);

  template <class EmitFn>
  inline void bulkContain_events_subset_qidx_marked(
      uint64_t value, const std::vector<uint32_t> &binMarks,
      uint32_t activeMark, EmitFn &&emit);

  template <std::ranges::range value_range_t, class CounterMatrix>
  inline void bulkCount_sparse(
      value_range_t &&values, CounterMatrix &result,
      std::vector<std::pair<uint32_t, uint16_t>> *touched = nullptr) {
    using PairCode = uint64_t;
    static thread_local std::vector<PairCode> uniq;

    using CounterRow = typename CounterMatrix::value_type;
    using Counter = typename CounterRow::value_type;
    static_assert(std::is_integral_v<Counter>,
                  "IMCF counter type must be integral");
    for (auto value : values) {
      uniq.clear();
      bulkContain_events(value, [&](uint32_t bin, uint16_t sp) {
        uniq.push_back((PairCode(bin) << 16) | PairCode(sp));
      });
      if (uniq.empty()) {
        continue;
      }
      std::sort(uniq.begin(), uniq.end());
      uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

      for (PairCode code : uniq) {
        size_t binIdx = static_cast<size_t>(code >> 16);
        size_t spIdx = static_cast<size_t>(code & 0xFFFFu);
        if (result.size() <= binIdx) {
          result.resize(binIdx + 1);
        }
        auto &row = result[binIdx];
        if (row.size() <= spIdx) {
          row.resize(spIdx + 1, Counter{0});
        }
        Counter &ref = row[spIdx];
        if (ref == Counter{0} && touched) {
          touched->emplace_back(static_cast<uint32_t>(binIdx),
                                static_cast<uint16_t>(spIdx));
        }
        if (ref < std::numeric_limits<Counter>::max()) {
          ++ref;
        }
      }
    }
  }

  template <std::ranges::range value_range_t, class CounterMatrix>
  inline void bulkCount_sparse_subset(
      value_range_t &&values, const std::vector<uint32_t> &binSubset,
      CounterMatrix &result,
      std::vector<std::pair<uint32_t, uint16_t>> *touched = nullptr) {
    using PairCode = uint64_t;
    static thread_local std::vector<PairCode> uniq;

    using CounterRow = typename CounterMatrix::value_type;
    using Counter = typename CounterRow::value_type;
    static_assert(std::is_integral_v<Counter>,
                  "IMCF counter type must be integral");
    for (auto value : values) {
      uniq.clear();
      bulkContain_events_subset(value, binSubset,
                                [&](uint32_t bin, uint16_t sp) {
                                  uniq.push_back((PairCode(bin) << 16) |
                                                PairCode(sp));
                                });
      if (uniq.empty()) {
        continue;
      }
      std::sort(uniq.begin(), uniq.end());
      uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

      for (PairCode code : uniq) {
        size_t binIdx = static_cast<size_t>(code >> 16);
        size_t spIdx = static_cast<size_t>(code & 0xFFFFu);
        if (result.size() <= binIdx) {
          result.resize(binIdx + 1);
        }
        auto &row = result[binIdx];
        if (row.size() <= spIdx) {
          row.resize(spIdx + 1, Counter{0});
        }
        Counter &ref = row[spIdx];
        if (ref == Counter{0} && touched) {
          touched->emplace_back(static_cast<uint32_t>(binIdx),
                                static_cast<uint16_t>(spIdx));
        }
        if (ref < std::numeric_limits<Counter>::max()) {
          ++ref;
        }
      }
    }
  }

  template <std::ranges::range value_range_t>
  inline void
  bulkCount_sparse_compat(value_range_t &&values,
                          std::vector<std::vector<size_t>> &result) {
    std::vector<std::vector<uint32_t>> tmp(result.size());
    for (size_t i = 0; i < result.size(); ++i) {
      tmp[i].resize(result[i].size(), 0);
    }
    bulkCount_sparse(std::forward<value_range_t>(values), tmp, nullptr);
    if (result.size() < tmp.size()) {
      result.resize(tmp.size());
    }
    for (size_t i = 0; i < tmp.size(); ++i) {
      if (result[i].size() < tmp[i].size()) {
        result[i].resize(tmp[i].size(), 0);
      }
      for (size_t j = 0; j < tmp[i].size(); ++j) {
        result[i][j] += static_cast<size_t>(tmp[i][j]);
      }
    }
  }
};

#include "interleaved-merged-cuckoo-filter-qimcf-impl.h"

} // namespace chimera::imcf
