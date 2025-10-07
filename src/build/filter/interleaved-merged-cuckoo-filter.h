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
#include <buildConfig.hpp>
#include <cassert>
#include <cereal/archives/binary.hpp>
#include <cereal/details/helpers.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <hashutil.h>
#include <iterator>
#include <kvec.h>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <ranges>
#include <robin_hood.h>
#include <sdsl/int_vector.hpp>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <xxhash.h>

#ifdef IMCF_MIRROR64
#include <sys/mman.h>
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

namespace chimera::imcf {
static inline uint64_t mix64(uint64_t value) {
  return XXH3_64bits(&value, sizeof(value));
}

static inline uint8_t ceil_log2(size_t value) {
  if (value <= 1) {
    return 0;
  }
#if defined(__cpp_lib_bitops)
  return static_cast<uint8_t>(std::bit_width(static_cast<unsigned long long>(value - 1)));
#else
  size_t v = value - 1;
  uint8_t result = 0;
  while (v > 0) {
    v >>= 1;
    ++result;
  }
  return result;
#endif
}

struct LayoutHeaderV2 {
  static constexpr uint16_t CurrentMajor = 3;
  static constexpr uint16_t CurrentMinor = 0;
  static constexpr uint32_t Magic = 0x32464D43u; // 'CMF2'

  uint16_t layoutMajor{CurrentMajor};
  uint16_t layoutMinor{CurrentMinor};
  uint8_t bucketEntries{4};
  uint8_t laneBits{16};
  uint8_t routeFingerprintBits{12};
  uint8_t fpHashKind{1};
  uint8_t stashMode{0};
  uint64_t fpSalt{ChimeraBuild::IMCFConfig::DefaultFingerprintSalt};
  uint64_t routeSalt{0x9E3779B97F4A7C15ull};
  std::vector<uint8_t> sBitsPerBin;
  std::vector<uint8_t> fingerprintBitsPerBin;

  template <class Archive> void serialize(Archive &ar) {
    ar(layoutMajor, layoutMinor, bucketEntries, laneBits, routeFingerprintBits,
       fpHashKind, stashMode, fpSalt, routeSalt, sBitsPerBin,
       fingerprintBitsPerBin);
  }
};

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
    int maxTaxidsPerGroup = 16) {
  if (hashCount.empty()) {
    return {};
  }
  if (maxTaxidsPerGroup <= 0) {
    throw std::invalid_argument(
        "IMCF partition: maxTaxidsPerGroup must be positive");
  }
  size_t maxTaxids = static_cast<size_t>(maxTaxidsPerGroup);

  // Step 1: Calculate the median of the hash counts and set a threshold
  std::vector<uint64_t> counts;
  counts.reserve(hashCount.size());
  for (const auto &[taxid, count] : hashCount) {
    counts.push_back(count);
  }

  std::sort(counts.begin(), counts.end());
  uint64_t median = counts[counts.size() / 2];

  uint64_t threshold = 0;
  if (median == 0) {
    threshold = 1;
  } else if (median > std::numeric_limits<uint64_t>::max() / 64ull) {
    threshold = std::numeric_limits<uint64_t>::max();
  } else {
    threshold = median * 64ull;
  }
  if (threshold == 0) {
    threshold = 1;
  }

  // Step 2: Create chunks for hash counts exceeding the threshold
  std::vector<HashChunk> hashChunks;
  for (const auto &[taxid, count] : hashCount) {
    if (count > threshold) {
      uint64_t numChunks = count / threshold;
      if (count % threshold != 0) {
        ++numChunks;
      }
      numChunks = std::max<uint64_t>(1, numChunks);
      uint64_t chunkSize = numChunks ? count / numChunks : count;
      for (uint64_t i = 0; i < numChunks; ++i) {
        uint64_t current = chunkSize;
        if (i == numChunks - 1) {
          current = count - chunkSize * (numChunks - 1);
        }
        hashChunks.push_back({taxid, current});
      }
    } else {
      hashChunks.push_back({taxid, count});
    }
  }

  // Step 3: Sort the chunks by hash count in descending order
  std::sort(hashChunks.begin(), hashChunks.end(),
            [](const HashChunk &a, const HashChunk &b) {
              return a.hashCount > b.hashCount;
            });

  // Step 4: Create groups and assign chunks using a priority queue
  size_t groupNum = (hashChunks.size() + maxTaxids - 1) / maxTaxids;
  groupNum = std::max<size_t>(1, groupNum);
  std::vector<Group> groups(groupNum);
  std::vector<robin_hood::unordered_flat_set<std::string>> used(groupNum);

  auto cmp = [&](const int a, const int b) -> bool {
    return groups[a].totalHash > groups[b].totalHash;
  };
  std::priority_queue<int, std::vector<int>, decltype(cmp)> minHeap(cmp);
  for (size_t i = 0; i < groupNum; ++i) {
    minHeap.push(static_cast<int>(i));
  }

  for (const auto &chunk : hashChunks) {
    std::vector<int> popped;
    int target = -1;
    while (!minHeap.empty()) {
      int cand = minHeap.top();
      minHeap.pop();
      if (groups[cand].taxids.size() >= maxTaxids ||
          used[cand].contains(chunk.taxid)) {
        popped.push_back(cand);
        continue;
      }
      target = cand;
      break;
    }

    if (target == -1) {
      for (int idx : popped) {
        minHeap.push(idx);
      }
      groups.emplace_back();
      used.emplace_back();
      target = static_cast<int>(groups.size() - 1);
    } else {
      for (int idx : popped) {
        minHeap.push(idx);
      }
    }

    if (groups[target].taxids.size() >= maxTaxids) {
      groups.emplace_back();
      used.emplace_back();
      target = static_cast<int>(groups.size() - 1);
    }

    groups[target].taxids.push_back(chunk.taxid);
    groups[target].assignedHashes.push_back(chunk.hashCount);
    groups[target].totalHash += chunk.hashCount;
    used[target].insert(chunk.taxid);
    minHeap.push(target);
  }

  return groups;
}

class InterleavedMergedCuckooFilter {
  typedef kvec_t(int) kvector;
  typedef kvec_t(bool) kvectorBool;
  LayoutHeaderV2 layoutHeader;
  sdsl::bit_vector data;
  size_t binNum;
  size_t binSize;
  size_t tagNum{4};
  int MaxCuckooCount{500};
  size_t hashSize{};
#ifdef IMCF_MIRROR64
  std::vector<uint32_t> bucketMirror;
  bool bitStorageReleased{false};
  size_t segmentShift{20};
  size_t segmentSizeBuckets{1ull << 20};
  size_t segmentCount{0};
  bool segmentAdviceEnabled{false};

  struct SegmentPrefetch {
    const InterleavedMergedCuckooFilter &parent;
    size_t touched[2]{};
    size_t count{0};

    SegmentPrefetch(const InterleavedMergedCuckooFilter &p, size_t h1,
                    size_t h2)
        : parent(p) {
      request(h1);
      if (h2 != h1) {
        request(h2);
      }
    }

    void request(size_t bucket) {
      if (!parent.hasSegmentAdvice()) {
        return;
      }
      size_t seg = parent.segmentIndex(bucket);
      for (size_t i = 0; i < count; ++i) {
        if (touched[i] == seg) {
          return;
        }
      }
      parent.adviseSegment(seg, MADV_WILLNEED);
      if (count < 2) {
        touched[count++] = seg;
      }
    }

    ~SegmentPrefetch() {
      if (!parent.hasSegmentAdvice()) {
        return;
      }
      for (size_t i = 0; i < count; ++i) {
        parent.adviseSegment(touched[i], MADV_DONTNEED);
      }
    }
  };
#endif
  std::vector<std::vector<uint32_t>> activeGroups;
  struct StashEntry {
    uint64_t bucket{0};
    uint32_t fingerprint{0};
    std::vector<uint32_t> species;

    template <class Archive> void serialize(Archive &ar) {
      ar(bucket, fingerprint, species);
    }
  };
  std::vector<std::vector<StashEntry>> stash;
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

  inline bool insertIntoStash(size_t binIndex, size_t bucket, uint32_t fingerprint,
                              size_t speciesIndex) {
    if (binIndex >= stash.size()) {
      return false;
    }
    auto &entries = stash[binIndex];
    for (auto &entry : entries) {
      if (entry.bucket == bucket && entry.fingerprint == fingerprint) {
        if (std::find(entry.species.begin(), entry.species.end(),
                      static_cast<uint32_t>(speciesIndex)) !=
            entry.species.end()) {
          return true;
        }
        entry.species.push_back(static_cast<uint32_t>(speciesIndex));
        return true;
      }
    }
    StashEntry newEntry;
    newEntry.bucket = bucket;
    newEntry.fingerprint = fingerprint;
    newEntry.species.push_back(static_cast<uint32_t>(speciesIndex));
    entries.push_back(std::move(newEntry));
    return true;
  }

  inline size_t bucketLinearIndex(size_t bucket, size_t bin) const {
    return bucket * binNum + bin;
  }

  inline size_t bucketBitOffset(size_t bucketLinear) const {
    return bucketLinear * tagNum * static_cast<size_t>(layoutHeader.laneBits);
  }

  inline size_t laneLinearIndex(size_t bucket, size_t bin, size_t lane) const {
    return (bucket * binNum + bin) * tagNum + lane;
  }

  inline size_t laneBitOffset(size_t bucket, size_t bin, size_t lane) const {
    return bucketBitOffset(bucketLinearIndex(bucket, bin)) +
           lane * static_cast<size_t>(layoutHeader.laneBits);
  }

  inline size_t laneCount() const { return tagNum; }

  inline uint32_t readLaneValue(size_t bucket, size_t bin, size_t lane) const {
#ifdef IMCF_MIRROR64
    if (!bucketMirror.empty()) {
      return bucketMirror[laneLinearIndex(bucket, bin, lane)];
    }
#endif
    size_t offset = laneBitOffset(bucket, bin, lane);
    return static_cast<uint32_t>(data.get_int(offset, layoutHeader.laneBits));
  }

  inline void writeLaneValue(size_t bucket, size_t bin, size_t lane,
                             uint32_t value) {
    size_t offset = laneBitOffset(bucket, bin, lane);
    data.set_int(offset, value, layoutHeader.laneBits);
#ifdef IMCF_MIRROR64
    if (!bucketMirror.empty()) {
      bucketMirror[laneLinearIndex(bucket, bin, lane)] = value;
    }
#endif
  }

  inline uint8_t speciesBits(size_t bin) const {
    if (bin >= layoutHeader.sBitsPerBin.size()) {
      return 0;
    }
    return layoutHeader.sBitsPerBin[bin];
  }

  inline uint8_t fingerprintBits(size_t bin) const {
    if (!layoutHeader.fingerprintBitsPerBin.empty() &&
        bin < layoutHeader.fingerprintBitsPerBin.size()) {
      return layoutHeader.fingerprintBitsPerBin[bin];
    }
    uint8_t sBits = speciesBits(bin);
    return static_cast<uint8_t>(layoutHeader.laneBits > sBits
                                    ? layoutHeader.laneBits - sBits
                                    : 0);
  }

  inline uint32_t fingerprintMask(size_t bin) const {
    uint8_t bits = fingerprintBits(bin);
    if (bits >= 32) {
      return 0xFFFFFFFFu;
    }
    if (bits == 0) {
      return 0u;
    }
    return (1u << bits) - 1u;
  }

  inline uint32_t speciesMask(size_t bin) const {
    uint8_t bits = speciesBits(bin);
    if (bits >= 32) {
      return 0xFFFFFFFFu;
    }
    if (bits == 0) {
      return 0u;
    }
    return (1u << bits) - 1u;
  }

  inline uint32_t routeMask() const {
    uint8_t bits = layoutHeader.routeFingerprintBits;
    if (bits >= 32) {
      return 0xFFFFFFFFu;
    }
    if (bits == 0) {
      return 0u;
    }
    return (1u << bits) - 1u;
  }

  inline uint32_t buildLaneValue(size_t bin, uint32_t fingerprint,
                                 uint32_t speciesIdx) const {
    uint32_t fpMask = fingerprintMask(bin);
    uint32_t spMask = speciesMask(bin);
    uint8_t fBits = fingerprintBits(bin);
    uint32_t packed = (speciesIdx & spMask) << fBits;
    packed |= (fingerprint & fpMask);
    return packed;
  }

  inline uint32_t laneFingerprint(uint32_t laneValue, size_t bin) const {
    return laneValue & fingerprintMask(bin);
  }

  inline uint32_t laneSpecies(uint32_t laneValue, size_t bin) const {
    uint8_t fBits = fingerprintBits(bin);
    uint32_t spMask = speciesMask(bin);
    return (laneValue >> fBits) & spMask;
  }

  struct FingerprintInfo {
    uint32_t full{0};
    uint32_t route{0};
  };

  inline FingerprintInfo makeFingerprint(uint64_t value, size_t bin, uint32_t presetRoute = 0u) const {
    uint64_t mixed = mix64(value ^ layoutHeader.fpSalt);
    uint32_t fpMask = fingerprintMask(bin);
    uint32_t full =
        fpMask == 0u ? static_cast<uint32_t>(mixed)
                     : (fpMask == 0xFFFFFFFFu ? static_cast<uint32_t>(mixed)
                                               : static_cast<uint32_t>(mixed & fpMask));
    if (full == 0u) {
      full = 1u;
    }
    uint32_t rMask = routeMask();
    uint32_t route = rMask == 0u ? 0u : (presetRoute & rMask);
    if (rMask != 0u && route == 0u) {
      route = static_cast<uint32_t>(mixed & rMask);
      if (route == 0u) {
        uint64_t remixed =
            mix64((static_cast<uint64_t>(full) ^ layoutHeader.routeSalt) +
                  value);
        route = static_cast<uint32_t>(remixed & rMask);
        if (route == 0u) {
          route = 1u;
        }
      }
    }
    if (rMask != 0u) {
      if ((full & rMask) != route) {
        full = (full & ~rMask) | route;
      }
    }
    return FingerprintInfo{full, route};
  }

  inline uint32_t computeRouteFingerprint(uint64_t value) const {
    uint32_t mask = routeMask();
    if (mask == 0u) {
      return 0u;
    }
    uint64_t mixed = mix64(value ^ layoutHeader.fpSalt);
    uint32_t route = static_cast<uint32_t>(mixed & mask);
    if (route == 0u) {
      uint64_t remixed =
          mix64((static_cast<uint64_t>(mixed) ^ layoutHeader.routeSalt) +
                value);
      route = static_cast<uint32_t>(remixed & mask);
      if (route == 0u) {
        route = 1u;
      }
    }
    return route;
  }

  inline void gatherCandidateBins(size_t bucket1, size_t bucket2,
                                  std::vector<uint32_t> &candidates) const {
    candidates.clear();
    if (activeGroups.empty()) {
      candidates.reserve(binNum);
      for (uint32_t bin = 0; bin < binNum; ++bin) {
        candidates.push_back(bin);
      }
      return;
    }

    const auto &list1 = activeGroups[bucket1];
    const auto &list2 = activeGroups[bucket2];
    size_t i = 0, j = 0;
    candidates.reserve(list1.size() + list2.size());
    while (i < list1.size() || j < list2.size()) {
      uint32_t valueOut;
      if (j >= list2.size() ||
          (i < list1.size() && list1[i] < list2[j])) {
        valueOut = list1[i++];
      } else if (i >= list1.size() || list2[j] < list1[i]) {
        valueOut = list2[j++];
      } else {
        valueOut = list1[i];
        ++i;
        ++j;
      }
      if (candidates.empty() || candidates.back() != valueOut) {
        candidates.push_back(valueOut);
      }
    }
  }

  inline void collectMatchesForBin(uint64_t value, uint32_t routeFingerprint,
                                   size_t bin, size_t bucket1, size_t bucket2,
                                   std::vector<uint32_t> &matches) const {
    matches.clear();
    FingerprintInfo fp = makeFingerprint(value, bin, routeFingerprint);

    auto probeBucket = [&](size_t bucket) {
      for (size_t lane = 0; lane < laneCount(); ++lane) {
        uint32_t laneValue = readLaneValue(bucket, bin, lane);
        if (laneValue == 0u) {
          continue;
        }
        if (laneFingerprint(laneValue, bin) == fp.full) {
          matches.push_back(laneSpecies(laneValue, bin));
        }
      }
      if (bin < stash.size()) {
        for (const auto &entry : stash[bin]) {
          if (entry.bucket != bucket || entry.fingerprint != fp.full) {
            continue;
          }
          matches.insert(matches.end(), entry.species.begin(),
                         entry.species.end());
        }
      }
    };

    probeBucket(bucket1);
    if (bucket2 != bucket1) {
      probeBucket(bucket2);
    }

    if (!matches.empty()) {
      std::sort(matches.begin(), matches.end());
      matches.erase(std::unique(matches.begin(), matches.end()),
                    matches.end());
    }
  }

#ifdef IMCF_MIRROR64
  inline size_t segmentIndex(size_t bucket) const {
    return bucket >> segmentShift;
  }

  inline void prepareSegmentAdvice() {
    segmentAdviceEnabled = false;
    segmentCount = 0;
  }

  inline void adviseSegment(size_t, int) const {}

  inline bool hasSegmentAdvice() const { return false; }

  inline void initializeMirrorStorage() {
    if (hashSize == 0 || binNum == 0) {
      bucketMirror.clear();
      bitStorageReleased = false;
      prepareSegmentAdvice();
      return;
    }
    bucketMirror.assign(hashSize * binNum * laneCount(), 0u);
    bitStorageReleased = false;
    prepareSegmentAdvice();
  }

  inline void rebuildMirrorFromBitVector() {
    if (hashSize == 0 || binNum == 0) {
      bucketMirror.clear();
      bitStorageReleased = false;
      prepareSegmentAdvice();
      return;
    }
    bucketMirror.assign(hashSize * binNum * laneCount(), 0u);
    for (size_t bucket = 0; bucket < hashSize; ++bucket) {
      for (size_t bin = 0; bin < binNum; ++bin) {
        for (size_t lane = 0; lane < laneCount(); ++lane) {
          size_t idx = laneLinearIndex(bucket, bin, lane);
          bucketMirror[idx] = static_cast<uint32_t>(
              data.get_int(laneBitOffset(bucket, bin, lane),
                           layoutHeader.laneBits));
        }
      }
    }
    bitStorageReleased = false;
    prepareSegmentAdvice();
  }
#endif

  inline void fetchActiveBins(size_t bucket, uint32_t routeFingerprint,
                              std::vector<uint32_t> &out) const {
    out.clear();
    if (bucket >= activeGroups.size()) {
      return;
    }
    const auto &candidates = activeGroups[bucket];
    if (candidates.empty()) {
      return;
    }
    out.reserve(candidates.size());
    for (uint32_t bin : candidates) {
      bool match = false;
      for (size_t lane = 0; lane < laneCount(); ++lane) {
        uint32_t laneValue = readLaneValue(bucket, bin, lane);
        if (laneValue == 0u) {
          continue;
        }
        if ((laneFingerprint(laneValue, bin) & routeMask()) == routeFingerprint) {
          match = true;
          break;
        }
      }
      if (!match) {
        const auto &stashEntries = stash[bin];
        for (const auto &entry : stashEntries) {
          if (entry.bucket == bucket &&
              (entry.fingerprint & routeMask()) == routeFingerprint) {
            match = true;
            break;
          }
        }
        if (!match) {
          continue;
        }
      }
      if (match) {
        out.push_back(bin);
      }
    }
  }

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

    layoutHeader = LayoutHeaderV2{};
    layoutHeader.bucketEntries = static_cast<uint8_t>(tagNum);
    layoutHeader.laneBits = 32;
    layoutHeader.routeFingerprintBits = layoutHeader.laneBits;
    layoutHeader.fpHashKind = 1;
    layoutHeader.stashMode = 0;
    layoutHeader.fpSalt = config.fpSalt;
    layoutHeader.routeSalt = 0x9E3779B97F4A7C15ull;
    layoutHeader.sBitsPerBin.assign(binNum, 0);
    layoutHeader.fingerprintBitsPerBin.assign(binNum, 0);
    uint8_t minFingerprintBits = layoutHeader.laneBits;
    for (size_t i = 0; i < binNum; ++i) {
      size_t speciesCount = groups[i].taxids.size();
      uint8_t sBits = ceil_log2(speciesCount);
      layoutHeader.sBitsPerBin[i] = sBits;
      uint8_t fBits = static_cast<uint8_t>(layoutHeader.laneBits > sBits
                                               ? layoutHeader.laneBits - sBits
                                               : 0);
      layoutHeader.fingerprintBitsPerBin[i] = fBits;
      if (fBits < minFingerprintBits) {
        minFingerprintBits = fBits;
      }
    }
    if (minFingerprintBits == 0) {
      minFingerprintBits = 1;
    }
    if (layoutHeader.routeFingerprintBits > minFingerprintBits) {
      layoutHeader.routeFingerprintBits = minFingerprintBits;
    }

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
    data = sdsl::bit_vector(static_cast<uint64_t>(binNum) *
                                static_cast<uint64_t>(binSize) *
                                static_cast<uint64_t>(tagNum) *
                                static_cast<uint64_t>(layoutHeader.laneBits),
                            0);

    config.binNum = binNum;
    config.binSize = binSize;
#ifdef IMCF_MIRROR64
    initializeMirrorStorage();
#endif
    stash.assign(binNum, {});
  }

  const LayoutHeaderV2 &layout() const { return layoutHeader; }
  const std::vector<uint8_t> &sBitsPerBin() const {
    return layoutHeader.sBitsPerBin;
  }

  template <class Archive> void serialize(Archive &ar) {
    ar(layoutHeader, data, binNum, binSize, tagNum, MaxCuckooCount, hashSize,
       stash);
#ifdef IMCF_MIRROR64
    if constexpr (Archive::is_loading::value) {
      rebuildMirrorFromBitVector();
    } else {
      if (bitStorageReleased) {
        throw std::runtime_error(
            "IMCF mirror: bit-vector released, cannot serialize");
      }
    }
#endif
    if constexpr (Archive::is_loading::value) {
      if (layoutHeader.layoutMajor != LayoutHeaderV2::CurrentMajor) {
        throw std::runtime_error(
            "IMCF 布局版本不兼容，请重新构建数据库 (layout_major=" +
            std::to_string(layoutHeader.layoutMajor) + ")");
      }
      if (layoutHeader.sBitsPerBin.size() != binNum) {
        throw std::runtime_error(
            "IMCF 布局信息损坏：sBits 数量与 bin 数不匹配");
      }
      if (layoutHeader.bucketEntries != static_cast<uint8_t>(tagNum)) {
        throw std::runtime_error(
            "IMCF 布局信息损坏：bucket_entries 与 tagNum 不一致");
      }
      if (layoutHeader.fingerprintBitsPerBin.size() != binNum) {
        if (layoutHeader.fingerprintBitsPerBin.empty()) {
          layoutHeader.fingerprintBitsPerBin.assign(binNum, 0);
          for (size_t i = 0; i < binNum; ++i) {
            uint8_t sBits = layoutHeader.sBitsPerBin[i];
            layoutHeader.fingerprintBitsPerBin[i] = static_cast<uint8_t>(
                layoutHeader.laneBits > sBits ? layoutHeader.laneBits - sBits
                                               : 0);
          }
        } else {
          throw std::runtime_error(
              "IMCF 布局信息损坏：fingerprintBits 数量与 bin 数不匹配");
        }
      }
      if (stash.size() != binNum) {
        stash.resize(binNum);
      }
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
inline size_t altHash(size_t bucket, uint32_t routeFingerprint) const {
  size_t mask = hashSize - 1;
  uint64_t mixed = mix64((static_cast<uint64_t>(routeFingerprint) ^ layoutHeader.routeSalt));
  size_t delta = (static_cast<size_t>(mixed) & mask) | 1ull;
  return (bucket ^ delta) & mask;
}

  inline void route(uint64_t value, std::vector<uint32_t> &bins) const {
    bins.clear();

    uint32_t routeFp = computeRouteFingerprint(value);
    size_t hash1 = hashIndex(value);
    size_t hash2 = altHash(hash1, routeFp);

    static thread_local std::vector<uint32_t> routeBuf1;
    static thread_local std::vector<uint32_t> routeBuf2;

    fetchActiveBins(hash1, routeFp, routeBuf1);
    if (hash2 != hash1) {
      fetchActiveBins(hash2, routeFp, routeBuf2);
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

  /**
   * @brief 插入指纹。
   *
   * 构造标签后先尝试两个桶，若均告满载则进入踢出流程。
   */
	inline bool insertTag(size_t binIndex, uint64_t value, size_t speciesIndex) {
    uint8_t sBits = speciesBits(binIndex);
    size_t speciesCapacity = sBits == 0 ? 1ull : (1ull << sBits);
    if (speciesIndex >= speciesCapacity) {
      return false;
    }

    uint32_t routeFp = computeRouteFingerprint(value);
    FingerprintInfo fp = makeFingerprint(value, binIndex, routeFp);
    uint32_t lanePacked =
        buildLaneValue(binIndex, fp.full, static_cast<uint32_t>(speciesIndex));

    size_t bucket1 = hashIndex(value);
    size_t bucket2 = altHash(bucket1, fp.route);

    auto tryInsertBucket = [&](size_t bucket) -> bool {
      for (size_t lane = 0; lane < laneCount(); ++lane) {
        uint32_t current = readLaneValue(bucket, binIndex, lane);
        if (current == 0u) {
          writeLaneValue(bucket, binIndex, lane, lanePacked);
          return true;
        }
        if (laneFingerprint(current, binIndex) == fp.full &&
            laneSpecies(current, binIndex) == speciesIndex) {
          return true;
        }
      }
      return false;
    };

    if (tryInsertBucket(bucket1)) {
      return true;
    }
    if (tryInsertBucket(bucket2)) {
      return true;
    }

    auto countMatching = [&](size_t bucket) -> int {
      int matches = 0;
      for (size_t lane = 0; lane < laneCount(); ++lane) {
        uint32_t current = readLaneValue(bucket, binIndex, lane);
        if (current != 0u && laneFingerprint(current, binIndex) == fp.full) {
          ++matches;
        }
      }
      return matches;
    };

    int sameFingerprint = countMatching(bucket1) +
                          (bucket1 != bucket2 ? countMatching(bucket2) : 0);
    if (sameFingerprint >= static_cast<int>(laneCount() * 2)) {
      bool stored = insertIntoStash(binIndex, bucket1, fp.full, speciesIndex);
      if (!stored && bucket2 != bucket1) {
        stored = insertIntoStash(binIndex, bucket2, fp.full, speciesIndex);
      }
      if (stored) {
        return true;
      }
      insertFailureTotal.fetch_add(1, std::memory_order_relaxed);
      insertFailureSaturated.fetch_add(1, std::memory_order_relaxed);
      return false;
    }

    size_t startBucket =
        (mix64(value ^ layoutHeader.routeSalt) & 1ull) ? bucket2 : bucket1;
    bool kicked = kickOut(binIndex, startBucket, fp, lanePacked);
    if (!kicked) {
      insertFailureTotal.fetch_add(1, std::memory_order_relaxed);
      insertFailureSaturated.fetch_add(1, std::memory_order_relaxed);
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
      uint64_t load = 0;
      for (size_t lane = 0; lane < laneCount(); ++lane) {
        if (readLaneValue(bucket, binIndex, lane) != 0u) {
          ++load;
        }
      }
      if (load < stats.occupancy.size()) {
        stats.occupancy[load]++;
      } else {
        stats.occupancy.back()++;
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
      size_t fullIndex =
          std::min<size_t>(stats.occupancy.size() - 1, laneCount());
      stats.percentFull = static_cast<double>(stats.occupancy[fullIndex]) /
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
          size_t count = entry.species.size();
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

  /**
   * @brief 踢出流程，只在两桶之间往返。
   *
   * 随机逐出当前桶一项换入新项，再借助可逆哈希跳往另一桶继续尝试，最多迭代
   * `MaxCuckooCount` 次。
   */
  inline bool kickOut(size_t binIndex, size_t startBucket, FingerprintInfo fp,
                      uint32_t laneValue) {
    size_t bucket = startBucket;

    static thread_local std::mt19937_64 gen{std::random_device{}()};
    std::uniform_int_distribution<size_t> dis(0, laneCount() - 1);

    for (int cnt = 0; cnt < MaxCuckooCount; ++cnt) {
      size_t lane = dis(gen);
      uint32_t victimValue = readLaneValue(bucket, binIndex, lane);
      writeLaneValue(bucket, binIndex, lane, laneValue);
      if (victimValue == 0u || victimValue == laneValue) {
        return true;
      }

      uint32_t victimFingerprint = laneFingerprint(victimValue, binIndex);
      uint32_t victimRoute = victimFingerprint & routeMask();
      if (victimRoute == 0u) {
        uint64_t remixed =
            mix64((static_cast<uint64_t>(victimFingerprint) ^
                   layoutHeader.routeSalt) + bucket);
        victimRoute = static_cast<uint32_t>(remixed & routeMask());
        if (victimRoute == 0u) {
          victimRoute = 1u;
        }
        victimFingerprint = (victimFingerprint & ~routeMask()) | victimRoute;
        writeLaneValue(bucket, binIndex, lane, victimValue);
      }

      laneValue = victimValue;
      fp.full = victimFingerprint;
      fp.route = victimRoute;
      bucket = altHash(bucket, fp.route);
    }
    return false;
  }

  template <typename Emit>
  inline void forEachStashMatch(size_t binIndex, size_t bucket,
                                uint32_t fingerprint, Emit &&emit) const {
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
      for (uint32_t species : entry.species) {
        emit(static_cast<uint16_t>(species));
      }
      return;
    }
  }

  /**
   * @brief Checks for the presence of a given minimizer hash in the Interleaved
   * Merged Cuckoo Filter and collects matching species indices per bin.
   *
   * The function computes the route fingerprint, probes both candidate buckets
   * for each relevant bin, and gathers the species indices whose fingerprints
   * match the query. Results are stored as a vector of species lists, where each
   * index corresponds to a bin in the filter.
   *
   * @param value The input minimizer hash value to be checked.
   * @param result Output container sized to the number of bins; each entry
   * contains the species indices matched within that bin.
   */
	inline void bulkContain(uint64_t value, std::vector<std::vector<uint32_t>> &result) {
    if (binNum == 0) {
      result.clear();
      return;
    }
    result.assign(binNum, {});

    uint32_t routeFp = computeRouteFingerprint(value);
    size_t hash1 = hashIndex(value);
    size_t hash2 = altHash(hash1, routeFp);

    static thread_local std::vector<uint32_t> candidates;
    static thread_local std::vector<uint32_t> matches;
    gatherCandidateBins(hash1, hash2, candidates);

    for (uint32_t bin : candidates) {
      collectMatchesForBin(value, routeFp, bin, hash1, hash2, matches);
      if (!matches.empty()) {
        result[bin] = matches;
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
   * The function uses an intermediate vector produced by `bulkContain` for each
   * minimizer and updates the per-bin counters accordingly:
   * - For every species index returned in `tmpResult`, the corresponding entry
   *   in `result` is incremented.
   * - If a sub-vector in `result` does not have enough space for a new tag,
   *   it is resized and initialized to zero.
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
    if (binNum == 0) {
      return;
    }

    uint32_t routeFp = computeRouteFingerprint(value);
    size_t hash1 = hashIndex(value);
    size_t hash2 = altHash(hash1, routeFp);

    static thread_local std::vector<uint32_t> candidates;
    static thread_local std::vector<uint32_t> matches;
    gatherCandidateBins(hash1, hash2, candidates);

    for (uint32_t bin : candidates) {
      collectMatchesForBin(value, routeFp, bin, hash1, hash2, matches);
      for (uint32_t species : matches) {
        emit(bin, static_cast<uint16_t>(species));
      }
    }
  }


  template <class EmitFn>
	inline void bulkContain_events_subset(uint64_t value,
                                        const std::vector<uint32_t> &binSubset,
                                        EmitFn &&emit) {
    if (binNum == 0 || binSubset.empty()) {
      return;
    }
    uint32_t routeFp = computeRouteFingerprint(value);
    size_t hash1 = hashIndex(value);
    size_t hash2 = altHash(hash1, routeFp);

    static thread_local std::vector<uint32_t> matches;
    for (uint32_t bin : binSubset) {
      if (bin >= binNum) {
        continue;
      }
      collectMatchesForBin(value, routeFp, bin, hash1, hash2, matches);
      for (uint32_t species : matches) {
        emit(bin, static_cast<uint16_t>(species));
      }
    }
  }


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

  inline void clearActiveGroups() { activeGroups.clear(); }

  inline bool hasActiveIndex() const {
    return activeGroups.size() == hashSize && !activeGroups.empty();
  }

#ifdef IMCF_MIRROR64
  inline void releaseBitStorage() {
    sdsl::bit_vector empty;
    empty.swap(data);
    bitStorageReleased = true;
  }
#endif

  inline void buildActiveGroups() {
    if (hashSize == 0 || binNum == 0) {
      activeGroups.clear();
      return;
    }

    activeGroups.assign(hashSize, {});

#pragma omp parallel for schedule(dynamic, 32)
    for (size_t bucket = 0; bucket < hashSize; ++bucket) {
      auto &groupList = activeGroups[bucket];
      groupList.clear();
      groupList.reserve(16);

      for (uint32_t bin = 0; bin < binNum; ++bin) {
        bool occupied = false;
        for (size_t lane = 0; lane < laneCount(); ++lane) {
          if (readLaneValue(bucket, bin, lane) != 0u) {
            occupied = true;
            break;
          }
        }
        if (occupied) {
          groupList.push_back(static_cast<uint32_t>(bin));
        }
      }
    }
  }

  inline bool saveActiveIndex(const std::string &path) const {
    if (hashSize == 0 || binNum == 0) {
      return false;
    }
    if (activeGroups.empty()) {
      return false;
    }

    struct IndexHeader {
      uint32_t magic{0x494D4349u}; // 'IMCI'
      uint16_t version{1};
      uint16_t reserved{0};
      uint64_t binCount{0};
      uint64_t hashCount{0};
    };

    std::ofstream out(path, std::ios::binary);
    if (!out.is_open()) {
      return false;
    }

    IndexHeader header;
    header.binCount = static_cast<uint64_t>(binNum);
    header.hashCount = static_cast<uint64_t>(hashSize);
    out.write(reinterpret_cast<const char *>(&header), sizeof(header));
    if (!out) {
      return false;
    }

    std::vector<uint64_t> offsets(hashSize + 1, 0);
    std::vector<uint8_t> payload;
    payload.reserve(activeGroups.size() * 4);

    for (size_t bucket = 0; bucket < hashSize; ++bucket) {
      offsets[bucket] = payload.size();
      const auto &groupList = activeGroups[bucket];
      uint32_t prev = 0;
      bool first = true;
      for (uint32_t bin : groupList) {
        uint32_t delta = first ? bin : (bin - prev);
        first = false;
        prev = bin;
        encodeVarint(payload, delta);
      }
    }
    offsets[hashSize] = payload.size();

    out.write(reinterpret_cast<const char *>(offsets.data()),
              offsets.size() * sizeof(uint64_t));
    if (!out) {
      return false;
    }
    if (!payload.empty()) {
      out.write(reinterpret_cast<const char *>(payload.data()), payload.size());
      if (!out) {
        return false;
      }
    }
    return true;
  }

  inline bool loadActiveIndex(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) {
      return false;
    }

    struct IndexHeader {
      uint32_t magic;
      uint16_t version;
      uint16_t reserved;
      uint64_t binCount;
      uint64_t hashCount;
    };

    IndexHeader header{};
    in.read(reinterpret_cast<char *>(&header), sizeof(header));
    if (!in) {
      return false;
    }
    if (header.magic != 0x494D4349u || header.version != 1) {
      return false;
    }
    if (header.binCount != static_cast<uint64_t>(binNum) ||
        header.hashCount != static_cast<uint64_t>(hashSize)) {
      return false;
    }

    std::vector<uint64_t> offsets(header.hashCount + 1, 0);
    in.read(reinterpret_cast<char *>(offsets.data()),
            offsets.size() * sizeof(uint64_t));
    if (!in) {
      return false;
    }

    std::vector<uint8_t> payload((std::istreambuf_iterator<char>(in)),
                                 std::istreambuf_iterator<char>());
    if (payload.empty() && offsets.back() != 0) {
      return false;
    }
    if (!payload.empty() && offsets.back() > payload.size()) {
      return false;
    }

    activeGroups.assign(hashSize, {});
    for (size_t bucket = 0; bucket < hashSize; ++bucket) {
      uint64_t start = offsets[bucket];
      uint64_t end = offsets[bucket + 1];
      if (start > end || end > payload.size()) {
        activeGroups.clear();
        return false;
      }
      size_t offset = static_cast<size_t>(start);
      uint32_t prev = 0;
      while (offset < end) {
        uint32_t delta = 0;
        if (!decodeVarint(payload, offset, delta)) {
          activeGroups.clear();
          return false;
        }
        uint32_t bin = prev + delta;
        prev = bin;
        activeGroups[bucket].push_back(bin);
      }
    }
    return true;
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
} // namespace chimera::imcf
