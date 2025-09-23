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
 *	which contains the basic operations of the Interleaved Merged Cuckoo Filter
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#pragma once

#include <sdsl/int_vector.hpp>
#include <hashutil.h>
#include <kvec.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/details/helpers.hpp>
#include <xxhash.h>
#include <simde/x86/avx2.h>
#include <simde/x86/sse2.h>
#include <random>
#include <robin_hood.h>
#include <queue>
#include <buildConfig.hpp>
#include <bitset>
#include <cassert>
#include <algorithm>
#include <ranges>
#include <utility>
#include <fstream>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <limits>
#include <type_traits>
#ifdef IMCF_MIRROR64
#include <sys/mman.h>
#endif

static inline size_t ceil_div_u64(uint64_t a, uint64_t b) {
	return (size_t)((a + b - 1) / b);
}

static inline size_t next_pow2(size_t v) {
	if (v <= 1) return 1;
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
	struct Group {
		std::vector<std::string> taxids;
		uint64_t totalHash{ 0 };
	};

	struct HashChunk {
		std::string taxid;
		uint64_t hashCount;
	};


	/**
	 * @brief Partitions hash counts into groups using a greedy algorithm for the Interleaved Merged Cuckoo Filter (IMCF).
	 *
	 * This function partitions the provided hash counts into groups to balance the total hash count in each group
	 * as evenly as possible, while adhering to a maximum number of taxids per group. The goal is to optimize the
	 * distribution for efficient processing within the IMCF.
	 *
	 * @param hashCount A map containing taxid strings as keys and their corresponding hash counts as values.
	 * @param mode A string specifying the partitioning mode. Currently supports "normal" mode for partitioning logic.
	 * @param maxTaxidsPerGroup The maximum number of taxids allowed in a single group. Default is 16.
	 *
	 * @return A vector of Group objects, where each Group contains a collection of taxids and the total hash count for that group.
	 *
	 * @details
	 * The function first calculates the median of the hash counts and sets a threshold as `median * 64`. This threshold
	 * helps determine if a taxid's hash count should be split into smaller chunks to distribute the load more evenly.
	 *
	 * - Step 1: Calculate the median of hash counts and establish a threshold.
	 * - Step 2: Split any hash count exceeding the threshold into smaller chunks and add them to `hashChunks`.
	 * - Step 3: Sort the `hashChunks` vector in descending order of hash counts.
	 * - Step 4: Use a priority queue to greedily assign chunks to groups, ensuring the group with the current lowest
	 *           total hash count receives the next chunk.
	 *
	 * The function ensures that each group receives chunks in a balanced manner, maintaining the load across all groups.
	 *
	 * @note
	 * - The function assumes that `HashChunk` and `Group` types are defined. `HashChunk` should have `taxid` and `hashCount` fields,
	 *   while `Group` should contain a `taxids` vector and `totalHash` field.
	 * - The `minHeap` priority queue helps maintain the group with the smallest total hash count for efficient chunk assignment.
	 */
	inline std::vector<Group> partitionHashCount(const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount, std::string mode, int maxTaxidsPerGroup = 16)
	{
		if (hashCount.empty()) {
			return {};
		}
		if (maxTaxidsPerGroup <= 0) {
			throw std::invalid_argument("IMCF partition: maxTaxidsPerGroup must be positive");
		}
		size_t maxTaxids = static_cast<size_t>(maxTaxidsPerGroup);

		if (mode == "normal")
		{
			// Step 1: Calculate the median of the hash counts and set a threshold
			std::vector<uint64_t> counts;
			counts.reserve(hashCount.size());
			for (const auto& [taxid, count] : hashCount) {
				counts.push_back(count);
			}

			std::sort(counts.begin(), counts.end());
			uint64_t median = counts[counts.size() / 2];

			uint64_t threshold = 0;
			if (median == 0) {
				threshold = 1;
			}
			else if (median > std::numeric_limits<uint64_t>::max() / 64ull) {
				threshold = std::numeric_limits<uint64_t>::max();
			}
			else {
				threshold = median * 64ull;
			}
			if (threshold == 0) {
				threshold = 1;
			}

			// Step 2: Create chunks for hash counts exceeding the threshold
			std::vector<HashChunk> hashChunks;
			for (const auto& [taxid, count] : hashCount) {
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
						hashChunks.push_back({ taxid, current });
					}
				}
				else {
					hashChunks.push_back({ taxid, count });
				}
			}

			// Step 3: Sort the chunks by hash count in descending order
			std::sort(hashChunks.begin(), hashChunks.end(),
				[](const HashChunk& a, const HashChunk& b) {
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

			for (const auto& chunk : hashChunks) {
				std::vector<int> popped;
				int target = -1;
				while (!minHeap.empty()) {
					int cand = minHeap.top();
					minHeap.pop();
				if (groups[cand].taxids.size() >= maxTaxids || used[cand].contains(chunk.taxid)) {
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
				}
				else {
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
				groups[target].totalHash += chunk.hashCount;
				used[target].insert(chunk.taxid);
				minHeap.push(target);
			}

			return groups;
		}
		else if (mode == "fast")
		{
			// Step 1: Sort taxids based on hash counts in descending order
			std::vector<std::pair<std::string, uint64_t>> sortedHashCounts;
			sortedHashCounts.reserve(hashCount.size());
			for (const auto& item : hashCount) {
				sortedHashCounts.emplace_back(item.first, item.second);
			}
			std::sort(sortedHashCounts.begin(), sortedHashCounts.end(),
				[](const std::pair<std::string, uint64_t>& a, const std::pair<std::string, uint64_t>& b) {
					return a.second > b.second;
				});
			// Calculate median and number of groups
			size_t taxidNum = sortedHashCounts.size();
			if (taxidNum == 0) {
				return {};
			}
			uint64_t median = sortedHashCounts[taxidNum / 2].second;
			size_t groupNum = (taxidNum + maxTaxids - 1) / maxTaxids;
			groupNum = std::max<size_t>(1, groupNum);
			// Initialize groups
			std::vector<Group> groups;
			groups.reserve(groupNum);
			groups.resize(groupNum);

			// Step 2: Create a priority queue to manage groups by total hash count
			auto cmp = [&](const int a, const int b) -> bool {
				return groups[a].totalHash > groups[b].totalHash;
				};
			std::priority_queue<int, std::vector<int>, decltype(cmp)> minHeap(cmp);
			for (size_t i = 0; i < groupNum; ++i) {
				minHeap.push(static_cast<int>(i));
			}
			// Step 3: Assign each taxid to the group with the smallest total hash count
			for (const auto& [taxid, count] : sortedHashCounts) {
				int smallestGroup = minHeap.top();
				minHeap.pop();
				// Check if the group can accommodate more taxids
				if (groups[smallestGroup].taxids.size() < maxTaxidsPerGroup) {
					groups[smallestGroup].taxids.push_back(taxid);
					groups[smallestGroup].totalHash += count;
					minHeap.push(smallestGroup);
				}
				else
				{
					// Handle case where all current groups are full
					bool assigned = false;
					std::vector<int> tempGroups;
					// Try to find a group that can accommodate the current taxid
					while (!minHeap.empty()) {
						int tempGroup = minHeap.top();
						minHeap.pop();
						if (groups[tempGroup].taxids.size() < maxTaxidsPerGroup) {
							groups[tempGroup].taxids.push_back(taxid);
							groups[tempGroup].totalHash += count;
							minHeap.push(tempGroup);
							assigned = true;
							break;
						}
						tempGroups.push_back(tempGroup);
					}
					// Reinsert groups back into the priority queue
					for (int tempGroup : tempGroups) {
						minHeap.push(tempGroup);
					}
					// Create a new group if no existing group can take the taxid
					if (!assigned) {
						groups.emplace_back(Group());
						groups.back().taxids.push_back(taxid);
						groups.back().totalHash += count;
						minHeap.push(static_cast<int>(groups.size() - 1));
					}
				}
			}


			//uint64_t maxTotalHash = 0;
			//uint64_t minTotalHash = UINT64_MAX;
			//for (const auto& group : groups) {
			//	maxTotalHash = std::max(maxTotalHash, group.totalHash);
			//	minTotalHash = std::min(minTotalHash, group.totalHash);
			//}

			//std::cout << "Difference between max and min totalHash: " << (maxTotalHash - minTotalHash) << std::endl;
			//std::cout << "Filter size: " << (maxTotalHash * groups.size() * 16.0 / 0.95) / 8 / 1024 / 1024 / 1024 << " GB" << std::endl;
			//std::cout << "Groups size: " << groups.size() << std::endl;
			//std::cout << "Before Split Groups size: " << (hashCount.size() + 16) / 16 << std::endl;

			return groups;
		}
		else
		{
			std::cerr << "Invalid mode: " << mode << std::endl;
			throw std::runtime_error("Invalid mode");
		}
		
	}


	class InterleavedMergedCuckooFilter {
		typedef kvec_t(int) kvector;
		typedef kvec_t(bool) kvectorBool;
sdsl::bit_vector data;
size_t binNum;
size_t binSize;
size_t tagNum{ 4 };
int MaxCuckooCount{ 500 };
size_t hashSize{};
#ifdef IMCF_MIRROR64
std::vector<uint64_t> bucketMirror;
bool bitStorageReleased{ false };
size_t segmentShift{ 20 };
size_t segmentSizeBuckets{ 1ull << 20 };
size_t segmentCount{ 0 };
bool segmentAdviceEnabled{ false };

struct SegmentPrefetch {
    const InterleavedMergedCuckooFilter& parent;
    size_t touched[2]{};
    size_t count{ 0 };

    SegmentPrefetch(const InterleavedMergedCuckooFilter& p, size_t h1, size_t h2)
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
std::vector<uint64_t> routerBucketOffsets;
std::vector<uint16_t> routerEntryFingerprints;
std::vector<uint64_t> routerEntryOffsets;
	std::vector<uint8_t> routerPayload;

	static inline void decodeVarintSequence(const std::vector<uint8_t>& dataBytes,
		size_t start, size_t end, std::vector<uint32_t>& out) {
		out.clear();
		if (start >= end) {
			return;
		}
		size_t offset = start;
		uint32_t prev = 0;
		while (offset < end) {
			uint32_t delta = 0;
			if (!decodeVarint(dataBytes, offset, delta)) {
				out.clear();
				return;
			}
			uint32_t value = prev + delta;
			out.push_back(value);
			prev = value;
		}
	}

	inline size_t bucketLinearIndex(size_t bucket, size_t bin) const {
		return bucket * binNum + bin;
	}

	inline size_t bucketBitOffset(size_t bucketLinear) const {
		return bucketLinear * tagNum * 16;
	}

	inline uint64_t readBucketRaw(size_t bucketLinear) const {
		return data.get_int(bucketBitOffset(bucketLinear), 64);
	}

	inline void writeBucketRaw(size_t bucketLinear, uint64_t value) {
		data.set_int(bucketBitOffset(bucketLinear), value, 64);
	}

inline uint64_t readBucket64(size_t bucket, size_t bin) const {
	size_t linear = bucketLinearIndex(bucket, bin);
#ifdef IMCF_MIRROR64
	return bucketMirror[linear];
#else
	return readBucketRaw(linear);
#endif
}

inline void writeBucket64(size_t bucket, size_t bin, uint64_t value) {
	size_t linear = bucketLinearIndex(bucket, bin);
#ifdef IMCF_MIRROR64
	bucketMirror[linear] = value;
	if (!bitStorageReleased) {
		writeBucketRaw(linear, value);
	}
#else
	writeBucketRaw(linear, value);
#endif
}

#ifdef IMCF_MIRROR64
	inline size_t segmentIndex(size_t bucket) const {
		return bucket >> segmentShift;
	}

	inline void prepareSegmentAdvice() {
		if (segmentShift >= sizeof(size_t) * 8 - 1) {
			segmentShift = static_cast<size_t>(sizeof(size_t) * 8 - 2);
		}
		segmentSizeBuckets = static_cast<size_t>(1) << segmentShift;
		if (segmentSizeBuckets == 0) {
			segmentSizeBuckets = 1ULL << 20;
		}
		if (bucketMirror.empty() || hashSize == 0) {
			segmentAdviceEnabled = false;
			segmentCount = 0;
			return;
		}
		segmentCount = (hashSize + segmentSizeBuckets - 1) / segmentSizeBuckets;
		segmentAdviceEnabled = segmentCount > 0;
		if (segmentAdviceEnabled) {
			void* ptr = const_cast<uint64_t*>(bucketMirror.data());
			size_t bytes = bucketMirror.size() * sizeof(uint64_t);
			::madvise(ptr, bytes, MADV_SEQUENTIAL);
		}
	}

	inline void adviseSegment(size_t segment, int advice) const {
		if (!segmentAdviceEnabled || segment >= segmentCount) {
			return;
		}
		size_t startBucket = segment * segmentSizeBuckets;
		if (startBucket >= hashSize) {
			return;
		}
		size_t buckets = std::min(segmentSizeBuckets, hashSize - startBucket);
		if (buckets == 0) {
			return;
		}
		size_t elements = buckets * binNum;
		void* ptr = const_cast<uint64_t*>(bucketMirror.data() + startBucket * binNum);
		size_t bytes = elements * sizeof(uint64_t);
		::madvise(ptr, bytes, advice);
	}

	inline bool hasSegmentAdvice() const {
		return segmentAdviceEnabled;
	}

	inline void initializeMirrorStorage() {
		bucketMirror.clear();
		if (hashSize == 0 || binNum == 0) {
			bitStorageReleased = false;
			prepareSegmentAdvice();
			return;
		}
		bucketMirror.assign(hashSize * binNum, 0);
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
		bucketMirror.resize(hashSize * binNum);
		for (size_t bucket = 0; bucket < hashSize; ++bucket) {
			for (size_t bin = 0; bin < binNum; ++bin) {
				size_t linear = bucketLinearIndex(bucket, bin);
				bucketMirror[linear] = readBucketRaw(linear);
			}
		}
		bitStorageReleased = false;
		prepareSegmentAdvice();
	}
#endif

	inline bool locateRouterEntry(size_t bucket, uint16_t fingerprint, size_t& entryIndex) const {
		if (routerBucketOffsets.size() != hashSize + 1) {
			return false;
		}
		if (bucket >= hashSize) {
			return false;
		}
		size_t begin = static_cast<size_t>(routerBucketOffsets[bucket]);
		size_t end = static_cast<size_t>(routerBucketOffsets[bucket + 1]);
		if (begin >= end) {
			return false;
		}
		size_t left = begin;
		size_t right = end;
		while (left < right) {
			size_t mid = left + ((right - left) >> 1);
			uint16_t midFp = routerEntryFingerprints[mid];
			if (midFp < fingerprint) {
				left = mid + 1;
			}
			else {
				right = mid;
			}
		}
		if (left >= end) {
			return false;
		}
		if (routerEntryFingerprints[left] != fingerprint) {
			return false;
		}
		entryIndex = left;
		return true;
	}

	inline void fetchRouterBins(size_t bucket, uint16_t fingerprint, std::vector<uint32_t>& out) const {
		size_t entryIndex = 0;
		if (!locateRouterEntry(bucket, fingerprint, entryIndex)) {
			out.clear();
			return;
		}
		if (entryIndex + 1 >= routerEntryOffsets.size()) {
			out.clear();
			return;
		}
		size_t start = static_cast<size_t>(routerEntryOffsets[entryIndex]);
		size_t end = static_cast<size_t>(routerEntryOffsets[entryIndex + 1]);
		decodeVarintSequence(routerPayload, start, end, out);
	}

	static inline void encodeVarint(std::vector<uint8_t>& out, uint32_t value) {
		while (value >= 0x80) {
			out.push_back(static_cast<uint8_t>((value & 0x7F) | 0x80));
			value >>= 7;
		}
		out.push_back(static_cast<uint8_t>(value));
	}

	static inline bool decodeVarint(const std::vector<uint8_t>& dataBytes, size_t& offset, uint32_t& value) {
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


		friend std::ostream& operator<<(std::ostream& os, const InterleavedMergedCuckooFilter& filter)
		{
			os << "Interleaved Merged Cuckoo Filter" << std::endl;
			os << "Bin number: " << filter.binNum << std::endl;
			os << "Bin size: " << filter.binSize << std::endl;

			return os;
		}

	public:
		InterleavedMergedCuckooFilter() = default;
		InterleavedMergedCuckooFilter(std::vector<Group>& groups, ChimeraBuild::IMCFConfig& config)
		{
			uint64_t maxTotalHash = 0;
			for (const auto& group : groups) {
				if (group.totalHash > maxTotalHash) {
					maxTotalHash = group.totalHash;
				}
			}
			binNum = groups.size();

			// 关键：向上取整 + 2 的幂，便于用按位掩码实现可逆备桶
			size_t needBinSize = ceil_div_u64(maxTotalHash, (uint64_t)(config.loadFactor * (double)tagNum));
			hashSize = next_pow2(std::max<size_t>(1, needBinSize));
			binSize = hashSize;

			// 保持旧布局：每桶 4×16bit，总位数 = binNum * binSize * tagNum * 16
			data = sdsl::bit_vector((uint64_t)binNum * (uint64_t)binSize * (uint64_t)tagNum * 16ull, 0);

			config.binNum = binNum;
			config.binSize = binSize;
#ifdef IMCF_MIRROR64
			initializeMirrorStorage();
#endif
		}

		template <class Archive>
		void serialize(Archive& ar) {
			ar(data, binNum, binSize, tagNum, MaxCuckooCount, hashSize);
#ifdef IMCF_MIRROR64
			if constexpr (Archive::is_loading::value) {
				rebuildMirrorFromBitVector();
			}
			else {
				if (bitStorageReleased) {
					throw std::runtime_error("IMCF mirror: bit-vector released, cannot serialize");
				}
			}
#endif
		}

		/**
		 * @brief 轻量主哈希索引。
		 *
		 * 使用若干位运算完成快速搅拌，并利用 `hashSize` 为 2 的幂这一前提直接按位掩码。
		 * 相比重型哈希大幅减少运算量，同时仍保持良好的分布特性。
		 */
		inline size_t hashIndex(uint64_t value) const {
			// 轻量搅拌 + 2^k 掩码（要求 hashSize 为 2 的幂；已在构造器中保证）
			uint64_t x = value ^ (value >> 33) ^ (value >> 17) ^ (value >> 9);
			return (size_t)x & (hashSize - 1);
		}


		/**
		 * @brief 可逆备桶哈希。
		 *
		 * 基于 fingerprint 乘黄金常数得到的扰动值与主桶索引 XOR，并依赖 2^k 大小做掩码，
		 * 保证 `altHash(altHash(b, fp), fp) == b`，便于在踢出链中双向定位。
		 */
		inline size_t altHash(size_t b, uint16_t fingerprint) const {
			// 对 fp 做一个短 hash（与 b 无关），再 XOR；最后做按位掩码
			uint64_t h = (uint64_t)fingerprint * 0x9E3779B185EBCA87ull;
			return (b ^ (size_t)h) & (hashSize - 1);
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
		inline uint16_t reduceTo12bitAndAddIndex(size_t value, size_t index) const
		{
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
		inline uint16_t reduceTo12bit(size_t value) const
		{
			uint64_t x = value ^ (value >> 33) ^ (value >> 17) ^ (value >> 9);
			uint16_t v = static_cast<uint16_t>(x & 0x0FFFu);
			return v ? v : 1;
		}

		inline bool hasRouterIndex() const {
			return routerBucketOffsets.size() == hashSize + 1 && routerEntryOffsets.size() == routerEntryFingerprints.size() + 1;
		}

		inline void clearRouterIndex() {
			routerBucketOffsets.clear();
			routerEntryFingerprints.clear();
			routerEntryOffsets.clear();
			routerPayload.clear();
		}

		inline void route(size_t value, std::vector<uint32_t>& bins) const {
			bins.clear();
			if (!hasRouterIndex()) {
				return;
			}

			uint16_t fingerprint = reduceTo12bit(value);
			size_t hash1 = hashIndex(value);
			size_t hash2 = altHash(hash1, fingerprint);

			static thread_local std::vector<uint32_t> routeBuf1;
			static thread_local std::vector<uint32_t> routeBuf2;

			fetchRouterBins(hash1, fingerprint, routeBuf1);
			if (hash2 != hash1) {
				fetchRouterBins(hash2, fingerprint, routeBuf2);
			}
			else {
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
				if (j >= routeBuf2.size() || (i < routeBuf1.size() && routeBuf1[i] < routeBuf2[j])) {
					valueOut = routeBuf1[i++];
				}
				else if (i >= routeBuf1.size() || routeBuf2[j] < routeBuf1[i]) {
					valueOut = routeBuf2[j++];
				}
				else {
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
		 * 构造标签后先命中主桶，再跳到可逆备桶，两桶任意存在空位或已有指纹即视为成功；
		 * 若均告满载则进入踢出流程。
		 */
		inline bool insertTag(size_t binIndex, size_t value, size_t index)
		{
			if (index >= 16) {
				assert(false && "IMCF group index overflow (taxids per group must be <=16)");
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
				if (chunk == tag) return true;
			}
			return false;
		};

			if (try_insert_bucket(b1)) return true;
			if (try_insert_bucket(b2)) return true;

			return kickOut(binIndex, value, tag);
		}


		/**
		 * @brief 踢出流程，只在两桶之间往返。
		 *
		 * 随机逐出当前桶一项换入新项，再借助可逆哈希跳往另一桶继续尝试，最多迭代 `MaxCuckooCount` 次。
		 */
		inline bool kickOut(size_t binIndex, size_t value, uint16_t tag)
		{
			uint16_t cur = tag;
			uint16_t fp = (uint16_t)(cur & 0x0FFFu);
			size_t b = hashIndex(value);

			static thread_local std::mt19937_64 gen{ std::random_device{}() };
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
				if (chunk == cur) return true;
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


		/**
		 * @brief Checks for the presence of a given minimizer hash in the Interleaved Merged Cuckoo Filter and updates the result vector.
		 *
		 * This function uses the SIMDe library to accelerate the process of checking whether a given minimizer hash exists
		 * in the Interleaved Merged Cuckoo Filter (IMCF). The hash value is reduced to 12 bits, and two hash positions are
		 * checked for matches. If a match is found, the function extracts the position index from the upper 4 bits and updates
		 * the result vector.
		 *
		 * @param value The input minimizer hash value to be checked.
		 * @param result A vector of std::bitset<16> where each bitset corresponds to a bin index. Each bit represents whether
		 *               a particular position index (species index) has been found in the filter.
		 *
		 * @details
		 * The function first reduces the input value to a 12-bit fingerprint using the reduceTo12bit function. It then computes
		 * two hash positions (hash1 and hash2) to identify potential bucket locations in the filter. The function iterates over
		 * each bin index, checks the data in these bucket positions, and loads them into SIMD vectors for comparison.
		 *
		 * - SIMDe Vectorization:
		 *   The function uses SIMD vectors to perform parallel comparisons:
		 *   - fingerprint_vec: Contains the 12-bit fingerprint replicated across all vector positions.
		 *   - fingerprint_mask: A mask used to extract the 12-bit fingerprint from each 16-bit entry.
		 *   - species_shift: The number of bits to shift to extract the species index from the upper 4 bits of each entry.
		 *
		 * The function updates the result vector by setting bits corresponding to matching species indices.
		 *
		 * @note Ensure the data object provides a method to retrieve 64-bit data (e.g., data.get_int). The function requires
		 * the SIMDe library for SIMD operations.
		 */
		inline void bulkContain(size_t value, std::vector<std::bitset<16>>& result)
		{
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
			for (size_t binIndex = 0; binIndex < binNum; binIndex++)
			{
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
				simde__m128i bucket_vec = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(entries));

				// Step 5: Extract fingerprints and species indices
				simde__m128i fingerprints = simde_mm_and_si128(bucket_vec, fingerprint_mask);
				simde__m128i species_indices = simde_mm_srli_epi16(bucket_vec, species_shift);

				// Compare extracted fingerprints with the input fingerprint
				simde__m128i cmp = simde_mm_cmpeq_epi16(fingerprints, fingerprint_vec);

				// Step 6: Update the result vector with matching species indices
				// Store comparison results and species indices to arrays for variable indexing
				uint16_t cmp_vals[8];
				uint16_t species_vals[8];
				simde_mm_storeu_si128(reinterpret_cast<simde__m128i*>(cmp_vals), cmp);
				simde_mm_storeu_si128(reinterpret_cast<simde__m128i*>(species_vals), species_indices);

				for (int i = 0; i < 8; ++i)
				{
					if (cmp_vals[i] == 0xFFFF)
					{
						uint16_t speciesIndex = species_vals[i];
						result[binIndex].set(speciesIndex);
					}
				}
			}
		}


		/**
		 * @brief Performs bulk counting of input minimizers and updates the result based on their presence in the filter.
		 *
		 * This function iterates over a range of input minimizers and checks their presence in the Interleaved Merged
		 * Cuckoo Filter (IMCF). For each minimizer, it performs a bulk containment query and updates the `result`
		 * vector with the counts based on the corresponding positions in the filter.
		 *
		 * @tparam value_range_t A range type representing the collection of minimizers to be checked.
		 * @param values The input range of minimizers to query in the IMCF.
		 * @param result A 2D vector that stores the count results. Each sub-vector corresponds to a specific position index,
		 *        and each element in the sub-vector represents the count for a particular tag at that index.
		 *
		 * @details
		 * The function uses a temporary vector of `std::bitset<16>` to store the intermediate results of the containment
		 * check for each minimizer. For each input value:
		 * - The `bitset` vector is reset to ensure no residual data from previous iterations.
		 * - The `bulkContain` function is called to determine the presence of the minimizer across multiple positions.
		 * - The results in `tmpResult` are analyzed, and the corresponding positions in the `result` vector are updated:
		 *   - If a bit is set in `tmpResult`, the respective count in `result` is incremented.
		 *   - If a sub-vector in `result` does not have enough space for a new tag, it is resized and initialized to zero.
		 *
		 * @note
		 * Ensure the `result` vector is pre-initialized with sub-vectors corresponding to the positions to be checked.
		 * The `bulkContain` function should be defined to handle the minimizer containment query.
		 */
		template <std::ranges::range value_range_t>
		inline void bulkCount(value_range_t&& values, std::vector<std::vector<size_t>>& result)
		{
			bulkCount_sparse_compat(std::forward<value_range_t>(values), result);
		}

		template <class EmitFn>
	inline void bulkContain_events(size_t value, EmitFn&& emit) {
		uint16_t fingerprint = reduceTo12bit(value);
		size_t hash1 = hashIndex(value);
		size_t hash2 = altHash(hash1, fingerprint);

#ifdef IMCF_MIRROR64
		[[maybe_unused]] SegmentPrefetch segmentPrefetch(*this, hash1, hash2);
#endif

		const uint64_t fpBroadcast = static_cast<uint64_t>(fingerprint) * 0x0001000100010001ULL;
		const uint64_t fpMask = 0x0FFF0FFF0FFF0FFFULL;
		const uint64_t idxMask = 0x000F000F000F000FULL;

			constexpr size_t BLOCK = 16;
			constexpr size_t SUB_BATCH = 4;

			auto process64 = [&](uint64_t q, uint32_t binIndex) {
				uint64_t diff = (q ^ fpBroadcast) & fpMask;
				uint64_t idxBits = (q >> 12) & idxMask;
				for (int lane = 0; lane < 4; ++lane) {
					if (((diff >> (lane * 16)) & 0x0FFFULL) == 0ULL) {
						uint16_t speciesIndex = static_cast<uint16_t>((idxBits >> (lane * 16)) & 0xF);
						emit(static_cast<uint32_t>(binIndex), speciesIndex);
					}
				}
			};

			if (activeGroups.empty()) {
				for (size_t base = 0; base < binNum; base += BLOCK) {
					size_t end = std::min(binNum, base + BLOCK);
					for (size_t batchStart = base; batchStart < end; batchStart += SUB_BATCH) {
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
						}
					}
				}
				return;
			}

			const auto& list1 = activeGroups[hash1];
			const auto& list2 = activeGroups[hash2];
			if (list1.empty() && list2.empty()) {
				return;
			}

			size_t i = 0, j = 0;
			uint32_t batchBins[SUB_BATCH];
			bool active1[SUB_BATCH];
			bool active2[SUB_BATCH];
			size_t count = 0;

			auto flush = [&]() {
				for (size_t idx = 0; idx < count; ++idx) {
					uint32_t bin = batchBins[idx];
					uint64_t qPrimary = 0;
					uint64_t qAlt = 0;
					if (active1[idx]) {
						qPrimary = readBucket64(hash1, bin);
					}
					if (active2[idx]) {
						qAlt = readBucket64(hash2, bin);
					}
					if (qPrimary != 0) {
						process64(qPrimary, bin);
					}
					if (qAlt != 0) {
						process64(qAlt, bin);
					}
				}
				count = 0;
			};

			while (i < list1.size() || j < list2.size()) {
				uint32_t bin;
				bool in1 = false;
				bool in2 = false;
				if (j >= list2.size() || (i < list1.size() && list1[i] <= list2[j])) {
					bin = list1[i++];
					in1 = true;
					if (j < list2.size() && list2[j] == bin) {
						in2 = true;
						++j;
					}
				} else {
					bin = list2[j++];
					in2 = true;
				}
				batchBins[count] = bin;
				active1[count] = in1;
				active2[count] = in2;
				++count;
				if (count == SUB_BATCH) {
					flush();
				}
			}
			if (count) {
				flush();
			}
		}

		template <class EmitFn>
	inline void bulkContain_events_subset(size_t value,
		const std::vector<uint32_t>& binSubset,
		EmitFn&& emit) {
		uint16_t fingerprint = reduceTo12bit(value);
		size_t hash1 = hashIndex(value);
		size_t hash2 = altHash(hash1, fingerprint);

#ifdef IMCF_MIRROR64
		[[maybe_unused]] SegmentPrefetch segmentPrefetch(*this, hash1, hash2);
#endif

			const uint64_t fpBroadcast = static_cast<uint64_t>(fingerprint) * 0x0001000100010001ULL;
			const uint64_t fpMask = 0x0FFF0FFF0FFF0FFFULL;
			const uint64_t idxMask = 0x000F000F000F000FULL;

			constexpr size_t SUB_BATCH = 4;

			auto process64 = [&](uint64_t q, uint32_t binIndex) {
				uint64_t diff = (q ^ fpBroadcast) & fpMask;
				uint64_t idxBits = (q >> 12) & idxMask;
				for (int lane = 0; lane < 4; ++lane) {
					if (((diff >> (lane * 16)) & 0x0FFFULL) == 0ULL) {
						uint16_t speciesIndex = static_cast<uint16_t>((idxBits >> (lane * 16)) & 0xF);
						emit(static_cast<uint32_t>(binIndex), speciesIndex);
					}
				}
			};

			if (hasRouterIndex()) {
				if (binSubset.empty()) {
					return;
				}
				assert(std::is_sorted(binSubset.begin(), binSubset.end()));
				static thread_local std::vector<uint32_t> routeBuf1;
				static thread_local std::vector<uint32_t> routeBuf2;
				static thread_local std::vector<uint32_t> mergedBins;
				static thread_local std::vector<uint8_t> mergedMask;

				fetchRouterBins(hash1, fingerprint, routeBuf1);
				if (hash2 != hash1) {
					fetchRouterBins(hash2, fingerprint, routeBuf2);
				}
				else {
					routeBuf2.clear();
				}

				if (routeBuf1.empty() && routeBuf2.empty()) {
					return;
				}

				mergedBins.clear();
				mergedMask.clear();
				mergedBins.reserve(routeBuf1.size() + routeBuf2.size());
				mergedMask.reserve(routeBuf1.size() + routeBuf2.size());

				size_t i = 0;
				size_t j = 0;
				while (i < routeBuf1.size() || j < routeBuf2.size()) {
					uint32_t candidate;
					uint8_t maskBits = 0;
					if (j >= routeBuf2.size() || (i < routeBuf1.size() && routeBuf1[i] < routeBuf2[j])) {
						candidate = routeBuf1[i++];
						maskBits |= 0x1;
					}
					else if (i >= routeBuf1.size() || routeBuf2[j] < routeBuf1[i]) {
						candidate = routeBuf2[j++];
						maskBits |= 0x2;
					}
					else {
						candidate = routeBuf1[i];
						maskBits |= 0x3;
						++i;
						++j;
					}
					if (!mergedBins.empty() && mergedBins.back() == candidate) {
						mergedMask.back() |= maskBits;
					} else {
						mergedBins.push_back(candidate);
						mergedMask.push_back(maskBits);
					}
				}

				for (size_t idx = 0; idx < mergedBins.size(); ++idx) {
					uint32_t binIndex = mergedBins[idx];
					if (!std::binary_search(binSubset.begin(), binSubset.end(), binIndex)) {
						continue;
					}
					if (binIndex >= binNum) {
						continue;
					}
					uint8_t maskBits = mergedMask[idx];
					if (maskBits & 0x1) {
						uint64_t q = readBucket64(hash1, binIndex);
						if (q != 0ULL) {
							process64(q, binIndex);
						}
					}
					if ((maskBits & 0x2) && hash2 != hash1) {
						uint64_t q = readBucket64(hash2, binIndex);
						if (q != 0ULL) {
							process64(q, binIndex);
						}
		}
	}
				return;
			}

			if (activeGroups.empty()) {
				for (size_t offset = 0; offset < binSubset.size(); offset += SUB_BATCH) {
					size_t batchEnd = std::min(binSubset.size(), offset + SUB_BATCH);
					uint64_t q1[SUB_BATCH];
					uint64_t q2[SUB_BATCH];
					uint32_t bins[SUB_BATCH];
					size_t cnt = 0;
					for (size_t i = offset; i < batchEnd; ++i) {
						uint32_t binIndex = binSubset[i];
						if (binIndex >= binNum) continue;
						bins[cnt] = binIndex;
						q1[cnt] = readBucket64(hash1, binIndex);
						q2[cnt] = readBucket64(hash2, binIndex);
						++cnt;
					}
					for (size_t i = 0; i < cnt; ++i) {
						if (q1[i] != 0) process64(q1[i], bins[i]);
						if (q2[i] != 0) process64(q2[i], bins[i]);
					}
				}
				return;
			}

			const auto& list1 = activeGroups[hash1];
			const auto& list2 = activeGroups[hash2];

			auto contains = [](const std::vector<uint32_t>& vec, uint32_t value) {
				return !vec.empty() && std::binary_search(vec.begin(), vec.end(), value);
			};

			for (size_t offset = 0; offset < binSubset.size(); offset += SUB_BATCH) {
				size_t batchEnd = std::min(binSubset.size(), offset + SUB_BATCH);
				size_t cnt = 0;
				uint32_t bins[SUB_BATCH];
				bool use1[SUB_BATCH];
				bool use2[SUB_BATCH];
				for (size_t i = offset; i < batchEnd; ++i) {
					uint32_t binIndex = binSubset[i];
					if (binIndex >= binNum) continue;
					bool in1 = contains(list1, binIndex);
					bool in2 = contains(list2, binIndex);
					if (!in1 && !in2) continue;
					bins[cnt] = binIndex;
					use1[cnt] = in1;
					use2[cnt] = in2;
					++cnt;
				}
				if (cnt == 0) continue;
				for (size_t i = 0; i < cnt; ++i) {
					uint32_t binIndex = bins[i];
					if (use1[i]) {
						uint64_t q = readBucket64(hash1, binIndex);
						if (q != 0) process64(q, binIndex);
					}
					if (use2[i]) {
						uint64_t q = readBucket64(hash2, binIndex);
						if (q != 0) process64(q, binIndex);
					}
				}
			}
		}

		template <std::ranges::range value_range_t, class CounterMatrix>
		inline void bulkCount_sparse(value_range_t&& values,
			CounterMatrix& result,
			std::vector<std::pair<uint32_t, uint16_t>>* touched = nullptr)
		{
			using CounterRow = typename CounterMatrix::value_type;
			using Counter = typename CounterRow::value_type;
			static_assert(std::is_integral_v<Counter>, "IMCF counter type must be integral");
			for (auto value : values) {
				bulkContain_events(value, [&](uint32_t bin, uint16_t sp) {
					size_t binIdx = static_cast<size_t>(bin);
					size_t spIdx = static_cast<size_t>(sp);
					if (result.size() <= binIdx) {
						result.resize(binIdx + 1);
					}
					auto& row = result[binIdx];
					if (row.size() <= spIdx) {
						row.resize(spIdx + 1, Counter{ 0 });
					}
					Counter& ref = row[spIdx];
					if (ref == Counter{ 0 } && touched) {
						touched->emplace_back(bin, sp);
					}
					if (ref < std::numeric_limits<Counter>::max()) {
						++ref;
					}
				});
			}
		}

		template <std::ranges::range value_range_t, class CounterMatrix>
		inline void bulkCount_sparse_subset(value_range_t&& values,
			const std::vector<uint32_t>& binSubset,
			CounterMatrix& result,
			std::vector<std::pair<uint32_t, uint16_t>>* touched = nullptr)
		{
			using CounterRow = typename CounterMatrix::value_type;
			using Counter = typename CounterRow::value_type;
			static_assert(std::is_integral_v<Counter>, "IMCF counter type must be integral");
			for (auto value : values) {
				bulkContain_events_subset(value, binSubset, [&](uint32_t bin, uint16_t sp) {
					size_t binIdx = static_cast<size_t>(bin);
					size_t spIdx = static_cast<size_t>(sp);
					if (result.size() <= binIdx) {
						result.resize(binIdx + 1);
					}
					auto& row = result[binIdx];
					if (row.size() <= spIdx) {
						row.resize(spIdx + 1, Counter{ 0 });
					}
					Counter& ref = row[spIdx];
					if (ref == Counter{ 0 } && touched) {
						touched->emplace_back(bin, sp);
					}
					if (ref < std::numeric_limits<Counter>::max()) {
						++ref;
					}
				});
			}
		}

		inline void clearActiveGroups() {
			activeGroups.clear();
		}

		inline bool hasActiveIndex() const {
			return !activeGroups.empty();
		}

#ifdef IMCF_MIRROR64
		inline void releaseBitStorage() {
			sdsl::bit_vector empty;
			empty.swap(data);
			bitStorageReleased = true;
		}
#endif

		inline void buildActiveGroups() {
			activeGroups.assign(hashSize, {});
			for (size_t bucket = 0; bucket < hashSize; ++bucket) {
				auto& groupList = activeGroups[bucket];
				for (size_t bin = 0; bin < binNum; ++bin) {
					uint64_t q = readBucket64(bucket, bin);
					if (q != 0ULL) {
						groupList.push_back(static_cast<uint32_t>(bin));
					}
				}
			}
		}

		inline void buildRouterIndex() {
			clearRouterIndex();
			if (hashSize == 0 || binNum == 0) {
				return;
			}

			routerBucketOffsets.assign(hashSize + 1, 0);
			routerEntryOffsets.clear();
			routerEntryOffsets.reserve(hashSize * 2);
			routerEntryOffsets.push_back(0);
			routerEntryFingerprints.clear();
			routerEntryFingerprints.reserve(hashSize * 2);
			routerPayload.clear();

			for (size_t bucket = 0; bucket < hashSize; ++bucket) {
				robin_hood::unordered_flat_map<uint16_t, std::vector<uint32_t>> perFp;
				perFp.reserve(32);
				for (size_t bin = 0; bin < binNum; ++bin) {
					uint64_t q = readBucket64(bucket, bin);
					if (q == 0ULL) {
						continue;
					}
					for (int lane = 0; lane < static_cast<int>(tagNum); ++lane) {
						uint16_t tag = static_cast<uint16_t>((q >> (lane * 16)) & 0xFFFFu);
						if (tag == 0u) {
							continue;
						}
						uint16_t fp = static_cast<uint16_t>(tag & 0x0FFFu);
						perFp[fp].push_back(static_cast<uint32_t>(bin));
					}
				}

				if (perFp.empty()) {
					routerBucketOffsets[bucket + 1] = routerBucketOffsets[bucket];
					continue;
				}

				std::vector<std::pair<uint16_t, std::vector<uint32_t>>> entries;
				entries.reserve(perFp.size());
				for (auto& kv : perFp) {
					auto& bins = kv.second;
					std::sort(bins.begin(), bins.end());
					bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
					if (!bins.empty()) {
						entries.emplace_back(kv.first, std::move(bins));
					}
				}

				if (entries.empty()) {
					routerBucketOffsets[bucket + 1] = routerBucketOffsets[bucket];
					continue;
				}

				std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
					return a.first < b.first;
				});

				for (auto& [fp, bins] : entries) {
					routerEntryFingerprints.push_back(fp);
					uint32_t prev = 0;
					bool first = true;
					for (uint32_t bin : bins) {
						uint32_t delta = first ? bin : (bin - prev);
						first = false;
						prev = bin;
						encodeVarint(routerPayload, delta);
					}
					routerEntryOffsets.push_back(routerPayload.size());
				}

				routerBucketOffsets[bucket + 1] = routerEntryFingerprints.size();
			}

			if (routerEntryOffsets.size() != routerEntryFingerprints.size() + 1) {
				clearRouterIndex();
			}
		}

		inline bool saveActiveIndex(const std::string& path) const {
			if (hashSize == 0 || binNum == 0) {
				return false;
			}
			if (activeGroups.empty()) {
				return false;
			}

			struct IndexHeader {
				uint32_t magic{ 0x494D4349u };// 'IMCI'
				uint16_t version{ 1 };
				uint16_t reserved{ 0 };
				uint64_t binCount{ 0 };
				uint64_t hashCount{ 0 };
			};

			std::ofstream out(path, std::ios::binary);
			if (!out.is_open()) {
				return false;
			}

			IndexHeader header;
			header.binCount = static_cast<uint64_t>(binNum);
			header.hashCount = static_cast<uint64_t>(hashSize);
			out.write(reinterpret_cast<const char*>(&header), sizeof(header));
			if (!out) {
				return false;
			}

			std::vector<uint64_t> offsets(hashSize + 1, 0);
			std::vector<uint8_t> payload;
			payload.reserve(activeGroups.size() * 4);

			for (size_t bucket = 0; bucket < hashSize; ++bucket) {
				offsets[bucket] = payload.size();
				const auto& groupList = activeGroups[bucket];
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

			out.write(reinterpret_cast<const char*>(offsets.data()), offsets.size() * sizeof(uint64_t));
			if (!out) {
				return false;
			}
			if (!payload.empty()) {
				out.write(reinterpret_cast<const char*>(payload.data()), payload.size());
				if (!out) {
					return false;
				}
			}
			return true;
		}

		inline bool loadActiveIndex(const std::string& path) {
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
			in.read(reinterpret_cast<char*>(&header), sizeof(header));
			if (!in) {
				return false;
			}
			if (header.magic != 0x494D4349u || header.version != 1) {
				return false;
			}
			if (header.binCount != static_cast<uint64_t>(binNum) || header.hashCount != static_cast<uint64_t>(hashSize)) {
				return false;
			}

			std::vector<uint64_t> offsets(header.hashCount + 1, 0);
			in.read(reinterpret_cast<char*>(offsets.data()), offsets.size() * sizeof(uint64_t));
			if (!in) {
				return false;
			}

			std::vector<uint8_t> payload((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
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

		inline bool saveRouterIndex(const std::string& path) const {
			if (!hasRouterIndex()) {
				return false;
			}

			struct RouterHeader {
				uint32_t magic{ 0x524D4349u }; // 'IMCR'
				uint16_t version{ 1 };
				uint16_t reserved{ 0 };
				uint64_t binCount{ 0 };
				uint64_t hashCount{ 0 };
				uint64_t entryCount{ 0 };
				uint64_t payloadBytes{ 0 };
			};

			std::ofstream out(path, std::ios::binary);
			if (!out.is_open()) {
				return false;
			}

			RouterHeader header;
			header.binCount = static_cast<uint64_t>(binNum);
			header.hashCount = static_cast<uint64_t>(hashSize);
			header.entryCount = static_cast<uint64_t>(routerEntryFingerprints.size());
			header.payloadBytes = static_cast<uint64_t>(routerPayload.size());
			out.write(reinterpret_cast<const char*>(&header), sizeof(header));
			if (!out) {
				return false;
			}

			out.write(reinterpret_cast<const char*>(routerBucketOffsets.data()), routerBucketOffsets.size() * sizeof(uint64_t));
			if (!out) {
				return false;
			}

			out.write(reinterpret_cast<const char*>(routerEntryOffsets.data()), routerEntryOffsets.size() * sizeof(uint64_t));
			if (!out) {
				return false;
			}

			out.write(reinterpret_cast<const char*>(routerEntryFingerprints.data()), routerEntryFingerprints.size() * sizeof(uint16_t));
			if (!out) {
				return false;
			}

			if (!routerPayload.empty()) {
				out.write(reinterpret_cast<const char*>(routerPayload.data()), routerPayload.size());
				if (!out) {
					return false;
				}
			}

			return true;
		}

		inline bool loadRouterIndex(const std::string& path) {
			std::ifstream in(path, std::ios::binary);
			if (!in.is_open()) {
				return false;
			}

			struct RouterHeader {
				uint32_t magic;
				uint16_t version;
				uint16_t reserved;
				uint64_t binCount;
				uint64_t hashCount;
				uint64_t entryCount;
				uint64_t payloadBytes;
			};

			RouterHeader header{};
			in.read(reinterpret_cast<char*>(&header), sizeof(header));
			if (!in) {
				return false;
			}
			if (header.magic != 0x524D4349u || header.version != 1) {
				return false;
			}
			if (header.binCount != static_cast<uint64_t>(binNum) || header.hashCount != static_cast<uint64_t>(hashSize)) {
				return false;
			}

			routerBucketOffsets.resize(static_cast<size_t>(header.hashCount) + 1);
			in.read(reinterpret_cast<char*>(routerBucketOffsets.data()), routerBucketOffsets.size() * sizeof(uint64_t));
			if (!in) {
				clearRouterIndex();
				return false;
			}

			routerEntryOffsets.resize(static_cast<size_t>(header.entryCount) + 1);
			in.read(reinterpret_cast<char*>(routerEntryOffsets.data()), routerEntryOffsets.size() * sizeof(uint64_t));
			if (!in) {
				clearRouterIndex();
				return false;
			}

			routerEntryFingerprints.resize(static_cast<size_t>(header.entryCount));
			in.read(reinterpret_cast<char*>(routerEntryFingerprints.data()), routerEntryFingerprints.size() * sizeof(uint16_t));
			if (!in) {
				clearRouterIndex();
				return false;
			}

			routerPayload.resize(static_cast<size_t>(header.payloadBytes));
			if (!routerPayload.empty()) {
				in.read(reinterpret_cast<char*>(routerPayload.data()), routerPayload.size());
				if (!in) {
					clearRouterIndex();
					return false;
				}
			}

			if (routerBucketOffsets.size() != hashSize + 1 || routerEntryOffsets.size() != routerEntryFingerprints.size() + 1) {
				clearRouterIndex();
				return false;
			}

			return true;
		}

		template <std::ranges::range value_range_t>
		inline void bulkCount_sparse_compat(value_range_t&& values,
			std::vector<std::vector<size_t>>& result)
		{
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
}
