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
#include <xxhash.h>
#include <simde/x86/avx2.h>
#include <simde/x86/sse2.h>
#include <random>
#include <robin_hood.h>
#include <queue>
#include <buildConfig.hpp>
#include <bitset>

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

			uint64_t threshold = median * 64;

			// Step 2: Create chunks for hash counts exceeding the threshold
			std::vector<HashChunk> hashChunks;
			for (const auto& [taxid, count] : hashCount) {
				if (count > threshold) {
					int numChunks = static_cast<int>(std::ceil(static_cast<double>(count) / threshold));
					uint64_t chunkSize = count / numChunks;
					for (int i = 0; i < numChunks; ++i) {
						uint64_t currentChunkSize = (i == numChunks - 1) ? (count - chunkSize * (numChunks - 1)) : chunkSize;
						hashChunks.push_back({ taxid, currentChunkSize });
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
			size_t groupNum = (hashChunks.size() + maxTaxidsPerGroup - 1) / maxTaxidsPerGroup;
			std::vector<Group> groups(groupNum);

			auto cmp = [&](const int a, const int b) -> bool {
				return groups[a].totalHash > groups[b].totalHash;
				};
			std::priority_queue<int, std::vector<int>, decltype(cmp)> minHeap(cmp);
			for (size_t i = 0; i < groupNum; ++i) {
				minHeap.push(static_cast<int>(i));
			}

			// Assign each chunk to the group with the lowest total hash count
			for (const auto& chunk : hashChunks) {
				int groupIndex = minHeap.top();
				minHeap.pop();

				groups[groupIndex].taxids.push_back(chunk.taxid);
				groups[groupIndex].totalHash += chunk.hashCount;

				minHeap.push(groupIndex);
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
			//std::cout << "Before Split Groups size: " << hashCount.size() / 16 << std::endl;


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
			size_t median = sortedHashCounts[taxidNum / 2].second;
			size_t groupNum = (taxidNum + maxTaxidsPerGroup) / maxTaxidsPerGroup;
			uint64_t threshold = median * maxTaxidsPerGroup;
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
		}

		template <class Archive>
		void serialize(Archive& ar) {
			ar(data, binNum, binSize, tagNum, MaxCuckooCount, hashSize);
		}

		/**
		 * @brief 轻量主哈希索引。
		 *
		 * 使用若干位运算完成快速搅拌，并利用 `hashSize` 为 2 的幂这一前提直接按位掩码。
		 * 相比重型哈希大幅减少运算量，同时仍保持良好的分布特性。
		 */
		inline size_t hashIndex(uint64_t value) {
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
		inline size_t altHash(size_t b, uint16_t fingerprint) {
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
		inline uint16_t reduceTo12bitAndAddIndex(size_t value, size_t index)
		{
			uint16_t reduced_value = reduceTo12bit(value);
			uint16_t combined = (index << 12) | reduced_value;
			return combined;
		}

		/**
		 * @brief 轻量 12 bit 指纹化。
		 *
		 * 复用主哈希的位混合策略，只截取低 12 bit，并保证不返回 0。
		 */
		inline uint16_t reduceTo12bit(size_t value)
		{
			uint64_t x = value ^ (value >> 33) ^ (value >> 17) ^ (value >> 9);
			uint16_t v = static_cast<uint16_t>(x & 0x0FFFu);
			return v ? v : 1;
		}


		/**
		 * @brief 插入指纹。
		 *
		 * 构造标签后先命中主桶，再跳到可逆备桶，两桶任意存在空位或已有指纹即视为成功；
		 * 若均告满载则进入踢出流程。
		 */
		inline bool insertTag(size_t binIndex, size_t value, size_t index)
		{
			uint16_t tag = reduceTo12bitAndAddIndex(value, index);
			uint16_t fp = (uint16_t)(tag & 0x0FFFu);

			size_t b1 = hashIndex(value);
			size_t b2 = altHash(b1, fp);

			auto try_insert_bucket = [&](size_t b) -> bool {
				size_t bucketPosition = b * binNum + binIndex;
				size_t indexStart = bucketPosition * tagNum * 16;
				uint64_t q = data.get_int(indexStart, 64);
				for (int i = 0; i < 4; ++i) {
					uint16_t chunk = (uint16_t)((q >> (i * 16)) & 0xFFFFu);
					if (chunk == 0u) {
						q &= ~((uint64_t)0xFFFFu << (i * 16));
						q |= ((uint64_t)tag << (i * 16));
						data.set_int(indexStart, q, 64);
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
				size_t bucketPosition = b * binNum + binIndex;
				size_t indexStart = bucketPosition * tagNum * 16;
				uint64_t q = data.get_int(indexStart, 64);

				for (int i = 0; i < 4; ++i) {
					uint16_t chunk = (uint16_t)((q >> (i * 16)) & 0xFFFFu);
					if (chunk == 0u) {
						q &= ~((uint64_t)0xFFFFu << (i * 16));
						q |= ((uint64_t)cur << (i * 16));
						data.set_int(indexStart, q, 64);
						return true;
					}
					if (chunk == cur) return true;
				}

				int rp = dis(gen);
				uint16_t victim = (uint16_t)((q >> (rp * 16)) & 0xFFFFu);
				q &= ~((uint64_t)0xFFFFu << (rp * 16));
				q |= ((uint64_t)cur << (rp * 16));
				data.set_int(indexStart, q, 64);

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
				// Retrieve data for the first hash position
				size_t bucketPosition1 = hash1 * binNum + binIndex;
				size_t indexStart1 = bucketPosition1 * tagNum * 16;
				uint64_t bucketData1 = data.get_int(indexStart1, 64);

				// Retrieve data for the second hash position
				size_t bucketPosition2 = hash2 * binNum + binIndex;
				size_t indexStart2 = bucketPosition2 * tagNum * 16;
				uint64_t bucketData2 = data.get_int(indexStart2, 64);

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
			// Temporary storage for query results, one bitset per position in the result vector
			std::vector<std::bitset<16>> tmpResult(result.size());

			// Temporary storage for query results, one bitset per position in the result vector
			for (auto value : values)
			{
				// Reset all bitsets for the current minimizer
				for (auto& bitset : tmpResult)
				{
					bitset.reset();
				}

				// Perform the bulk containment query for the current minimizer
				bulkContain(value, tmpResult);

				// Update the result vector based on the positions where the minimizer is present
				for (size_t i = 0; i < result.size(); ++i)
				{
					// Check if the j-th bit is set in the i-th bitset
					for (size_t j = 0; j < 16; ++j)
					{
						// Resize the sub-vector if necessary to accommodate the current tag position
						if (tmpResult[i].test(j))
						{
							if (result[i].size() <= j)
							{
								result[i].resize(j + 1, 0);
							}
							// Increment the count for the current tag at the current position
							result[i][j]++;
						}
					}
				}
			}
		}
	};
}
