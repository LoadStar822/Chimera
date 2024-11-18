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
#include <robin_hood.h>
#include <queue>
#include <buildConfig.hpp>
#include <bitset>

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
			binSize = maxTotalHash / (config.loadFactor * tagNum);
			data = sdsl::bit_vector(binNum * binSize * tagNum * 16, 0);
			hashSize = binSize;
			config.binNum = binNum;
			config.binSize = binSize;
		}

		template <class Archive>
		void serialize(Archive& ar) {
			ar(data, binNum, binSize, tagNum, MaxCuckooCount, hashSize);
		}

		/**
		 * @brief Computes the hash index for a given 64-bit value.
		 *
		 * This function generates a 64-bit hash of the input value using the XXH64 hashing algorithm
		 * and then maps the hash to an index within a predefined hash table size (`hashSize`) by
		 * taking the modulus of the hash value.
		 *
		 * @param value The 64-bit unsigned integer value to be hashed.
		 *
		 * @return The computed hash index as a `size_t`, representing the position in the hash table.
		 *
		 * @details
		 * - The function uses the XXH64 hashing algorithm, known for its speed and good distribution properties.
		 * - The hash is generated with a seed value of 0, which can be adjusted if needed for different hash distributions.
		 * - The final index is obtained by taking the modulus of the hash with `hashSize`, ensuring the index falls within
		 *   the bounds of the hash table.
		 *
		 * @note
		 * - Ensure that `hashSize` is properly defined and represents the size of the hash table to avoid out-of-bounds errors.
		 * - The choice of seed value in XXH64 can be modified based on specific requirements or to achieve different hashing behaviors.
		 */
		inline size_t hashIndex(uint64_t value) {
			uint64_t hash = XXH64(&value, sizeof(value), 0);
			return hash % hashSize;
		}


		/**
		 * @brief Computes an alternative hash index based on a given hash value and fingerprint.
		 *
		 * The primary function of this method is to enable two hash values to mutually compute each other under the same tag.
		 * This is particularly useful for calculating backup positions during kick-out operations in Interleaved Merged Cuckoo Filters.
		 *
		 * This function combines a base hash value with a fingerprint to produce a new hash index. The fingerprint
		 * is used to perturb the original hash value, providing a way to generate multiple hash indices from a single
		 * base hash. The resulting value is then hashed using XXH64 and mapped to an index within the hash table size (`hashSize`)
		 * by taking the modulus.
		 *
		 * @param hashValue The original hash value (typically obtained from a primary hash function).
		 * @param fingerprint A 16-bit unsigned integer used to modify the original hash value, allowing for alternative hashing.
		 *
		 * @return The computed alternative hash index as a `size_t`, representing a secondary position in the hash table.
		 *
		 * @details
		 * - The function first combines `hashValue` with the `fingerprint` by XOR-ing `hashValue` with `fingerprint` multiplied by a
		 *   constant (0x5bd1e995). This mixing step ensures that different fingerprints result in distinct hash modifications.
		 * - The mixed value (`data`) is then hashed using the XXH64 algorithm with a seed value of 0 to generate a new 64-bit hash.
		 * - The final index is obtained by taking the modulus of the new hash with `hashSize`, ensuring the index falls within
		 *   the bounds of the hash table.
		 * - By enabling two hash values to compute each other under the same tag, `altHash` facilitates the calculation of
		 *   backup positions during kick-out operations in Interleaved Merged Cuckoo Filters. This ensures efficient collision
		 *   resolution and optimal placement of entries within the filter.
		 *
		 * @note
		 * - Ensure that `hashSize` is properly defined and represents the size of the hash table to avoid out-of-bounds errors.
		 * - The multiplication by the constant `0x5bd1e995` is intended to introduce a pseudo-random transformation to the fingerprint,
		 *   enhancing the distribution of the resulting hash indices.
		 * - This function is useful in scenarios where multiple hash indices are needed for a single key, such as in cuckoo hashing.
		 * - Specifically, `altHash` enables two hash values to mutually compute each other under the same tag, facilitating the calculation of backup positions during kick-out operations in Interleaved Merged Cuckoo Filters.
		 */
		inline size_t altHash(size_t hashValue, uint16_t fingerprint) {
			uint64_t data = hashValue ^ (fingerprint * 0x5bd1e995);
			uint64_t hash = XXH64(&data, sizeof(data), 0);
			return hash % hashSize;
		}


		/**
		 * @brief Reduces a hashed value to 12 bits and embeds a species index into the upper 4 bits.
		 *
		 * This function performs the following steps:
		 * 1. Hashes the input `value` using the XXH64 hashing algorithm to produce a 64-bit hash.
		 * 2. Reduces the 64-bit hash to a 12-bit value by applying a bitmask.
		 * 3. Ensures that the reduced value is not zero by setting it to `1` if it is.
		 * 4. Embeds the provided `index` into the upper 4 bits of the resulting 16-bit value.
		 *
		 * The final 16-bit `combined` value encodes both the species index and the reduced hash value, where:
		 * - The upper 4 bits represent the `index` (allowing for up to 16 species).
		 * - The lower 12 bits represent the reduced hash value (ranging from 1 to 4095).
		 *
		 * @param value The input value to be hashed and reduced. Typically, this could be a k-mer or similar entity.
		 * @param index The species index to be embedded into the upper 4 bits of the resulting value. Must be less than 16.
		 *
		 * @return A 16-bit unsigned integer combining the species index and the reduced hash value.
		 *
		 * @details
		 * - Hashing: Utilizes the XXH64 algorithm for hashing, known for its speed and effective distribution.
		 * - Reduction to 12 Bits:
		 *   - Applies a bitmask (`0x0FFF`) to retain only the lower 12 bits of the hash.
		 *   - If the reduced value is `0`, it is set to `1` to ensure a non-zero hash value.
		 * - Embedding Species Index:
		 *   - Shifts the `index` left by 12 bits to occupy the upper 4 bits of the 16-bit `combined` value.
		 *   - Combines the shifted `index` with the reduced hash value using the bitwise OR operation.
		 *
		 * @note
		 * - The `index` should be within the range [0, 15] to fit into 4 bits.
		 * - Ensure that `hashSize` is appropriately defined elsewhere in the code to prevent overflow or unintended behavior.
		 */
		inline uint16_t reduceTo12bitAndAddIndex(size_t value, size_t index)
		{
			uint64_t hash = XXH64(&value, sizeof(value), 0);
			uint16_t reduced_value = static_cast<uint16_t>(hash & 0x0FFF);
			reduced_value = reduced_value == 0 ? 1 : reduced_value;
			uint16_t combined = (index << 12) | reduced_value;
			return combined;
		}

		/**
		 * @brief Reduces a hashed value to 12 bits.
		 *
		 * This function performs the following steps:
		 * 1. Hashes the input `value` using the XXH64 hashing algorithm to produce a 64-bit hash.
		 * 2. Reduces the 64-bit hash to a 12-bit value by applying a bitmask.
		 * 3. Ensures that the reduced value is not zero by setting it to `1` if it is.
		 *
		 * The resulting 12-bit `reduced_value` ranges from 1 to 4095.
		 *
		 * @param value The input value to be hashed and reduced. Typically, this could be a k-mer or similar entity.
		 *
		 * @return A 16-bit unsigned integer containing the reduced 12-bit hash value.
		 *
		 * @details
		 * - Hashing: Utilizes the XXH64 algorithm for hashing, known for its speed and effective distribution.
		 * - Reduction to 12 Bits:
		 *   - Applies a bitmask (`0x0FFF`) to retain only the lower 12 bits of the hash.
		 *   - If the reduced value is `0`, it is set to `1` to ensure a non-zero hash value.
		 *
		 * @note
		 * - Ensure that `hashSize` is appropriately defined elsewhere in the code to prevent overflow or unintended behavior.
		 * - This function is useful when a compact representation of the hash value is required, such as in memory-constrained environments.
		 */
		inline uint16_t reduceTo12bit(size_t value)
		{
			uint64_t hash = XXH64(&value, sizeof(value), 0);
			uint16_t reduced_value = static_cast<uint16_t>(hash & 0x0FFF);
			reduced_value = reduced_value == 0 ? 1 : reduced_value;
			return reduced_value;
		}


		/**
		 * @brief Inserts a fingerprint into the Interleaved Merged Cuckoo Filter.
		 *
		 * This function handles the insertion of a fingerprint into the Interleaved Merged Cuckoo Filter (IMCF).
		 * The process involves:
		 * 1. Reducing the computed minimizer hash to 12 bits and embedding the species index into the upper 4 bits
		 *    to form a complete 16-bit fingerprint.
		 * 2. Calculating the hash index to determine the appropriate position within the IMCF.
		 * 3. Checking each of the four slots in the targeted bucket:
		 *    - If a slot is empty, the fingerprint is inserted into that slot.
		 *    - If a slot already contains the same fingerprint, the function acknowledges that the fingerprint
		 *      is already present.
		 * 4. If all four slots are occupied and none contain the fingerprint, the function initiates a kick-out
		 *    operation to relocate existing fingerprints and make space for the new one.
		 *
		 * @param binIndex The index of the bin within the bucket where the fingerprint is to be inserted.
		 * @param value The value (e.g., a k-mer) to be hashed and inserted as a fingerprint.
		 * @param index The species index associated with the fingerprint, which will be embedded into the upper 4 bits
		 *        of the fingerprint.
		 *
		 * @return Returns true if the fingerprint is successfully inserted or already exists in the filter.
		 *         Returns false if the insertion fails due to the filter being full.
		 *
		 * @throws std::runtime_error Throws an exception if the filter is full and the fingerprint cannot be
		 *         inserted after attempting a kick-out operation.
		 *
		 * @details
		 * - Fingerprint Construction:
		 *   The function first calls `reduceTo12bitAndAddIndex` with `value` and `index` to generate a 16-bit fingerprint (tag).
		 *   The lower 12 bits of `tag` represent the reduced minimizer hash, while the upper 4 bits represent the species index.
		 *
		 * - Hash Calculation:
		 *   The primary hash index (b) is computed using `hashIndex(value)`, and the exact bucket position is determined
		 *   using `b * binNum + binIndex`.
		 *
		 * - Slot Inspection and Insertion:
		 *   The function retrieves the current state of the four slots in the bucket and checks for an empty slot or
		 *   a matching fingerprint. If an empty slot is found, the fingerprint is inserted; if the fingerprint is already
		 *   present, the function returns true.
		 *
		 * - Kick-Out Operation:
		 *   If all slots are full and none contain the fingerprint, a kick-out operation is performed to make space.
		 *   If the kick-out fails, an exception is thrown.
		 *
		 * @note Ensure `index` is less than 16 to fit within the 4-bit limit. Proper thread synchronization should be
		 *       maintained if the function is used concurrently.
		 */
		inline bool insertTag(size_t binIndex, size_t value, size_t index)
		{
			// Step 1: Reduce the minimizer hash to 12 bits and add the species index to form a complete 16-bit fingerprint
			uint16_t tag = reduceTo12bitAndAddIndex(value, index);
			// Step 2: Compute the primary hash index for the given value
			size_t b = hashIndex(value); 
			// Step 3: Determine the exact bucket position within the IMCF
			size_t bucketPosition = b * binNum + binIndex; 
			// Step 4: Calculate the starting index for the four slots in the bucket
			size_t indexStart = bucketPosition * tagNum * 16; 
			// Step 5: Retrieve the current state of the four slots in the bucket
			uint64_t query = data.get_int(indexStart, 64); 
			// Step 6: Iterate through each of the four slots to find an empty slot or a matching fingerprint
			for (int i = 0; i < 4; i++)
			{
				// Extract the 16-bit chunk corresponding to the current slot
				uint16_t chunk = (query >> (i << 4)) & 0xFFFF;
				if (chunk == 0)
				{
					// If the slot is empty, insert the fingerprint by updating the corresponding bits
					query &= ~((uint64_t)0xFFFF << (i << 4));
					query |= ((uint64_t)tag << (i << 4));
					data.set_int(indexStart, query, 64);
					return true;
				}
				else if (chunk == tag)
				{
					// If the fingerprint already exists in the slot, acknowledge its presence
					return true;
				}
			}
			// Step 7: If all slots are occupied and none contain the fingerprint, attempt a kick-out operation
			if (!kickOut(binIndex, value, tag))
			{
				// If the kick-out operation fails, throw an exception indicating the filter is full
				throw std::runtime_error("Filter is full. Cannot insert more tags.");
				return false;
			}
			// If the kick-out operation succeeds, return true
			return true;
		}


		/**
		 * @brief Attempts to insert a fingerprint into the Interleaved Merged Cuckoo Filter by kicking out existing entries.
		 *
		 * This function handles cases where all slots in the target bucket are occupied during an insertion. It uses
		 * an alternative hash to find a backup position and checks for an empty slot. If no empty slot is found,
		 * a random slot is selected, and its fingerprint is evicted to make space for the new one. The evicted fingerprint
		 * is then inserted in subsequent iterations, repeating this process up to `MaxCuckooCount` times. If an empty
		 * slot is not found after these attempts, the function returns false to indicate that the filter is full.
		 *
		 * @param binIndex The index of the bin within the bucket where the operation starts.
		 * @param value The original value (e.g., a k-mer) whose hash determines the starting position.
		 * @param tag The 16-bit fingerprint to be inserted or relocated.
		 *
		 * @return Returns true if the fingerprint is successfully inserted or found. Returns false if the filter is
		 *         full after the maximum number of kick-out attempts (`MaxCuckooCount`).
		 *
		 * @details
		 * - The function calculates an alternative hash using `altHash` to find a backup bucket position.
		 * - The function iterates through each slot in the backup bucket to check for an empty slot or matching fingerprint.
		 * - If a slot is full, it randomly selects a slot to evict and inserts the new fingerprint. The evicted fingerprint
		 *   is then used in the next iteration.
		 * - If the function reaches `MaxCuckooCount` attempts without finding an empty slot, it returns false.
		 *
		 * @note Ensure `hashSize` is appropriately defined to map hash values correctly. Use proper synchronization if
		 *       used in a multi-threaded environment.
		 */
		inline bool kickOut(size_t binIndex, size_t value, uint16_t tag)
		{
			// Step 1: Compute the primary hash index for the given value
			size_t b = hashIndex(value);
			uint16_t oldTag = tag;
			uint64_t query;
			size_t bucketPosition;
			size_t indexStart;

			// Thread-local random number generator for selecting a random slot
			static thread_local std::mt19937 gen(std::random_device{}());
			std::uniform_int_distribution<> dis(0, tagNum - 1);

			// Step 2: Attempt to insert the fingerprint up to MaxCuckooCount times
			for (int count = 0; count < MaxCuckooCount; count++)
			{
				// Compute an alternative hash to find a backup bucket position
				bucketPosition = b * binNum + binIndex;
				indexStart = bucketPosition * tagNum * 16;
				query = data.get_int(indexStart, 64);

				// Step 3: Check each slot for an empty space or a matching fingerprint
				for (size_t i = 0; i < tagNum; i++)
				{
					uint16_t chunk = (query >> (i * 16)) & 0xFFFF;
					if (chunk == 0)
					{
						// Insert the fingerprint into the empty slot and return true
						query &= ~((uint64_t)0xFFFF << (i * 16));
						query |= ((uint64_t)oldTag << (i * 16));
						data.set_int(indexStart, query, 64);
						return true;
					}
					else if (chunk == oldTag)
					{
						// Fingerprint is already present, return true
						return true;
					}
				}

				// Step 4: Randomly evict a slot and insert the new fingerprint
				size_t randPos = dis(gen);
				uint16_t existingTag = (query >> (randPos * 16)) & 0xFFFF;

				query &= ~((uint64_t)0xFFFF << (randPos * 16));
				query |= ((uint64_t)oldTag << (randPos * 16));
				data.set_int(indexStart, query, 64);
				// Use the evicted fingerprint in the next iteration
				oldTag = existingTag;
				b = altHash(b, oldTag & 0x0FFF) % hashSize;
			}
			// Return false if insertion fails after MaxCuckooCount attempts
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

			// Step 3: Prepare SIMD vectors for comparison
			simde__m256i fingerprint_vec = simde_mm256_set1_epi16(fingerprint);

			simde__m256i fingerprint_mask = simde_mm256_set1_epi16(0x0FFF);

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

				// Load data into SIMD vectors
				simde__m256i bucket_vec = simde_mm256_loadu_si256(reinterpret_cast<const simde__m256i*>(entries));

				// Step 5: Extract fingerprints and species indices
				simde__m256i fingerprints = simde_mm256_and_si256(bucket_vec, fingerprint_mask);
				simde__m256i species_indices = simde_mm256_srli_epi16(bucket_vec, species_shift);

				// Compare extracted fingerprints with the input fingerprint
				simde__m256i cmp = simde_mm256_cmpeq_epi16(fingerprints, fingerprint_vec);

				// Step 6: Update the result vector with matching species indices
				for (int i = 0; i < 8; ++i)
				{
					uint16_t cmp_result = simde_mm256_extract_epi16(cmp, i);
					if (cmp_result == 0xFFFF)	// If a match is found
					{
						uint16_t speciesIndex = simde_mm256_extract_epi16(species_indices, i);
						result[binIndex].set(speciesIndex);	// Mark the matching species index
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