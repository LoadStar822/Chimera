/*
 * -----------------------------------------------------------------------------
 * Filename:      interleaved-merged-cuckoo-filter-nosimd.h
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-11-12
 *
 * Last Modified: 2024-11-15
 *
 * Description:
 *  This is the header file of the Interleaved Merged Cuckoo Filter without SIMD,
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

namespace chimera::imcf {
	struct Group {
		std::vector<std::string> taxids;
		uint64_t totalHash{ 0 };
	};

	struct HashChunk {
		std::string taxid;
		uint64_t hashCount;
	};

	inline std::vector<Group> partitionHashCount(const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount, std::string mode, int maxTaxidsPerGroup = 16)
	{
		if (mode == "normal")
		{
			std::vector<uint64_t> counts;
			counts.reserve(hashCount.size());
			for (const auto& [taxid, count] : hashCount) {
				counts.push_back(count);
			}

			std::sort(counts.begin(), counts.end());
			uint64_t median = counts[counts.size() / 2];

			uint64_t threshold = median * 64;

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

			std::sort(hashChunks.begin(), hashChunks.end(),
				[](const HashChunk& a, const HashChunk& b) {
					return a.hashCount > b.hashCount;
				});

			size_t groupNum = (hashChunks.size() + maxTaxidsPerGroup - 1) / maxTaxidsPerGroup;
			std::vector<Group> groups(groupNum);

			auto cmp = [&](const int a, const int b) -> bool {
				return groups[a].totalHash > groups[b].totalHash;
				};
			std::priority_queue<int, std::vector<int>, decltype(cmp)> minHeap(cmp);
			for (size_t i = 0; i < groupNum; ++i) {
				minHeap.push(static_cast<int>(i));
			}

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
			std::vector<std::pair<std::string, uint64_t>> sortedHashCounts;
			sortedHashCounts.reserve(hashCount.size());
			for (const auto& item : hashCount) {
				sortedHashCounts.emplace_back(item.first, item.second);
			}
			std::sort(sortedHashCounts.begin(), sortedHashCounts.end(),
				[](const std::pair<std::string, uint64_t>& a, const std::pair<std::string, uint64_t>& b) {
					return a.second > b.second;
				});
			size_t taxidNum = sortedHashCounts.size();
			size_t median = sortedHashCounts[taxidNum / 2].second;
			size_t groupNum = (taxidNum + maxTaxidsPerGroup) / maxTaxidsPerGroup;
			uint64_t threshold = median * maxTaxidsPerGroup;
			std::vector<Group> groups;
			groups.reserve(groupNum);
			groups.resize(groupNum);

			auto cmp = [&](const int a, const int b) -> bool {
				return groups[a].totalHash > groups[b].totalHash;
				};
			std::priority_queue<int, std::vector<int>, decltype(cmp)> minHeap(cmp);
			for (size_t i = 0; i < groupNum; ++i) {
				minHeap.push(static_cast<int>(i));
			}
			for (const auto& [taxid, count] : sortedHashCounts) {
				int smallestGroup = minHeap.top();
				minHeap.pop();
				if (groups[smallestGroup].taxids.size() < maxTaxidsPerGroup) {
					groups[smallestGroup].taxids.push_back(taxid);
					groups[smallestGroup].totalHash += count;
					minHeap.push(smallestGroup);
				}
				else
				{
					bool assigned = false;
					std::vector<int> tempGroups;
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
					for (int tempGroup : tempGroups) {
						minHeap.push(tempGroup);
					}
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

		inline size_t hashIndex(uint64_t value) {
			uint64_t hash = XXH64(&value, sizeof(value), 0);
			return hash % hashSize;
		}

		inline size_t altHash(size_t hashValue, uint16_t fingerprint) {
			uint64_t data = hashValue ^ (fingerprint * 0x5bd1e995);
			uint64_t hash = XXH64(&data, sizeof(data), 0);
			return hash % hashSize;
		}

		inline uint16_t reduceTo12bitAndAddIndex(size_t value, size_t index)
		{
			uint64_t hash = XXH64(&value, sizeof(value), 0);
			uint16_t reduced_value = static_cast<uint16_t>(hash & 0x0FFF);
			reduced_value = reduced_value == 0 ? 1 : reduced_value;
			uint16_t combined = (index << 12) | reduced_value;
			return combined;
		}

		inline uint16_t reduceTo12bit(size_t value)
		{
			uint64_t hash = XXH64(&value, sizeof(value), 0);
			uint16_t reduced_value = static_cast<uint16_t>(hash & 0x0FFF);
			reduced_value = reduced_value == 0 ? 1 : reduced_value;
			return reduced_value;
		}

		inline bool insertTag(size_t binIndex, size_t value, size_t index)
		{
			uint16_t tag = reduceTo12bitAndAddIndex(value, index);
			size_t b = hashIndex(value);
			size_t bucketPosition = b * binNum + binIndex;
			size_t indexStart = bucketPosition * tagNum * 16;
			uint64_t query = data.get_int(indexStart, 64);
			for (int i = 0; i < 4; i++)
			{
				uint16_t chunk = (query >> (i << 4)) & 0xFFFF;
				if (chunk == 0)
				{
					query &= ~((uint64_t)0xFFFF << (i << 4));
					query |= ((uint64_t)tag << (i << 4));
					data.set_int(indexStart, query, 64);
					return true;
				}
				else if (chunk == tag)
				{
					return true;
				}
			}
			if (!kickOut(binIndex, value, tag))
			{
				throw std::runtime_error("Filter is full. Cannot insert more tags.");
				return false;
			}
			return true;
		}

		inline bool kickOut(size_t binIndex, size_t value, uint16_t tag)
		{
			size_t b = hashIndex(value);
			uint16_t oldTag = tag;
			uint64_t query;
			size_t bucketPosition;
			size_t indexStart;

			static thread_local std::mt19937 gen(std::random_device{}());
			std::uniform_int_distribution<> dis(0, tagNum - 1);

			for (int count = 0; count < MaxCuckooCount; count++)
			{
				bucketPosition = b * binNum + binIndex;
				indexStart = bucketPosition * tagNum * 16;
				query = data.get_int(indexStart, 64);

				for (size_t i = 0; i < tagNum; i++)
				{
					uint16_t chunk = (query >> (i * 16)) & 0xFFFF;
					if (chunk == 0)
					{
						query &= ~((uint64_t)0xFFFF << (i * 16));
						query |= ((uint64_t)oldTag << (i * 16));
						data.set_int(indexStart, query, 64);
						return true;
					}
					else if (chunk == oldTag)
					{
						return true;
					}
				}

				size_t randPos = dis(gen);
				uint16_t existingTag = (query >> (randPos * 16)) & 0xFFFF;

				query &= ~((uint64_t)0xFFFF << (randPos * 16));
				query |= ((uint64_t)oldTag << (randPos * 16));
				data.set_int(indexStart, query, 64);

				oldTag = existingTag;
				b = altHash(b, oldTag & 0x0FFF) % hashSize;
			}
			return false;
		}

		inline void bulkContain(size_t& value, std::vector<std::vector<bool>>& result)
		{
			uint16_t fingerprint = reduceTo12bit(value);
			size_t hash1 = hashIndex(value);
			size_t hash2 = altHash(hash1, fingerprint);

			for (size_t binIndex = 0; binIndex < binNum; binIndex++)
			{
				size_t bucketPosition1 = hash1 * binNum + binIndex;
				size_t indexStart1 = bucketPosition1 * tagNum * 16;
				uint64_t bucketData1 = data.get_int(indexStart1, 64);

				size_t bucketPosition2 = hash2 * binNum + binIndex;
				size_t indexStart2 = bucketPosition2 * tagNum * 16;
				uint64_t bucketData2 = data.get_int(indexStart2, 64);

				for (size_t i = 0; i < tagNum; i++)
				{
					uint16_t entry1 = (bucketData1 >> (i * 16)) & 0xFFFF;
					if (entry1 != 0)
					{
						uint16_t storedFingerprint = entry1 & 0x0FFF;
						uint16_t speciesIndex = entry1 >> 12;
						if (storedFingerprint == fingerprint)
						{
							result[binIndex][speciesIndex] = true;
						}
					}

					uint16_t entry2 = (bucketData2 >> (i * 16)) & 0xFFFF;
					if (entry2 != 0)
					{
						uint16_t storedFingerprint = entry2 & 0x0FFF;
						uint16_t speciesIndex = entry2 >> 12;
						if (storedFingerprint == fingerprint)
						{
							result[binIndex][speciesIndex] = true;
						}
					}
				}
			}
		}


		template <std::ranges::range value_range_t>
		inline void bulkCount(value_range_t&& values, std::vector<std::vector<size_t>>& result)
		{
			std::vector<std::vector<bool>> tmpResult(result.size());
			for (size_t i = 0; i < result.size(); ++i)
			{
				tmpResult[i].resize(result[i].size(), false);
			}

			for (auto value : values)
			{
				for (auto& row : tmpResult)
				{
					std::fill(row.begin(), row.end(), false);
				}

				bulkContain(value, tmpResult);

				for (int i = 0; i < result.size(); i++)
				{
					for (int j = 0; j < result[i].size(); j++)
					{
						if (tmpResult[i][j])
						{
							result[i][j]++;
						}
					}
				}
			}
		}
	};
}