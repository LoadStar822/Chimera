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
#include <atomic>

namespace chimera::imcf {
	struct Group {
		std::vector<std::string> taxids;
		uint64_t totalHash{ 0 };
	};

	struct HashChunk {
		std::string taxid;
		uint64_t hashCount;
	};

	inline std::vector<Group> partitionHashCount(const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount, int maxTaxidsPerGroup = 16)
	{
		if (hashCount.empty()) {
			return {};
		}
		if (maxTaxidsPerGroup <= 0) {
			throw std::invalid_argument("IMCF partition: maxTaxidsPerGroup must be positive");
		}

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

		std::sort(hashChunks.begin(), hashChunks.end(),
			[](const HashChunk& a, const HashChunk& b) {
				return a.hashCount > b.hashCount;
			});

		size_t maxTaxids = static_cast<size_t>(maxTaxidsPerGroup);
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

	class InterleavedMergedCuckooFilter {
		typedef kvec_t(int) kvector;
		typedef kvec_t(bool) kvectorBool;
	sdsl::bit_vector data;
	size_t binNum;
	size_t binSize;
	size_t tagNum{ 4 };
	int MaxCuckooCount{ 500 };
	size_t hashSize{};
	std::atomic<uint64_t> insertFailureTotal{ 0 };
	std::atomic<uint64_t> insertFailureSaturated{ 0 };


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
				uint16_t fp = tag & 0x0FFF;
				size_t b1 = hashIndex(value);
				size_t b2 = altHash(b1, fp);

				auto tryInsert = [&](size_t bucket) -> bool {
					size_t bucketPosition = bucket * binNum + binIndex;
					size_t indexStart = bucketPosition * tagNum * 16;
					uint64_t query = data.get_int(indexStart, 64);
					for (int i = 0; i < 4; ++i) {
						uint16_t chunk = (uint16_t)((query >> (i * 16)) & 0xFFFF);
						if (chunk == 0u) {
							query &= ~((uint64_t)0xFFFFu << (i * 16));
							query |= ((uint64_t)tag << (i * 16));
							data.set_int(indexStart, query, 64);
							return true;
						}
						if (chunk == tag) {
							return true;
						}
					}
					return false;
				};

				if (tryInsert(b1)) return true;
				if (tryInsert(b2)) return true;

				auto countMatching = [&](size_t bucket) -> int {
					size_t bucketPosition = bucket * binNum + binIndex;
					size_t indexStart = bucketPosition * tagNum * 16;
					uint64_t query = data.get_int(indexStart, 64);
					int matches = 0;
					for (int i = 0; i < 4; ++i) {
						uint16_t chunk = (uint16_t)((query >> (i * 16)) & 0xFFFF);
						if ((chunk & 0x0FFFu) == fp) {
							++matches;
						}
					}
					return matches;
				};

				int sameFingerprint = countMatching(b1) + countMatching(b2);
				if (sameFingerprint >= 8) {
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
