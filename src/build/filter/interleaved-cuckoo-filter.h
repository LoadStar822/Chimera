/*
 * -----------------------------------------------------------------------------
 * Filename:      interleaved-cuckoo-filter.h
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-07-09
 *
 * Last Modified: 2024-07-26
 *
 * Description:
 *  This is the header file of the Interleaved Cuckoo Filter,
 *	which contains the basic operations of the Interleaved Cuckoo Filter
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#ifndef INTERLEAVED_CUCKOO_FILTER_H_
#define INTERLEAVED_CUCKOO_FILTER_H_

#include <cstddef>
#include <stdexcept>
#include <sdsl/int_vector.hpp>
#include <hashutil.h>
#include <kvec.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>

namespace chimera {
	class InterleavedCuckooFilter {
		typedef kvec_t(int) kvector;
		size_t bins{}; // number of user bins
		size_t tc_bins{}; // number of technical bins
		size_t bin_size{}; // size of each bin
		size_t bin_words{}; // number of 64-bit integers required to store the bins, where each 64-bit integer can store 64 bins
		sdsl::bit_vector data; // data structure to store the bins
		int MaxCuckooCount{ 500 }; // maximum number of cuckoo kicks
		size_t TagNum{ 4 }; // number of tags per bin
		size_t hashSize{}; // size of the filter

	public:
		InterleavedCuckooFilter() = default;
		InterleavedCuckooFilter(size_t bins, size_t bin_size)
		{
			if (bins <= 0) {
				throw std::invalid_argument("Invalid number of bins. Must be greater than 0.");
			}
			if (bin_size <= 0) {
				throw std::invalid_argument("Invalid bin size. Must be greater than 0.");
			}

			this->bins = bins;
			this->bin_size = bin_size;

			// Calculate tc_bins and bin_words using bitwise operations
			this->bin_words = (bins + 1) >> 1; // Equivalent to ceil(bins / 2)
			this->tc_bins = this->bin_words << 1; // Equivalent to tc_bins * 64

			//// Calculate the load factor
			//double loadFactor = static_cast<double>(bins * bin_size) / (tc_bins * bin_size * TagNum);
			//if (loadFactor > 0.95) this->tc_bins = this->tc_bins << 1; // Equivalent to tc_bins * 2

			// Resize the data structure to store the bins
			this->data = sdsl::bit_vector(tc_bins * bin_size * TagNum * 8, 0);
			hashSize = tc_bins * bin_size * TagNum - tc_bins * (TagNum + 1);
		}

		// 使用掩码批量插入八位数到 bit_vector
		void batch_insert_to_bit_vector(uint8_t value, size_t position) {
			uint64_t mask = static_cast<uint64_t>(value) << (position % 64);  // 将8位数移到正确位置
			size_t idx = position / 64;
			data.data()[idx] |= mask;  // 批量写入
			if ((position % 64) > 56) {
				// 处理跨越64位边界的情况
				data.data()[idx + 1] |= (value >> (64 - (position % 64)));
			}
		}

		// 查询 bit_vector 中的值
		uint8_t query_bit_vector(size_t position) {
			uint64_t mask = 0xFFULL << (position % 64);
			size_t idx = position / 64;
			uint64_t chunk = (data.data()[idx] & mask) >> (position % 64);
			if ((position % 64) > 56) {
				chunk |= (data.data()[idx + 1] & 0xFFULL) << (64 - (position % 64));
			}
			return static_cast<uint8_t>(chunk);
		}

		size_t hashIndex(uint64_t value) {
			return value % hashSize;
		}

		size_t altHash(size_t pos, uint8_t tag) {
			return (pos ^ (tag * 0x5bd1e995)) % hashSize;
		}

		inline uint8_t reduce_to_8bit(uint64_t value)
		{
			uint8_t reduced_value = static_cast<uint8_t>(((value * 2654435761U) >> 24) & 0xFF);
			return reduced_value == 0 ? 1 : reduced_value; // 确保返回值大于0
		}

		/**
		 * Insert a tag into the specified bin of the Interleaved Cuckoo Filter.
		 *
		 * @param binIndex The index of the bin to insert the tag into.
		 * @param value The value of the tag to insert.
		 * @return True if the tag was successfully inserted, false otherwise.
		 */
		bool insertTag(size_t binIndex, size_t value)
		{
			assert(binIndex < bins);
			size_t indexStart, index;
			uint8_t query;
			auto tag = reduce_to_8bit(value);
			// Calculate the starting index for the current hash function
			indexStart = hashIndex(value) + binIndex * TagNum;
			for (size_t j = 0; j < TagNum; j++)
			{
				// Calculate the current index
				index = (indexStart + j) << 3;
				query = query_bit_vector(index);
				if (query == 0 || query == tag)
				{
					batch_insert_to_bit_vector(tag, index);

					return true;
				}
			}
			// Kick out a tag if all bins are occupied
			if (!kickOut(binIndex, value, tag))
			{
				throw std::runtime_error("Filter is full. Cannot insert more tags.");
				return false;
			}
			return true;
		}

		bool kickOut(size_t binIndex, size_t value, uint8_t tag)
		{
			size_t oldIndex{ hashIndex(value) }, oldTag{ tag }, newIndex, newTag, index;
			uint8_t query;
			for (int count = 0; count < MaxCuckooCount; count++)
			{
				oldIndex = altHash(oldIndex, oldTag) + binIndex * TagNum;
				for (size_t j = 0; j < TagNum; j++)
				{
					index = (oldIndex + j) << 3;
					query = query_bit_vector(index);
					if (query == 0 || query == oldTag)
					{
						batch_insert_to_bit_vector(oldTag, index);
						return true;
					}
				}
				newIndex = (oldIndex + count % TagNum) << 3;
				newTag = query_bit_vector(newIndex);
				batch_insert_to_bit_vector(oldTag, newIndex);
				oldTag = newTag;
			}
			return false;
		}

		/**
		 * Serialize the Interleaved Cuckoo Filter to the specified output stream.
		 *
		 * @param os The output stream to serialize the Interleaved Cuckoo Filter to.
		 */
		template <class Archive>
		void serialize(Archive& ar) {
			ar(bins);
			ar(tc_bins);
			ar(bin_size);
			ar(bin_words);
			ar(data);
			ar(TagNum);
			ar(MaxCuckooCount);
			ar(hashSize);
		}

		/**
		 * Get the data structure that stores the bins.
		 *
		 * @return The data structure that stores the bins.
		 */
		sdsl::bit_vector getdata() {
			return data;
		}

		typedef kvec_t(bool) kvector_bool;
		kvector_bool bulk_contain(size_t value)
		{
			kvector_bool result;
			kv_init(result);
			kv_resize(bool, result, bins);
			std::memset(result.a, 0, sizeof(bool) * bins);
			kv_size(result) = bins;
			uint8_t tag = reduce_to_8bit(value);
			size_t hash1 = hashIndex(value);
			size_t hash2 = altHash(hash1, tag) << 3;
			hash1 = hash1 << 3;
			size_t tmp1{ 0 }, tmp2{ 0 };
			for (size_t batch = 0; batch < bin_words; batch++)
			{
				tmp1 = data.get_int(hash1);
				tmp2 = data.get_int(hash2);
				hash1 += 64;
				hash2 += 64;
				for (int bit = 0; bit < 64; bit += 8)
				{
					int bin_index = (batch * 2) + (bit / 32);
					uint8_t tag_bits1 = (tmp1 >> bit) & 0xFF;
					uint8_t tag_bits2 = (tmp2 >> bit) & 0xFF;
					if (tag_bits1 == tag || tag_bits2 == tag)
					{
						kv_A(result, bin_index) = true;
					}
				}
			}
			return result;
		}

		template <std::ranges::range value_range_t>
		kvector bulk_count(value_range_t&& values)
		{
			kvector result;
			kv_init(result);
			kv_resize(int, result, bins);
			std::memset(result.a, 0, sizeof(int) * bins);
			kv_size(result) = bins;
			for (auto value : values)
			{
				kvector_bool tmp = bulk_contain(value);
				for (size_t i = 0; i < tmp.n; i++)
				{
					kv_A(result, i) += kv_A(tmp, i);
				}
				kv_destroy(tmp);
			}
			return result;
		}
	};
} // namespace chimera

#endif // INTERLEAVED_CUCKOO_FILTER_H_
