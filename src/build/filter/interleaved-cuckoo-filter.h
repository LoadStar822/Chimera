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

namespace chimera {
	class InterleavedCuckooFilter {
		typedef kvec_t(int) kvector;
		size_t bins{}; // number of user bins
		size_t tc_bins{}; // number of technical bins
		size_t bin_size{}; // size of each bin
		size_t bin_words{}; // number of 64-bit integers required to store the bins, where each 64-bit integer can store 64 bins
		sdsl::int_vector<8> data; // data structure to store the bins
		int MaxCuckooCount{ 500 }; // maximum number of cuckoo kicks
		size_t TagNum{ 4 }; // number of tags per bin
		size_t hashSize{}; // size of the filter
		kvector result;

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
			this->data = sdsl::int_vector<8>(tc_bins * bin_size * TagNum, 0);
			hashSize = tc_bins * bin_size * TagNum - tc_bins * (TagNum + 1);
		}

		size_t hashIndex(uint64_t value) {
			return value % hashSize;
		}

		size_t altHash(size_t pos, uint8_t tag) {
			return (pos ^ (tag * 0x5bd1e995)) % hashSize;
		}

		inline uint8_t reduce_to_8bit(uint64_t value)
		{
			return static_cast<uint8_t>(value & 0xFF);
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
			auto tag = reduce_to_8bit(value);
			// Calculate the starting index for the current hash function
			indexStart = hashIndex(value) + binIndex * TagNum;
			for (size_t j = 0; j < TagNum; j++)
			{
				// Calculate the current index
				index = indexStart + j;
				if (data[index] == 0 || data[index] == tag)
				{
					// Insert the tag if the current slot is empty
					data[index] = tag;
					for (int i = (index << 3) - 32; i < (index << 3) + 64; i++)
					{
						size_t y = data.get_int(i, 64);
					}

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
			size_t oldIndex{ hashIndex(value) }, oldTag{ tag }, newTag;
			for (int count = 0; count < MaxCuckooCount; count++)
			{
				oldIndex = altHash(oldIndex, oldTag) + binIndex * TagNum;
				for (size_t j = 0; j < TagNum; j++)
				{
					oldIndex++;
					if (data[oldIndex] == 0 || data[oldIndex] == oldTag)
					{
						data[oldIndex] = oldTag;
						return true;
					}
				}
				newTag = data[oldIndex];
				data[oldIndex] = oldTag;
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
		}

		/**
		 * Get the data structure that stores the bins.
		 *
		 * @return The data structure that stores the bins.
		 */
		sdsl::int_vector<8> getdata() {
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
			size_t hash1 = hashIndex(value) << 3;
			size_t hash2 = altHash(hash1, tag) << 3;
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
