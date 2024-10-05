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
 * Last Modified: 2024-10-05
 *
 * Description:
 *  This is the header file of the Interleaved Cuckoo Filter,
 *	which contains the basic operations of the Interleaved Cuckoo Filter
 *
 * Version:
 *  1.3
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
#include <xxhash.h>
#include <simde/x86/avx2.h>
#include <vector>
#include <cstdint>

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
		size_t bitNum{ 16 };

		friend std::ostream& operator<<(std::ostream& os, const InterleavedCuckooFilter& filter) {
			os << "InterleavedCuckooFilter Details:\n";
			os << "  User Bins: " << filter.bins << "\n";
			os << "  Technical Bins: " << filter.tc_bins << "\n";
			os << "  Bin Size: " << filter.bin_size << " bits\n";
			os << "  Bin Words: " << filter.bin_words << " 64-bit words\n";
			os << "  Max Cuckoo Count: " << filter.MaxCuckooCount << "\n";
			os << "  Tags per Bin: " << filter.TagNum << "\n";
			os << "  Hash Size: " << filter.hashSize << " bits\n";
			os << "  Data Size: " << filter.data.size() << " bits\n";

			return os;
		}

	public:
		InterleavedCuckooFilter() = default;
		InterleavedCuckooFilter(size_t bins, size_t bin_size, int bitNum)
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
			this->bitNum = bitNum;
			this->data = sdsl::bit_vector(tc_bins * bin_size * TagNum * bitNum, 0);
			hashSize = tc_bins * bin_size * TagNum - tc_bins * (TagNum + 1);
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: batch_insert_to_bit_vector
		 *
		 * Description:
		 * Inserts an 16-bit value into the bit_vector at the specified position using a mask.
		 *c
		 * Parameters:
		 * - value: The 16-bit value to be inserted.
		 * - position: The position in the bit_vector where the value should be inserted.
		 *
		 * Returns: void
		 * -----------------------------------------------------------------------------------------------
		 */
		inline void batch_insert_to_bit_vector_16bit(uint16_t value, size_t position) {
			size_t idx = position / 64;
			size_t offset = position % 64;

			data.data()[idx] |= static_cast<uint64_t>(value) << offset;

			if (offset > 64 - bitNum) {
				size_t bits_in_first_word = 64 - offset;
				uint16_t overflow_value = value >> bits_in_first_word;
				data.data()[idx + 1] |= static_cast<uint64_t>(overflow_value);
			}
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: batch_insert_to_bit_vector
		 *
		 * Description:
		 * Inserts an 8-bit value into the bit_vector at the specified position using a mask.
		 *c
		 * Parameters:
		 * - value: The 8-bit value to be inserted.
		 * - position: The position in the bit_vector where the value should be inserted.
		 *
		 * Returns: void
		 * -----------------------------------------------------------------------------------------------
		 */
		void batch_insert_to_bit_vector_8bit(uint8_t value, size_t position) {
			uint64_t mask = static_cast<uint64_t>(value) << (position % 64);  // Move the 8-bit value to the correct position
			size_t idx = position / 64;
			data.data()[idx] |= mask;  // Batch write the value
			if ((position % 64) > 56) {
				// Handle the case where the value crosses a 64-bit boundary
				data.data()[idx + 1] |= (value >> (64 - (position % 64)));
			}
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: query_bit_vector
		 *
		 * Description:
		 * Queries the value in the bit_vector at the specified position.
		 *
		 * Parameters:
		 * - position: The position in the bit_vector to query.
		 *
		 * Returns: The 16-bit value at the specified position.
		 * -----------------------------------------------------------------------------------------------
		 */
		inline uint16_t query_bit_vector_16bit(size_t position) {
			size_t idx = position / 64;
			size_t offset = position % 64;

			uint64_t chunk = data.data()[idx] >> offset;

			if (offset > 64 - bitNum) {
				size_t bits_in_first_word = 64 - offset;
				uint64_t next_chunk = data.data()[idx + 1];
				chunk |= next_chunk << bits_in_first_word;
			}

			uint16_t mask = (1ULL << bitNum) - 1;  // 0xFFFF
			return static_cast<uint16_t>(chunk & mask);
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: query_bit_vector
		 *
		 * Description:
		 * Queries the value in the bit_vector at the specified position.
		 *
		 * Parameters:
		 * - position: The position in the bit_vector to query.
		 *
		 * Returns: The 8-bit value at the specified position.
		 * -----------------------------------------------------------------------------------------------
		 */
		inline uint8_t query_bit_vector_8bit(size_t position) {
			uint64_t mask = 0xFFULL << (position % 64);
			size_t idx = position / 64;
			uint64_t chunk = (data.data()[idx] & mask) >> (position % 64);
			if ((position % 64) > 56) {
				chunk |= (data.data()[idx + 1] & 0xFFULL) << (64 - (position % 64));
			}
			return static_cast<uint8_t>(chunk);
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: hashIndex
		 *
		 * Description:
		 * Calculates the hash index for a given value.
		 *
		 * Parameters:
		 * - value: The value for which to calculate the hash index.
		 *
		 * Returns: The hash index.
		 * -----------------------------------------------------------------------------------------------
		 */
		inline size_t hashIndex(uint64_t value) {
			return value % hashSize;
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: altHash
		 *
		 * Description:
		 * Calculates the alternative hash for a given position and tag.
		 *
		 * Parameters:
		 * - pos: The position for which to calculate the alternative hash.
		 * - tag: The tag for which to calculate the alternative hash.
		 *
		 * Returns: The alternative hash.
		 * -----------------------------------------------------------------------------------------------
		 */
		inline size_t altHash(size_t pos, uint8_t tag) {
			return (pos ^ (tag * 0x5bd1e995)) % hashSize;
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: reduce_to_16bit
		 *
		 * Description:
		 * Reduces a 64-bit value to an 16-bit value.
		 *
		 * Parameters:
		 * - value: The 64-bit value to be reduced.
		 *
		 * Returns: The reduced 16-bit value.
		 * -----------------------------------------------------------------------------------------------
		 */
		inline uint16_t reduce_to_16bit(uint64_t value)
		{
			uint64_t hash = XXH64(&value, sizeof(value), 0);
			uint16_t reduced_value = static_cast<uint16_t>(hash & 0xFFFF);
			return reduced_value == 0 ? 1 : reduced_value;
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: reduce_to_8bit
		 *
		 * Description:
		 * Reduces a 64-bit value to an 8-bit value.
		 *
		 * Parameters:
		 * - value: The 64-bit value to be reduced.
		 *
		 * Returns: The reduced 8-bit value.
		 * -----------------------------------------------------------------------------------------------
		 */
		inline uint8_t reduce_to_8bit(uint64_t value)
		{
			uint64_t hash = XXH64(&value, sizeof(value), 0);
			uint8_t reduced_value = static_cast<uint8_t>(hash & 0xFF);
			return reduced_value == 0 ? 1 : reduced_value;
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: insertTag
		 *
		 * Description:
		 * Inserts a tag into the specified bin of the Interleaved Cuckoo Filter.
		 *
		 * Parameters:
		 * - binIndex: The index of the bin to insert the tag into.
		 * - value: The value of the tag to insert.
		 *
		 * Returns: True if the tag was successfully inserted, false otherwise.
		 * -----------------------------------------------------------------------------------------------
		 */
		bool insertTag(size_t binIndex, size_t value)
		{
			assert(binIndex < bins);
			if (bitNum == 16)
			{
				size_t indexStart, position;
				uint16_t query;
				auto tag = reduce_to_16bit(value);
				indexStart = (hashIndex(value) + binIndex * TagNum) * bitNum;
				for (size_t j = 0; j < TagNum; j++)
				{
					position = indexStart + j * bitNum;
					query = query_bit_vector_16bit(position);
					if (query == 0 || query == tag)
					{
						batch_insert_to_bit_vector_16bit(tag, position);
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
			else if (bitNum == 8)
			{
				size_t indexStart, index;
				uint8_t query;
				auto tag = reduce_to_8bit(value);
				// Calculate the starting index for the current hash function
				indexStart = hashIndex(value) + binIndex * TagNum;
				for (size_t j = 0; j < TagNum; j++)
				{
					// Calculate the current index
					index = (indexStart + j) << 3;
					query = query_bit_vector_8bit(index);
					if (query == 0 || query == tag)
					{
						batch_insert_to_bit_vector_8bit(tag, index);

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
			return false;
		}

		/*
		 * -----------------------------------------------------------------------------------------------
		 * Function: kickOut
		 *
		 * Description:
		 * Kicks out a tag from the specified bin and inserts a new tag in its place.
		 *
		 * Parameters:
		 * - binIndex: The index of the bin to kick out the tag from.
		 * - value: The value of the tag to insert.
		 * - tag: The tag to insert.
		 *
		 * Returns: True if the tag was successfully kicked out and inserted, false otherwise.
		 * -----------------------------------------------------------------------------------------------
		 */
		bool kickOut(size_t binIndex, size_t value, uint16_t tag)
		{
			if (bitNum == 16)
			{
				size_t oldIndex{ hashIndex(value) }, oldTag{ tag }, newTag, position;
				uint16_t query;
				for (int count = 0; count < MaxCuckooCount; count++)
				{
					oldIndex = altHash(oldIndex, oldTag) + binIndex * TagNum;
					size_t indexStart = oldIndex * bitNum;
					for (size_t j = 0; j < TagNum; j++)
					{
						position = indexStart + j * bitNum;
						query = query_bit_vector_16bit(position);
						if (query == 0 || query == oldTag)
						{
							batch_insert_to_bit_vector_16bit(oldTag, position);
							return true;
						}
					}
					size_t randPos = rand() % TagNum;
					position = indexStart + randPos * bitNum;
					newTag = query_bit_vector_16bit(position);
					batch_insert_to_bit_vector_16bit(oldTag, position);
					oldTag = newTag;
				}
				return false;
			}
			else if (bitNum == 8)
			{
				size_t oldIndex{ hashIndex(value) }, oldTag{ tag }, newIndex, newTag, index;
				uint8_t query;
				for (int count = 0; count < MaxCuckooCount; count++)
				{
					oldIndex = altHash(oldIndex, oldTag) + binIndex * TagNum;
					for (size_t j = 0; j < TagNum; j++)
					{
						index = (oldIndex + j) << 3;
						query = query_bit_vector_8bit(index);
						if (query == 0 || query == oldTag)
						{
							batch_insert_to_bit_vector_8bit(oldTag, index);
							return true;
						}
					}
					newIndex = (oldIndex + count % TagNum) << 3;
					newTag = query_bit_vector_8bit(newIndex);
					batch_insert_to_bit_vector_8bit(oldTag, newIndex);
					oldTag = newTag;
				}
				return false;
			}
			else
			{
				throw std::runtime_error("Invalid bitNum.");
				return false;
			}
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
			ar(bitNum);
		}

		/**
		* -----------------------------------------------------------------------------------------------
		* Function: bulk_contain
		*
		* Description:
		* Checks if the given value is present in the Interleaved Cuckoo Filter.
		*
		* Parameters:
		* - value: The value to check for.
		*
		* Returns: A kvector_bool indicating the presence of the value in each bin.
		* -----------------------------------------------------------------------------------------------
		*/
		typedef kvec_t(bool) kvector_bool;
		inline void bulk_contain(size_t value, kvector_bool& result)
		{
			if (bitNum == 16)
			{
				uint16_t tag = reduce_to_16bit(value);
				size_t hash1 = hashIndex(value);
				size_t hash2 = altHash(hash1, tag);
				size_t position1 = hash1 << 4;
				size_t position2 = hash2 << 4;

				const size_t stride = (TagNum) << 4;
				size_t binOffset = 0;

				const size_t N = 2;
				size_t bin = 0;

				simde__m256i tag_vec = simde_mm256_set1_epi16((int16_t)tag);

				for (; bin + N <= bins; bin += N)
				{
					uint64_t tag_bits1_0 = data.get_int(position1 + binOffset + 0 * stride);
					uint64_t tag_bits1_1 = data.get_int(position1 + binOffset + 1 * stride);
					uint64_t tag_bits2_0 = data.get_int(position2 + binOffset + 0 * stride);
					uint64_t tag_bits2_1 = data.get_int(position2 + binOffset + 1 * stride);

					simde__m256i tag_bits_vec = simde_mm256_set_epi64x(
						tag_bits2_1, tag_bits1_1, tag_bits2_0, tag_bits1_0);

					simde__m256i cmp_mask = simde_mm256_cmpeq_epi16(tag_bits_vec, tag_vec);

					uint32_t cmp_result = simde_mm256_movemask_epi8(cmp_mask);

					if (cmp_result & 0xFFFF)
					{
						kv_A(result, bin) = true;
					}
					if (cmp_result & 0xFFFF0000)
					{
						kv_A(result, bin + 1) = true;
					}

					binOffset += N * stride;
				}
				for (; bin < bins; ++bin)
				{
					size_t pos1 = position1 + binOffset;
					size_t pos2 = position2 + binOffset;
					uint64_t tag_bits1 = data.get_int(pos1);
					uint64_t tag_bits2 = data.get_int(pos2);

					simde__m128i tag_bits_vec = simde_mm_set_epi64x(tag_bits2, tag_bits1);
					simde__m128i tag_vec_128 = simde_mm_set1_epi16((int16_t)tag);
					simde__m128i cmp_mask = simde_mm_cmpeq_epi16(tag_bits_vec, tag_vec_128);
					int cmp_result = simde_mm_movemask_epi8(cmp_mask);

					if (cmp_result != 0)
					{
						kv_A(result, bin) = true;
					}

					binOffset += stride;
				}
			}
			else if (bitNum == 8)
			{
				uint8_t tag = reduce_to_8bit(value);
				size_t hash1 = hashIndex(value);
				size_t hash2 = altHash(hash1, tag) << 3;
				hash1 = hash1 << 3;
				size_t tmp1{ 0 }, tmp2{ 0 };

				simde__m128i tag_vec = simde_mm_set1_epi8(static_cast<char>(tag));

				for (size_t batch = 0; batch < bin_words; batch++)
				{
					tmp1 = data.get_int(hash1);
					tmp2 = data.get_int(hash2);
					hash1 += 64;
					hash2 += 64;

					simde__m128i data_vec = simde_mm_set_epi64x(static_cast<long long>(tmp2), static_cast<long long>(tmp1));

					simde__m128i cmp = simde_mm_cmpeq_epi8(data_vec, tag_vec);

					int mask = simde_mm_movemask_epi8(cmp);

					bool match0 = (mask & 0x0F0F) != 0;
					bool match1 = (mask & 0xF0F0) != 0;

					kv_A(result, batch * 2 + 0) = match0;
					kv_A(result, batch * 2 + 1) = match1;
				}
			}
		}

		/**
		* -----------------------------------------------------------------------------------------------
		* Function: bulk_count
		*
		* Description:
		* Counts the occurrences of values in the Interleaved Cuckoo Filter.
		*
		* Parameters:
		* - values: A range of values to count.
		*
		* Returns: A kvector indicating the count of each value in each bin.
		* -----------------------------------------------------------------------------------------------
		*/
		template <std::ranges::range value_range_t>
		inline kvector bulk_count(value_range_t&& values)
		{
			kvector result;
			kv_init(result);
			kv_resize(int, result, bins);
			std::memset(result.a, 0, sizeof(int) * bins);
			kv_size(result) = bins;

			kvector_bool tmp_result;
			kv_init(tmp_result);
			kv_resize(bool, tmp_result, bins);
			kv_size(tmp_result) = bins;

			const size_t simd_width = 8; // 256 bits / 32 bits per int
			size_t i = 0;

			for (auto value : values)
			{
				std::fill(tmp_result.a, tmp_result.a + bins, false);
				bulk_contain(value, tmp_result);

				for (i = 0; i + simd_width <= bins; i += simd_width)
				{
					simde__m256i mask = simde_mm256_loadu_si256(reinterpret_cast<const simde__m256i*>(&tmp_result.a[i]));

					simde__m256i ones = simde_mm256_set1_epi8(1);

					simde__m256i increment = simde_mm256_and_si256(mask, ones);

					simde__m256i increment_int = simde_mm256_cvtepi8_epi32(simde_mm256_castsi256_si128(increment));

					simde__m256i current = simde_mm256_loadu_si256(reinterpret_cast<const simde__m256i*>(&result.a[i]));

					simde__m256i updated = simde_mm256_add_epi32(current, increment_int);

					simde_mm256_storeu_si256(reinterpret_cast<simde__m256i*>(&result.a[i]), updated);
				}

				for (; i < bins; ++i)
				{
					result.a[i] += tmp_result.a[i] ? 1 : 0;
				}
			}
			kv_destroy(tmp_result);
			return result;
		}
	};
} // namespace chimera

#endif // INTERLEAVED_CUCKOO_FILTER_H_
