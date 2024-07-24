#ifndef INTERLEAVED_CUCKOO_FILTER_H_
#define INTERLEAVED_CUCKOO_FILTER_H_

#include <cstddef>
#include <stdexcept>
#include <sdsl/int_vector.hpp>
#include <hashutil.h>

namespace interleaved_cuckoo_filter {
	class InterleavedCuckooFilter {
		TwoIndependentMultiplyShift hasher;
		size_t bins{}; // number of user bins
		size_t tc_bins{}; // number of technical bins
		size_t bin_size{}; // size of each bin
		size_t bin_words{}; // number of 64-bit integers required to store the bins, where each 64-bit integer can store 64 bins
		sdsl::int_vector<1> data; // data structure to store the bins
		sdsl::int_vector<64> tags; // data structure to store the tags
		int MaxCuckooCount{ 500 }; // maximum number of cuckoo kicks
		size_t TagNum{ 4 }; // number of tags per bin
		size_t hash_shift{}; // number of bits to shift the hash value
		size_t hash_func_count{ 2 }; // number of hash functions to use
		static constexpr std::array<size_t, 5> hash_seeds{ 13572355802537770549ULL, // 2**64 / (e/2)
			   13043817825332782213ULL, // 2**64 / sqrt(2)
			   10650232656628343401ULL, // 2**64 / sqrt(3)
			   16499269484942379435ULL, // 2**64 / (sqrt(5)/2)
			   4893150838803335377ULL }; // 2**64 / (3*pi/5)
		struct VictimCache {
			size_t index; // index of the victim cache
			size_t tag; // tag value of the victim cache
			bool used; // flag indicating if the victim cache is used
		};
		VictimCache victim;

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
			this->hash_shift = std::countl_zero(bin_size);

			// Calculate tc_bins and bin_words using bitwise operations
			this->bin_words = (bins + 63) >> 6; // Equivalent to ceil(bins / 64)
			this->tc_bins = this->bin_words << 6; // Equivalent to tc_bins * 64

			// Calculate the load factor
			double loadFactor = static_cast<double>(bins * bin_size) / (tc_bins * bin_size * TagNum);
			if (loadFactor > 0.95) this->tc_bins = this->tc_bins << 1; // Equivalent to tc_bins * 2

			// Resize the data structure to store the bins
			this->data = sdsl::int_vector<1>(tc_bins * bin_size * TagNum, 0);
			this->tags = sdsl::int_vector<64>(tc_bins * bin_size * TagNum, 0);
		}

		/**
		* Calculate the hash value for the given index using the specified seed.
		*
		* @param i The index to calculate the hash value for.
		* @param seed The seed value to use for the hash calculation.
		* @return The calculated hash value.
		*/
		size_t hash(size_t i, size_t seed) const
		{
			// Multiply i by the seed
			i *= seed;
			// Ensure that hash_shift is less than 64
			assert(hash_shift < 64);
			// XOR i with the right-shifted value of i by hash_shift
			i ^= i >> hash_shift;
			// Multiply i by a constant value
			i *= 11400714819323198485ULL;
			// Take the modulo of i with bin_size
			i %= bin_size;
			// Multiply i by tc_bins
			i *= tc_bins;
			// Return the calculated hash value
			return i;
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
			for (size_t i = 0; i < hash_func_count; i++)
			{
				// Calculate the starting index for the current hash function
				indexStart = hash(value, hash_seeds[i]) + binIndex * TagNum;
				for (size_t j = 0; j < TagNum; j++)
				{
					// Calculate the current index
					index = indexStart + j;
					if (data[index] == 0 || tags[index] == value)
					{
						// Insert the tag if the current slot is empty
						data[index] = 1;
						tags[index] = value;
						return true;
					}

					else
					{
						// Continue to the next slot if the current bin is occupied
						continue;
					}
				}
			}
			// Kick out a tag if all bins are occupied
			if (!kickOut(binIndex, value))
			{
				throw std::runtime_error("Filter is full. Cannot insert more tags.");
				return false;
			}
			return true;
		}

		/**
		 * Kick out a tag to an empty slot in the specified bin of the Interleaved Cuckoo Filter.
		 *
		 * @param binIndex The index of the bin to kick out the tag from.
		 * @param value The value of the tag to kick out.
		 * @return True if the tag was successfully kicked out to an empty slot, false otherwise.
		 */
		bool kickOut(size_t binIndex, size_t value)
		{
			size_t oldIndex, oldTag{ value }, newTag;
			for (int count = 0; count < MaxCuckooCount; count++)
			{
				// Calculate the index for the current cuckoo kick
				oldIndex = hash(oldTag, hash_seeds[count % hash_func_count]) + binIndex * TagNum + count % TagNum;
				// Check if the current slot is occupied
				if (data[oldIndex] == 1)
				{
					// Perform a cuckoo kick by swapping the tags
					newTag = tags[oldIndex];
					tags[oldIndex] = oldTag;
					oldTag = newTag;
					continue;
				}
				else
				{
					// Insert the tag into the current slot
					tags[oldIndex] = oldTag;
					data[oldIndex] = 1;
					return true;
				}
			}
			return false;
		}

		bool lookupTag(size_t binIndex, size_t value)
		{
			assert(binIndex < bins);
			size_t indexStart, index;
			for (size_t i = 0; i < hash_func_count; i++)
			{
				// Calculate the starting index for the current hash function
				indexStart = hash(value, hash_seeds[i]) + binIndex * TagNum;
				for (size_t j = 0; j < TagNum; j++)
				{
					// Calculate the current index
					index = indexStart + j;
					if (data[index] == 1 && tags[index] == value)
					{
						return true;
					}
				}
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
			ar(hash_shift);
			ar(hash_func_count);
		}

		/**
		 * Get the data structure that stores the bins.
		 *
		 * @return The data structure that stores the bins.
		 */
		sdsl::int_vector<1> getdata() {
			return data;
		}

		void generateHash(const size_t& value, size_t* index, uint32_t* tag)
		{
			const uint64_t hash = hasher(value);
		}
	};

	class TwoIndependentMultiplyShift {
		//unsigned __int128 multiply_, add_;

		int64_t multiply_, add_;

	public:
		TwoIndependentMultiplyShift() {
			::std::random_device random;
			for (auto v : { &multiply_, &add_ }) {
				*v = random();
				for (int i = 1; i <= 4; ++i) {
					*v = *v << 32;
					*v |= random();
				}
			}
		}

		uint64_t operator()(uint64_t key) const {
			return (add_ + multiply_ * static_cast<decltype(multiply_)>(key)) >> 64;
		}
	};
} // namespace interleaved_cuckoo_filter

#endif // INTERLEAVED_CUCKOO_FILTER_H_
