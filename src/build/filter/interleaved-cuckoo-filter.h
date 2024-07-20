#ifndef INTERLEAVED_CUCKOO_FILTER_H_
#define INTERLEAVED_CUCKOO_FILTER_H_

#include <cstddef>
#include <stdexcept>
#include <sdsl/int_vector.hpp>
#include <hashutil.h>

namespace interleaved_cuckoo_filter {

class InterleavedCuckooFilter {
	size_t bins{}; // number of user bins
	size_t tc_bins{}; // number of technical bins
	size_t bin_size{}; // size of each bin
	size_t bin_words{}; // number of 64-bit integers required to store the bins, where each 64-bit integer can store 64 bins
	sdsl::int_vector<1> data; // data structure to store the bins
	size_t kMaxCuckooCount = 500;

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
		this->bin_words = (bins + 63) >> 6; // Equivalent to ceil(bins / 64)
		this->tc_bins = this->bin_words << 6; // Equivalent to tc_bins * 64
		// Resize the data structure to store the bins
		this->data = sdsl::int_vector<1>(tc_bins * bin_size, 0);
	}

    void insert(size_t binIndex, size_t value) 
	{
        assert(binIndex < bins);
            
        size_t binStart = binIndex * bin_size;
        size_t binEnd = binStart + bin_size;
            
        for (size_t i = binStart; i < binEnd; i++) {
            data[i] = value;
        }
    }
};

} // namespace interleaved_cuckoo_filter

#endif // INTERLEAVED_CUCKOO_FILTER_H_
