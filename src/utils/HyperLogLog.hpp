/*
 * -----------------------------------------------------------------------------
 * Filename:      HyperLogLog.hpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-10-28
 *
 * Last Modified: 2024-10-28
 *
 * Description:
 * This file defines the HyperLogLog class for Chimera
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>
#include <simde/x86/avx2.h>
#include <xxhash.h>

class HyperLogLog
{
public:
	explicit HyperLogLog(uint8_t bits = 5)
		: bits(bits), size(1ULL << bits), registers(size, 0)
	{
		if (bits < 5 || bits > 32)
			throw std::invalid_argument("Bit width must be in the range [5,32].");

		switch (size)
		{
		case 16:
			alpha = 0.673;
			break;
		case 32:
			alpha = 0.697;
			break;
		case 64:
			alpha = 0.709;
			break;
		default:
			alpha = 0.7213 / (1 + 1.079 / size);
			break;
		}
		correction_factor = alpha * size * size;
	}

	void add(uint64_t value)
	{
		uint64_t hash = XXH64(&value, sizeof(value), 0);
		uint64_t index = hash >> (64 - bits);
		uint8_t rank = leadingZeros((hash << bits) | ((1ULL << bits) - 1)) + 1;
		registers[index] = std::max(registers[index], rank);
	}

	double estimate() const
	{
		double sum = 0.0;
		size_t zero_count = 0;

		// using SIMD to accelerate
		size_t vec_size = registers.size() / 32 * 32; // 32 registers per __m256d
		simde__m256d vec_sum = simde_mm256_setzero_pd();

		const double* lookup = getLookupTable();

		size_t i = 0;
		for (; i < vec_size; i += 32)
		{
			for (int k = 0; k < 32; k += 4)
			{
				simde__m256d v = simde_mm256_set_pd(
					lookup[registers[i + k + 0]],
					lookup[registers[i + k + 1]],
					lookup[registers[i + k + 2]],
					lookup[registers[i + k + 3]]);
				vec_sum = simde_mm256_add_pd(vec_sum, v);
			}

			for (size_t j = i; j < i + 32; ++j)
				if (registers[j] == 0)
					++zero_count;
		}

		// process the remaining elements
		for (; i < size; ++i)
		{
			sum += lookup[registers[i]];
			if (registers[i] == 0)
				++zero_count;
		}

		// accumulate the result of SIMD part
		alignas(32) double temp[4];
		simde_mm256_store_pd(temp, vec_sum);
		sum += temp[0] + temp[1] + temp[2] + temp[3];

		double estimate = correction_factor / sum;

		// small cardinality correction
		if (estimate <= 2.5 * size)
		{
			if (zero_count != 0)
				estimate = size * std::log(static_cast<double>(size) / zero_count);
		}

		return estimate;
	}

	void merge(const HyperLogLog& other)
	{
		assert(bits == other.bits);
		for (size_t i = 0; i < size; ++i)
		{
			registers[i] = std::max(registers[i], other.registers[i]);
		}
	}

	void reset()
	{
		std::fill(registers.begin(), registers.end(), 0);
	}

	// serialize using cereal
	template<class Archive>
	void serialize(Archive& ar)
	{
		ar(bits, size, alpha, correction_factor, registers);
	}

private:
	uint8_t bits;
	uint64_t size;
	double alpha;
	double correction_factor;
	std::vector<uint8_t> registers;

	static uint8_t leadingZeros(uint64_t x)
	{
		return static_cast<uint8_t>(std::countl_zero(x));
	}

	static const double* getLookupTable()
	{
		static double lookup[65] = { 0 };
		static bool initialized = false;
		if (!initialized)
		{
			for (int i = 0; i <= 64; ++i)
				lookup[i] = std::pow(2.0, -static_cast<double>(i));
			initialized = true;
		}
		return lookup;
	}
};
