/*
 * -----------------------------------------------------------------------------
 * Filename:      buildConfig.hpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-07-30
 *
 * Last Modified: 2024-11-18
 *
 * Description:
 *  this file defines the configuration for the build module
 *
 * Version:
 *  1.2
 * -----------------------------------------------------------------------------
 */
#pragma once
#ifndef BUILDCONFIG_HPP
#define BUILDCONFIG_HPP
#include <iostream>
#include <iomanip>
#include <cstdint>

namespace ChimeraBuild {
	struct BuildConfig {
		std::string input_file;
		std::string output_file;
		std::string mode;
		std::string filter{ "hicf" };
		uint8_t kmer_size;
		uint16_t window_size;
		uint64_t min_length;
		uint16_t threads;
		bool verbose = true;
		double load_factor{ 0.95 };
		size_t max_hashes_per_taxid = 0;
		uint8_t fixed_cutoff = 0;
		double alpha{ 1.2 };
		double relaxedLoadFactor{ 0.95 };
	};

	inline std::ostream& operator<<(std::ostream& os, const BuildConfig& config) {
		os << std::string(50, '=') << std::endl;
		os << " Build Configuration " << std::endl;
		os << std::string(50, '=') << std::endl;

		os << std::left
			<< std::setw(25) << "Input file:" << config.input_file << std::endl
			<< std::setw(25) << "Output file:" << config.output_file << std::endl
			<< std::setw(25) << "Filter:" << config.filter << std::endl
			<< std::setw(25) << "Mode:" << config.mode << std::endl
			<< std::setw(25) << "Kmer size:" << (int)config.kmer_size << std::endl
			<< std::setw(25) << "Window size:" << config.window_size << std::endl
			<< std::setw(25) << "Minimum length:" << config.min_length << std::endl
			<< std::setw(25) << "Threads:" << config.threads << std::endl
			<< std::setw(25) << "Load factor:" << config.load_factor << std::endl
			<< std::setw(25) << "Relaxed load factor:" << config.relaxedLoadFactor << std::endl
			<< std::setw(25) << "Alpha:" << config.alpha << std::endl
			<< std::setw(25) << "Max hashes per taxid:" << config.max_hashes_per_taxid << std::endl
			<< std::setw(25) << "Fixed cutoff:" << (int)config.fixed_cutoff << std::endl
			<< std::setw(25) << "Verbose:" << config.verbose << std::endl;

		os << std::string(50, '=') << std::endl;

		return os;
	}

	struct FileInfo {
		size_t fileNum = 0;
		size_t invalidNum = 0;
		size_t sequenceNum = 0;
		size_t skippedNum = 0;
		size_t skippedSeqNum = 0;
		size_t bpLength = 0;

		void operator+=(const FileInfo& other) {
			fileNum += other.fileNum;
			invalidNum += other.invalidNum;
			sequenceNum += other.sequenceNum;
			skippedNum += other.skippedNum;
			bpLength += other.bpLength;
		}
	};

	struct ICFConfig {
		uint8_t kmer_size;
		uint16_t window_size;
		size_t bins;
		size_t bin_size;
		int bitNum = 16;

		template <class Archive>
		void serialize(Archive& archive) {
			archive(kmer_size, window_size, bins, bin_size, bitNum);
		}
	};

	struct HICFConfig {
		size_t userBinsNum;
		size_t technicalBinsNum;
		size_t technicalBinsMaxNum;
		uint8_t tagNum{ 4 };
		uint8_t tagBits{ 16 };
		int MaxCuckooCount{ 500 };
		double loadFactor{ 0.95 };
		double relaxedLoadFactor{ 0.95 };
		uint8_t kmerSize;
		uint16_t windowSize;
		size_t totalHashes;


		template <class Archive>
		void serialize(Archive& archive) {
			archive(userBinsNum, technicalBinsNum, tagNum, tagBits, MaxCuckooCount, loadFactor);
			archive(relaxedLoadFactor, kmerSize, windowSize, technicalBinsMaxNum, totalHashes);
		}
	};

	struct IMCFConfig {
		size_t binNum;
		size_t binSize;
		uint8_t kmerSize;
		uint16_t windowSize;
		int MaxCuckooCount{ 500 };
		double loadFactor{ 0.95 };

		template <class Archive>
		void serialize(Archive& archive) {
			archive(binNum, binSize, MaxCuckooCount, loadFactor);
			archive(kmerSize, windowSize);
		}
	};
}

#endif // BUILDCONFIG_HPP