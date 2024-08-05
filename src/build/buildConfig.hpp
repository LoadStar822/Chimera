#pragma once
#ifndef BUILDCONFIG_HPP
#define BUILDCONFIG_HPP
#include <iostream>
#include <iomanip>
namespace ChimeraBuild {
	struct BuildConfig {
		std::string input_file;
		std::string output_file;
		std::string mode;
		uint8_t kmer_size;
		uint16_t window_size;
		uint64_t min_length;
		uint16_t threads;
		bool verbose;
		double filter_size;
	};

	inline std::ostream& operator<<(std::ostream& os, const BuildConfig& config) {
		os << std::string(40, '=') << std::endl;
		os << " Build Configuration " << std::endl;
		os << std::string(40, '=') << std::endl;

		os << std::left
			<< std::setw(20) << "Input file:" << config.input_file << std::endl
			<< std::setw(20) << "Output file:" << config.output_file << std::endl
			<< std::setw(20) << "Mode:" << config.mode << std::endl
			<< std::setw(20) << "Kmer size:" << (int)config.kmer_size << std::endl
			<< std::setw(20) << "Window size:" << config.window_size << std::endl
			<< std::setw(20) << "Minimum length:" << config.min_length << std::endl
			<< std::setw(20) << "Threads:" << config.threads << std::endl
			<< std::setw(20) << "Filter size:" << config.filter_size << std::endl;

		os << std::string(40, '=') << std::endl;

		return os;
	}

	struct FileInfo {
		size_t fileNum = 0;
		size_t invalidNum = 0;
		size_t sequenceNum = 0;
		size_t skippedNum = 0;
	};

	struct ICFConfig {
		uint8_t kmer_size;
		uint16_t window_size;
		size_t bins;
	};
}

#endif // BUILDCONFIG_HPP