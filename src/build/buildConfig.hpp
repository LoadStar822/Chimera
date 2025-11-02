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
#include <string>
#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>

namespace ChimeraBuild {
	inline constexpr uint64_t adjust_seed(uint8_t const kmer_size,
		uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
	{
		unsigned double_k = static_cast<unsigned>(kmer_size) * 2u;
		unsigned shift = double_k >= 64u ? 0u : (64u - double_k);
		return seed >> shift;
	}

	struct BuildConfig {
		std::string taxonomy_kind{ "auto" };
		std::string taxonomy_version{ "auto" };
		std::string input_file;
		std::string output_file;
		std::string filter{ "imcf" };
		std::string feature{ "strobemer" }; // syncmer | strobemer | auto
		uint8_t strobemer_k{ 28 };
		uint8_t strobemer_order{ 2 };
		uint16_t strobemer_w_min{ 12 };
		uint16_t strobemer_w_max{ 32 };
		uint8_t kmer_size{ 31 };
		uint16_t smer_size{ 16 };
		uint16_t syncmer_position{ 7 };
		uint64_t min_length{ 0 };
		uint16_t threads;
		bool verbose = true;
		double load_factor{ 0.95 };
		size_t max_hashes_per_taxid = 0;
		bool adaptive_cutoff = false;
	};

	inline std::ostream& operator<<(std::ostream& os, const BuildConfig& config) {
		os << std::string(50, '=') << std::endl;
		os << " Build Configuration " << std::endl;
		os << std::string(50, '=') << std::endl;

		os << std::left
			<< std::setw(25) << "Input file:" << config.input_file << std::endl
			<< std::setw(25) << "Output file:" << config.output_file << std::endl
			<< std::setw(25) << "Filter:" << config.filter << std::endl
			<< std::setw(25) << "Feature method:" << config.feature << std::endl
			<< std::setw(25) << "Strobemer k:" << static_cast<int>(config.strobemer_k) << std::endl
			<< std::setw(25) << "Strobemer order:" << static_cast<int>(config.strobemer_order) << std::endl
			<< std::setw(25) << "Strobemer w_min:" << config.strobemer_w_min << std::endl
			<< std::setw(25) << "Strobemer w_max:" << config.strobemer_w_max << std::endl
			<< std::setw(25) << "Taxonomy kind:" << config.taxonomy_kind << std::endl
			<< std::setw(25) << "Taxonomy version:" << config.taxonomy_version << std::endl
			<< std::setw(25) << "Kmer size:" << (int)config.kmer_size << std::endl
			<< std::setw(25) << "Syncmer s-mer size:" << config.smer_size << std::endl
			<< std::setw(25) << "Syncmer offset:" << config.syncmer_position << std::endl
			<< std::setw(25) << "Minimum length:" << config.min_length << std::endl
			<< std::setw(25) << "Threads:" << config.threads << std::endl
			<< std::setw(25) << "Load factor:" << config.load_factor << std::endl
			<< std::setw(25) << "Max hashes per taxid:" << config.max_hashes_per_taxid << std::endl
			<< std::setw(25) << "Adaptive cutoff:" << std::boolalpha << config.adaptive_cutoff << std::noboolalpha << std::endl
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

	struct IMCFConfig {
		inline static constexpr uint64_t DefaultFingerprintSalt = 0xD1B54A32D192ED03ULL;
		inline static constexpr uint8_t CurrentHashVersion = 2;

		size_t binNum{};
		size_t binSize{};
		uint8_t kmerSize{};
		uint16_t smerSize{};
		uint16_t syncmerPosition{ 0 };
		int MaxCuckooCount{ 500 };
		double loadFactor{ 0.95 };
		uint64_t seed64{ 0 };
		uint64_t fpSalt{ DefaultFingerprintSalt };
		uint8_t hashVersion{ 0 };
		std::string taxonomyKind{ "ncbi" };
		std::string taxonomyVersion{ "ncbi-taxdump" };
		uint8_t featureMethod{ 0 }; // 0=syncmer, 1=strobemer
	uint8_t strobeOrder{ 0 };
	uint16_t strobeWmin{ 0 };
	uint16_t strobeWmax{ 0 };
	uint8_t strobeK{ 0 };

	template <class Archive>
	void serialize(Archive& archive) {
		archive(binNum,
		       binSize,
		       MaxCuckooCount,
		       loadFactor,
		       kmerSize,
		       smerSize,
		       syncmerPosition,
		       seed64,
		       fpSalt,
		       hashVersion,
		       taxonomyKind,
		       taxonomyVersion,
		       featureMethod,
		       strobeOrder,
		       strobeWmin,
		       strobeWmax,
		       strobeK);
	}
	};
}

#endif // BUILDCONFIG_HPP
