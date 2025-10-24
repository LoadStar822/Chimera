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
		std::string input_file;
		std::string output_file;
		std::string filter{ "imcf" };
		std::string taxonomy_kind{ "auto" };
		std::string taxonomy_version{ "auto" };
		uint8_t kmer_size{ 31 };
		uint16_t smer_size{ 16 };
		uint16_t syncmer_position{ 7 };
		uint64_t min_length;
		uint16_t threads;
		bool verbose = true;
		double load_factor{ 0.95 };
		size_t max_hashes_per_taxid = 0;
		bool adaptive_cutoff = false;
		double toxic_quantile{ 0.999 };
		uint64_t toxic_top_n{ 0 };
		double toxic_min_fraction{ 1e-4 };
		double toxic_safety_fraction{ 0.1 };
		uint64_t toxic_safety_min{ 1024 };
		bool enable_strobemers{ true };
		bool strobemer_auto{ true };
		uint16_t strobemer_w_min{ 0 };
		uint16_t strobemer_w_max{ 0 };
		uint16_t strobemer_q{ 0 };
		uint32_t strobemer_max_dist{ 0 };
		uint16_t strobemer_aux_len{ 15 };
		double strobemer_weight{ 0.0 };
		double strobemer_ratio{ 0.6 };
		uint64_t strobemer_seed{ 0x6F37B5E1C943A7DDULL };
	};

	inline std::ostream& operator<<(std::ostream& os, const BuildConfig& config) {
		os << std::string(50, '=') << std::endl;
		os << " Build Configuration " << std::endl;
		os << std::string(50, '=') << std::endl;

		os << std::left
			<< std::setw(25) << "Input file:" << config.input_file << std::endl
			<< std::setw(25) << "Output file:" << config.output_file << std::endl
			<< std::setw(25) << "Filter:" << config.filter << std::endl
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
			<< std::setw(25) << "Toxic quantile:" << config.toxic_quantile << std::endl
			<< std::setw(25) << "Toxic top-N:" << config.toxic_top_n << std::endl
			<< std::setw(25) << "Toxic min fraction:" << config.toxic_min_fraction << std::endl
			<< std::setw(25) << "Toxic safety fraction:" << config.toxic_safety_fraction << std::endl
			<< std::setw(25) << "Toxic safety min:" << config.toxic_safety_min << std::endl
			<< std::setw(25) << "Enable strobemers:" << config.enable_strobemers << std::endl
			<< std::setw(25) << "Strobemer auto:" << config.strobemer_auto << std::endl
			<< std::setw(25) << "Strobemer w_min:" << config.strobemer_w_min << std::endl
			<< std::setw(25) << "Strobemer w_max:" << config.strobemer_w_max << std::endl
			<< std::setw(25) << "Strobemer q:" << config.strobemer_q << std::endl
			<< std::setw(25) << "Strobemer max dist:" << config.strobemer_max_dist << std::endl
			<< std::setw(25) << "Strobemer aux len:" << config.strobemer_aux_len << std::endl
			<< std::setw(25) << "Strobemer weight:" << config.strobemer_weight << std::endl
			<< std::setw(25) << "Strobemer ratio:" << config.strobemer_ratio << std::endl
			<< std::setw(25) << "Strobemer seed:" << config.strobemer_seed << std::endl
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
		inline static constexpr uint8_t CurrentHashVersion = 3;

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
		bool enableStrobemers{ false };
		bool strobemerAuto{ true };
		uint16_t strobemerWMin{ 0 };
		uint16_t strobemerWMax{ 0 };
		uint16_t strobemerQ{ 0 };
		uint32_t strobemerMaxDist{ 0 };
		uint16_t strobemerAuxLen{ 15 };
		double strobemerWeight{ 0.0 };
		double strobemerRatio{ 0.6 };
		uint64_t strobemerSeed{ 0 };

		template <class Archive>
		void save(Archive& archive, const std::uint32_t version) const {
			archive(binNum, binSize, MaxCuckooCount, loadFactor,
				kmerSize, smerSize, syncmerPosition, seed64, fpSalt, hashVersion);
			if (version >= 1) {
				archive(taxonomyKind, taxonomyVersion);
			}
			if (version >= 2) {
				archive(enableStrobemers, strobemerAuto, strobemerWMin, strobemerWMax,
				        strobemerQ, strobemerMaxDist, strobemerAuxLen,
				        strobemerWeight, strobemerRatio, strobemerSeed);
			}
		}

		template <class Archive>
		void load(Archive& archive, const std::uint32_t version) {
			archive(binNum, binSize, MaxCuckooCount, loadFactor,
				kmerSize, smerSize, syncmerPosition, seed64, fpSalt, hashVersion);
			if (version >= 1) {
				archive(taxonomyKind, taxonomyVersion);
			}
			else {
				taxonomyKind = "ncbi";
				taxonomyVersion = "ncbi-taxdump";
			}
			if (version >= 2) {
				archive(enableStrobemers, strobemerAuto, strobemerWMin, strobemerWMax,
				        strobemerQ, strobemerMaxDist, strobemerAuxLen,
				        strobemerWeight, strobemerRatio, strobemerSeed);
			} else {
				enableStrobemers = false;
				strobemerAuto = true;
				strobemerWMin = strobemerWMax = strobemerQ = 0;
				strobemerMaxDist = 0;
				strobemerAuxLen = 15;
				strobemerWeight = 0.0;
				strobemerRatio = 0.6;
				strobemerSeed = 0;
			}
		}
	};
}

CEREAL_CLASS_VERSION(ChimeraBuild::IMCFConfig, 2);

#endif // BUILDCONFIG_HPP
