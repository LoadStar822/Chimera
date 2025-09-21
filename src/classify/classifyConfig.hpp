/*
 * -----------------------------------------------------------------------------
 * Filename:      classifyConfig.hpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-08-10
 *
 * Last Modified: 2024-11-18
 *
 * Description:
 *  Classify configuration for Chimera
 *
 * Version:
 *  1.2
 * -----------------------------------------------------------------------------
 */
#pragma once
#ifndef CLASSIFYCONFIG_HPP
#define CLASSIFYCONFIG_HPP
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdint>
#include <atomic>
#include "robin_hood.h"
#include <optional>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace ChimeraClassify {
	struct ClassifyConfig {
		std::vector<std::string> singleFiles;
		std::vector<std::string> pairedFiles;
		std::string outputFile;
		std::string dbFile;
		std::string taxFile;
		std::string filter{ "hicf" };
		double shotThreshold;
		uint16_t threads;
		std::string mode;
		bool verbose = true;
		size_t batchSize;
		bool lca = false;
		bool em = false;
		bool vem = false;
		double emThreshold;
		size_t emIter;
		double post_thres = 0.9;
		double post_margin = 0.2;
		double post_ratio = std::numeric_limits<double>::quiet_NaN();
		double post_pi_min = 1e-4;
		bool lca_fallback = false;
		bool output_posterior = true;
		bool skip_post_filter = true;
	};

	inline std::ostream& operator<<(std::ostream& os, const ClassifyConfig& config) {
		os << std::string(40, '=') << std::endl;
		os << " Classify Configuration " << std::endl;
		os << std::string(40, '=') << std::endl;

		os << std::left
			<< std::setw(20) << "Single files:" << std::endl;
		for (const auto& file : config.singleFiles) {
			os << std::setw(20) << "" << file << std::endl;
		}
		os << std::setw(20) << "Paired files:" << std::endl;
		for (const auto& file : config.pairedFiles) {
			os << std::setw(20) << "" << file << std::endl;
		}
		os << std::setw(20) << "Output file:" << config.outputFile << std::endl
			<< std::setw(20) << "Database file:" << config.dbFile << std::endl
			<< std::setw(20) << "Shot threshold:" << config.shotThreshold << std::endl
			<< std::setw(20) << "Mode:" << config.mode << std::endl
			<< std::setw(20) << "Filter:" << config.filter << std::endl
			<< std::setw(20) << "Batch size:" << config.batchSize << std::endl
			<< std::setw(20) << "LCA:" << config.lca << std::endl
			<< std::setw(20) << "EM:" << config.em << std::endl
			<< std::setw(20) << "VEM:" << config.vem << std::endl
			<< std::setw(20) << "Threads:" << config.threads << std::endl
			<< std::setw(20) << "Verbose:" << config.verbose << std::endl
			<< std::setw(20) << "Posterior thres:" << config.post_thres << std::endl
			<< std::setw(20) << "Posterior margin:" << config.post_margin << std::endl
			<< std::setw(20) << "Posterior ratio:" << (std::isnan(config.post_ratio) ? std::string("nan") : std::to_string(config.post_ratio)) << std::endl
			<< std::setw(20) << "Posterior pi min:" << config.post_pi_min << std::endl
			<< std::setw(20) << "LCA fallback:" << config.lca_fallback << std::endl
			<< std::setw(20) << "Output posterior:" << config.output_posterior << std::endl
			<< std::setw(20) << "Skip post filter:" << config.skip_post_filter << std::endl;

		os << std::string(40, '=') << std::endl;

		return os;
	}

	struct FileInfo {
		size_t fileNum{ 0 };
		size_t sequenceNum{ 0 };
		size_t classifiedNum{ 0 };
		size_t unclassifiedNum{ 0 };
		size_t lcaNum = 0;
		size_t minLen = 0;
		size_t maxLen = 0;
		size_t avgLen = 0;
		size_t bpLength = 0;
		std::unordered_set<std::string> uniqueTaxids;
		robin_hood::unordered_flat_map<std::string, size_t> taxidTotalMatches;
		robin_hood::unordered_flat_map<std::string, size_t> taxidUniqueMatches;
	};

	struct batchReads {
		std::vector< std::string >                 ids;
		std::vector< std::vector< seqan3::dna4 > > seqs;
		std::vector< std::vector< seqan3::dna4 > > seqs2{};
	};

	struct classifyResult {
		std::string id;
		std::vector<std::pair<std::string, size_t>> taxidCount;
		std::vector<std::pair<std::string, double>> posteriors;
	};

	struct DecisionConfig {
		double posterior_threshold = 0.9;
		double margin_delta = 0.2;
		double margin_ratio = std::numeric_limits<double>::quiet_NaN();
		double min_class_weight = 1e-4;
		bool use_lca_fallback = false;
	};
}
#endif // !CLASSIFYCONFIG_HPP
