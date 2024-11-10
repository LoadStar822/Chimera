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
 * Last Modified: 2024-09-20
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
#include <string>
#include <iomanip>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <cstdint>
#include <atomic>

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
			<< std::setw(20) << "Verbose:" << config.verbose << std::endl;

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
		std::unordered_map<std::string, size_t> taxidTotalMatches;
		std::unordered_map<std::string, size_t> taxidUniqueMatches;
	};

	struct batchReads {
		std::vector< std::string >                 ids;
		std::vector< std::vector< seqan3::dna4 > > seqs;
		std::vector< std::vector< seqan3::dna4 > > seqs2{};
	};

	struct classifyResult {
		std::string id;
		std::vector<std::pair<std::string, size_t>> taxidCount;
	};
}
#endif // !CLASSIFYCONFIG_HPP