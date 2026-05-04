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
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace ChimeraClassify {
		struct ClassifyConfig {
			std::vector<std::string> singleFiles;
			std::vector<std::string> pairedFiles;
			std::string outputFile;
			std::string dbFile;
		double shotThreshold = 0.70;
		double firstFilterBeta = 0.8;
		bool firstFilterBeta_user = false; // set when user or auto-override explicitly chooses beta
			double presence_pi = 1e-3;
			double presence_tau = 4.6;
			uint32_t presence_breadth_bits = 2048; // breadth sketch bits (power of 2 suggested)
			uint16_t threads;
			size_t batchSize;
		size_t emIter;
		double em_prune_ratio = 2e-4;   // relative to max_expected in EM sparsity
		double em_conf_power = 2.0;     // confidence exponent for EM M-step (0 disables)
			double post_pi_min = 5e-4;
	// Continuous sample-level community dispersion in [0,1].
	double community_dispersion_u = 1.0;
	double community_dispersion_s = 1.0;
		size_t hash_sample_min = 16;
		size_t hash_sample_max = 96;
			double idf_max = 8.0;
			};

	struct FileInfo {
		size_t fileNum{ 0 };
		size_t sequenceNum{ 0 };
		size_t classifiedNum{ 0 };
		size_t unclassifiedNum{ 0 };
		size_t minLen = 0;
		size_t maxLen = 0;
		size_t avgLen = 0;
		size_t bpLength = 0;
		std::unordered_set<std::string> uniqueTaxids;
	};

	struct batchReads {
		std::vector< std::string >                 ids;
		std::vector< std::vector< seqan3::dna4 > > seqs;
		std::vector< std::vector< seqan3::dna4 > > seqs2{};
	};

			struct classifyResult {
				std::string id;
				std::vector<std::pair<std::string, double>> taxidCount;
				std::vector<std::pair<std::string, double>> posteriors;
				double evaluated{ 0.0 }; // 实际参与判别的 feature 数，用于归一化
				std::string reject_reason; // 为空表示未拒绝或接受
				std::string best_taxid_hint; // 最佳候选 taxid（即使未被接受）
				std::vector<std::pair<std::string, double>> abundanceCount;
				std::vector<std::pair<std::string, double>> sampleMixturePosteriors;
				std::vector<std::pair<std::string, double>> sampleMixtureLocalScores;
				double sampleMixtureTopScore{ 0.0 };
			};

			struct DecisionConfig {
					double min_class_weight = 1e-4;
				};
		}
	#endif // !CLASSIFYCONFIG_HPP
