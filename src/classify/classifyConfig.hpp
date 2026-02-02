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
			std::string weight_map_file; // optional contig/read weight map (CAMI mapping.tsv or id<TAB>weight)
			std::string outputFile;
			std::string dbFile;
		std::string taxonomyKind{ "auto" };
		std::string taxonomyVersion{ "auto" };
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
	// postEmDecision posterior tail controls (affect taxidCount sparsity; not POST_TOPK dump)
	bool low_div_active = false;    // internal: set when low-div branch is enabled
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
			double evaluated{ 0.0 }; // 实际参与判别的 syncmer 数，用于归一化
			double sample_weight{ 0.0 }; // 可选样本权重（如 CAMI number_reads）；0 表示未提供
			std::string reject_reason; // 为空表示未拒绝或接受
			std::string best_taxid_hint; // 最佳候选 taxid（即使未被接受）
		};

			struct DecisionConfig {
				double min_class_weight = 1e-4;
				bool allow_fallback_on_reject = false; // 高多样性：放宽拒绝，减少 unclassified
				double fallback_gap_min = 0.10; // fallback 仅在 top1-top2 gap 足够大时触发（抑制 FP）
			};
		}
	#endif // !CLASSIFYCONFIG_HPP
