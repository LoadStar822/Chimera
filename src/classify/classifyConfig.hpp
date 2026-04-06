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
	// Unified tail-risk control in [0,1], where 0=head-heavy and 1=tail-rich.
	double tail_risk_u = 1.0;
	double tail_risk_s = 1.0;
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

	enum class ReadRegime {
		ShortLike,
		LongLike,
	};

	struct CandidatePolicy {
		bool enable_taxpool{ true };
	};

	struct PostDecisionPolicy {
		bool disable_fallback{ false };
		bool enable_selective_reject{ false };
	};

	struct AutoClassifyPolicy {
		ReadRegime regime{ ReadRegime::LongLike };
		CandidatePolicy candidate{};
		PostDecisionPolicy post{};
	};

	inline AutoClassifyPolicy derive_auto_policy(const FileInfo &fi,
	                                             const ClassifyConfig &cfg) {
		(void)cfg;
		AutoClassifyPolicy policy{};
		const bool short_like = fi.avgLen < 1000;
		policy.regime = short_like ? ReadRegime::ShortLike : ReadRegime::LongLike;
		policy.candidate.enable_taxpool = true;
		policy.post.disable_fallback = short_like;
		policy.post.enable_selective_reject = !short_like;
		return policy;
	}

	inline CandidatePolicy derive_candidate_policy(const ClassifyConfig &cfg) {
		(void)cfg;
		CandidatePolicy policy{};
		policy.enable_taxpool = true;
		return policy;
	}

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
			std::string reject_reason; // 为空表示未拒绝或接受
			std::string best_taxid_hint; // 最佳候选 taxid（即使未被接受）
		};

			struct DecisionConfig {
					double min_class_weight = 1e-4;
					double fallback_strength = 0.0; // 连续 fallback 强度（0=禁用，1=最积极）
					double selective_reject_strength = 0.0; // 连续 selective reject 强度
					double fallback_gap_min = 0.10; // fallback 仅在 top1-top2 gap 足够大时触发（抑制 FP）
				};
		}
	#endif // !CLASSIFYCONFIG_HPP
