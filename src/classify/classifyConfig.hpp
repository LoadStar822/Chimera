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
#include <array>
#include <unordered_set>
#include <string>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdint>
#include <atomic>
#include <memory>
#include "robin_hood.h"
#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace ChimeraClassify {
		struct ClassifyConfig {
			std::vector<std::string> singleFiles;
			std::vector<std::string> pairedFiles;
			std::string weight_map_file; // optional contig/read weight map (CAMI mapping.tsv or id<TAB>weight)
			std::string outputFile;
			std::string dbFile;
		std::string feature{ "auto" };
		uint8_t strobemer_k{ 0 };      // inherit from DB unless >0
		uint8_t strobemer_order{ 0 };
		uint16_t strobemer_w_min{ 0 };
		uint16_t strobemer_w_max{ 0 };
		std::string taxonomyKind{ "auto" };
		std::string taxonomyVersion{ "auto" };
		double shotThreshold = 0.70;
		double firstFilterBeta = 0.8;
		bool firstFilterBeta_user = false; // set when user or auto-override explicitly chooses beta
		// High-div / EM only: fixed-budget keepalive (replace tail, no pad expansion).
		// Protect a strong per-read hint candidate from being pruned out of the
		// pre-EM topK, so the correct branch can enter EM/posterior lists.
		double preem_keepalive_min_ratio = 0.20;      // require hint_score >= top1_score * ratio
			double preem_keepalive_replace_ratio = 1.00;  // require hint_score >= tail_score * ratio
		double preem_keepalive_abs_min = 2.0;         // absolute hint score floor
		// High-div / EM only: experimental beta gate relaxation (disabled by default).
		bool preem_beta_relax = false;
		// NCBI-only experimental knobs for strain/assembly saturation:
		// - collapse_strain_hits: collapse per-hash hit lists to 1 representative taxid per species
		//   (affects scoring/deg on the hot path; can change behavior).
		// - collapse_strain_candidates: collapse the pre-EM candidate list to 1 representative per species
		//   before topK truncation (lighter; mainly affects candidate diversity).
	bool collapse_strain_hits = true;
		bool collapse_strain_candidates = true;
		// NCBI-only: compute deg/exclusivity using species groups instead of raw taxid
		// multiplicity. This helps avoid "strain saturation" (many strain taxids) from
		// over-penalizing within-genus shared evidence, while keeping output taxids
		// unchanged (still DB taxids).
	bool deg_by_species = false;
			double presence_pi = 1e-3;
			double presence_tau = 4.6;
			uint32_t presence_breadth_bits = 2048; // breadth sketch bits (power of 2 suggested)
			uint16_t threads;
			bool verbose = true;
			size_t batchSize;
	bool em = true;
		double emThreshold;
		size_t emIter;
	double em_prune_ratio = 2e-4;   // relative to max_expected in EM sparsity
	double em_prior_strength = 1.0; // Dirichlet mass; 0 uses alpha only
	double em_conf_power = 2.0;     // confidence exponent for EM M-step (0 disables)
	double post_thres = 0.56;
	double post_pi_min = 5e-4;
	// postEmDecision posterior tail controls (affect taxidCount sparsity; not POST_TOPK dump)
	double post_min_fraction = 0.01; // per-read posterior fraction cutoff (0 disables)
	uint32_t dump_post_topk = 256; // 输出 POST_TOPK=...（用于 profile 侧属内纠错）
	bool low_div_auto = true;       // internal: auto low-diversity branch
	bool low_div_active = false;    // internal: set when low-div branch is enabled
	uint32_t low_div_probe_reads = 200000; // internal: probe reads for low-div detect
	size_t max_reads = 0;           // internal: stop after N reads (0 disables)
	size_t hash_sample_min = 16;
	size_t hash_sample_max = 96;
		double idf_max = 8.0;
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
			os << std::setw(20) << "Weight map:"
				<< (config.weight_map_file.empty() ? "none" : config.weight_map_file) << std::endl;
			os << std::setw(20) << "Output file:" << config.outputFile << std::endl
				<< std::setw(20) << "Database file:" << config.dbFile << std::endl
				<< std::setw(20) << "Feature method:" << config.feature << std::endl
			<< std::setw(20) << "Strobemer k:" << static_cast<int>(config.strobemer_k) << std::endl
			<< std::setw(20) << "Strobemer order:" << static_cast<int>(config.strobemer_order) << std::endl
			<< std::setw(20) << "Strobemer w_min:" << config.strobemer_w_min << std::endl
			<< std::setw(20) << "Strobemer w_max:" << config.strobemer_w_max << std::endl
			<< std::setw(20) << "Taxonomy kind:" << config.taxonomyKind << std::endl
			<< std::setw(20) << "Taxonomy version:" << config.taxonomyVersion << std::endl
			<< std::setw(20) << "Shot threshold:" << config.shotThreshold << std::endl
			<< std::setw(20) << "First filter beta:" << config.firstFilterBeta << std::endl
			<< std::setw(20) << "Pre-EM keepalive min ratio:" << config.preem_keepalive_min_ratio << std::endl
			<< std::setw(20) << "Pre-EM keepalive repl ratio:" << config.preem_keepalive_replace_ratio << std::endl
			<< std::setw(20) << "Pre-EM keepalive abs min:" << config.preem_keepalive_abs_min << std::endl
			<< std::setw(20) << "Pre-EM beta relax:" << config.preem_beta_relax << std::endl
			<< std::setw(20) << "Collapse hits:" << config.collapse_strain_hits << std::endl
			<< std::setw(20) << "Collapse cands:" << config.collapse_strain_candidates << std::endl
			<< std::setw(20) << "Deg by species:" << config.deg_by_species << std::endl
			<< std::setw(20) << "Presence pi:" << config.presence_pi << std::endl
				<< std::setw(20) << "Presence tau:" << config.presence_tau << std::endl
				<< std::setw(20) << "Breadth bits:" << config.presence_breadth_bits << std::endl
			<< std::setw(20) << "Batch size:" << config.batchSize << std::endl
			<< std::setw(20) << "EM:" << config.em << std::endl
			<< std::setw(20) << "Threads:" << config.threads << std::endl
			<< std::setw(20) << "Verbose:" << config.verbose << std::endl
			<< std::setw(20) << "EM prune ratio:" << config.em_prune_ratio << std::endl
			<< std::setw(20) << "EM prior strength:" << config.em_prior_strength << std::endl
			<< std::setw(20) << "EM conf power:" << config.em_conf_power << std::endl
			<< std::setw(20) << "Hash sample min:" << config.hash_sample_min << std::endl
			<< std::setw(20) << "Hash sample max:" << config.hash_sample_max << std::endl
			<< std::setw(20) << "IDF max:" << config.idf_max << std::endl
			<< std::setw(20) << "Posterior thres:" << config.post_thres << std::endl
			<< std::setw(20) << "Posterior pi min:" << config.post_pi_min << std::endl
			<< std::setw(20) << "Post min frac:" << config.post_min_fraction << std::endl
			<< std::setw(20) << "Dump POST_TOPK:" << config.dump_post_topk << std::endl;
		os << std::string(40, '=') << std::endl;

		return os;
	}

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
		robin_hood::unordered_flat_map<std::string, size_t> taxidTotalMatches;
		robin_hood::unordered_flat_map<std::string, size_t> taxidUniqueMatches;
		robin_hood::unordered_flat_map<std::string, size_t> rejectReasons;
		robin_hood::unordered_flat_map<std::string,
			robin_hood::unordered_flat_map<std::string, size_t>> rejectByTaxid;
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
			bool presence_passed{ false }; // 是否通过 presence 层
		};

			struct DecisionConfig {
				double posterior_threshold = 0.56;
				double min_class_weight = 1e-4;
				double posterior_min_fraction = 0.01; // 软分配时的最小 posterior 占比阈值
				double posterior_head_mass = 0.95;   // 每条 read 只保留 posterior^alpha 的头部质量（抑制长尾）
				uint32_t posterior_max_taxa = 8;     // 每条 read 最多输出的 taxon 数（抑制长尾）
				bool allow_fallback_on_reject = false; // 高多样性：放宽拒绝，减少 unclassified
				double fallback_gap_min = 0.10; // fallback 仅在 top1-top2 gap 足够大时触发（抑制 FP）
			};
		}
	#endif // !CLASSIFYCONFIG_HPP
