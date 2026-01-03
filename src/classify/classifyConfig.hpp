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
#include <seqan3/core/debug_stream.hpp>

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
		bool adaptive_shot = true;
		double firstFilterBeta = 0.8;
		bool firstFilterBeta_user = false; // set when user or auto-override explicitly chooses beta
		size_t preEmTopK = 16;
		// High-div / EM only: when pre-EM candidate list collapses to a singleton
		// species (after strain collapse), keep a tiny floor of additional species
		// candidates alive to avoid over-pruning sister species. This is NOT a
		// pad-to-32 expansion.
		uint32_t preem_floor_target = 4;           // target #non-unclassified species
		double preem_floor_min_ratio = 0.20;       // trigger gate: top2/top1 >= ratio
		double preem_floor_add_min_ratio = 0.10;   // add gate: cand_score >= top1*ratio
		// High-div / EM only: fixed-budget keepalive (replace tail, no pad expansion).
		// Protect a strong per-read hint candidate from being pruned out of the
		// pre-EM topK, so the correct branch can enter EM/posterior lists.
		double preem_keepalive_min_ratio = 0.20;      // require hint_score >= top1_score * ratio
			double preem_keepalive_replace_ratio = 1.00;  // require hint_score >= tail_score * ratio
		double preem_keepalive_abs_min = 2.0;         // absolute hint score floor
		// High-div / EM only: underfull fixed-budget fill (no pad expansion).
		// Disabled by default (can be enabled for experiments).
		bool preem_underfull_fill = false;
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
		std::string decoy_mode{ "imcf-edge-shuffle" };
		uint32_t decoy_reps = 3;
		double exclusive_gamma = 1.2;
			double presence_pi = 1e-3;
			double presence_tau = 4.6;
			double presence_noise = 0.0; // <=0 表示自动估计
			uint32_t presence_u_min = 1;
			bool presence_pre_filter = false; // allow pre-EM presence filter to hard-prune candidates
			uint32_t presence_breadth_bits = 2048; // breadth sketch bits (power of 2 suggested)
			double presence_breadth_min_ratio = 0.0; // minimum breadth ratio (0 disables gate)
			uint32_t presence_breadth_min_obs = 0; // minimum observed uniques (0 disables gate)
			double presence_breadth_penalty = 0.0; // subtract from logPosterior if breadth low (0 disables)
			uint16_t threads;
			bool verbose = true;
			size_t batchSize;
	bool em = true;
		double emThreshold;
		size_t emIter;
	double em_prune_ratio = 2e-4;   // relative to max_expected in EM sparsity
	double em_prior_strength = 1.0; // Dirichlet mass; 0 uses alpha only
		double em_coexist_penalty = 0.6; // penalty for near-tied taxa in EM softmax
	double em_conf_power = 2.0;     // confidence exponent for EM M-step (0 disables)
	double post_thres = 0.56;
	double post_pi_min = 5e-4;
	// postEmDecision posterior tail controls (affect taxidCount sparsity; not POST_TOPK dump)
	double post_min_fraction = 0.01; // per-read posterior fraction cutoff (0 disables)
	double post_power = 1.5;         // posterior^alpha sharpening (>=1 recommended)
	double post_head_mass = 0.0;     // 0 => auto (based on pi tuning); else (0,1]
	uint32_t post_max_taxa = 0;      // 0 => auto (based on pi tuning); else >=1
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
			<< std::setw(20) << "Adaptive shot:" << config.adaptive_shot << std::endl
			<< std::setw(20) << "First filter beta:" << config.firstFilterBeta << std::endl
			<< std::setw(20) << "Pre-EM topK:" << config.preEmTopK << std::endl
			<< std::setw(20) << "Pre-EM floor:" << config.preem_floor_target << std::endl
			<< std::setw(20) << "Pre-EM floor min ratio:" << config.preem_floor_min_ratio << std::endl
			<< std::setw(20) << "Pre-EM floor add ratio:" << config.preem_floor_add_min_ratio << std::endl
			<< std::setw(20) << "Pre-EM keepalive min ratio:" << config.preem_keepalive_min_ratio << std::endl
			<< std::setw(20) << "Pre-EM keepalive repl ratio:" << config.preem_keepalive_replace_ratio << std::endl
			<< std::setw(20) << "Pre-EM keepalive abs min:" << config.preem_keepalive_abs_min << std::endl
			<< std::setw(20) << "Pre-EM underfull fill:" << config.preem_underfull_fill << std::endl
			<< std::setw(20) << "Pre-EM beta relax:" << config.preem_beta_relax << std::endl
			<< std::setw(20) << "Collapse hits:" << config.collapse_strain_hits << std::endl
			<< std::setw(20) << "Collapse cands:" << config.collapse_strain_candidates << std::endl
			<< std::setw(20) << "Deg by species:" << config.deg_by_species << std::endl
			<< std::setw(20) << "Presence pi:" << config.presence_pi << std::endl
				<< std::setw(20) << "Presence tau:" << config.presence_tau << std::endl
				<< std::setw(20) << "Presence noise:" << config.presence_noise << std::endl
				<< std::setw(20) << "Presence U min:" << config.presence_u_min << std::endl
				<< std::setw(20) << "Presence pre:" << config.presence_pre_filter << std::endl
				<< std::setw(20) << "Breadth bits:" << config.presence_breadth_bits << std::endl
				<< std::setw(20) << "Breadth min r:" << config.presence_breadth_min_ratio << std::endl
				<< std::setw(20) << "Breadth min u:" << config.presence_breadth_min_obs << std::endl
				<< std::setw(20) << "Breadth penalty:" << config.presence_breadth_penalty << std::endl
				<< std::setw(20) << "Decoy mode:" << config.decoy_mode << std::endl
				<< std::setw(20) << "Decoy reps:" << config.decoy_reps << std::endl
			<< std::setw(20) << "Exclusive gamma:" << config.exclusive_gamma << std::endl
			<< std::setw(20) << "Batch size:" << config.batchSize << std::endl
			<< std::setw(20) << "EM:" << config.em << std::endl
			<< std::setw(20) << "Threads:" << config.threads << std::endl
			<< std::setw(20) << "Verbose:" << config.verbose << std::endl
			<< std::setw(20) << "EM prune ratio:" << config.em_prune_ratio << std::endl
			<< std::setw(20) << "EM prior strength:" << config.em_prior_strength << std::endl
			<< std::setw(20) << "EM coexist penalty:" << config.em_coexist_penalty << std::endl
			<< std::setw(20) << "EM conf power:" << config.em_conf_power << std::endl
			<< std::setw(20) << "Hash sample min:" << config.hash_sample_min << std::endl
			<< std::setw(20) << "Hash sample max:" << config.hash_sample_max << std::endl
			<< std::setw(20) << "IDF max:" << config.idf_max << std::endl
			<< std::setw(20) << "Posterior thres:" << config.post_thres << std::endl
			<< std::setw(20) << "Posterior pi min:" << config.post_pi_min << std::endl
			<< std::setw(20) << "Post min frac:" << config.post_min_fraction << std::endl
			<< std::setw(20) << "Post power:" << config.post_power << std::endl
			<< std::setw(20) << "Post head mass:" << (config.post_head_mass > 0.0 ? std::to_string(config.post_head_mass) : std::string("auto")) << std::endl
			<< std::setw(20) << "Post max taxa:" << (config.post_max_taxa > 0 ? std::to_string(config.post_max_taxa) : std::string("auto")) << std::endl
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
		size_t preem_route_reads{ 0 };
		size_t preem_route_rare_reads{ 0 };
		size_t preem_cap_checks{ 0 };
		size_t preem_cap_expanded{ 0 };
		size_t preem_floor_checks{ 0 };
		size_t preem_floor_applied{ 0 };
		size_t preem_floor_added{ 0 };
		size_t preem_floor_skipped_dominant{ 0 };
		size_t preem_floor_filtered_weak{ 0 };
		size_t preem_keepalive_attempt{ 0 };
		size_t preem_keepalive_applied{ 0 };
		size_t preem_keepalive_blocked_low_ratio{ 0 };
		size_t preem_keepalive_blocked_low_gain{ 0 };
		size_t preem_keepalive_blocked_low_abs{ 0 };
		// Pre-EM beta relax / debugging counters (high-div only).
		size_t preem_beta_relax_seen{ 0 };
		size_t preem_beta_relax_dom_beta{ 0 };
		size_t preem_beta_relax_strict_le_halfk{ 0 };
		size_t preem_beta_relax_eff_lt_min{ 0 };
		size_t preem_beta_relax_beta_user{ 0 };
		size_t preem_beta_relax_checks{ 0 };
		size_t preem_beta_relax_applied{ 0 };
		size_t preem_beta_relax_thr_drop_sum{ 0 };
		size_t preem_beta_relax_suppressed_overflow{ 0 };
		size_t preem_dynamic_topk_96{ 0 };
		// High-div / EM only: fixed-budget underfull fill (no pad expansion).
		size_t preem_underfull_fill_checks{ 0 };
		size_t preem_underfull_fill_applied{ 0 };
		size_t preem_underfull_fill_added{ 0 };
		size_t preem_underfull_fill_stage2_added{ 0 };
		// High-div / EM only: nearTie/binOverflow breakdown and finalK histogram.
		size_t preem_neartie_gap{ 0 };
		size_t preem_neartie_ratio{ 0 };
		size_t preem_overflow_topbins{ 0 };
		size_t preem_overflow_size{ 0 };
		size_t preem_finalk_le16{ 0 };
		size_t preem_finalk_17_32{ 0 };
		size_t preem_finalk_33_64{ 0 };
		size_t preem_finalk_65_96{ 0 };
		size_t preem_finalk_eq96{ 0 };
		// Hit-level IDF clamp audit (minimizer scoring).
			size_t hit_idf_total{ 0 };
			size_t hit_idf_raw_lt0p5{ 0 };
			std::array<uint64_t, 4> hit_idf_raw_bins{};
			double hit_idf_contrib_sum_old{ 0.0 };
				double hit_idf_contrib_sum_new{ 0.0 };
					std::array<uint64_t, 64> hit_idf_raw_hist{};
					std::array<uint64_t, 64> hit_idf_eff_hist{};
					std::array<uint64_t, 64> hit_idf_power_hist{};
					// TF saturation audit (shared minimizer de-dup / pile-up suppression).
					size_t tf_sat_enabled_reads{ 0 };
					size_t tf_sat_shared_hits{ 0 };
					size_t tf_sat_damped_hits{ 0 };
					double tf_sat_base_sum{ 0.0 };
					double tf_sat_drop_sum{ 0.0 };
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
				double posterior_power = 1.5;        // posterior^alpha 压尖，>1 越尖锐
				double posterior_head_mass = 0.95;   // 每条 read 只保留 posterior^alpha 的头部质量（抑制长尾）
				uint32_t posterior_max_taxa = 8;     // 每条 read 最多输出的 taxon 数（抑制长尾）
				bool allow_fallback_on_reject = false; // 高多样性：放宽拒绝，减少 unclassified
				double fallback_gap_min = 0.10; // fallback 仅在 top1-top2 gap 足够大时触发（抑制 FP）
			};
		}
	#endif // !CLASSIFYCONFIG_HPP
