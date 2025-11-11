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
#include <memory>
#include "robin_hood.h"
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace ChimeraClassify {
namespace dbg {
  struct TraceRecord;
}
	struct ClassifyConfig {
		std::vector<std::string> singleFiles;
		std::vector<std::string> pairedFiles;
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
		size_t preEmTopK = 16;
		bool adaptive_fdr = true;
		double fdr_z = 3.0;
		size_t min_eval_count = 0;
		bool evidence_override = false;
		std::string presence_caller{ "tdfdr" };
		double presence_q = 0.05;
		std::string decoy_mode{ "imcf-edge-shuffle" };
		uint32_t decoy_reps = 3;
		double exclusive_gamma = 1.2;
		uint32_t min_unique_evidence = 5;
		bool presence_report_only = false;
		double presence_abundance_prior = 0.30;
		bool auto_q_tune = true;
		uint16_t threads;
		bool verbose = true;
		bool progress = true;
		size_t progressStep = 5000;
		double progressInterval = 1.0;
		size_t batchSize;
		bool em = false;
		bool vem = false;
		double emThreshold;
		size_t emIter;
		double post_thres = 0.56;
		double post_margin = 0.03;
		double post_ratio = 1.30;
		double post_relax_abs = 0.50;
		double post_relax_ratio = 1.50;
		double post_relax_delta = 0.08;
		double post_relax_delta_abs = 0.52;
		double post_pi_min = 5e-4;
		bool output_posterior = true;
		double em_temp = 1.10;
		double em_prior_strength = 0.25;
		double em_coexist_penalty = 0.20;
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
			<< std::setw(20) << "Presence caller:" << config.presence_caller << std::endl
			<< std::setw(20) << "Presence q:" << config.presence_q << std::endl
			<< std::setw(20) << "Presence report only:" << config.presence_report_only << std::endl
			<< std::setw(20) << "Presence abundance prior:" << config.presence_abundance_prior << std::endl
			<< std::setw(20) << "Auto q tune:" << config.auto_q_tune << std::endl
			<< std::setw(20) << "Decoy mode:" << config.decoy_mode << std::endl
			<< std::setw(20) << "Decoy reps:" << config.decoy_reps << std::endl
			<< std::setw(20) << "Exclusive gamma:" << config.exclusive_gamma << std::endl
			<< std::setw(20) << "Min unique evidence:" << config.min_unique_evidence << std::endl
			<< std::setw(20) << "Adaptive FDR:" << config.adaptive_fdr << std::endl
			<< std::setw(20) << "FDR Z:" << config.fdr_z << std::endl
			<< std::setw(20) << "Min eval count:" << config.min_eval_count << std::endl
			<< std::setw(20) << "Evidence override:" << config.evidence_override << std::endl
			<< std::setw(20) << "Batch size:" << config.batchSize << std::endl
			<< std::setw(20) << "EM:" << config.em << std::endl
			<< std::setw(20) << "VEM:" << config.vem << std::endl
			<< std::setw(20) << "Threads:" << config.threads << std::endl
			<< std::setw(20) << "Verbose:" << config.verbose << std::endl
			<< std::setw(20) << "Progress:" << config.progress << std::endl
			<< std::setw(20) << "Progress step:" << config.progressStep << std::endl
			<< std::setw(20) << "Progress interval:" << config.progressInterval << std::endl
			<< std::setw(20) << "Posterior thres:" << config.post_thres << std::endl
			<< std::setw(20) << "Posterior margin:" << config.post_margin << std::endl
			<< std::setw(20) << "Posterior ratio:" << (std::isnan(config.post_ratio) ? std::string("nan") : std::to_string(config.post_ratio)) << std::endl
			<< std::setw(20) << "Post relax abs:" << config.post_relax_abs << std::endl
			<< std::setw(20) << "Post relax ratio:" << config.post_relax_ratio << std::endl
			<< std::setw(20) << "Post relax delta:" << config.post_relax_delta << std::endl
			<< std::setw(20) << "Post relax delta abs:" << config.post_relax_delta_abs << std::endl
			<< std::setw(20) << "Posterior pi min:" << config.post_pi_min << std::endl
			<< std::setw(20) << "EM temp:" << config.em_temp << std::endl
			<< std::setw(20) << "EM prior strength:" << config.em_prior_strength << std::endl
			<< std::setw(20) << "EM coexist penalty:" << config.em_coexist_penalty << std::endl
			<< std::setw(20) << "Output posterior:" << config.output_posterior << std::endl;

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
		double evaluated{ 0.0 }; // 实际参与判别的 syncmer 数，用于归一化
		std::shared_ptr<dbg::TraceRecord> trace;
	};

	struct DecisionConfig {
		double posterior_threshold = 0.56;
		double margin_delta = 0.03;
		double margin_ratio = 1.30;
		double relax_abs = 0.50;
		double relax_ratio = 1.50;
		double relax_delta = 0.08;
		double relax_delta_abs = 0.52;
		double min_class_weight = 1e-4;
		bool evidence_override = true;
	};
}
#endif // !CLASSIFYCONFIG_HPP
