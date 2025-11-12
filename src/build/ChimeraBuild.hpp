/*
 * -----------------------------------------------------------------------------
 * Filename:      Chimera.hpp
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
 * this file defines the configuration for the build module
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#pragma once
#include <buildConfig.hpp>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>
#include <robin_hood.h>
#include <filesystem>
#include <string_view>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <dna4_traits.hpp>
#include <thread>
#include <vector>
#include <mutex>
#include <cstdint>
#include <interleaved-merged-cuckoo-filter.h>

namespace fs = std::filesystem;

namespace ChimeraBuild {
	std::vector<std::vector<std::string>> buildIMCF(
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		const std::vector<chimera::imcf::Group>& groups,
		const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount,
		std::string_view featureSuffix);

	// 写 ANI 草图（FracMinHash 边车）
	bool writeAniSketch(const std::string &db_prefix,
	                   const IMCFConfig &imcfConfig,
	                   const std::vector<std::vector<std::string>> &indexToTaxid,
	                   uint32_t scale,
	                   std::string_view featureSuffix,
	                   bool verbose);
	void run(BuildConfig config);
}
