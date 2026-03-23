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
#include <string_view>
#include <string>
#include <omp.h>
#include <robin_hood.h>
#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <dna4_traits.hpp>
#include <thread>
#include <vector>
#include <mutex>
#include <cstdint>
#include <interleaved-merged-cuckoo-filter.h>

namespace fs = std::filesystem;

namespace ChimeraBuild {
	struct HashFrequencyContext;
	struct FeatureBuildLayout;
}

namespace chimera::presence {
	struct CoverageMeta;
}

namespace ChimeraBuild {
	std::vector<std::vector<std::string>> buildIMCF(
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		const std::vector<chimera::imcf::Group>& groups,
		const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount,
		const HashFrequencyContext* hashFreqContext,
		const FeatureBuildLayout* featureLayout,
		uint16_t effectiveSpan,
		uint16_t refReadLen,
		uint32_t uniqueDegThreshold,
		chimera::presence::CoverageMeta* coverageMeta);
	void run(BuildConfig config);
}
