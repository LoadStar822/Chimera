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
#include <interleaved-cuckoo-filter.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <robin_hood.h>
#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <chrono>
#include <dna4_traits.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/utility.hpp>
#include <thread>
#include <vector>
#include <mutex>
#include <cstdint>
#include <execution>
#include <HyperLogLog.hpp>
#include <toolbox.hpp>
#include <hierarchical-interleaved-cuckoo-filter.h>
#include <interleaved-merged-cuckoo-filter.h>

namespace fs = std::filesystem;

namespace ChimeraBuild {
	void run(BuildConfig config);
	size_t computeLayoutHICF(HICFConfig& hicfConfig,
		BuildConfig& config,
		chimera::hicf::DataStore& data);
}