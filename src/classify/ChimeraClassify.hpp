/*
 * -----------------------------------------------------------------------------
 * Filename:      ChimeraClassify.hpp
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
 *  This is ChimeraClassify module for Chimera
 *
 * Version:
 *  1.3
 * -----------------------------------------------------------------------------
 */
#include "classifyConfig.hpp"
#include "kthread.h"
#include <EM.hpp>
#include <VEM.hpp>
#include <buildConfig.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/utility.hpp>
#include <chrono>
#include <cmath>
#include <concurrentqueue.h>
#include <cstdint>
#include <dna4_traits.hpp>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>
#include <interleaved-merged-cuckoo-filter.h>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <queue>
#include <robin_hood.h>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/contrib/stream/bgzf.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

namespace ChimeraClassify {
void run(ClassifyConfig config);
void postEmDecision(
    std::vector<classifyResult> &results, const DecisionConfig &decisionConfig,
    const std::unordered_map<std::string, double> &classWeights);
} // namespace ChimeraClassify
