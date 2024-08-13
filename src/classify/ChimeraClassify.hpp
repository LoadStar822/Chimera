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
 * Last Modified: 2024-08-10
 *
 * Description:
 *  This is a simple C++ program that outputs "Hello, World!".
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#include <classifyConfig.hpp>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <seqan3/contrib/stream/bgzf.hpp>
#include <dna4_traits.hpp>
#include <omp.h>
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
#include <seqan3/utility/views/chunk.hpp>
#include <chrono>
#include <dna4_traits.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/utility.hpp>
#include <thread>
#include <vector>
#include <mutex>
#include <queue>
#include <buildConfig.hpp>
namespace ChimeraClassify {
	void run(ClassifyConfig config);
}