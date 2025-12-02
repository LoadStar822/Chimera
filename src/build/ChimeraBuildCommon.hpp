// Internal helpers shared by Chimera build modules.
#pragma once

#include "ChimeraBuild.hpp"

#include <atomic>
#include <memory>
#include <string_view>

#include <utils/CountMinSketch.hpp>
#include <utils/FeatureHasher.hpp>
#include <utils/PresenceModel.hpp>

namespace ChimeraBuild {

uint64_t mix64(uint64_t value);

struct TaxidShardPlan {
  size_t groupIndex;
  size_t slotIndex;
  uint64_t expectedCount;
};

struct HashFrequencyContext {
  double quantile{0.999};
  std::unique_ptr<CountMinSketch> sketch;
  chimera::presence::HashFreqStats stats{};
  std::atomic<uint64_t> passA_total_hashes{0};
  std::atomic<uint64_t> passB_total_hashes{0};
  std::atomic<uint64_t> passB_filtered_hashes{0};

  bool enabled() const { return static_cast<bool>(sketch); }
};

void print_build_time(long long milliseconds);

chimera::feature::Method resolve_feature_method(const BuildConfig &config);
chimera::feature::Params make_feature_params(const BuildConfig &config,
                                             chimera::feature::Method &selected,
                                             uint64_t &seed_out);

void build_hash_frequency_sketch(
    const BuildConfig &config,
    const robin_hood::unordered_flat_map<std::string,
                                         std::vector<std::string>> &inputFiles,
    HashFrequencyContext &context);

void parseInputFile(const std::string &filePath,
                    robin_hood::unordered_flat_map<
                        std::string, std::vector<std::string>> &inputFiles,
                    robin_hood::unordered_flat_map<std::string, uint64_t>
                        &hashCount,
                    FileInfo &fileInfo);

void createOrResetDirectory(const std::string &dir, const BuildConfig &config);

void syncmer_count(
    BuildConfig &config,
    robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    FileInfo &fileInfo, HashFrequencyContext *hashFreqContext,
    robin_hood::unordered_flat_map<std::string, uint64_t> *bpCount);

std::string tmp_hash_path(const std::string &taxid,
                          std::string_view suffix);
std::string tmp_hash_thread_path(const std::string &taxid,
                                 std::string_view suffix, int thread_id);

void saveIMCF(chimera::imcf::InterleavedMergedCuckooFilter &imcf,
              const std::string &output_file,
              std::vector<std::vector<std::string>> &indexToTaxid,
              IMCFConfig &imcfConfig, bool needIndexStructures,
              const chimera::presence::CoverageMeta *presenceMeta);

std::vector<std::vector<std::string>> buildIMCF(
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    const std::vector<chimera::imcf::Group> &groups,
    const robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    std::string_view featureSuffix,
    const HashFrequencyContext *hashFreqContext,
    const robin_hood::unordered_flat_map<std::string, uint64_t> *bpCount,
    uint16_t effectiveSpan, uint16_t refReadLen, uint32_t uniqueDegThreshold,
    chimera::presence::CoverageMeta *coverageMeta);

chimera::presence::CoverageMeta compute_presence_meta(
    const robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    std::string_view featureSuffix,
    const HashFrequencyContext *hashFreqContext,
    const robin_hood::unordered_flat_map<std::string, uint64_t> *bpCount,
    uint16_t effectiveSpan, uint16_t refReadLen,
    uint32_t uniqueDegThreshold, uint16_t threads);

} // namespace ChimeraBuild
