// Internal helpers shared by Chimera build modules.
#pragma once

#include "ChimeraBuild.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <utils/CountMinSketch.hpp>
#include <utils/FeatureHasher.hpp>
#include <utils/PresenceModel.hpp>

namespace ChimeraBuild {

uint64_t mix64(uint64_t value);

struct TaxidShardPlan {
  size_t groupIndex;
  size_t slotIndex;
  uint64_t count;           // exact number of hashes assigned to this shard
  uint64_t fileOffsetStart; // starting hash index inside the taxid virtual stream
};

struct FeatureChunkRef {
  uint16_t threadId{0};
  uint64_t spoolOffset{0};
  uint64_t count{0};
};

struct TaxidFeatureLayout {
  std::vector<FeatureChunkRef> chunkRefs;
  uint64_t totalSignatures{0};
  uint64_t uniqueSignatures{0};
  uint64_t genomeBp{0};
};

struct FeatureBuildLayout {
  std::vector<std::string> taxids;
  std::vector<TaxidFeatureLayout> perTaxid;
  std::vector<std::string> threadSpoolPaths;
};

struct CmsThresholdBitmask {
  uint32_t depth{0};
  uint32_t width{0};
  size_t wordsPerRow{0};
  std::vector<uint64_t> words;

  void clear() {
    depth = 0;
    width = 0;
    wordsPerRow = 0;
    words.clear();
  }

  void reset(uint32_t depthIn, uint32_t widthIn) {
    depth = depthIn;
    width = widthIn;
    wordsPerRow = (static_cast<size_t>(widthIn) + 63u) / 64u;
    words.assign(static_cast<size_t>(depthIn) * wordsPerRow, 0ull);
  }

  bool empty() const { return words.empty(); }

  void set(uint32_t row, uint32_t column) {
    const size_t wordIndex =
        static_cast<size_t>(row) * wordsPerRow + (column >> 6);
    words[wordIndex] |= (1ull << (column & 63u));
  }

  bool test(uint32_t row, uint32_t column) const {
    const size_t wordIndex =
        static_cast<size_t>(row) * wordsPerRow + (column >> 6);
    return (words[wordIndex] >> (column & 63u)) & 1ull;
  }
};

struct HashFrequencyContext {
  struct HashDecision {
    bool filtered{false};
    bool unique{true};
  };

  double quantile{0.999};
  std::unique_ptr<CountMinSketch> sketch;
  chimera::presence::HashFreqStats stats{};
  uint32_t unique_deg_threshold{1};
  CmsThresholdBitmask high_df_mask;
  CmsThresholdBitmask gt_unique_mask;
  std::atomic<uint64_t> passA_total_hashes{0};
  std::atomic<uint64_t> passB_total_hashes{0};
  std::atomic<uint64_t> passB_filtered_hashes{0};

  bool enabled() const { return static_cast<bool>(sketch); }

  bool mask_all_rows_set(const CmsThresholdBitmask &mask,
                         uint64_t hash) const {
    if (mask.depth == 4) {
      return mask.test(0, CountMinSketch::column_index(0, hash, mask.width)) &&
             mask.test(1, CountMinSketch::column_index(1, hash, mask.width)) &&
             mask.test(2, CountMinSketch::column_index(2, hash, mask.width)) &&
             mask.test(3, CountMinSketch::column_index(3, hash, mask.width));
    }
    for (uint32_t row = 0; row < mask.depth; ++row) {
      if (!mask.test(row,
                     CountMinSketch::column_index(row, hash, mask.width))) {
        return false;
      }
    }
    return true;
  }

  bool should_filter_hash(uint64_t hash) const {
    return !high_df_mask.empty() && mask_all_rows_set(high_df_mask, hash);
  }

  bool is_unique_signature(uint64_t hash) const {
    return gt_unique_mask.empty() || !mask_all_rows_set(gt_unique_mask, hash);
  }

  HashDecision decide_hash(uint64_t hash) const {
    if (high_df_mask.empty() && gt_unique_mask.empty()) {
      return {false, true};
    }
    if (!high_df_mask.empty() && !gt_unique_mask.empty() &&
        high_df_mask.depth == gt_unique_mask.depth &&
        high_df_mask.width == gt_unique_mask.width) {
      bool highAll = true;
      bool uniqueAll = true;
      for (uint32_t row = 0; row < high_df_mask.depth; ++row) {
        const uint64_t column =
            CountMinSketch::column_index(row, hash, high_df_mask.width);
        highAll = highAll && high_df_mask.test(row, column);
        uniqueAll = uniqueAll && gt_unique_mask.test(row, column);
      }
      return {highAll, !uniqueAll};
    }
    const bool filtered = should_filter_hash(hash);
    return {filtered, filtered ? false : is_unique_signature(hash)};
  }
};

void print_build_time(long long milliseconds);

chimera::feature::Params make_feature_params(const BuildConfig &config,
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

void feature_count(
    BuildConfig &config,
    robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    FileInfo &fileInfo, HashFrequencyContext *hashFreqContext,
    robin_hood::unordered_flat_map<std::string, uint64_t> *bpCount,
    FeatureBuildLayout *featureLayout);

std::string tmp_hash_spool_path(std::string_view suffix, int thread_id);

void set_tmp_work_dir(const std::filesystem::path &dir);
const std::filesystem::path &tmp_work_dir();

void saveIMCF(chimera::imcf::InterleavedMergedCuckooFilter &imcf,
              const std::string &output_file,
              std::vector<std::vector<std::string>> &indexToTaxid,
              IMCFConfig &imcfConfig,
              const chimera::presence::CoverageMeta *presenceMeta);

std::vector<std::vector<std::string>> buildIMCF(
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    const std::vector<chimera::imcf::Group> &groups,
    const robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    const HashFrequencyContext *hashFreqContext,
    const FeatureBuildLayout *featureLayout,
    uint16_t effectiveSpan, uint16_t refReadLen, uint32_t uniqueDegThreshold,
    chimera::presence::CoverageMeta *coverageMeta,
    size_t groupIndexOffset,
    bool verifyShardTotals,
    const robin_hood::unordered_flat_map<std::string, uint64_t>
        *shardOffsetBase,
    size_t groupRangeBegin,
    size_t groupRangeEnd);

void populateCoverageMeta(const HashFrequencyContext *hashFreqContext,
                          const FeatureBuildLayout *featureLayout,
                          uint16_t effectiveSpan, uint16_t refReadLen,
                          uint32_t uniqueDegThreshold,
                          chimera::presence::CoverageMeta &coverageMeta);

} // namespace ChimeraBuild
