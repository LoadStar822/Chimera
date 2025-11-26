// Internal helpers shared across Chimera classify modules.
#pragma once

#include "ChimeraClassify.hpp"

#include <array>
#include <atomic>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <concurrentqueue.h>
#include <interleaved-merged-cuckoo-filter.h>
#include <robin_hood.h>

#include <utils/CountMinSketch.hpp>
#include <utils/FeatureHasher.hpp>
#include <utils/PresenceModel.hpp>

namespace ChimeraClassify {

using FeatureMethod = chimera::feature::Method;
inline constexpr size_t kInvalidLength = std::numeric_limits<size_t>::max();

FeatureMethod parse_feature_method_string(std::string feature);
std::string feature_method_to_string(FeatureMethod method);

struct MarginDecision {
  bool accept{false};
  size_t need{0};
  double margin{0.0};
  double ratio{0.0};
};

MarginDecision decide_high_conf(size_t best, size_t second, double eff_eval);

struct WeightingContext {
  const CountMinSketch *freqSketch{nullptr};
  chimera::presence::HashFreqStats freqStats{};
  double freqQuantile{0.0};

  bool enabled() const { return freqSketch != nullptr; }
};

struct ReadStats {
  size_t count{0};
  size_t total_len{0};
  size_t min_len{kInvalidLength};
  size_t max_len{0};

  void update(size_t len);
};

ReadStats sample_read_stats(const ClassifyConfig &config,
                            size_t max_reads = 20000);

chimera::feature::Params prepare_feature_params_for_classify(
    const ChimeraBuild::IMCFConfig &imcfConfig, ClassifyConfig &config,
    FeatureMethod method, size_t &feature_min_len);

// --- fast taxid dictionary ---
struct TaxDict {
  std::vector<std::vector<uint32_t>> idx2id;  // [bin][species] -> tid_id
  std::vector<std::string> id2str;            // tid_id -> taxid string
  std::vector<std::vector<uint32_t>> tid2bin; // tid_id -> 所在 bin 列表
  robin_hood::unordered_flat_map<std::string, uint32_t> str2id; // taxid -> tid_id
};

TaxDict build_tax_dict(const std::vector<std::vector<std::string>> &idx2tax);

void print_classify_time(long long milliseconds);

void parseReads(moodycamel::ConcurrentQueue<batchReads> &readQueue,
                ClassifyConfig config, FileInfo &fileInfo);

struct IMCFIndexStatus {
  bool builtActive{false};
  long long activeMs{0};
};

IMCFIndexStatus loadFilter(const std::string &input_file,
                           chimera::imcf::InterleavedMergedCuckooFilter &imcf,
                           ChimeraBuild::IMCFConfig &imcfConfig,
                           std::vector<std::vector<std::string>> &indexToTaxid,
                           chimera::presence::CoverageMeta *coverageMeta =
                               nullptr);

void saveResult(std::vector<classifyResult> classifyResults,
                ClassifyConfig config);

struct GroupHeat {
  std::vector<uint32_t> score;
  uint32_t decay_shift = 5;   // divide by 32
  uint32_t decay_period = 64; // decay every 64 sequences
  uint32_t counter = 0;

  void ensure(size_t bins);
  void decay_if_needed();
  void boost(uint32_t bin, uint32_t delta);
};

struct PresenceStats {
  double score{0.0};
  double uniqueScore{0.0};
  double hits{0.0};
  uint64_t uniqueHits{0};
  std::vector<double> decoys;
};

struct PresenceAccumulator {
  size_t decoyReps{0};
  robin_hood::unordered_flat_map<uint32_t, PresenceStats> stats;

  explicit PresenceAccumulator(size_t reps = 0);
  PresenceStats &touch(uint32_t tid);
  void add_target(uint32_t tid, double weight, bool uniqueEdge);
  void add_decoy(size_t rep, uint32_t tid, double weight);
};

struct PresenceSummary {
  size_t decoyReps{0};
  robin_hood::unordered_flat_map<uint32_t, PresenceStats> stats;

  explicit PresenceSummary(size_t reps = 0);
  void merge(const PresenceAccumulator &acc);
};

struct PresenceDecision {
  std::unordered_set<uint32_t> accepted;
  robin_hood::unordered_flat_map<uint32_t, double> qValues;
  robin_hood::unordered_flat_map<uint32_t, double> posteriors;
  robin_hood::unordered_flat_map<uint32_t, double> logPosteriors;
  robin_hood::unordered_flat_map<uint32_t, double> lambdaHats;
  double threshold{1.0};
  double noiseMu{0.0};
  double priorPi{0.0};
  size_t tested{0};
  size_t acceptedCount{0};
  size_t decoyPositives{0};
};

PresenceDecision evaluate_presence_coverage(
    const PresenceSummary &summary, const TaxDict &tax,
    const ClassifyConfig &config, const chimera::presence::CoverageMeta &meta,
    size_t totalReads, size_t meanReadLen);

struct PresenceFilterStats {
  size_t trimmedAssignments{0};
  size_t forcedUnclassified{0};
};

PresenceFilterStats apply_presence_filter(
    const PresenceDecision &decision, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo);

void processSequence(
    const std::vector<uint64_t> &hashs1, size_t readLen,
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    ClassifyConfig &config, const WeightingContext &weightCtx, GroupHeat &heat,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const std::string &id,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    PresenceAccumulator *presenceAcc, uint64_t decoySeed);

void processBatch(
    batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    ClassifyConfig &config, chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<classifyResult> &classifyResults,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    FileInfo &fileInfo, GroupHeat &heat, const WeightingContext &weightCtx,
    PresenceAccumulator *presenceAcc, uint64_t decoySeed);

void classify_streaming(
    ChimeraBuild::IMCFConfig &imcfConfig,
    moodycamel::ConcurrentQueue<batchReads> &readQueue, ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    std::atomic<bool> &producer_done,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, PresenceSummary *presenceSummary,
    uint64_t decoySeed);

void classify(ChimeraBuild::IMCFConfig &imcfConfig,
              moodycamel::ConcurrentQueue<batchReads> &readQueue,
              ClassifyConfig &config,
              chimera::imcf::InterleavedMergedCuckooFilter &imcf,
              std::vector<std::vector<std::string>> &indexToTaxid,
              const TaxDict &tax, std::vector<classifyResult> &classifyResults,
              FileInfo &fileInfo,
              const chimera::feature::Params &feature_params,
              size_t feature_min_len, const WeightingContext &weightCtx);

} // namespace ChimeraClassify
