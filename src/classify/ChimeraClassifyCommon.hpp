// Internal helpers shared across Chimera classify modules.
#pragma once

#include "ChimeraClassify.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
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

// Optional NCBI taxonomy helper: map a taxid (including strain/subspecies) to its
// species-level ancestor taxid by walking nodes.dmp parent pointers.
struct NcbiTaxdump {
  std::vector<uint32_t> parent;
  std::vector<uint8_t> is_species; // 1 if rank == species

  bool enabled() const {
    return !parent.empty() && parent.size() == is_species.size();
  }

  uint32_t to_species(uint32_t tid) const {
    if (!enabled() || tid == 0 || tid >= is_species.size()) {
      return tid;
    }
    if (is_species[tid]) {
      return tid;
    }
    uint32_t cur = tid;
    // Walk up until we hit species or root; keep a small step cap to avoid cycles.
    for (int steps = 0; steps < 128; ++steps) {
      if (cur == 0 || cur >= parent.size()) {
        break;
      }
      uint32_t p = parent[cur];
      if (p == 0 || p == cur) {
        break;
      }
      cur = p;
      if (cur < is_species.size() && is_species[cur]) {
        return cur;
      }
    }
    return tid;
  }
};

struct WeightingContext {
  const CountMinSketch *freqSketch{nullptr};
  chimera::presence::HashFreqStats freqStats{};
  double freqQuantile{0.0};
  const std::unordered_map<std::string, double> *sampleWeights{nullptr};
  bool freq_trusted{true};
  // Optional NCBI taxonomy helper (when available in current environment).
  const NcbiTaxdump *ncbiTaxdump{nullptr};
  // Optional NCBI-only: map internal tid_id -> representative tid_id of its
  // species (collapsing strain/subspecies taxids). This avoids per-hit taxdump
  // lookups in hot loops.
  const std::vector<uint32_t> *tid2speciesRep{nullptr};
  // Optional NCBI-only: map internal tid_id -> species group id (numeric NCBI
  // species taxid, or a stable synthetic id for non-numeric taxids). Used for
  // computing deg/exclusivity at the species level without changing output
  // taxids or breaking presence coverage meta.
  const std::vector<uint32_t> *tid2speciesGroup{nullptr};

  bool enabled() const { return freqSketch != nullptr; }
  bool has_sample_weights() const { return sampleWeights != nullptr; }
};

inline bool allow_unique_edge(std::size_t deg_effective, std::size_t df_bins,
                              bool has_freq, bool freq_trusted,
                              uint32_t df_est, uint32_t unique_deg_threshold) {
  if (deg_effective != 1 || df_bins != 1) {
    return false;
  }
  if (!has_freq || !freq_trusted) {
    return true;
  }
  const uint32_t threshold = std::max<uint32_t>(1u, unique_deg_threshold);
  return df_est <= threshold;
}

inline bool allow_low_df_boost(std::size_t df_bins, bool has_freq,
                               bool freq_trusted, uint32_t df_est,
                               uint32_t unique_deg_threshold) {
  if (df_bins > 2) {
    return false;
  }
  if (!has_freq || !freq_trusted) {
    return true;
  }
  const uint32_t threshold = std::max<uint32_t>(1u, unique_deg_threshold);
  const uint32_t df_gate = static_cast<uint32_t>(2u * threshold + 1u);
  return df_est <= df_gate;
}

inline std::vector<uint64_t> select_rare_route_values(
    const std::vector<std::pair<uint32_t, uint64_t>> &scored,
    std::size_t budget) {
  std::vector<uint64_t> out;
  if (budget == 0 || scored.empty()) {
    return out;
  }
  std::vector<std::pair<uint32_t, uint64_t>> ordered = scored;
  std::sort(ordered.begin(), ordered.end(),
            [](const auto &a, const auto &b) {
              if (a.first != b.first) {
                return a.first < b.first;
              }
              return a.second < b.second;
            });
  robin_hood::unordered_flat_set<uint64_t> seen;
  seen.reserve(std::min<std::size_t>(budget * 2, ordered.size()));
  for (const auto &kv : ordered) {
    if (!seen.insert(kv.second).second) {
      continue;
    }
    out.push_back(kv.second);
    if (out.size() >= budget) {
      break;
    }
  }
  return out;
}

inline double compute_head_mass(
    const std::vector<std::pair<uint32_t, uint32_t>> &ranked,
    std::size_t topN, uint64_t total) {
  if (total == 0 || ranked.empty() || topN == 0) {
    return 0.0;
  }
  const std::size_t limit = std::min<std::size_t>(topN, ranked.size());
  uint64_t sum = 0;
  for (std::size_t i = 0; i < limit; ++i) {
    sum += ranked[i].second;
  }
  return static_cast<double>(sum) / static_cast<double>(total);
}

inline std::size_t compute_candidate_cap(std::size_t baseCap,
                                         std::size_t maxCap,
                                         double headMass,
                                         double headMassThresh,
                                         std::size_t rankedSize) {
  if (rankedSize == 0) {
    return 0;
  }
  if (baseCap == 0) {
    baseCap = 1;
  }
  if (maxCap < baseCap) {
    maxCap = baseCap;
  }
  std::size_t cap = std::min<std::size_t>(baseCap, rankedSize);
  if (rankedSize > baseCap && headMass < headMassThresh) {
    cap = std::min<std::size_t>(maxCap, rankedSize);
  }
  return cap;
}

struct LowDivStats {
  double shannon{0.0};
  double top_mass{0.0};
  double eff_species{0.0};
  std::size_t total{0};
  std::size_t unclassified{0};
};

inline LowDivStats compute_low_div_stats(const std::vector<double> &counts,
                                         std::size_t unclassified,
                                         std::size_t topk) {
  LowDivStats stats;
  stats.unclassified = unclassified;
  double total = 0.0;
  for (double c : counts) {
    if (c > 0.0) {
      total += c;
    }
  }
  stats.total = static_cast<std::size_t>(std::llround(total + unclassified));
  if (!(total > 0.0)) {
    return stats;
  }

  std::vector<double> rels;
  rels.reserve(counts.size());
  for (double c : counts) {
    if (c <= 0.0) {
      continue;
    }
    double p = c / total;
    stats.shannon -= p * std::log(p);
    rels.push_back(p);
  }

  std::sort(rels.begin(), rels.end(), std::greater<double>());
  if (topk > 0 && !rels.empty()) {
    const std::size_t limit = std::min(topk, rels.size());
    for (std::size_t i = 0; i < limit; ++i) {
      stats.top_mass += rels[i];
    }
  }
  stats.eff_species = (stats.shannon > 0.0) ? std::exp(stats.shannon) : 0.0;
  return stats;
}

inline bool is_low_diversity(const LowDivStats &stats,
                             double top_mass_thresh,
                             double eff_species_max) {
  return (stats.top_mass >= top_mass_thresh) &&
         (stats.eff_species <= eff_species_max);
}

inline void apply_low_div_overrides(ClassifyConfig &config) {
  config.firstFilterBeta = 0.5;
  config.firstFilterBeta_user = true;
  config.preEmTopK = 256;
  config.exclusive_gamma = 0.0;
  config.em_conf_power = 1.0;
  config.em_coexist_penalty = 0.0;
  // Low-diversity: still dump a small POST_TOPK so downstream profile can
  // leverage a second evidence stream without exploding I/O.
  config.dump_post_topk = 16;
  config.low_div_active = true;
}

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

void parseReads(std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
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
  uint64_t score{0};
  uint64_t uniqueScore{0};
  uint64_t hits{0};
  uint64_t uniqueHits{0};
  uint64_t readHits{0};
  uint64_t uniqueReads{0};
  std::vector<uint64_t> decoys;
  std::vector<uint64_t> uniqueSketch;
  std::vector<uint64_t> breadthSketch;
};

struct PresenceAccumulator {
  size_t decoyReps{0};
  uint32_t sketchBits{0};
  size_t sketchWords{0};
  robin_hood::unordered_flat_map<uint32_t, PresenceStats> stats;

  explicit PresenceAccumulator(size_t reps = 0, uint32_t breadthBits = 0);
  PresenceStats &touch(uint32_t tid);
  void add_target(uint32_t tid, double hit_weight, double score_weight,
                  bool uniqueEdge, bool localUniqueEdge,
                  uint32_t unique_bucket,
                  uint32_t breadth_bucket);
  void add_read_support(uint32_t tid, bool uniqueRead);
  void add_decoy(size_t rep, uint32_t tid, double weight);
  bool sketches_enabled() const { return sketchWords > 0; }
};

struct PresenceSummary {
  size_t decoyReps{0};
  uint32_t sketchBits{0};
  size_t sketchWords{0};
  robin_hood::unordered_flat_map<uint32_t, PresenceStats> stats;

  explicit PresenceSummary(size_t reps = 0, uint32_t breadthBits = 0);
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

PresenceDecision evaluate_presence_coverage_unique(
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
    std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
    ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    std::atomic<bool> &producer_done,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, PresenceSummary *presenceSummary,
    uint64_t decoySeed);

void classify(ChimeraBuild::IMCFConfig &imcfConfig,
              std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
              ClassifyConfig &config,
              chimera::imcf::InterleavedMergedCuckooFilter &imcf,
              std::vector<std::vector<std::string>> &indexToTaxid,
              const TaxDict &tax, std::vector<classifyResult> &classifyResults,
              FileInfo &fileInfo,
              const chimera::feature::Params &feature_params,
              size_t feature_min_len, const WeightingContext &weightCtx);

} // namespace ChimeraClassify
