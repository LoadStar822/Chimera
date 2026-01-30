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
  std::vector<uint8_t> is_genus;   // 1 if rank == genus

  bool enabled() const {
    return !parent.empty() && parent.size() == is_species.size() &&
           parent.size() == is_genus.size();
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

  uint32_t to_genus(uint32_t tid) const {
    if (!enabled() || tid == 0 || tid >= is_genus.size()) {
      return 0;
    }
    if (is_genus[tid]) {
      return tid;
    }
    uint32_t cur = tid;
    for (int steps = 0; steps < 128; ++steps) {
      if (cur == 0 || cur >= parent.size()) {
        break;
      }
      uint32_t p = parent[cur];
      if (p == 0 || p == cur) {
        break;
      }
      cur = p;
      if (cur < is_genus.size() && is_genus[cur]) {
        return cur;
      }
    }
    return 0;
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

// Treat a minimizer hit as "species-unique" even when the underlying species
// spans multiple bins (df_bins>1). This corrects bin-level fragmentation that
// would otherwise downweight IDF and disable unique-edge boosts.
//
// This is strictly a high-div behavior: low-div branch keeps the original
// df_bins to avoid destabilizing low-div scoring.
inline std::size_t effective_df_bins(std::size_t deg_effective,
                                     std::size_t df_bins,
                                     bool low_div_active) {
  if (low_div_active) {
    return df_bins;
  }
  if (deg_effective == 1 && df_bins > 1) {
    return 1;
  }
  return df_bins;
}

// Choose which DF to use for IDF computation.
//
// Empirically, using df_eff for IDF helps high-div "short contigs" avoid
// posterior_weight rejects, while using df_bins for long reads reduces FP
// sensitivity (especially under subset/topBins).
//
// Policy:
// - low-div branch: always df_bins (keep ATCC-like behavior stable)
// - high-div: if (df_eff < df_bins) and readLen <= max_len, use df_eff;
//             otherwise use df_bins
inline std::size_t df_for_idf(std::size_t df_bins, std::size_t df_eff,
                              bool low_div_active, std::size_t readLen,
                              std::size_t max_len = 4096) {
  if (low_div_active) {
    return df_bins;
  }
  if (df_eff < df_bins && readLen <= max_len) {
    return df_eff;
  }
  return df_bins;
}

// Scale down the unique-edge bonus when the "uniqueness" comes from a
// fragmented (multi-bin) species, to avoid over-boosting.
// - df_bins==1  => keep full base_bonus (e.g. 3.0)
// - df_bins>1   => smoothly decay towards 1.0 as df_bins grows
inline double unique_edge_bonus(double base_bonus, std::size_t df_bins) {
  if (df_bins <= 1 || !(base_bonus > 1.0)) {
    return base_bonus;
  }
  const double denom = 1.0 + std::log2(static_cast<double>(df_bins));
  if (!(denom > 0.0)) {
    return base_bonus;
  }
  const double scaled = 1.0 + (base_bonus - 1.0) / denom;
  return std::clamp(scaled, 1.0, base_bonus);
}

// Compute raw IDF from a DF count. Callers typically pass df_bins, but may
// choose df_eff under conservative guards (see df_for_idf()).
inline double idf_raw_from_df_bins(double totalBins, std::size_t df_bins) {
  const double denom = static_cast<double>(df_bins) + 1.0;
  return std::log2((totalBins + 1.0) / denom);
}

// Local-unique edge is defined strictly at the bin level (df_bins==1). We do
// NOT treat fragmented species (df_bins>1 but df_eff==1) as locally unique,
// because this would change presence sketch semantics.
inline bool is_local_unique_edge(std::size_t deg_effective,
                                 std::size_t df_bins) {
  return deg_effective == 1 && df_bins == 1;
}

inline double clamp_idf(double idf_raw, bool low_div_active, double idf_max,
                        double idf_power) {
  const double idf_min = low_div_active ? 0.5 : 0.0;
  const double idf_max_eff = std::max(idf_min, idf_max);
  const double idf0 = std::clamp(idf_raw, idf_min, idf_max_eff);
  if (low_div_active) {
    return idf0;
  }
  if (idf_max_eff <= 0.0) {
    return 0.0;
  }
  const double p = std::clamp(idf_power, 1.0, 2.0);
  if (p <= 1.0) {
    return idf0;
  }
  if (p >= 2.0) {
    return (idf0 * idf0) / idf_max_eff;
  }
  // Evidence-adaptive downweighting for high-DF (low-IDF) minimizers while
  // keeping the same upper bound (idf_max). We exponentiate on the normalized
  // [0,1] scale:
  //   idf_eff = idf_max * (idf0/idf_max)^p
  const double x = (idf_max_eff > 0.0) ? (idf0 / idf_max_eff) : 0.0;
  return idf_max_eff * std::pow(x, p);
}

struct AutoPostPiMinTune {
  double tuned{0.0};
  double pi_hi{0.0};
  double pi_lo{0.0};
  double t{0.0}; // interpolation factor in [0,1]
  bool applied{false};
};

// High-div only: choose a more permissive (smaller) post_pi_min when the sample
// has short average read length (weak evidence), which helps reduce
// posterior_weight/em_post rejects without touching posterior thresholds.
//
// We use a smooth log-space interpolation to avoid hard length boundaries:
//   avgLen >= L0  => pi=pi_hi
//   avgLen <= L1  => pi=pi_lo
//   else          => log10(pi) = (1-t)log10(pi_hi) + tlog10(pi_lo)
inline AutoPostPiMinTune tune_post_pi_min_by_avg_len(
    double pi_hi, double pi_lo, std::size_t avgLen, bool low_div_active,
    std::size_t L0 = 2000, std::size_t L1 = 800) {
  AutoPostPiMinTune out;
  out.pi_hi = pi_hi;
  out.pi_lo = pi_lo;

  if (low_div_active || !(pi_hi > 0.0) || !(pi_lo > 0.0) || avgLen == 0 ||
      L0 <= L1 || pi_lo >= pi_hi) {
    out.tuned = pi_hi;
    return out;
  }

  if (avgLen >= L0) {
    out.tuned = pi_hi;
    return out;
  }

  if (avgLen <= L1) {
    out.tuned = pi_lo;
    out.t = 1.0;
    out.applied = true;
    return out;
  }

  const double denom =
      static_cast<double>(L0) - static_cast<double>(L1);
  double t = (static_cast<double>(L0) - static_cast<double>(avgLen)) / denom;
  t = std::clamp(t, 0.0, 1.0);
  out.t = t;

  const double log_hi = std::log10(pi_hi);
  const double log_lo = std::log10(pi_lo);
  const double log_pi = (1.0 - t) * log_hi + t * log_lo;
  double tuned = std::pow(10.0, log_pi);
  if (!(tuned > 0.0)) {
    tuned = pi_lo;
  }
  out.tuned = std::clamp(tuned, pi_lo, pi_hi);
  out.applied = (out.tuned < pi_hi);
  return out;
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
  config.em_conf_power = 1.0;
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
    const ChimeraBuild::IMCFConfig &imcfConfig,
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

void loadFilter(const std::string &input_file,
                chimera::imcf::InterleavedMergedCuckooFilter &imcf,
                ChimeraBuild::IMCFConfig &imcfConfig,
                std::vector<std::vector<std::string>> &indexToTaxid,
                chimera::presence::CoverageMeta *coverageMeta = nullptr);

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
  std::vector<uint64_t> uniqueSketch;
  std::vector<uint64_t> breadthSketch;
};

struct PresenceAccumulator {
  uint32_t sketchBits{0};
  size_t sketchWords{0};
  robin_hood::unordered_flat_map<uint32_t, PresenceStats> stats;

  explicit PresenceAccumulator(uint32_t breadthBits = 0);
  PresenceStats &touch(uint32_t tid);
  void add_target(uint32_t tid, double hit_weight, double score_weight,
                  bool uniqueEdge, bool localUniqueEdge,
                  uint32_t unique_bucket,
                  uint32_t breadth_bucket);
  void add_read_support(uint32_t tid, bool uniqueRead);
  bool sketches_enabled() const { return sketchWords > 0; }
};

struct PresenceSummary {
  uint32_t sketchBits{0};
  size_t sketchWords{0};
  robin_hood::unordered_flat_map<uint32_t, PresenceStats> stats;

  explicit PresenceSummary(uint32_t breadthBits = 0);
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
    PresenceAccumulator *presenceAcc);

void processBatch(
    batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    ClassifyConfig &config, chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<classifyResult> &classifyResults,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    FileInfo &fileInfo, GroupHeat &heat, const WeightingContext &weightCtx,
    PresenceAccumulator *presenceAcc);

void classify_streaming(
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
    ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    std::atomic<bool> &producer_done,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, PresenceSummary *presenceSummary);

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
