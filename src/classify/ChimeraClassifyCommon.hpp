// Internal helpers shared across Chimera classify modules.
#pragma once

#include "ChimeraClassify.hpp"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cmath>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <memory>
#include <mutex>
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

inline constexpr size_t kInvalidLength = std::numeric_limits<size_t>::max();
inline constexpr size_t kClassifyQueueMaxPendingBytes =
    64ull * 1024ull * 1024ull;

struct QueueThrottle {
  std::mutex mutex;
  std::condition_variable cv;
  size_t pending_bytes{0};
  size_t peak_pending_bytes{0};
  size_t max_pending_bytes{kClassifyQueueMaxPendingBytes};
};

inline size_t estimate_batch_bytes(const batchReads &batch) {
  size_t total = sizeof(batchReads);
  for (const auto &id : batch.ids) {
    total += id.size();
  }
  for (const auto &seq : batch.seqs) {
    total += seq.size() * sizeof(seqan3::dna4);
  }
  for (const auto &seq : batch.seqs2) {
    total += seq.size() * sizeof(seqan3::dna4);
  }
  return total;
}

inline void acquire_queue_slot(QueueThrottle *throttle, size_t bytes) {
  if (throttle == nullptr) {
    return;
  }
  bytes = std::max<size_t>(bytes, 1);
  std::unique_lock<std::mutex> lock(throttle->mutex);
  throttle->cv.wait(lock, [&]() {
    return throttle->pending_bytes == 0 ||
           throttle->pending_bytes + bytes <= throttle->max_pending_bytes;
  });
  throttle->pending_bytes += bytes;
  throttle->peak_pending_bytes =
      std::max(throttle->peak_pending_bytes, throttle->pending_bytes);
}

inline void release_queue_slot(QueueThrottle *throttle, size_t bytes) {
  if (throttle == nullptr) {
    return;
  }
  bytes = std::max<size_t>(bytes, 1);
  {
    std::lock_guard<std::mutex> lock(throttle->mutex);
    throttle->pending_bytes =
        (bytes >= throttle->pending_bytes) ? 0
                                           : throttle->pending_bytes - bytes;
  }
  throttle->cv.notify_one();
}

inline double clamp01(double x) { return std::clamp(x, 0.0, 1.0); }

inline double lerp(double lo, double hi, double t) {
  const double tt = clamp01(t);
  return lo + (hi - lo) * tt;
}

inline double smoothstep(double x, double edge0, double edge1) {
  if (!(edge1 > edge0)) {
    return 0.0;
  }
  const double t = clamp01((x - edge0) / (edge1 - edge0));
  return t * t * (3.0 - 2.0 * t);
}

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
  bool freq_trusted{true};
  // Optional NCBI taxonomy helper (when available in current environment).
  const NcbiTaxdump *ncbiTaxdump{nullptr};
  // Optional NCBI-only: map internal tid_id -> representative tid_id of its
  // species (collapsing strain/subspecies taxids). This avoids per-hit taxdump
  // lookups in hot loops.
  const std::vector<uint32_t> *tid2speciesRep{nullptr};
  // Optional NCBI-only: map internal tid_id -> genus taxid for fast
  // same-genus grouping in coarse candidate routing.
  const std::vector<uint32_t> *tid2genus{nullptr};
  // Optional NCBI-only: map internal tid_id -> species group id (numeric NCBI
  // species taxid, or a stable synthetic id for non-numeric taxids). Used for
  // computing deg/exclusivity at the species level without changing output
  // taxids or breaking presence coverage meta.

  bool enabled() const { return freqSketch != nullptr; }
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
inline std::size_t effective_df_bins(std::size_t deg_effective,
                                     std::size_t df_bins) {
  if (deg_effective == 1 && df_bins > 1) {
    return 1;
  }
  return df_bins;
}

// Choose which DF to use for IDF computation.
//
// Empirically, using df_eff for IDF helps tail-rich "short contigs" avoid
// posterior_weight rejects, while using df_bins for long reads reduces FP
// sensitivity (especially under subset/topBins).
//
// Policy:
// - readLen<=min_len: allow full fragmented-species correction (w=tail_s)
// - readLen>=max_len: disable correction (w=0)
// - middle zone: smoothstep transition
inline double df_for_idf(std::size_t df_bins, std::size_t df_eff,
                         std::size_t readLen, double tail_risk_s,
                         std::size_t max_len = 4096,
                         std::size_t min_len = 1024) {
  const double w_s = clamp01(tail_risk_s);
  double w_len = 0.0;
  if (readLen <= min_len) {
    w_len = 1.0;
  } else if (readLen >= max_len) {
    w_len = 0.0;
  } else {
    w_len = 1.0 - smoothstep(static_cast<double>(readLen),
                             static_cast<double>(min_len),
                             static_cast<double>(max_len));
  }
  const double w = clamp01(w_s * w_len);
  return ((1.0 - w) * static_cast<double>(df_bins)) +
         (w * static_cast<double>(df_eff));
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
inline double idf_raw_from_df_bins(double totalBins, double df_bins) {
  const double denom = std::max(0.0, df_bins) + 1.0;
  return std::log2((totalBins + 1.0) / denom);
}

// Local-unique edge is defined strictly at the bin level (df_bins==1). We do
// NOT treat fragmented species (df_bins>1 but df_eff==1) as locally unique,
// because this would change presence sketch semantics.
inline bool is_local_unique_edge(std::size_t deg_effective,
                                 std::size_t df_bins) {
  return deg_effective == 1 && df_bins == 1;
}

inline double clamp_idf(double idf_raw, double tail_risk_s, double idf_max,
                        double idf_power) {
  const double s = clamp01(tail_risk_s);
  const double idf_min = 0.5 * (1.0 - s);
  const double idf_max_eff = std::max(idf_min, idf_max);
  const double idf0 = std::clamp(idf_raw, idf_min, idf_max_eff);
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
    double pi_hi, double pi_lo, std::size_t avgLen, double relax_strength,
    std::size_t L0 = 2000, std::size_t L1 = 800) {
  AutoPostPiMinTune out;
  out.pi_hi = pi_hi;
  out.pi_lo = pi_lo;

  const double relax = clamp01(relax_strength);
  if (!(pi_hi > 0.0) || !(pi_lo > 0.0) || avgLen == 0 || relax <= 0.0 ||
      L0 <= L1 || pi_lo >= pi_hi) {
    out.tuned = pi_hi;
    return out;
  }

  double t_len = 0.0;
  if (avgLen >= L0) {
    t_len = 0.0;
  } else if (avgLen <= L1) {
    t_len = 1.0;
  } else {
    const double denom = static_cast<double>(L0) - static_cast<double>(L1);
    t_len = (static_cast<double>(L0) - static_cast<double>(avgLen)) / denom;
    t_len = std::clamp(t_len, 0.0, 1.0);
  }
  double t = clamp01(t_len * relax);
  out.t = t;

  if (t <= 0.0) {
    out.tuned = pi_hi;
    return out;
  }

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

struct TailRiskProbeStats {
  double shannon{0.0};
  double top_mass{0.0};
  double eff_species{0.0};
  std::size_t total{0};
  std::size_t unclassified{0};
};

inline TailRiskProbeStats
compute_tail_risk_probe_stats(const std::vector<double> &counts,
                              std::size_t unclassified, std::size_t topk) {
  TailRiskProbeStats stats;
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

inline double compute_tail_risk_u(double tail_risk_r, double eff_species,
                                  double r_anchor, double e_anchor) {
  const double r = clamp01(tail_risk_r);
  const double r0 = std::max(r_anchor, 1e-9);
  const double u_r = clamp01(r / (r + r0));

  const double e0 = std::max(e_anchor, 1e-9);
  const double e = std::max(0.0, eff_species);
  const double u_e = clamp01((e - e0) / (e + e0));
  return clamp01(1.0 - ((1.0 - u_r) * (1.0 - u_e)));
}

inline double compute_tail_risk_s(double tail_risk_u, int beta) {
  const double u = clamp01(tail_risk_u);
  const int b = std::max(1, beta);
  const double up = std::pow(u, static_cast<double>(b));
  const double down = std::pow(1.0 - u, static_cast<double>(b));
  return clamp01(up / (up + down + 1e-12));
}

chimera::feature::Params prepare_feature_params_for_classify(
    const ChimeraBuild::IMCFConfig &imcfConfig, size_t &feature_min_len);

// --- fast taxid dictionary ---
struct TaxDict {
  std::vector<std::vector<uint32_t>> idx2id;  // [bin][species] -> tid_id
  std::vector<std::string> id2str;            // tid_id -> taxid string
  std::vector<std::vector<uint32_t>> tid2bin; // tid_id -> 所在 bin 列表
  robin_hood::unordered_flat_map<std::string, uint32_t> str2id; // taxid -> tid_id
};

TaxDict build_tax_dict(const std::vector<std::vector<std::string>> &idx2tax);

struct SpoolCandidate {
  uint32_t tid{0};
  double score{0.0};
};

struct SpoolReadRecord {
  std::string id;
  double evaluated{0.0};
  uint32_t best_taxid_hint{0};
  std::string reject_reason;
  std::vector<SpoolCandidate> candidates;
};

inline constexpr uint32_t kSpoolUnclassifiedTid =
    std::numeric_limits<uint32_t>::max();

struct CompactClassifyResult {
  std::string id;
  double evaluated{0.0};
  uint32_t best_taxid_hint{kSpoolUnclassifiedTid};
  std::string reject_reason;
  std::vector<SpoolCandidate> candidates;
};

void write_spool_header(std::ostream &os);
void read_spool_header(std::istream &is, const std::string &path);
void write_spool_record(std::ostream &os, const SpoolReadRecord &record);
void write_spool_record(std::ostream &os,
                        const CompactClassifyResult &record);
bool read_spool_record(std::istream &is, SpoolReadRecord &record);

void parseReads(std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
                ClassifyConfig config, FileInfo &fileInfo,
                size_t max_reads = 0,
                std::vector<QueueThrottle> *queueThrottles = nullptr);

void loadFilter(const std::string &input_file,
                chimera::imcf::InterleavedMergedCuckooFilter &imcf,
                ChimeraBuild::IMCFConfig &imcfConfig,
                std::vector<std::vector<std::string>> &indexToTaxid,
                chimera::presence::CoverageMeta *coverageMeta = nullptr);

void saveResult(const std::vector<classifyResult> &classifyResults,
                const ClassifyConfig &config);

void writeResultRecord(std::ostream &os, const classifyResult &result,
                       std::ostringstream &postTopkOss);

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
                  double unique_strength, bool localUniqueEdge,
                  uint32_t unique_bucket,
                  uint32_t breadth_bucket);
  void add_targets(const std::vector<uint32_t> &tids, double hit_weight,
                   double score_weight, double unique_strength,
                   bool localUniqueEdge, uint32_t unique_bucket,
                   uint32_t breadth_bucket);
  void add_read_support(uint32_t tid, bool uniqueRead);
};

struct PresenceSummary {
  uint32_t sketchBits{0};
  size_t sketchWords{0};
  robin_hood::unordered_flat_map<uint32_t, PresenceStats> stats;

  explicit PresenceSummary(uint32_t breadthBits = 0);
  void merge(const PresenceAccumulator &acc);
};

struct PresenceDecision {
  robin_hood::unordered_flat_map<uint32_t, double> posteriors;
  robin_hood::unordered_flat_map<uint32_t, double> logPosteriors;
  double threshold{1.0};
};

struct TaxpoolGenusPick {
  uint32_t genus_id{0u};
  uint32_t rep_id{0u};
  uint32_t rare_support{0u};
  uint64_t coarse_score{0u};
};

struct ProcessScratch {
  std::vector<uint64_t> hashs1;
  std::vector<uint64_t> sampleVals;
  std::vector<uint64_t> routeVals;
  std::vector<std::pair<uint32_t, uint64_t>> sampleScored;
  std::vector<std::pair<uint32_t, uint64_t>> routeScored;
  std::vector<std::vector<uint32_t>> sampleCount;
  std::vector<std::pair<uint32_t, uint16_t>> touchedS;
  std::vector<std::pair<uint32_t, uint32_t>> rankedBins;
  std::vector<size_t> deferredEval;
  std::vector<std::pair<uint32_t, uint64_t>> weighted;
  std::vector<uint32_t> topBins;
  std::vector<std::pair<uint32_t, uint64_t>> repRanked;
  std::vector<std::pair<uint32_t, uint64_t>> genusRanked;
  std::vector<uint32_t> repPool;
  std::vector<uint32_t> dominantCandidates;
  std::vector<TaxpoolGenusPick> nonDominantGenusPicks;
  std::vector<std::pair<uint32_t, double>> rankedTidScores;
  std::vector<uint32_t> minimizerTids;
  std::vector<uint32_t> minimizerBins;
  std::vector<uint32_t> minimizerTidEpoch;
  uint32_t minimizerTidEpochValue{0};
  std::vector<uint32_t> minimizerBinEpoch;
  uint32_t minimizerBinEpochValue{0};
  std::vector<double> tidScoreDense;
  std::vector<uint32_t> tidScoreEpoch;
  std::vector<uint32_t> activeTidScores;
  uint32_t tidScoreEpochValue{0};
  std::vector<double> uniqueHitsDense;
  std::vector<uint32_t> uniqueHitsEpoch;
  std::vector<uint32_t> activeUniqueHits;
  uint32_t uniqueHitsEpochValue{0};
  std::vector<uint32_t> routed;
  std::vector<uint32_t> rareRepHits;
  std::vector<uint32_t> rareGenusHits;
};

PresenceDecision evaluate_presence_coverage(
    const PresenceSummary &summary, const TaxDict &tax,
    const ClassifyConfig &config, const chimera::presence::CoverageMeta &meta,
    size_t totalReads, size_t meanReadLen);

void processSequence(
    const std::vector<uint64_t> &hashs1, size_t readLen,
    ChimeraBuild::IMCFConfig &imcfConfig, const TaxDict &tax,
    ClassifyConfig &config, const WeightingContext &weightCtx,
    const AutoClassifyPolicy &autoPolicy, GroupHeat &heat,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const std::string &id,
    std::vector<classifyResult> *classifyResults,
    std::vector<CompactClassifyResult> *compactResults, FileInfo &fileInfo,
    PresenceAccumulator *presenceAcc, ProcessScratch &scratch);

void processBatch(
    batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
    const TaxDict &tax, ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<classifyResult> &classifyResults,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    FileInfo &fileInfo, GroupHeat &heat, const WeightingContext &weightCtx,
    const AutoClassifyPolicy &autoPolicy, PresenceAccumulator *presenceAcc,
    ProcessScratch &scratch);

void processBatchCompact(
    batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
    const TaxDict &tax, ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<CompactClassifyResult> &classifyResults,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    FileInfo &fileInfo, GroupHeat &heat, const WeightingContext &weightCtx,
    const AutoClassifyPolicy &autoPolicy, PresenceAccumulator *presenceAcc,
    ProcessScratch &scratch);

void classify_streaming(
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
    ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    std::atomic<bool> &producer_done,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, const AutoClassifyPolicy &autoPolicy,
    PresenceSummary *presenceSummary,
    std::vector<QueueThrottle> *queueThrottles = nullptr);

void classify_streaming_spool(
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
    ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const TaxDict &tax,
    const std::vector<std::string> &spoolPaths, FileInfo &fileInfo,
    std::atomic<bool> &producer_done,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, const AutoClassifyPolicy &autoPolicy,
    PresenceSummary *presenceSummary,
    std::vector<QueueThrottle> *queueThrottles = nullptr);

void classify(ChimeraBuild::IMCFConfig &imcfConfig,
              std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
              ClassifyConfig &config,
              chimera::imcf::InterleavedMergedCuckooFilter &imcf,
              const TaxDict &tax, std::vector<classifyResult> &classifyResults,
              FileInfo &fileInfo,
              const chimera::feature::Params &feature_params,
              size_t feature_min_len, const WeightingContext &weightCtx,
              const AutoClassifyPolicy &autoPolicy);

} // namespace ChimeraClassify
