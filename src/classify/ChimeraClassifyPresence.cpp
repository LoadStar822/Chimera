#include "ChimeraClassifyCommon.hpp"

#include <utils/Parse.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

namespace ChimeraClassify {

namespace {
enum class PostEmPresenceLevel {
  kUnknown = 0,
  kAccepted = 1,
  kRejected = 2,
};

struct PostEmPruneAudit {
  bool pruned_any{false};
};

struct PostEmPruneResult {
  std::vector<std::pair<std::string, double>> posterior;
  PostEmPruneAudit audit;
};

struct PostEmPruneTop1Override {
  std::string taxid;
  double local_pi_min_floor{0.0};
  bool active{false};
};

template <class PresenceLevelFn>
inline PostEmPruneResult prune_post_em_posterior(
    const std::vector<std::pair<std::string, double>> &posterior_full_sorted,
    double min_class_weight,
    const std::unordered_map<std::string, double> &classWeights,
    PresenceLevelFn presence_level, double presence_pi_floor,
    double reject_factor,
    const PostEmPruneTop1Override *top1_override = nullptr) {
  PostEmPruneResult out;
  out.posterior = posterior_full_sorted;

  const double pi_prune = std::min(min_class_weight, 1e-4);
  const bool pruning_active =
      (pi_prune > 0.0 && min_class_weight > 0.0 && !classWeights.empty());
  if (!pruning_active || posterior_full_sorted.empty()) {
    return out;
  }

  std::vector<std::pair<std::string, double>> pruned;
  pruned.reserve(posterior_full_sorted.size());
  double sum = 0.0;
  for (const auto &kv : posterior_full_sorted) {
    const PostEmPresenceLevel pres = presence_level(kv.first);
    double local_pi_min = pi_prune;
    if (pres == PostEmPresenceLevel::kAccepted) {
      local_pi_min = std::min(pi_prune, presence_pi_floor);
    } else if (pres == PostEmPresenceLevel::kRejected) {
      local_pi_min = pi_prune * reject_factor;
    }
    if (top1_override && top1_override->active && !top1_override->taxid.empty() &&
        kv.first == top1_override->taxid && pres == PostEmPresenceLevel::kRejected &&
        top1_override->local_pi_min_floor > 0.0) {
      local_pi_min = std::min(local_pi_min, top1_override->local_pi_min_floor);
    }
    auto it = classWeights.find(kv.first);
    const double w = (it != classWeights.end()) ? it->second : 0.0;
    if (w >= local_pi_min) {
      pruned.push_back(kv);
      sum += kv.second;
    }
  }

  // Soft behavior: if nothing survives, keep the original (full) posterior.
  if (pruned.empty()) {
    return out;
  }

  out.audit.pruned_any = (pruned.size() < posterior_full_sorted.size());
  if (sum > 0.0) {
    for (auto &kv : pruned) {
      kv.second /= sum;
    }
  }
  out.posterior = std::move(pruned);
  return out;
}

struct TopkLookup {
  bool found{false};
  std::size_t rank{0};
  double prob{0.0};
};

inline TopkLookup lookup_taxid_in_topk(
    const std::vector<std::pair<std::string, double>> &posterior_sorted,
    const std::string &taxid, std::size_t topk) {
  TopkLookup out;
  if (taxid.empty() || topk == 0 || posterior_sorted.empty()) {
    return out;
  }
  const std::size_t k = std::min<std::size_t>(topk, posterior_sorted.size());
  for (std::size_t i = 0; i < k; ++i) {
    if (posterior_sorted[i].first == taxid) {
      out.found = true;
      out.rank = i;
      out.prob = posterior_sorted[i].second;
      return out;
    }
  }
  return out;
}

constexpr double kPresenceFixedScale = 1048576.0; // 2^20
constexpr uint32_t kInvalidBucket = std::numeric_limits<uint32_t>::max();

inline uint64_t to_fixed(double value) {
  if (!(value > 0.0)) {
    return 0;
  }
  long double scaled = static_cast<long double>(value) *
                       static_cast<long double>(kPresenceFixedScale);
  if (!(scaled > 0.0)) {
    return 0;
  }
  long double rounded = std::llround(scaled);
  if (!(rounded > 0.0)) {
    return 0;
  }
  long double max_u64 = static_cast<long double>(
      std::numeric_limits<uint64_t>::max());
  if (rounded >= max_u64) {
    return std::numeric_limits<uint64_t>::max();
  }
  return static_cast<uint64_t>(rounded);
}

inline double from_fixed(uint64_t value) {
  return static_cast<double>(value) / kPresenceFixedScale;
}

inline size_t sketch_words(uint32_t bits) {
  if (bits == 0) {
    return 0;
  }
  return static_cast<size_t>((bits + 63u) / 64u);
}

inline void ensure_sketch(PresenceStats &entry, size_t words) {
  if (words == 0) {
    return;
  }
  if (entry.uniqueSketch.size() != words) {
    entry.uniqueSketch.assign(words, 0);
  }
  if (entry.breadthSketch.size() != words) {
    entry.breadthSketch.assign(words, 0);
  }
}

inline void set_sketch_bit(std::vector<uint64_t> &sketch, size_t words,
                           uint32_t bucket) {
  if (words == 0 || bucket == kInvalidBucket) {
    return;
  }
  size_t idx = static_cast<size_t>(bucket >> 6);
  if (idx >= words) {
    return;
  }
  sketch[idx] |= (1ULL << (bucket & 63u));
}

inline uint32_t popcount_vec(const std::vector<uint64_t> &sketch) {
  uint32_t total = 0;
  for (uint64_t word : sketch) {
    total += static_cast<uint32_t>(__builtin_popcountll(word));
  }
  return total;
}
} // namespace

PresenceAccumulator::PresenceAccumulator(uint32_t breadthBits)
    : sketchBits(breadthBits), sketchWords(sketch_words(breadthBits)) {}

PresenceStats &PresenceAccumulator::touch(uint32_t tid) {
  auto [it, inserted] = stats.try_emplace(tid);
  (void)inserted;
  return it->second;
}

void PresenceAccumulator::add_target(uint32_t tid, double hit_weight,
                                     double score_weight, double unique_strength,
                                     bool localUniqueEdge,
                                     uint32_t unique_bucket,
                                     uint32_t breadth_bucket) {
  PresenceStats &entry = touch(tid);
  const uint64_t hit_inc = to_fixed(hit_weight);
  const uint64_t score_inc = to_fixed(score_weight);
  entry.score += score_inc;
  entry.hits += hit_inc;
  const double u = std::clamp(unique_strength, 0.0, 1.0);
  if (u > 0.0) {
    entry.uniqueHits += to_fixed(hit_weight * u);
    entry.uniqueScore += to_fixed(score_weight * u);
  }
  if (localUniqueEdge && sketchWords > 0) {
    ensure_sketch(entry, sketchWords);
    set_sketch_bit(entry.uniqueSketch, sketchWords, unique_bucket);
    set_sketch_bit(entry.breadthSketch, sketchWords, breadth_bucket);
  }
}

void PresenceAccumulator::add_read_support(uint32_t tid, bool uniqueRead) {
  PresenceStats &entry = touch(tid);
  entry.readHits += 1;
  if (uniqueRead) {
    entry.uniqueReads += 1;
  }
}

PresenceSummary::PresenceSummary(uint32_t breadthBits)
    : sketchBits(breadthBits), sketchWords(sketch_words(breadthBits)) {}

void PresenceSummary::merge(const PresenceAccumulator &acc) {
  for (const auto &[tid, entry] : acc.stats) {
    auto &dst = stats[tid];
    dst.score += entry.score;
    dst.uniqueScore += entry.uniqueScore;
    dst.hits += entry.hits;
    dst.uniqueHits += entry.uniqueHits;
    dst.readHits += entry.readHits;
    dst.uniqueReads += entry.uniqueReads;
    if (sketchWords > 0) {
      if (!entry.uniqueSketch.empty()) {
        ensure_sketch(dst, sketchWords);
        size_t copy = std::min(sketchWords, entry.uniqueSketch.size());
        for (size_t i = 0; i < copy; ++i) {
          dst.uniqueSketch[i] |= entry.uniqueSketch[i];
        }
      }
      if (!entry.breadthSketch.empty()) {
        ensure_sketch(dst, sketchWords);
        size_t copy = std::min(sketchWords, entry.breadthSketch.size());
        for (size_t i = 0; i < copy; ++i) {
          dst.breadthSketch[i] |= entry.breadthSketch[i];
        }
      }
    }
  }
}

static PresenceDecision evaluate_presence_coverage_impl(
    const PresenceSummary &summary, const TaxDict &tax,
    const ClassifyConfig &config, const chimera::presence::CoverageMeta &meta,
    size_t totalReads, size_t meanReadLen, bool unique_only) {
  PresenceDecision decision;
  decision.threshold = config.presence_tau;
  if (summary.stats.empty()) {
    return decision;
  }

  robin_hood::unordered_flat_map<std::string, uint64_t> uniqueMap;
  uniqueMap.reserve(meta.entries.size());
  robin_hood::unordered_flat_map<std::string, double> densityMap;
  robin_hood::unordered_flat_map<std::string, double> expectedRefMap;
  robin_hood::unordered_flat_map<std::string, uint64_t> totalMap;
  robin_hood::unordered_flat_map<std::string, uint64_t> genomeMap;
  for (const auto &entry : meta.entries) {
    uniqueMap.emplace(entry.taxid, entry.unique_signatures);
    densityMap.emplace(entry.taxid, entry.unique_density);
    expectedRefMap.emplace(entry.taxid, entry.expected_unique_per_ref_read);
    totalMap.emplace(entry.taxid, entry.total_signatures);
    genomeMap.emplace(entry.taxid, entry.genome_length);
  }

  std::vector<uint64_t> uniqueCounts(tax.id2str.size(), 0);
  std::vector<uint64_t> totalCounts(tax.id2str.size(), 0);
  std::vector<uint64_t> genomeCounts(tax.id2str.size(), 0);
  std::vector<double> densityVec(tax.id2str.size(), 0.0);
  std::vector<double> expectedRefVec(tax.id2str.size(), 0.0);
  for (size_t i = 0; i < tax.id2str.size(); ++i) {
    auto it = uniqueMap.find(tax.id2str[i]);
    if (it != uniqueMap.end()) {
      uniqueCounts[i] = it->second;
    }
    if (auto tit = totalMap.find(tax.id2str[i]); tit != totalMap.end()) {
      totalCounts[i] = tit->second;
    }
    if (auto git = genomeMap.find(tax.id2str[i]); git != genomeMap.end()) {
      genomeCounts[i] = git->second;
    }
    if (auto dit = densityMap.find(tax.id2str[i]); dit != densityMap.end()) {
      densityVec[i] = dit->second;
    }
    if (auto eit = expectedRefMap.find(tax.id2str[i]);
        eit != expectedRefMap.end()) {
      expectedRefVec[i] = eit->second;
    }
  }

  auto resolve_unique = [&](uint32_t tid) -> double {
    if (tid < uniqueCounts.size() && uniqueCounts[tid] > 0) {
      return static_cast<double>(uniqueCounts[tid]);
    }
    return 0.0;
  };
  auto resolve_density = [&](uint32_t tid) -> double {
    if (tid < densityVec.size()) {
      return densityVec[tid];
    }
    return 0.0;
  };
  auto resolve_expected_ref = [&](uint32_t tid) -> double {
    if (tid < expectedRefVec.size()) {
      return expectedRefVec[tid];
    }
    return 0.0;
  };
  auto resolve_total = [&](uint32_t tid) -> double {
    if (tid < totalCounts.size()) {
      return static_cast<double>(totalCounts[tid]);
    }
    return 0.0;
  };
  auto resolve_genome = [&](uint32_t tid) -> double {
    if (tid < genomeCounts.size()) {
      return static_cast<double>(genomeCounts[tid]);
    }
    return 0.0;
  };

  const uint16_t span_used =
      (meta.effective_span > 0) ? meta.effective_span
                                : static_cast<uint16_t>(1);
  const size_t read_len_used =
      (meanReadLen > 0) ? meanReadLen
                        : static_cast<size_t>(meta.ref_read_length);
  const double window_current =
      std::max<int64_t>(1, static_cast<int64_t>(read_len_used) -
                               static_cast<int64_t>(span_used) + 1);
  const double window_ref =
      std::max<int64_t>(1, static_cast<int64_t>(meta.ref_read_length) -
                               static_cast<int64_t>(span_used) + 1);
  auto resolve_exposure = [&](uint32_t tid) -> double {
    double density = resolve_density(tid);
    double expected_per_read = 0.0;
    if (density > 0.0) {
      expected_per_read = density * window_current;
    } else {
      double expected_ref = resolve_expected_ref(tid);
      if (expected_ref > 0.0 && window_ref > 0.0) {
        expected_per_read = expected_ref * (window_current / window_ref);
      } else {
        double total = resolve_total(tid);
        double genome = resolve_genome(tid);
        if (total > 0.0 && genome > 0.0) {
          double density_total = total / genome;
          expected_per_read = density_total * window_current;
        }
      }
    }
    double exposure = expected_per_read * static_cast<double>(totalReads);
    if (exposure <= 0.0) {
      exposure = resolve_unique(tid);
    }
    return exposure;
  };

  double mu = 1e-4;
  double derived_u_min = 0.0;
  {
    std::vector<double> exposures;
    exposures.reserve(summary.stats.size());
    for (const auto &[tid, stats] : summary.stats) {
      double e = resolve_exposure(tid);
      if (e > 0.0) {
        exposures.push_back(e);
      }
    }
    if (!exposures.empty()) {
      std::sort(exposures.begin(), exposures.end());
      size_t idx = static_cast<size_t>(
          std::floor(0.10 * static_cast<double>(exposures.size())));
      if (idx >= exposures.size())
        idx = exposures.size() - 1;
      derived_u_min = exposures[idx];
    }
  }
  double u_min_effective = derived_u_min;

  auto percentile = [](std::vector<double> values, double q) -> double {
    if (values.empty()) {
      return 0.0;
    }
    if (q <= 0.0) {
      return *std::min_element(values.begin(), values.end());
    }
    if (q >= 1.0) {
      return *std::max_element(values.begin(), values.end());
    }
    std::sort(values.begin(), values.end());
    size_t idx = static_cast<size_t>(
        std::floor(q * static_cast<double>(values.size() - 1)));
    if (idx >= values.size()) {
      idx = values.size() - 1;
    }
    return values[idx];
  };

  std::vector<double> score_vals;
  score_vals.reserve(summary.stats.size());
  for (const auto &[tid, stats] : summary.stats) {
    (void)tid;
    score_vals.push_back(from_fixed(stats.score));
  }
  const double score_cut = percentile(score_vals, 0.30);
  const double score_high = percentile(score_vals, 0.90);

  std::vector<double> noise_unique;
  std::vector<double> noise_breadth;
  std::vector<double> noise_rate;
  std::vector<double> noise_read_ratio;
  noise_unique.reserve(summary.stats.size() / 3 + 1);
  noise_breadth.reserve(summary.stats.size() / 3 + 1);
  noise_rate.reserve(summary.stats.size() / 3 + 1);
  noise_read_ratio.reserve(summary.stats.size() / 3 + 1);
  for (const auto &[tid, stats] : summary.stats) {
    (void)tid;
    double score = from_fixed(stats.score);
    if (score > score_cut) {
      continue;
    }
    double unique_obs = 0.0;
    double breadth_obs = 0.0;
    if (!stats.uniqueSketch.empty()) {
      unique_obs = static_cast<double>(popcount_vec(stats.uniqueSketch));
      noise_unique.push_back(unique_obs);
    } else {
      noise_unique.push_back(0.0);
    }
    if (!stats.breadthSketch.empty()) {
      breadth_obs = static_cast<double>(popcount_vec(stats.breadthSketch));
      noise_breadth.push_back(breadth_obs);
    } else {
      noise_breadth.push_back(0.0);
    }
    double breadth_ratio = 0.0;
    if (unique_obs > 0.0) {
      if (breadth_obs > 0.0) {
        breadth_ratio = std::min(1.0, breadth_obs / unique_obs);
      } else {
        breadth_ratio = 1.0;
      }
    }
    double local_strength = unique_obs * breadth_ratio;
    double rate = local_strength / std::max(1.0, score);
    noise_rate.push_back(rate);
    double read_ratio = 0.0;
    if (stats.readHits > 0) {
      read_ratio = static_cast<double>(stats.uniqueReads) /
                   static_cast<double>(stats.readHits);
    }
    noise_read_ratio.push_back(read_ratio);
  }
  const double unique_cut = percentile(noise_unique, 0.95);
  const double breadth_cut = percentile(noise_breadth, 0.95);
  const double rate_cut = percentile(noise_rate, 0.95);
  const double read_ratio_cut = percentile(noise_read_ratio, 0.95);

  if (!(mu > 0.0)) {
    mu = 1e-4;
  }
  mu = std::max(mu, 1e-8);

  double pi = std::clamp(config.presence_pi, 1e-9, 1.0 - 1e-6);
  double logPriorOdds = std::log(pi) - std::log1p(-pi);

  for (const auto &[tid, stats] : summary.stats) {
    double exposure = resolve_exposure(tid);
    if (exposure <= 0.0) {
      exposure = resolve_unique(tid);
    }
    if (exposure <= 0.0) {
      exposure = std::max(1.0, from_fixed(stats.hits));
    }
    double u_eff = std::max<double>(exposure, u_min_effective);
    double unique_hits = from_fixed(stats.uniqueHits);
    const double unique_score = from_fixed(stats.uniqueScore);
    const double score = from_fixed(stats.score);

    double unique_obs = 0.0;
    double breadth_obs = 0.0;
    if (!stats.uniqueSketch.empty()) {
      unique_obs = static_cast<double>(popcount_vec(stats.uniqueSketch));
    }
    if (!stats.breadthSketch.empty()) {
      breadth_obs = static_cast<double>(popcount_vec(stats.breadthSketch));
    }
    double breadth_ratio = 0.0;
    if (unique_obs > 0.0) {
      if (breadth_obs > 0.0) {
        breadth_ratio = std::min(1.0, breadth_obs / unique_obs);
      } else {
        breadth_ratio = 1.0;
      }
    }
    double unique_effective = unique_obs * breadth_ratio;
    if (unique_effective > 0.0) {
      double unique_total = resolve_unique(tid);
      double avg_mult = 1.0;
      if (unique_total > 0.0) {
        avg_mult = exposure / unique_total;
        if (!(avg_mult > 1.0)) {
          avg_mult = 1.0;
        }
      }
      double unique_cap = unique_effective * avg_mult;
      if (unique_cap > 0.0) {
        if (unique_hits > 0.0) {
          unique_hits = std::min(unique_hits, unique_cap);
        } else {
          unique_hits = unique_cap;
        }
      }
    }

    double C = 0.0;
    if (unique_only) {
      C = (unique_hits > 0.0) ? unique_hits
                              : std::max(0.0, unique_score);
    } else {
      C = (unique_hits > 0.0)
              ? unique_hits
              : ((unique_score > 0.0) ? unique_score : score);
    }
    double local_factor = 1.0;
    if (score_high > 0.0 && score < score_high) {
      if (unique_cut > 0.0) {
        local_factor =
            std::min(local_factor, unique_obs / std::max(1e-9, unique_cut));
      }
      if (breadth_cut > 0.0) {
        local_factor = std::min(
            local_factor, breadth_obs / std::max(1e-9, breadth_cut));
      }
      if (rate_cut > 0.0) {
        double rate = unique_effective / std::max(1.0, score);
        local_factor = std::min(local_factor,
                                rate / std::max(1e-12, rate_cut));
      }
      if (read_ratio_cut > 0.0) {
        double read_ratio = 0.0;
        if (stats.readHits > 0) {
          read_ratio = static_cast<double>(stats.uniqueReads) /
                       static_cast<double>(stats.readHits);
        }
        local_factor = std::min(
            local_factor, read_ratio / std::max(1e-12, read_ratio_cut));
      }
      local_factor = std::clamp(local_factor, 0.0, 1.0);
      if (local_factor < 1.0) {
        C *= local_factor;
      }
    }
    double lambda_hat = std::max(0.0, (C / u_eff) - mu);
    double logBF = 0.0;
    if (mu > 0.0) {
      double ratio = (lambda_hat + mu) / mu;
      if (ratio > 0.0 && C > 0.0) {
        logBF = C * std::log(ratio) - lambda_hat * u_eff;
      } else {
        logBF = -lambda_hat * u_eff;
      }
    }
    double logPosterior = logBF + logPriorOdds;
    decision.logPosteriors[tid] = logPosterior;
    double posteriorProb = 0.5;
    if (logPosterior >= 0.0) {
      posteriorProb = 1.0 / (1.0 + std::exp(-logPosterior));
    } else {
      double e = std::exp(logPosterior);
      posteriorProb = e / (1.0 + e);
    }
    decision.posteriors[tid] = posteriorProb;
  }
  return decision;
}

PresenceDecision evaluate_presence_coverage(
    const PresenceSummary &summary, const TaxDict &tax,
    const ClassifyConfig &config, const chimera::presence::CoverageMeta &meta,
    size_t totalReads, size_t meanReadLen) {
  return evaluate_presence_coverage_impl(summary, tax, config, meta, totalReads,
                                        meanReadLen, false);
}

void postEmDecision(
    std::vector<classifyResult> &results, const DecisionConfig &decisionConfig,
    const std::unordered_map<std::string, double> &classWeights,
    const TaxDict &tax, const PresenceDecision *presenceDecision,
    const NcbiTaxdump *ncbiTaxdump, const PostDecisionPolicy &postPolicy) {
  constexpr const char *kUnclassified = "unclassified";

		  using PresenceLevel = PostEmPresenceLevel;

  constexpr double kPresencePiFloor = 1e-6;
  constexpr double kRejectFactor = 2.0;
  const double fallback_strength =
      std::clamp(decisionConfig.fallback_strength, 0.0, 1.0);
  const double selective_strength =
      std::clamp(decisionConfig.selective_reject_strength, 0.0, 1.0);
  const bool selective_reject_fulln_enabled =
      selective_strength > 1e-9 && postPolicy.enable_selective_reject;
  const size_t selective_fulln_min = static_cast<size_t>(
      std::llround(lerp(18.0, 10.0, selective_strength)));
  const size_t selective_pruned_max = static_cast<size_t>(
      std::llround(lerp(4.0, 6.0, selective_strength)));

  double presence_tau = std::numeric_limits<double>::infinity();
  double presence_strict_tau = std::numeric_limits<double>::infinity();
  if (presenceDecision && !presenceDecision->logPosteriors.empty()) {
    presence_tau = presenceDecision->threshold;
    presence_strict_tau = presence_tau;
    constexpr size_t kPresenceCap = 2048;
    std::vector<double> vals;
    vals.reserve(presenceDecision->logPosteriors.size());
    for (const auto &kv : presenceDecision->logPosteriors) {
      vals.push_back(kv.second);
    }
    if (vals.size() > kPresenceCap) {
      auto nth = vals.end() - static_cast<std::ptrdiff_t>(kPresenceCap);
      std::nth_element(vals.begin(), nth, vals.end());
      presence_strict_tau = std::max(presence_strict_tau, *nth);
    }
  }

	  auto presence_level = [&](const std::string &taxid) -> PresenceLevel {
    if (!presenceDecision) {
      return PresenceLevel::kUnknown;
    }
    if (presenceDecision->logPosteriors.empty()) {
      return PresenceLevel::kUnknown;
    }
    if (taxid == kUnclassified) {
      return PresenceLevel::kUnknown;
    }
    auto it = tax.str2id.find(taxid);
    if (it == tax.str2id.end()) {
      return PresenceLevel::kUnknown;
    }
    uint32_t tid = it->second;
    auto it_lp = presenceDecision->logPosteriors.find(tid);
    if (it_lp == presenceDecision->logPosteriors.end()) {
      return PresenceLevel::kUnknown;
    }
	    const double lp = it_lp->second;
	    if (lp >= presence_strict_tau) {
	      return PresenceLevel::kAccepted;
	    }
    // Use symmetric strong-negative rejection to avoid over-penalizing
    // mid-confidence taxa (especially in tail-rich samples).
	    // Accept: lp >= strict_tau (possibly raised by cap).
	    // Reject:  lp <= -tau (strong evidence of absence).
	    // Else:    Unknown (neutral).
		    if (lp <= -presence_tau) {
		      return PresenceLevel::kRejected;
		    }
		    return PresenceLevel::kUnknown;
			  };

			  constexpr size_t kFallbackHintTopK = 16;
			  constexpr double kFallbackGenusConflictGapMax = 0.20;

		  for (auto &result : results) {
		    // Keep the full posterior list (sorted) for dumping/analysis (POST_TOPK),
		    // but use a pruned+renormalized view for final decisions and taxidCount
		    // to avoid exploding long-tail allocations.
		    auto posterior_full = std::move(result.posteriors);
		    if (posterior_full.empty()) {
		      result.taxidCount.clear();
		      result.taxidCount.emplace_back(kUnclassified, 1.0);
		      continue;
		    }
	    std::sort(posterior_full.begin(), posterior_full.end(),
	              [](const auto &a, const auto &b) { return a.second > b.second; });

	    auto pruned = prune_post_em_posterior(
	        posterior_full, decisionConfig.min_class_weight, classWeights,
	        presence_level, kPresencePiFloor, kRejectFactor);
	    std::vector<std::pair<std::string, double>> posterior =
	        std::move(pruned.posterior);

	

	    // Ensure POST_TOPK uses the same (pruned/renormalized) posterior view as the
    // final decision logic. This reduces long-tail noise in head-heavy
	    // profile aggregation (e.g., avoids peaky-rescue picking hitchhiking taxa).
    result.posteriors = posterior;

	    const auto &top = posterior.front();
	    double top_score = top.second;
	    double second_score = (posterior.size() > 1) ? posterior[1].second : 0.0;
	    // Use the *full* (pre-prune) posterior probabilities for thresholding when
	    // possible. The post-EM prune step can renormalize the distribution and
	    // inflate top1/top2 scores (especially when low-weight taxa are dropped),
	    // which may cause overconfident wrong calls. We keep pruning for output
	    // stability, but gate using the unpruned probabilities to stay conservative.
	    auto find_full_prob = [&](const std::string &taxid,
	                              double fallback) -> double {
	      if (taxid.empty()) {
	        return fallback;
	      }
	      // posterior_full is already sorted by prob desc.
	      for (const auto &kv : posterior_full) {
	        if (kv.first == taxid) {
	          return kv.second;
	        }
	      }
	      return fallback;
	    };
	    const double full_top_score = find_full_prob(top.first, top_score);
	    const double full_second_score =
	        (posterior.size() > 1)
	            ? find_full_prob(posterior[1].first, second_score)
	            : 0.0;
	    const double gate_top_score = std::min(top_score, full_top_score);
	    const double gate_second_score = std::min(second_score, full_second_score);

	    PresenceLevel top_presence = presence_level(top.first);

    double class_weight = 0.0;
    bool weight_ok = true;
    if (!classWeights.empty()) {
      auto weight_it = classWeights.find(top.first);
      if (weight_it != classWeights.end()) {
        class_weight = weight_it->second;
        double pi_min = decisionConfig.min_class_weight;
        if (top_presence == PresenceLevel::kAccepted) {
          pi_min = std::min(pi_min, kPresencePiFloor);
	        } else if (top_presence == PresenceLevel::kRejected) {
	          pi_min = std::min(1.0, pi_min * kRejectFactor);
	        }
        weight_ok = (class_weight >= pi_min);
      }
    }
    bool pass = weight_ok;

			    if (pass) {
			      const size_t full_n = posterior_full.size();
			      const size_t pruned_n = posterior.size();
          if (selective_reject_fulln_enabled && pruned.audit.pruned_any &&
              full_n >= selective_fulln_min &&
              pruned_n <= selective_pruned_max) {
			        bool rescued = false;
				        if (!result.best_taxid_hint.empty() &&
				            result.best_taxid_hint != kUnclassified) {
				          const auto hint_hit = lookup_taxid_in_topk(
				              posterior, result.best_taxid_hint, posterior.size());
				          // Only rescue when the hint is the pruned POST_TOPK rank-1 (top2).
				          // Empirically this bucket has much higher precision than "hint==top1"
				          // and avoids leaking high-confidence wrong reads back into output.
				          if (hint_hit.found && hint_hit.rank == 1) {
				            double total_evidence =
				                (result.sample_weight > 0.0) ? result.sample_weight
				                                            : result.evaluated;
				            if (!(total_evidence > 0.0)) {
				              total_evidence = 1.0;
				            }
				            double fallback = static_cast<double>(
				                std::max<double>(1.0, std::llround(total_evidence)));
				            result.taxidCount.clear();
				            result.taxidCount.emplace_back(result.best_taxid_hint,
				                                           fallback);
				            result.reject_reason.clear();
				            rescued = true;
				          }
				        }
		        if (rescued) {
		          continue;
		        }
		        result.taxidCount.clear();
		        result.taxidCount.emplace_back(kUnclassified, 1.0);
		        result.reject_reason = "selective_cert";
		        continue;
		      }

	      result.taxidCount.clear();

      double total_evidence =
          (result.sample_weight > 0.0) ? result.sample_weight : result.evaluated;
      if (!(total_evidence > 0.0)) {
        total_evidence = 1.0;
      }

	      auto eligible = [&](const std::string &taxid, double post) -> bool {
	        if (!(post > 0.0)) {
	          return false;
	        }
	        if (taxid != top.first &&
	            presence_level(taxid) == PresenceLevel::kRejected) {
	          return false;
	        }
	        return true;
	      };
	
	      double sum_adj_total = 0.0;
	      for (const auto &kv : posterior) {
	        const std::string &taxid = kv.first;
	        double post = kv.second;
	        if (!eligible(taxid, post)) {
	          continue;
	        }
        sum_adj_total += post;
	      }
	
	      if (sum_adj_total <= 0.0) {
	        double fallback =
	            static_cast<double>(std::max<double>(1.0, std::llround(total_evidence)));
	        result.taxidCount.emplace_back(top.first, fallback);
	        result.reject_reason.clear();
	        continue;
	      }
	
      for (const auto &kv : posterior) {
        const std::string &taxid = kv.first;
        double post = kv.second;
        if (!eligible(taxid, post)) {
          continue;
        }
        double adj = post;
        double weight = adj / sum_adj_total;
        double frac = weight * total_evidence;
        double count_val = static_cast<double>(std::llround(frac));
        if (count_val > 0.0) {
          result.taxidCount.emplace_back(taxid, count_val);
        }
      }
	
	      if (result.taxidCount.empty()) {
	        double fallback =
	            static_cast<double>(std::max<double>(1.0, std::llround(total_evidence)));
        result.taxidCount.emplace_back(top.first, fallback);
      }

      result.reject_reason.clear();
      continue;
    }

    // keep_multi 逻辑会把丰度拆给多个物种，EM 后期希望赢家通吃，故跳过。
    // double gap = top_score - second_score;
    // bool keep_multi = (second_score > 0.0 && gap < 0.10 &&
    //                    top_score >= 0.35 && second_score >= 0.25);
	    // if (keep_multi) { ... }

        // Continuous fallback: strength=0 keeps reject, strength=1 behaves like
        // previous tail-rich fallback.
        if (fallback_strength > 1e-9) {
          const double gap = top_score - second_score;
          const double base_gap_min =
              std::max(0.0, std::min(1.0, decisionConfig.fallback_gap_min));
          const double gap_min = lerp(1.0, base_gap_min, fallback_strength);
          if (gap < gap_min) {
            // keep rejected
          } else {
			        bool genus_conflict_block = false;
			        if (ncbiTaxdump && ncbiTaxdump->enabled() &&
			            posterior.size() > 1 && gap < kFallbackGenusConflictGapMax) {
			          uint32_t tid1 = 0;
			          uint32_t tid2 = 0;
			          if (chimera::utils::try_parse_u32(top.first, tid1) &&
			              chimera::utils::try_parse_u32(posterior[1].first, tid2)) {
			            const uint32_t g1 = ncbiTaxdump->to_genus(tid1);
			            const uint32_t g2 = ncbiTaxdump->to_genus(tid2);
			            genus_conflict_block =
			                (g1 != 0 && g2 != 0 && g1 != g2);
			          }
			        }
				        if (genus_conflict_block) {
				          // keep rejected
			        } else {
			          std::string fallback_taxid = top.first;
			          bool using_hint = false;
			          if (!result.best_taxid_hint.empty() &&
			              result.best_taxid_hint != kUnclassified) {
			            fallback_taxid = result.best_taxid_hint;
			            using_hint = true;
			          }

		          if (using_hint) {
		            const size_t topk = std::min(
		                kFallbackHintTopK, static_cast<size_t>(posterior.size()));
		            bool hint_in_topk = false;
		            for (size_t i = 0; i < topk; ++i) {
		              if (posterior[i].first == fallback_taxid) {
		                hint_in_topk = true;
		                break;
		              }
		            }
			            if (!hint_in_topk) {
			              // keep rejected
			            } else {
			              double total_evidence =
		                  (result.sample_weight > 0.0) ? result.sample_weight
			                                              : result.evaluated;
			              if (!(total_evidence > 0.0)) {
			                total_evidence = 1.0;
			              }
			              double fallback = static_cast<double>(
			                  std::max<double>(1.0, std::llround(total_evidence)));
				              result.taxidCount.clear();
					              result.taxidCount.emplace_back(fallback_taxid, fallback);
					              result.reject_reason.clear();
					              continue;
				            }
				          } else {
			            double total_evidence = (result.sample_weight > 0.0)
			                                        ? result.sample_weight
			                                        : result.evaluated;
			            if (!(total_evidence > 0.0)) {
			              total_evidence = 1.0;
			            }
			            double fallback = static_cast<double>(
			                std::max<double>(1.0, std::llround(total_evidence)));
				            result.taxidCount.clear();
					            result.taxidCount.emplace_back(fallback_taxid, fallback);
					            result.reject_reason.clear();
					            continue;
				          }
            }
          }
        }

		    result.taxidCount.clear();
	    result.taxidCount.emplace_back(kUnclassified, 1.0);
		    if (result.reject_reason.empty()) {
		      result.reject_reason = weight_ok ? "em_post" : "posterior_weight";
		    }
			  }

				}

} // namespace ChimeraClassify
