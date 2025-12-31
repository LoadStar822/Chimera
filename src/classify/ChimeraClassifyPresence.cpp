#include "ChimeraClassifyCommon.hpp"
#include "post_em_prune_audit.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

namespace ChimeraClassify {

namespace {
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

PresenceAccumulator::PresenceAccumulator(size_t reps, uint32_t breadthBits)
    : decoyReps(reps), sketchBits(breadthBits),
      sketchWords(sketch_words(breadthBits)) {}

PresenceStats &PresenceAccumulator::touch(uint32_t tid) {
  auto [it, inserted] = stats.try_emplace(tid);
  if (decoyReps > 0 && it->second.decoys.size() != decoyReps) {
    it->second.decoys.resize(decoyReps, 0);
  }
  (void)inserted;
  return it->second;
}

void PresenceAccumulator::add_target(uint32_t tid, double hit_weight,
                                     double score_weight, bool uniqueEdge,
                                     bool localUniqueEdge,
                                     uint32_t unique_bucket,
                                     uint32_t breadth_bucket) {
  PresenceStats &entry = touch(tid);
  const uint64_t hit_inc = to_fixed(hit_weight);
  const uint64_t score_inc = to_fixed(score_weight);
  entry.score += score_inc;
  entry.hits += hit_inc;
  if (uniqueEdge) {
    entry.uniqueHits += hit_inc;
    entry.uniqueScore += score_inc;
  }
  if (localUniqueEdge && sketchWords > 0) {
    ensure_sketch(entry, sketchWords);
    set_sketch_bit(entry.uniqueSketch, sketchWords, unique_bucket);
    set_sketch_bit(entry.breadthSketch, sketchWords, breadth_bucket);
  }
}

void PresenceAccumulator::add_decoy(size_t rep, uint32_t tid, double weight) {
  if (decoyReps == 0 || rep >= decoyReps) {
    return;
  }
  PresenceStats &entry = touch(tid);
  entry.decoys[rep] += to_fixed(weight);
}

void PresenceAccumulator::add_read_support(uint32_t tid, bool uniqueRead) {
  PresenceStats &entry = touch(tid);
  entry.readHits += 1;
  if (uniqueRead) {
    entry.uniqueReads += 1;
  }
}

PresenceSummary::PresenceSummary(size_t reps, uint32_t breadthBits)
    : decoyReps(reps), sketchBits(breadthBits),
      sketchWords(sketch_words(breadthBits)) {}

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
    if (decoyReps > 0) {
      if (dst.decoys.size() != decoyReps) {
        dst.decoys.resize(decoyReps, 0);
      }
      size_t copy = std::min(decoyReps, entry.decoys.size());
      for (size_t i = 0; i < copy; ++i) {
        dst.decoys[i] += entry.decoys[i];
      }
    }
  }
}

static PresenceDecision evaluate_presence_tdFDR(const PresenceSummary &summary,
                                                const TaxDict &tax,
                                                const ClassifyConfig &config) {
  PresenceDecision decision;
  decision.threshold = config.presence_tau;
  decision.priorPi = config.presence_pi;
  decision.noiseMu = 0.0;
  decision.tested = summary.stats.size();
  if (summary.stats.empty()) {
    return decision;
  }

  std::vector<double> pooled;
  for (const auto &[tid, stats] : summary.stats) {
    double decoyMax = 0.0;
    if (!stats.decoys.empty()) {
      decoyMax = from_fixed(
          *std::max_element(stats.decoys.begin(), stats.decoys.end()));
    }
    pooled.push_back(decoyMax);
  }
  std::sort(pooled.begin(), pooled.end(), std::greater<double>());
  auto decoy_q = [&](double score) -> double {
    size_t positives = 0;
    for (double d : pooled) {
      if (d >= score) {
        ++positives;
      } else {
        break;
      }
    }
    if (pooled.empty()) {
      return 1.0;
    }
    double fdr = static_cast<double>(positives) /
                 static_cast<double>(std::max<size_t>(1, pooled.size()));
    return std::min(1.0, fdr);
  };

  auto prior_odds = std::log(config.presence_pi) -
                    std::log1p(-config.presence_pi);
  double tau = config.presence_tau;
  for (const auto &[tid, stats] : summary.stats) {
    double score = from_fixed(stats.score);
    if (score <= 0.0) {
      continue;
    }
    double q = decoy_q(score);
    decision.qValues[tid] = q;
    if (q <= 1e-9) {
      q = 1e-9;
    }
    if (q >= 1.0) {
      q = 1.0 - 1e-9;
    }
    double logPost = std::log(q) - std::log1p(-q) + prior_odds;
    decision.logPosteriors[tid] = logPost;
    double posterior = 0.5;
    if (logPost >= 0.0) {
      posterior = 1.0 / (1.0 + std::exp(-logPost));
    } else {
      double e = std::exp(logPost);
      posterior = e / (1.0 + e);
    }
    decision.posteriors[tid] = posterior;
    if (posterior >= tau) {
      decision.accepted.insert(tid);
    }
  }
  decision.acceptedCount = decision.accepted.size();
  return decision;
}

PresenceFilterStats apply_presence_filter(const PresenceDecision &decision,
                                          const TaxDict &tax,
                                          std::vector<classifyResult>
                                              &classifyResults,
                                          FileInfo &fileInfo) {
  PresenceFilterStats stats;
  if (decision.accepted.empty()) {
    return stats;
  }
  auto keepTaxid = [&](const std::string &taxid) {
    if (taxid == "unclassified") {
      return true;
    }
    auto it = tax.str2id.find(taxid);
    if (it == tax.str2id.end()) {
      return true;
    }
    return decision.accepted.find(it->second) != decision.accepted.end();
  };
  for (auto &result : classifyResults) {
    if (result.taxidCount.empty()) {
      continue;
    }
    auto before = result.taxidCount.size();
    result.taxidCount.erase(
        std::remove_if(result.taxidCount.begin(), result.taxidCount.end(),
                       [&](const auto &kv) { return !keepTaxid(kv.first); }),
        result.taxidCount.end());
    bool passed = false;
    for (auto &kv : result.taxidCount) {
      auto itid = tax.str2id.find(kv.first);
      if (itid != tax.str2id.end() &&
          decision.accepted.find(itid->second) != decision.accepted.end()) {
        passed = true;
        break;
      }
    }
    if (passed) {
      result.presence_passed = true;
    }
    stats.trimmedAssignments += before - result.taxidCount.size();
    if (result.taxidCount.empty()) {
      result.taxidCount.emplace_back("unclassified", 1);
      ++stats.forcedUnclassified;
      if (result.reject_reason.empty()) {
        result.reject_reason = "presence_filter";
      }
    }
    if (!result.posteriors.empty()) {
      result.posteriors.erase(
          std::remove_if(result.posteriors.begin(), result.posteriors.end(),
                         [&](const auto &kv) { return !keepTaxid(kv.first); }),
          result.posteriors.end());
    }
  }
  for (auto it = fileInfo.uniqueTaxids.begin();
       it != fileInfo.uniqueTaxids.end();) {
    auto mapIt = tax.str2id.find(*it);
    if (mapIt != tax.str2id.end() &&
        decision.accepted.find(mapIt->second) == decision.accepted.end()) {
      it = fileInfo.uniqueTaxids.erase(it);
    } else {
      ++it;
    }
  }
  return stats;
}

static PresenceDecision evaluate_presence_coverage_impl(
    const PresenceSummary &summary, const TaxDict &tax,
    const ClassifyConfig &config, const chimera::presence::CoverageMeta &meta,
    size_t totalReads, size_t meanReadLen, bool unique_only) {
  PresenceDecision decision;
  decision.threshold = config.presence_tau;
  decision.priorPi = config.presence_pi;
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

  double mu = config.presence_noise;
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
  double u_min_effective =
      std::max<double>(config.presence_u_min, derived_u_min);

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

  if (config.verbose) {
    std::cout << "[presence][auto] local_unique_gate="
              << ((unique_cut > 0.0 || breadth_cut > 0.0) ? "on" : "off")
              << " score_q30=" << std::scientific << score_cut
              << " score_q90=" << score_high
              << " unique_q95=" << unique_cut
              << " breadth_q95=" << breadth_cut
              << " rate_q95=" << rate_cut
              << " readratio_q95=" << read_ratio_cut
              << " noise_n=" << noise_unique.size()
              << " total=" << summary.stats.size() << std::defaultfloat
              << std::endl;
  }

  if (!(mu > 0.0)) {
    double muAccum = 0.0;
    double weightSum = 0.0;
    size_t usable_taxa = 0;
    for (const auto &[tid, stats] : summary.stats) {
      double exposure = resolve_exposure(tid);
      if (exposure <= 0.0) {
        exposure = resolve_unique(tid);
      }
      if (exposure <= 0.0) {
        exposure = std::max(1.0, from_fixed(stats.hits));
      }
      if (!stats.decoys.empty()) {
        double decoyMean = 0.0;
        for (uint64_t d : stats.decoys) {
          decoyMean += from_fixed(d);
        }
        decoyMean /= static_cast<double>(stats.decoys.size());
        double mu_j = decoyMean / exposure;
        muAccum += mu_j * exposure;
        weightSum += exposure;
        ++usable_taxa;
      }
    }
    if (config.verbose) {
      std::cout << "presence mu-est: usable_taxa=" << usable_taxa
                << ", weightSum=" << std::scientific << weightSum
                << std::defaultfloat << std::endl;
    }
    if (weightSum > 0.0) {
      mu = muAccum / weightSum;
    }
  }
  if (!(mu > 0.0)) {
    mu = 1e-4;
  }
  mu = std::max(mu, 1e-8);
  decision.noiseMu = mu;

  double pi = std::clamp(config.presence_pi, 1e-9, 1.0 - 1e-6);
  double logPriorOdds = std::log(pi) - std::log1p(-pi);

  decision.tested = summary.stats.size();
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
    decision.lambdaHats[tid] = lambda_hat;
    double posteriorProb = 0.5;
    if (logPosterior >= 0.0) {
      posteriorProb = 1.0 / (1.0 + std::exp(-logPosterior));
    } else {
      double e = std::exp(logPosterior);
      posteriorProb = e / (1.0 + e);
    }
    decision.qValues[tid] = posteriorProb;
    decision.posteriors[tid] = posteriorProb;
    if (logPosterior >= config.presence_tau) {
      decision.accepted.insert(tid);
    }
  }
  decision.acceptedCount = decision.accepted.size();
  return decision;
}

PresenceDecision evaluate_presence_coverage(
    const PresenceSummary &summary, const TaxDict &tax,
    const ClassifyConfig &config, const chimera::presence::CoverageMeta &meta,
    size_t totalReads, size_t meanReadLen) {
  return evaluate_presence_coverage_impl(summary, tax, config, meta, totalReads,
                                        meanReadLen, false);
}

PresenceDecision evaluate_presence_coverage_unique(
    const PresenceSummary &summary, const TaxDict &tax,
    const ClassifyConfig &config, const chimera::presence::CoverageMeta &meta,
    size_t totalReads, size_t meanReadLen) {
  return evaluate_presence_coverage_impl(summary, tax, config, meta, totalReads,
                                        meanReadLen, true);
}

void postEmDecision(
    std::vector<classifyResult> &results, const DecisionConfig &decisionConfig,
    const std::unordered_map<std::string, double> &classWeights,
    const TaxDict &tax, const PresenceDecision *presenceDecision,
    const NcbiTaxdump *ncbiTaxdump) {
  constexpr const char *kUnclassified = "unclassified";
  std::ostream *dump_post = nullptr;
  std::ofstream dump_file;
  const char *dump_path = std::getenv("CHIMERA_DUMP_POSTERIOR");
  if (dump_path && *dump_path) {
    dump_file.open(dump_path);
    if (dump_file.good()) {
      dump_post = &dump_file;
    }
  }

  auto format_val = [](double value) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(4) << value;
    return oss.str();
  };

  using PresenceLevel = PostEmPresenceLevel;

  constexpr double kPresencePiFloor = 1e-6;
  constexpr double kRejectFactor = 2.0;
  constexpr double kRejectDynPostBoost = 0.04;

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
	    // mid-confidence taxa (especially in high-diversity samples).
	    // Accept: lp >= strict_tau (possibly raised by cap).
	    // Reject:  lp <= -tau (strong evidence of absence).
	    // Else:    Unknown (neutral).
	    if (lp <= -presence_tau) {
	      return PresenceLevel::kRejected;
	    }
	    return PresenceLevel::kUnknown;
		  };

	  struct FallbackGateStats {
	    size_t rejected_total{0};
	    size_t rejected_em_post{0};
	    size_t rejected_posterior_weight{0};
	    size_t fallback_applied{0};
	    size_t fallback_blocked_hint_not_in_topk{0};
	    size_t fallback_blocked_genus_conflict{0};
	    size_t fallback_blocked_gap_em_post{0};
	    size_t fallback_blocked_gap_posterior_weight{0};
	    size_t fallback_used_hint{0};
	    std::vector<double> gaps_rejected_em_post;
	    std::vector<double> gaps_rejected_posterior_weight;
	  };

		  FallbackGateStats fb_stats;
		  constexpr size_t kFallbackHintTopK = 16;
		  constexpr double kFallbackGenusConflictGapMax = 0.20;

		  struct PruneAuditStats {
		    size_t reads_total{0};
		    size_t reads_empty_posterior{0};
		    size_t pruning_active{0};
		    size_t pruned_any{0};
		    size_t top1_removed{0};
		    size_t pruned_all_fallback_full{0};
		    size_t top1_removed_w_lt_1e_5{0};
		    std::vector<double> dropped_mass_vals;
		    std::vector<double> top1_removed_fullpost_vals;
		    std::vector<double> top1_removed_gap_vals;
		    double dropped_mass_sum{0.0};
		  };
		  PruneAuditStats prune_stats;

	  auto parse_u32 = [](const std::string &s, uint32_t &out) -> bool {
	    if (s.empty()) {
	      return false;
	    }
	    for (unsigned char c : s) {
	      if (!std::isdigit(c)) {
	        return false;
	      }
	    }
	    try {
	      unsigned long v = std::stoul(s);
	      if (v > std::numeric_limits<uint32_t>::max()) {
	        return false;
	      }
	      out = static_cast<uint32_t>(v);
	      return true;
	    } catch (...) {
	      return false;
	    }
	  };

	  for (auto &result : results) {
	    prune_stats.reads_total += 1;
	    // Keep the full posterior list (sorted) for dumping/analysis (POST_TOPK),
	    // but use a pruned+renormalized view for final decisions and taxidCount
	    // to avoid exploding long-tail allocations.
	    auto posterior_full = std::move(result.posteriors);
	    if (posterior_full.empty()) {
	      prune_stats.reads_empty_posterior += 1;
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

	    if (pruned.audit.pruning_active) {
	      prune_stats.pruning_active += 1;
	      if (pruned.audit.fallback_full) {
	        prune_stats.pruned_all_fallback_full += 1;
	      }
	      if (pruned.audit.pruned_any) {
	        prune_stats.pruned_any += 1;
	        prune_stats.dropped_mass_sum += pruned.audit.dropped_mass;
	        prune_stats.dropped_mass_vals.push_back(pruned.audit.dropped_mass);
	      }
	      if (pruned.audit.top1_removed) {
	        prune_stats.top1_removed += 1;
	        prune_stats.top1_removed_fullpost_vals.push_back(
	            pruned.audit.full_top1_post);
	        prune_stats.top1_removed_gap_vals.push_back(pruned.audit.full_gap);
	        if (pruned.audit.full_top1_weight < 1e-5) {
	          prune_stats.top1_removed_w_lt_1e_5 += 1;
	        }
	      }
	    }

	    // Ensure POST_TOPK uses the same (pruned/renormalized) posterior view as the
	    // final decision logic. This reduces long-tail noise in low-diversity
	    // profile aggregation (e.g., avoids peaky-rescue picking hitchhiking taxa).
    result.posteriors = posterior;

    if (dump_post) {
      (*dump_post) << result.id;
      for (const auto &kv : posterior) {
        (*dump_post) << '\t' << kv.first << ':' << kv.second;
      }
      (*dump_post) << '\n';
    }

    const auto &top = posterior.front();
    double top_score = top.second;
    double second_score = (posterior.size() > 1) ? posterior[1].second : 0.0;

    PresenceLevel top_presence = presence_level(top.first);
    if (top_presence == PresenceLevel::kAccepted) {
      result.presence_passed = true;
    }

    double dyn_post = decisionConfig.posterior_threshold;
    double evalWeight = result.evaluated;
    if (evalWeight < 24.0) {
      dyn_post = std::max(dyn_post, 0.60);
    } else if (evalWeight < 48.0) {
      dyn_post = std::max(dyn_post, 0.56);
    } else {
      dyn_post = std::max(dyn_post, 0.52);
    }
    if (top_presence == PresenceLevel::kRejected) {
      dyn_post = std::max(dyn_post,
                          decisionConfig.posterior_threshold + kRejectDynPostBoost);
    }
    if (result.presence_passed) {
      dyn_post = std::max(0.40, dyn_post - 0.10);
    } else {
      dyn_post = std::max(0.45, dyn_post);
    }
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
        // Soft gate: when the read-level posterior is very confident and well-separated,
        // allow a lower global pi threshold to avoid rejecting low-abundance true taxa.
        // This is intentionally conservative (needs both high posterior and a clear gap).
        if (top_presence != PresenceLevel::kRejected && pi_min > 0.0 &&
            top_score > dyn_post && dyn_post < 0.999) {
          auto clamp01 = [](double x) {
            return std::max(0.0, std::min(1.0, x));
          };
          double soft =
              clamp01((top_score - dyn_post) / (1.0 - dyn_post));
          double gap = top_score - second_score;
          double gap_factor = clamp01(gap / 0.20);
          soft *= gap_factor;
          // emphasize only very confident cases
          soft *= soft;
          // do not relax all the way down to presence floor unless presence already passed
          constexpr double kSoftPiFloor = 1e-5;
          const double floor = std::max(kPresencePiFloor, kSoftPiFloor);
          pi_min = (1.0 - soft) * pi_min + soft * floor;
        }
        weight_ok = (class_weight >= pi_min);
      }
    }
    bool pass = weight_ok && (top_score >= dyn_post);

    if (pass) {
      result.taxidCount.clear();

      double total_evidence =
          (result.sample_weight > 0.0) ? result.sample_weight : result.evaluated;
      if (!(total_evidence > 0.0)) {
        total_evidence = 1.0;
      }

	      double alpha = std::max(1.0, decisionConfig.posterior_power);
	      double head_mass = decisionConfig.posterior_head_mass;
	      if (!(head_mass > 0.0 && head_mass <= 1.0)) {
	        head_mass = 0.95;
	      }
	      size_t max_taxa = static_cast<size_t>(decisionConfig.posterior_max_taxa);
	      if (max_taxa < 1) {
	        max_taxa = 1;
	      }

	      // Winner-take-all for high-confidence reads to suppress long-tail spillover.
	      constexpr double kPosteriorGapWTA = 0.15;
	      if ((top_score - second_score) >= kPosteriorGapWTA) {
	        head_mass = 1.0;
	        max_taxa = 1;
	      }
	
	      auto eligible = [&](const std::string &taxid, double post) -> bool {
	        if (!(post > 0.0)) {
	          return false;
	        }
	        if (post < decisionConfig.posterior_min_fraction) {
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
	        sum_adj_total += std::pow(post, alpha);
	      }
	
	      if (sum_adj_total <= 0.0) {
	        double fallback =
	            static_cast<double>(std::max<double>(1.0, std::llround(total_evidence)));
	        result.taxidCount.emplace_back(top.first, fallback);
	        result.reject_reason.clear();
	        continue;
	      }
	
	      // Per-read posterior head truncation:
	      // keep only the top taxa covering `head_mass` of sum(post^alpha),
	      // and at most `max_taxa` taxa, to suppress long-tail non-zero outputs.
	      double sum_adj_kept = 0.0;
	      double cum_adj = 0.0;
	      size_t kept = 0;
	      for (const auto &kv : posterior) {
	        const std::string &taxid = kv.first;
	        double post = kv.second;
	        if (!eligible(taxid, post)) {
	          continue;
	        }
	        double adj = std::pow(post, alpha);
	        cum_adj += adj;
	        sum_adj_kept += adj;
	        kept += 1;
	        if (kept >= max_taxa || cum_adj >= head_mass * sum_adj_total) {
	          break;
	        }
	      }
	      if (sum_adj_kept <= 0.0) {
	        double fallback =
	            static_cast<double>(std::max<double>(1.0, std::llround(total_evidence)));
	        result.taxidCount.emplace_back(top.first, fallback);
	        result.reject_reason.clear();
	        continue;
	      }
	
	      cum_adj = 0.0;
	      kept = 0;
	      for (const auto &kv : posterior) {
	        const std::string &taxid = kv.first;
	        double post = kv.second;
	        if (!eligible(taxid, post)) {
	          continue;
	        }
	        double adj = std::pow(post, alpha);
	        cum_adj += adj;
	        kept += 1;
	
	        double weight = adj / sum_adj_kept;
	        double frac = weight * total_evidence;
	        double count_val = static_cast<double>(std::llround(frac));
	        if (count_val > 0.0) {
	          result.taxidCount.emplace_back(taxid, count_val);
	        }
	
	        if (kept >= max_taxa || cum_adj >= head_mass * sum_adj_total) {
	          break;
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

		    // High-diversity: allow a conservative fallback on reject to reduce
		    // unclassified reads, but only when the posterior is clearly separated
		    // (large top1-top2 gap). This avoids turning ambiguous reads into FP+FN.
		    if (decisionConfig.allow_fallback_on_reject) {
			      const double gap = top_score - second_score;
			      const double base_gap_min =
			          std::max(0.0, std::min(1.0, decisionConfig.fallback_gap_min));
			      const bool reject_by_weight = !weight_ok;
			      const bool reject_is_em_post = !reject_by_weight;
			      const double gap_min = base_gap_min;

		      fb_stats.rejected_total += 1;
		      if (reject_is_em_post) {
		        fb_stats.rejected_em_post += 1;
		        fb_stats.gaps_rejected_em_post.push_back(gap);
		      } else {
		        fb_stats.rejected_posterior_weight += 1;
		        fb_stats.gaps_rejected_posterior_weight.push_back(gap);
		      }

			      if (gap < gap_min) {
			        if (reject_is_em_post) {
			          fb_stats.fallback_blocked_gap_em_post += 1;
			        } else {
			          fb_stats.fallback_blocked_gap_posterior_weight += 1;
			        }
			        // keep rejected
			      } else {
			        bool genus_conflict_block = false;
			        if (ncbiTaxdump && ncbiTaxdump->enabled() &&
			            posterior.size() > 1 && gap < kFallbackGenusConflictGapMax) {
			          uint32_t tid1 = 0;
			          uint32_t tid2 = 0;
			          if (parse_u32(top.first, tid1) &&
			              parse_u32(posterior[1].first, tid2)) {
			            const uint32_t g1 = ncbiTaxdump->to_genus(tid1);
			            const uint32_t g2 = ncbiTaxdump->to_genus(tid2);
			            genus_conflict_block =
			                (g1 != 0 && g2 != 0 && g1 != g2);
			          }
			        }
			        if (genus_conflict_block) {
			          fb_stats.fallback_blocked_genus_conflict += 1;
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
			              fb_stats.fallback_blocked_hint_not_in_topk += 1;
			              // keep rejected
			            } else {
			              fb_stats.fallback_used_hint += 1;
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
			              fb_stats.fallback_applied += 1;
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
			            fb_stats.fallback_applied += 1;
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

	  if (decisionConfig.allow_fallback_on_reject && fb_stats.rejected_total > 0) {
    auto q = [](std::vector<double> vals, double p) -> double {
      if (vals.empty()) {
        return 0.0;
      }
      p = std::max(0.0, std::min(1.0, p));
      const size_t idx =
          static_cast<size_t>(std::floor(p * static_cast<double>(vals.size() - 1)));
      std::nth_element(vals.begin(),
                       vals.begin() + static_cast<std::ptrdiff_t>(idx), vals.end());
      return vals[idx];
    };

    std::cout << "PostEM fallback gate: rejected=" << fb_stats.rejected_total
              << " (em_post=" << fb_stats.rejected_em_post
              << ", posterior_weight=" << fb_stats.rejected_posterior_weight
              << "), applied=" << fb_stats.fallback_applied
	              << ", used_hint=" << fb_stats.fallback_used_hint
	              << ", blocked_hint_not_in_topk="
	              << fb_stats.fallback_blocked_hint_not_in_topk
	              << ", blocked_genus_conflict="
	              << fb_stats.fallback_blocked_genus_conflict
	              << ", blocked_gap_em_post=" << fb_stats.fallback_blocked_gap_em_post
	              << ", blocked_gap_posterior_weight="
	              << fb_stats.fallback_blocked_gap_posterior_weight
	              << ", gap_em_post(p50/p90)="
              << q(fb_stats.gaps_rejected_em_post, 0.50) << '/'
              << q(fb_stats.gaps_rejected_em_post, 0.90)
              << ", gap_posterior_weight(p50/p90)="
              << q(fb_stats.gaps_rejected_posterior_weight, 0.50) << '/'
              << q(fb_stats.gaps_rejected_posterior_weight, 0.90) << std::endl;
	  }

	  const bool prune_active_sample =
	      (decisionConfig.min_class_weight > 0.0 && !classWeights.empty());
	  if (prune_active_sample && prune_stats.reads_total > 0) {
	    auto q = [](std::vector<double> vals, double p) -> double {
	      if (vals.empty()) {
	        return 0.0;
	      }
	      p = std::max(0.0, std::min(1.0, p));
	      const size_t idx = static_cast<size_t>(
	          std::floor(p * static_cast<double>(vals.size() - 1)));
	      std::nth_element(vals.begin(),
	                       vals.begin() + static_cast<std::ptrdiff_t>(idx),
	                       vals.end());
	      return vals[idx];
	    };

	    double dropped_mean = 0.0;
	    if (prune_stats.pruned_any > 0) {
	      dropped_mean =
	          prune_stats.dropped_mass_sum / static_cast<double>(prune_stats.pruned_any);
	    }

	    std::cout << "PostEM prune audit: reads=" << prune_stats.reads_total
	              << ", empty=" << prune_stats.reads_empty_posterior
	              << ", pruned_any=" << prune_stats.pruned_any
	              << ", top1_removed=" << prune_stats.top1_removed
	              << ", pruned_all_fallback_full="
	              << prune_stats.pruned_all_fallback_full
	              << ", dropped_mass(mean/p50/p90)=" << format_val(dropped_mean)
	              << '/' << format_val(q(prune_stats.dropped_mass_vals, 0.50))
	              << '/' << format_val(q(prune_stats.dropped_mass_vals, 0.90))
	              << ", top1_removed_fullpost(p50/p90)="
	              << format_val(q(prune_stats.top1_removed_fullpost_vals, 0.50))
	              << '/' << format_val(q(prune_stats.top1_removed_fullpost_vals, 0.90))
	              << ", top1_removed_gap(p50/p90)="
	              << format_val(q(prune_stats.top1_removed_gap_vals, 0.50)) << '/'
	              << format_val(q(prune_stats.top1_removed_gap_vals, 0.90))
	              << ", top1_removed_w<1e-5=" << prune_stats.top1_removed_w_lt_1e_5
	              << std::endl;
	  }

				}

} // namespace ChimeraClassify
