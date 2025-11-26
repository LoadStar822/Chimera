#include "ChimeraClassifyCommon.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

namespace ChimeraClassify {

PresenceAccumulator::PresenceAccumulator(size_t reps) : decoyReps(reps) {}

PresenceStats &PresenceAccumulator::touch(uint32_t tid) {
  auto [it, inserted] = stats.try_emplace(tid);
  if (decoyReps > 0 && it->second.decoys.size() != decoyReps) {
    it->second.decoys.resize(decoyReps, 0.0);
  }
  (void)inserted;
  return it->second;
}

void PresenceAccumulator::add_target(uint32_t tid, double weight,
                                     bool uniqueEdge) {
  PresenceStats &entry = touch(tid);
  entry.score += weight;
  entry.hits += 1.0;
  if (uniqueEdge) {
    entry.uniqueHits += 1;
    entry.uniqueScore += weight;
  }
}

void PresenceAccumulator::add_decoy(size_t rep, uint32_t tid, double weight) {
  if (decoyReps == 0 || rep >= decoyReps) {
    return;
  }
  PresenceStats &entry = touch(tid);
  entry.decoys[rep] += weight;
}

PresenceSummary::PresenceSummary(size_t reps) : decoyReps(reps) {}

void PresenceSummary::merge(const PresenceAccumulator &acc) {
  for (const auto &[tid, entry] : acc.stats) {
    auto &dst = stats[tid];
    dst.score += entry.score;
    dst.uniqueScore += entry.uniqueScore;
    dst.hits += entry.hits;
    dst.uniqueHits += entry.uniqueHits;
    if (decoyReps > 0) {
      if (dst.decoys.size() != decoyReps) {
        dst.decoys.resize(decoyReps, 0.0);
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
      decoyMax =
          *std::max_element(stats.decoys.begin(), stats.decoys.end());
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
    double score = stats.score;
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

PresenceDecision evaluate_presence_coverage(
    const PresenceSummary &summary, const TaxDict &tax,
    const ClassifyConfig &config, const chimera::presence::CoverageMeta &meta,
    size_t totalReads, size_t meanReadLen) {
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
  for (const auto &entry : meta.entries) {
    uniqueMap.emplace(entry.taxid, entry.unique_signatures);
    densityMap.emplace(entry.taxid, entry.unique_density);
    expectedRefMap.emplace(entry.taxid, entry.expected_unique_per_ref_read);
  }

  std::vector<uint64_t> uniqueCounts(tax.id2str.size(), 0);
  std::vector<double> densityVec(tax.id2str.size(), 0.0);
  std::vector<double> expectedRefVec(tax.id2str.size(), 0.0);
  for (size_t i = 0; i < tax.id2str.size(); ++i) {
    auto it = uniqueMap.find(tax.id2str[i]);
    if (it != uniqueMap.end()) {
      uniqueCounts[i] = it->second;
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

  if (!(mu > 0.0)) {
    double muAccum = 0.0;
    double weightSum = 0.0;
    for (const auto &[tid, stats] : summary.stats) {
      double exposure = resolve_exposure(tid);
      if (exposure <= 0.0) {
        continue;
      }
      if (!stats.decoys.empty()) {
        double decoyMean = std::accumulate(stats.decoys.begin(),
                                           stats.decoys.end(), 0.0) /
                           static_cast<double>(stats.decoys.size());
        double mu_j = decoyMean / exposure;
        muAccum += mu_j * exposure;
        weightSum += exposure;
      }
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
      exposure = std::max(1.0, stats.hits);
    }
    double u_eff = std::max<double>(exposure, u_min_effective);
    double C = (stats.uniqueHits > 0)
                   ? static_cast<double>(stats.uniqueHits)
                   : ((stats.uniqueScore > 0.0) ? stats.uniqueScore
                                                : stats.score);
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

void postEmDecision(
    std::vector<classifyResult> &results, const DecisionConfig &decisionConfig,
    const std::unordered_map<std::string, double> &classWeights) {
  constexpr const char *kUnclassified = "unclassified";

  auto format_val = [](double value) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(4) << value;
    return oss.str();
  };

  auto prune_by_global_pi = [&](classifyResult &res, double pi_min) {
    if (pi_min <= 0.0 || res.posteriors.empty() || classWeights.empty()) {
      return;
    }
    std::vector<std::pair<std::string, double>> kept;
    kept.reserve(res.posteriors.size());
    double sum = 0.0;
    for (const auto &kv : res.posteriors) {
      auto weight_it = classWeights.find(kv.first);
      double w = (weight_it != classWeights.end()) ? weight_it->second : 0.0;
      if (w >= pi_min) {
        kept.push_back(kv);
        sum += kv.second;
      }
    }
    if (kept.empty()) {
      res.posteriors.clear();
      return;
    }
    if (sum > 0.0) {
      for (auto &kv : kept) {
        kv.second /= sum;
      }
    }
    res.posteriors.swap(kept);
  };

  for (auto &result : results) {
    if (decisionConfig.min_class_weight > 0.0) {
      prune_by_global_pi(result, decisionConfig.min_class_weight);
    }
    if (result.posteriors.empty()) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(kUnclassified, 1);
      continue;
    }

    auto posterior = result.posteriors;
    std::sort(posterior.begin(), posterior.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    const auto &top = posterior.front();
    double top_score = top.second;
    double second_score = (posterior.size() > 1) ? posterior[1].second : 0.0;

    double dyn_post = decisionConfig.posterior_threshold;
    double evalWeight = result.evaluated;
    if (evalWeight < 24.0) {
      dyn_post = std::max(dyn_post, 0.60);
    } else if (evalWeight < 48.0) {
      dyn_post = std::max(dyn_post, 0.56);
    } else {
      dyn_post = std::max(dyn_post, 0.52);
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
        weight_ok = (class_weight >= decisionConfig.min_class_weight);
      }
    }
    bool pass = weight_ok && (top_score >= dyn_post);

    result.posteriors = std::move(posterior);

    if (pass) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(top.first, 0);
      result.reject_reason.clear();
      continue;
    }

    double gap = top_score - second_score;
    bool keep_multi = (second_score > 0.0 && gap < 0.10 &&
                       top_score >= 0.35 && second_score >= 0.25);
    if (keep_multi) {
      result.taxidCount.clear();
      auto to_count = [&](double p) -> size_t {
        if (result.evaluated <= 0.0) {
          return static_cast<size_t>(std::ceil(p * 10.0));
        }
        double c = std::round(p * result.evaluated);
        return static_cast<size_t>(std::max<double>(1.0, c));
      };
      result.taxidCount.emplace_back(top.first, to_count(top_score));
      result.taxidCount.emplace_back(result.posteriors[1].first,
                                     to_count(result.posteriors[1].second));
      result.reject_reason.clear();
      continue;
    }

    result.taxidCount.clear();
    result.taxidCount.emplace_back(kUnclassified, 1);
    if (result.reject_reason.empty()) {
      result.reject_reason = weight_ok ? "em_post" : "posterior_weight";
    }
  }
}

} // namespace ChimeraClassify
