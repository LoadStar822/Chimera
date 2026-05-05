#include "ChimeraClassifyAutoPolicy.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace ChimeraClassify {

double compute_sample_evidence_strength(
    const FileInfo &fileInfo,
    const chimera::presence::CoverageMeta &coverageMeta,
    const CommunityDispersionStats &dispersionStats) {
  const size_t avg_len =
      (fileInfo.avgLen > 0)
          ? fileInfo.avgLen
          : static_cast<size_t>(coverageMeta.ref_read_length);
  const size_t span =
      (coverageMeta.effective_span > 0)
          ? static_cast<size_t>(coverageMeta.effective_span)
          : static_cast<size_t>(1);
  const size_t ref_len =
      static_cast<size_t>(std::max<uint16_t>(coverageMeta.ref_read_length, 1));
  const double sample_window = static_cast<double>(
      std::max<int64_t>(1, static_cast<int64_t>(avg_len) -
                               static_cast<int64_t>(span) + 1));
  const double ref_window = static_cast<double>(
      std::max<int64_t>(1, static_cast<int64_t>(ref_len) -
                               static_cast<int64_t>(span) + 1));
  double community_complexity = dispersionStats.simpson_species;
  if (!(community_complexity > 0.0)) {
    community_complexity = dispersionStats.eff_species;
  }
  if (!(community_complexity > 0.0)) {
    community_complexity = 1.0;
  }
  const double required_window = community_complexity * ref_window;
  if (!(sample_window > required_window)) {
    return 0.0;
  }
  const double excess = sample_window - required_window;
  const double denom = sample_window + required_window;
  if (!(denom > 0.0)) {
    return 1.0;
  }
  return clamp01(excess / denom);
}

CommunityAutoPolicy
derive_community_auto_policy(const ClassifyConfig &config,
                             const CommunityDispersionStats &dispersionStats) {
  CommunityAutoPolicy policy;
  policy.community_dispersion_u =
      compute_community_dispersion_u(dispersionStats.top_mass,
                                     dispersionStats.eff_species);
  policy.community_dispersion_s =
      compute_community_dispersion_s(dispersionStats.top_mass,
                                     dispersionStats.eff_species);
  policy.first_filter_beta =
      config.firstFilterBeta_user
          ? config.firstFilterBeta
          : lerp(0.5, 0.8, policy.community_dispersion_s);
  policy.em_conf_power = lerp(1.0, 2.0, policy.community_dispersion_s);
  return policy;
}

AbundanceAutoPolicy derive_abundance_auto_policy(const ClassifyConfig &config,
                                                 const EMOptions &options) {
  AbundanceAutoPolicy policy;
  policy.abundance_weight = clamp01(config.community_dispersion_s);
  policy.prune_ratio =
      std::pow(options.prune_ratio, 0.5 + (0.5 * policy.abundance_weight));
  policy.iterations = static_cast<size_t>(
      std::ceil(std::sqrt(static_cast<double>(config.emIter)) *
                policy.abundance_weight));
  return policy;
}

PostPiAutoPolicy derive_post_pi_auto_policy(
    const ClassifyConfig &config, const FileInfo &fileInfo,
    const chimera::presence::CoverageMeta &coverageMeta,
    const CommunityDispersionStats &dispersionStats,
    const std::unordered_map<std::string, double> &classWeights) {
  PostPiAutoPolicy policy;
  policy.tuned = config.post_pi_min;
  policy.weight_tuned = policy.tuned;
  policy.dominance = clamp01(1.0 - config.community_dispersion_s);
  policy.relax = clamp01(config.community_dispersion_s);
  policy.sample_evidence_strength =
      compute_sample_evidence_strength(fileInfo, coverageMeta,
                                       dispersionStats);

  if (policy.tuned > 0.0 && !classWeights.empty()) {
    std::vector<double> weights;
    weights.reserve(classWeights.size());
    double total_mass = 0.0;
    for (const auto &kv : classWeights) {
      const double w = kv.second;
      if (w > 0.0) {
        weights.push_back(w);
        total_mass += w;
      }
    }
    if (total_mass > 0.0 && weights.size() > 1) {
      std::sort(weights.begin(), weights.end(),
                [](double a, double b) { return a > b; });

      auto low_weight_mass = [&](double pi_min) -> double {
        if (!(pi_min > 0.0)) {
          return 0.0;
        }
        double mass = 0.0;
        for (size_t i = weights.size(); i-- > 0;) {
          const double w = weights[i];
          if (w >= pi_min) {
            break;
          }
          mass += w;
        }
        return mass / total_mass;
      };

      const double base = policy.tuned;
      policy.low_weight_mass = low_weight_mass(base);
      double square_sum = 0.0;
      const double inv_total_mass = 1.0 / total_mass;
      for (double w : weights) {
        const double p = w * inv_total_mass;
        square_sum += p * p;
      }
      const double simpson_eff =
          (square_sum > 0.0) ? (1.0 / square_sum) : 1.0;
      policy.low_weight_target = 1.0 / (1.0 + simpson_eff);
      if (policy.low_weight_mass > policy.low_weight_target &&
          policy.relax > 0.0) {
        const double target_mass = policy.low_weight_target * total_mass;
        double tuned_candidate = base;
        double pruned_mass = 0.0;
        for (size_t i = weights.size(); i-- > 0;) {
          const double w = weights[i];
          if (w >= base) {
            break;
          }
          if (pruned_mass + w > target_mass) {
            tuned_candidate = w;
            break;
          }
          pruned_mass += w;
        }
        if (tuned_candidate < base) {
          const double log_base = std::log10(base);
          const double log_tuned = std::log10(tuned_candidate);
          const double log_mix =
              ((1.0 - policy.relax) * log_base) +
              (policy.relax * log_tuned);
          double tuned = std::pow(10.0, log_mix);
          tuned = std::clamp(tuned, tuned_candidate, base);
          policy.tuned = tuned;
        }
      }
    }
  }
  policy.weight_tuned = policy.tuned;

  constexpr double kAutoPostPiMinLoDefault = 1e-4;
  auto evidence_tune = tune_post_pi_min_by_evidence_strength(
      policy.weight_tuned, kAutoPostPiMinLoDefault,
      policy.sample_evidence_strength, policy.relax);
  policy.tuned = evidence_tune.tuned;
  policy.evidence_t = evidence_tune.t;
  return policy;
}

} // namespace ChimeraClassify
