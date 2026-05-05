// Centralized automatic classify policy derivation.
#pragma once

#include "ChimeraClassifyCommon.hpp"

#include <string>
#include <unordered_map>

namespace ChimeraClassify {

struct CommunityAutoPolicy {
  double community_dispersion_u{1.0};
  double community_dispersion_s{1.0};
  double first_filter_beta{0.8};
  double em_conf_power{2.0};
};

struct AbundanceAutoPolicy {
  double prune_ratio{0.0};
  double abundance_weight{0.0};
  size_t iterations{0};
};

struct PostPiAutoPolicy {
  double tuned{5e-4};
  double weight_tuned{5e-4};
  double low_weight_mass{0.0};
  double low_weight_target{0.0};
  double dominance{0.0};
  double relax{0.0};
  double evidence_t{0.0};
  double sample_evidence_strength{1.0};
};

double compute_sample_evidence_strength(
    const FileInfo &fileInfo,
    const chimera::presence::CoverageMeta &coverageMeta,
    const CommunityDispersionStats &dispersionStats);

CommunityAutoPolicy
derive_community_auto_policy(const ClassifyConfig &config,
                             const CommunityDispersionStats &dispersionStats);

AbundanceAutoPolicy derive_abundance_auto_policy(const ClassifyConfig &config,
                                                 const EMOptions &options);

PostPiAutoPolicy derive_post_pi_auto_policy(
    const ClassifyConfig &config, const FileInfo &fileInfo,
    const chimera::presence::CoverageMeta &coverageMeta,
    const CommunityDispersionStats &dispersionStats,
    const std::unordered_map<std::string, double> &classWeights);

} // namespace ChimeraClassify
