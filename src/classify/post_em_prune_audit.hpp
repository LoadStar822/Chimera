// Post-EM posterior pruning + audit helpers (header-only).
#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ChimeraClassify {

enum class PostEmPresenceLevel {
  kUnknown = 0,
  kAccepted = 1,
  kRejected = 2,
};

struct PostEmPruneAudit {
  bool pruning_active{false};
  bool pruned_any{false};
  bool fallback_full{false};
  bool top1_removed{false};
  double dropped_mass{0.0};
  double full_top1_post{0.0};
  double full_gap{0.0};
  double full_top1_weight{0.0};
  double full_top1_local_pi_min{0.0};
  bool full_top1_would_survive{true};
};

struct PostEmPruneResult {
  std::vector<std::pair<std::string, double>> posterior;
  PostEmPruneAudit audit;
};

template <class PresenceLevelFn>
inline PostEmPruneResult prune_post_em_posterior(
    const std::vector<std::pair<std::string, double>> &posterior_full_sorted,
    double min_class_weight,
    const std::unordered_map<std::string, double> &classWeights,
    PresenceLevelFn presence_level, double presence_pi_floor,
    double reject_factor) {
  PostEmPruneResult out;
  out.posterior = posterior_full_sorted;

  if (!posterior_full_sorted.empty()) {
    out.audit.full_top1_post = posterior_full_sorted.front().second;
    out.audit.full_gap =
        posterior_full_sorted.front().second -
        ((posterior_full_sorted.size() > 1) ? posterior_full_sorted[1].second
                                            : 0.0);
  }

  const double pi_prune = std::min(min_class_weight, 1e-4);
  const bool pruning_active =
      (pi_prune > 0.0 && min_class_weight > 0.0 && !classWeights.empty());
  out.audit.pruning_active = pruning_active;
  if (!pruning_active || posterior_full_sorted.empty()) {
    return out;
  }

  // Audit full top1 survival under the prune gate (weight vs local_pi_min).
  {
    const auto &top = posterior_full_sorted.front();
    const PostEmPresenceLevel pres = presence_level(top.first);
    double local_pi_min = pi_prune;
    if (pres == PostEmPresenceLevel::kAccepted) {
      local_pi_min = std::min(pi_prune, presence_pi_floor);
    } else if (pres == PostEmPresenceLevel::kRejected) {
      local_pi_min = pi_prune * reject_factor;
    }
    out.audit.full_top1_local_pi_min = local_pi_min;
    auto it = classWeights.find(top.first);
    out.audit.full_top1_weight = (it != classWeights.end()) ? it->second : 0.0;
    out.audit.full_top1_would_survive =
        (out.audit.full_top1_weight >= local_pi_min);
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
    auto it = classWeights.find(kv.first);
    const double w = (it != classWeights.end()) ? it->second : 0.0;
    if (w >= local_pi_min) {
      pruned.push_back(kv);
      sum += kv.second;
    }
  }

  // Soft behavior: if nothing survives, keep the original (full) posterior.
  if (pruned.empty()) {
    out.audit.fallback_full = true;
    return out;
  }

  out.audit.pruned_any = (pruned.size() < posterior_full_sorted.size());
  out.audit.top1_removed = (pruned.front().first != posterior_full_sorted.front().first);
  if (out.audit.pruned_any) {
    out.audit.dropped_mass = std::max(0.0, 1.0 - sum);
  }
  if (sum > 0.0) {
    for (auto &kv : pruned) {
      kv.second /= sum;
    }
  }
  out.posterior = std::move(pruned);
  return out;
}

} // namespace ChimeraClassify
