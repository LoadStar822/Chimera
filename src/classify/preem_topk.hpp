#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <unordered_set>
#include <vector>

namespace ChimeraClassify {

struct PreemBetaRelaxDecision {
  bool applied{false};
  double beta_local{0.0};
  size_t thr_final_used{0};
};

inline PreemBetaRelaxDecision decide_preem_beta_relax(
    bool use_em, bool low_div_active, bool beta_user, double base_beta,
    double maxEvidence, double eff_eval, double best_ratio, double unique_ratio,
    size_t n_strict, size_t thr_beta, size_t thr_eval, size_t thr_min_eval,
    size_t thr_final_raw, size_t delta, double eff_eval_min) {
  PreemBetaRelaxDecision out;
  out.beta_local = base_beta;
  out.thr_final_used = thr_final_raw;

  if (!use_em) {
    return out;
  }
  if (low_div_active) {
    return out;
  }
  if (beta_user) {
    return out;
  }
  if (!(base_beta > 0.0)) {
    return out;
  }
  if (!(maxEvidence > 0.0)) {
    return out;
  }
  if (!(eff_eval >= eff_eval_min)) {
    return out;
  }
  (void)n_strict;

  const size_t thr_beta_eval_raw = std::min(thr_beta, thr_eval);
  if (thr_final_raw != thr_beta_eval_raw) {
    return out;
  }
  if (thr_beta_eval_raw < thr_min_eval + delta) {
    return out;
  }

  const double u_ratio = std::clamp((3.0 - best_ratio) / 2.0, 0.0, 1.0);
  const double u_uniq =
      std::clamp((0.12 - unique_ratio) / 0.12, 0.0, 1.0);
  const double u = std::max(u_ratio, u_uniq);

  const double beta_local = std::clamp(base_beta * (1.0 - 0.25 * u), 0.0,
                                       base_beta);
  const size_t thr_beta_relaxed = static_cast<size_t>(
      std::floor(beta_local * std::max(0.0, maxEvidence)));
  const size_t thr_beta_eval_relaxed = std::min(thr_beta_relaxed, thr_eval);
  const size_t thr_final_relaxed =
      std::max({thr_beta_eval_relaxed, static_cast<size_t>(1), thr_min_eval});
  const size_t thr_final_used =
      std::min<size_t>(thr_final_raw, thr_final_relaxed);
  if (thr_final_used < thr_final_raw) {
    out.applied = true;
    out.beta_local = beta_local;
    out.thr_final_used = thr_final_used;
  }
  return out;
}

inline void pad_preem_candidates(
    std::vector<std::pair<std::string, double>> &items,
    const std::vector<std::pair<std::string, double>> &ranked, size_t target) {
  if (target == 0 || items.size() >= target || ranked.empty()) {
    return;
  }

  std::unordered_set<std::string> seen;
  seen.reserve(items.size() + 8);
  for (const auto &kv : items) {
    seen.insert(kv.first);
  }

  for (const auto &kv : ranked) {
    if (items.size() >= target) {
      break;
    }
    if (kv.first == "unclassified") {
      continue;
    }
    if (!(kv.second > 0.0)) {
      continue;
    }
    if (!seen.insert(kv.first).second) {
      continue;
    }
    items.push_back(kv);
  }
}

inline void ensure_preem_floor_candidates(
    std::vector<std::pair<std::string, double>> &items,
    const std::vector<std::pair<std::string, double>> &ranked,
    size_t floor_k) {
  pad_preem_candidates(items, ranked, floor_k);
}

struct PreemUnderfullFillResult {
  bool applied{false};
  size_t stage1_added{0};
  size_t stage2_added{0};
};

inline PreemUnderfullFillResult fill_preem_candidates_underfull(
    std::vector<std::pair<std::string, double>> &items,
    const std::vector<std::pair<std::string, double>> &ranked_strict,
    const std::vector<std::pair<std::string, double>> &ranked_loose,
    size_t target, size_t stage2_cap) {
  PreemUnderfullFillResult out;
  if (target == 0 || items.size() >= target) {
    return out;
  }

  const size_t before = items.size();
  pad_preem_candidates(items, ranked_strict, target);
  const size_t after_stage1 = items.size();
  out.stage1_added = after_stage1 - before;

  if (items.size() < target && stage2_cap > 0 && !ranked_loose.empty()) {
    const size_t want = std::min(stage2_cap, ranked_loose.size());
    std::vector<std::pair<std::string, double>> limited;
    limited.reserve(want);
    for (size_t i = 0; i < want; ++i) {
      limited.push_back(ranked_loose[i]);
    }
    pad_preem_candidates(items, limited, target);
  }
  out.stage2_added = items.size() - after_stage1;
  out.applied = (items.size() > before);
  return out;
}

inline bool should_apply_preem_floor(double top1_score, double top2_score,
                                     double min_ratio) {
  if (!(top1_score > 0.0)) {
    return false;
  }
  if (!(top2_score > 0.0)) {
    return false;
  }
  if (!(min_ratio > 0.0)) {
    return false;
  }
  return (top2_score / top1_score) >= min_ratio;
}

enum class KeepaliveResult {
  kSkippedEmpty,
  kSkippedInvalidProtected,
  kSkippedAlreadyPresent,
  kSkippedNoNonUnclassified,
  kSkippedNoTail,
  kBlockedLowAbs,
  kBlockedLowRatio,
  kBlockedLowGain,
  kApplied,
};

inline KeepaliveResult preem_keepalive_replace_tail(
    std::vector<std::pair<std::string, double>> &items,
    const std::pair<std::string, double> &protected_item, double min_ratio,
    double replace_ratio, double abs_min) {
  if (items.empty()) {
    return KeepaliveResult::kSkippedEmpty;
  }
  const std::string &protected_taxid = protected_item.first;
  const double protected_score = protected_item.second;
  if (protected_taxid.empty() || protected_taxid == "unclassified" ||
      !(protected_score > 0.0)) {
    return KeepaliveResult::kSkippedInvalidProtected;
  }

  for (const auto &kv : items) {
    if (kv.first == protected_taxid) {
      return KeepaliveResult::kSkippedAlreadyPresent;
    }
  }

  double top1_score = 0.0;
  size_t top1_idx = items.size();
  bool has_top1 = false;
  for (size_t i = 0; i < items.size(); ++i) {
    const auto &kv = items[i];
    if (kv.first == "unclassified") {
      continue;
    }
    if (!(kv.second > 0.0)) {
      continue;
    }
    top1_score = kv.second;
    top1_idx = i;
    has_top1 = true;
    break;
  }
  if (!has_top1) {
    return KeepaliveResult::kSkippedNoNonUnclassified;
  }

  size_t tail_idx = items.size();
  for (size_t i = items.size(); i-- > 0;) {
    if (items[i].first == "unclassified") {
      continue;
    }
    tail_idx = i;
    break;
  }
  if (tail_idx >= items.size()) {
    return KeepaliveResult::kSkippedNoNonUnclassified;
  }
  if (tail_idx <= top1_idx) {
    return KeepaliveResult::kSkippedNoTail;
  }
  const double tail_score = items[tail_idx].second;

  if (abs_min > 0.0 && protected_score < abs_min) {
    return KeepaliveResult::kBlockedLowAbs;
  }
  if (min_ratio > 0.0 && protected_score < (top1_score * min_ratio)) {
    return KeepaliveResult::kBlockedLowRatio;
  }
  if (replace_ratio > 0.0 && protected_score < (tail_score * replace_ratio)) {
    return KeepaliveResult::kBlockedLowGain;
  }

  items[tail_idx] = protected_item;

  auto cmp = [](const auto &a, const auto &b) {
    if (a.second != b.second) {
      return a.second > b.second;
    }
    return a.first < b.first;
  };
  std::sort(items.begin(), items.end(), cmp);
  return KeepaliveResult::kApplied;
}

} // namespace ChimeraClassify
