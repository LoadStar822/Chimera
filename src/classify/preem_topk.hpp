#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
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

} // namespace ChimeraClassify
