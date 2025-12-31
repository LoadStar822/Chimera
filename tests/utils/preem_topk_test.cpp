#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "classify/preem_topk.hpp"

namespace {

bool expect_true(const std::string &name, bool got, std::string &message) {
  if (got) {
    return true;
  }
  message = name + " 期望 true, 实际 false";
  return false;
}

bool expect_eq_size_t(const std::string &name, size_t got, size_t want,
                      std::string &message) {
  if (got == want) {
    return true;
  }
  message = name + " 期望 " + std::to_string(want) + ", 实际 " +
            std::to_string(got);
  return false;
}

bool expect_near(const std::string &name, double got, double want, double eps,
                 std::string &message) {
  if (std::abs(got - want) <= eps) {
    return true;
  }
  message = name + " 期望 " + std::to_string(want) + "±" + std::to_string(eps) +
            ", 实际 " + std::to_string(got);
  return false;
}

} // namespace

int main() {
  int failures = 0;
  std::vector<std::string> failure_messages;

  {
    std::string message;
    std::vector<std::pair<std::string, double>> items = {
        {"a", 2}, {"b", 6}, {"c", 5}, {"d", 4}, {"e", 3}, {"f", 1}};

    ChimeraClassify::normalize_preem_topk(items, 5);

    bool ok = (items.size() == 5);
    if (!expect_true("truncate_to_k", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }

    double max_val = -1.0;
    for (const auto &kv : items) {
      max_val = std::max(max_val, kv.second);
    }
    ok = (!items.empty() && items.front().second == max_val);
    if (!expect_true("front_is_max_after_normalize", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }

    ok = std::is_sorted(items.begin(), items.end(),
                        [](const auto &x, const auto &y) {
                          if (x.second != y.second) {
                            return x.second > y.second;
                          }
                          return x.first < y.first;
                        });
    if (!expect_true("sorted_desc_then_taxid", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> items = {{"b", 1.0},
                                                         {"a", 1.0},
                                                         {"c", 2.0}};
    ChimeraClassify::normalize_preem_topk(items, 16);
    bool ok = (items.size() == 3);
    if (!expect_true("no_expand_when_k_larger", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
    ok = (items.size() >= 3 && items[0].first == "c" && items[1].first == "a" &&
          items[2].first == "b");
    if (!expect_true("stable_tie_break_by_taxid", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> items = {{"x", 1.0}};
    std::vector<std::pair<std::string, double>> ranked = {
        {"a", 10.0}, {"x", 9.0}, {"unclassified", 8.0}, {"b", 7.0}, {"c", 6.0}};

    ChimeraClassify::pad_preem_candidates(items, ranked, 3);
    ChimeraClassify::normalize_preem_topk(items, 3);

    bool ok = (items.size() == 3);
    if (!expect_true("pad_to_target", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
    ok = (items.size() >= 3 && items[0].first == "a" && items[1].first == "b" &&
          items[2].first == "x");
    if (!expect_true("pad_skips_duplicates_and_unclassified", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> items = {{"x", 9.0}};
    std::vector<std::pair<std::string, double>> ranked = {
        {"a", 10.0}, {"x", 9.0}, {"unclassified", 8.0}, {"b", 7.0}, {"c", 0.0}};

    ChimeraClassify::ensure_preem_floor_candidates(items, ranked, 3);
    ChimeraClassify::normalize_preem_topk(items, 3);

    bool ok = (items.size() == 3);
    if (!expect_true("ensure_floor_to_target", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
    ok = (items.size() >= 3 && items[0].first == "a" && items[1].first == "x" &&
          items[2].first == "b");
    if (!expect_true("ensure_floor_uses_ranked_and_keeps_existing", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    bool ok = !ChimeraClassify::should_apply_preem_floor(100.0, 10.0, 0.2);
    if (!expect_true("floor_gate_blocks_dominant_read", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }

    ok = ChimeraClassify::should_apply_preem_floor(100.0, 20.0, 0.2);
    if (!expect_true("floor_gate_allows_competitive_read", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> items = {
        {"a", 10.0}, {"b", 5.0}, {"c", 1.0}};

    auto res = ChimeraClassify::preem_keepalive_replace_tail(
        items, {"hint", 3.0}, /*min_ratio=*/0.2, /*replace_ratio=*/1.2,
        /*abs_min=*/2.0);

    bool ok = (res == ChimeraClassify::KeepaliveResult::kApplied);
    if (!expect_true("keepalive_applies_when_strong_enough", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }

    ok = (items.size() == 3 && items[0].first == "a" && items[1].first == "b" &&
          items[2].first == "hint");
    if (!expect_true("keepalive_replaces_tail_and_keeps_sort", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    // beta relax: should apply when EM/high-div, thr_final dominated by thr_beta_eval,
    // strict candidates are few, and evidence is not too weak.
    auto decision = ChimeraClassify::decide_preem_beta_relax(
        /*use_em=*/true, /*low_div_active=*/false, /*beta_user=*/false,
        /*base_beta=*/0.45, /*maxEvidence=*/100.0, /*eff_eval=*/100.0,
        /*best_ratio=*/1.0, /*unique_ratio=*/0.05,
        /*base_topk=*/16, /*n_strict=*/3,
        /*thr_beta=*/45, /*thr_eval=*/40, /*thr_min_eval=*/4,
        /*thr_final_raw=*/40, /*delta=*/8, /*eff_eval_min=*/48.0);

    bool ok = decision.applied;
    if (!expect_true("beta_relax_applies", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
    ok = expect_eq_size_t("beta_relax_thr_final_used", decision.thr_final_used,
                          33, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(message);
    }
    ok = expect_near("beta_relax_beta_local", decision.beta_local, 0.3375,
                     1e-6, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    // beta relax: should NOT apply in low-div branch.
    auto decision = ChimeraClassify::decide_preem_beta_relax(
        /*use_em=*/true, /*low_div_active=*/true, /*beta_user=*/false,
        /*base_beta=*/0.45, /*maxEvidence=*/100.0, /*eff_eval=*/100.0,
        /*best_ratio=*/1.0, /*unique_ratio=*/0.05,
        /*base_topk=*/16, /*n_strict=*/3,
        /*thr_beta=*/45, /*thr_eval=*/40, /*thr_min_eval=*/4,
        /*thr_final_raw=*/40, /*delta=*/8, /*eff_eval_min=*/48.0);
    bool ok = !decision.applied && decision.thr_final_used == 40;
    if (!expect_true("beta_relax_skips_lowdiv", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    // binOverflow size suppression: when preem_relax_applied, do not expand dynamicTopK due to size>baseTopK.
    bool overflow_no_suppress = ChimeraClassify::preem_bin_overflow(
        /*fallback_full=*/true, /*topBins_size=*/0, /*taxidCount_size=*/20,
        /*baseTopK=*/16, /*suppress_size_overflow=*/false);
    bool overflow_suppress = ChimeraClassify::preem_bin_overflow(
        /*fallback_full=*/true, /*topBins_size=*/0, /*taxidCount_size=*/20,
        /*baseTopK=*/16, /*suppress_size_overflow=*/true);
    bool ok = overflow_no_suppress && !overflow_suppress;
    if (!expect_true("preem_bin_overflow_suppression", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  if (failures == 0) {
    std::cout << "All tests passed." << std::endl;
    return 0;
  }

  std::cerr << "Failures: " << failures << std::endl;
  for (const auto &msg : failure_messages) {
    std::cerr << "- " << msg << std::endl;
  }
  return 1;
}
