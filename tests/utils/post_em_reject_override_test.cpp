#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "classify/ChimeraClassify.hpp"
#include "classify/ChimeraClassifyCommon.hpp"

namespace {

bool expect_true(const std::string &name, bool got, std::string &message) {
  if (got) {
    return true;
  }
  message = name + " 期望 true, 实际 false";
  return false;
}

} // namespace

int main() {
  int failures = 0;
  std::vector<std::string> failure_messages;

  {
    std::string message;
    ChimeraClassify::DecisionConfig dc;
    dc.posterior_threshold = 0.56;
    dc.min_class_weight = 2e-4; // pi_prune=min(2e-4,1e-4)=1e-4
    dc.allow_fallback_on_reject = true;
    dc.fallback_gap_min = 0.10;

    ChimeraClassify::classifyResult r;
    r.id = "r_reject_override";
    r.evaluated = 100.0;
    r.best_taxid_hint = "456"; // typical high-div fallback uses hint
    r.posteriors.emplace_back("123", 0.90);
    r.posteriors.emplace_back("456", 0.10);

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    // 123 will be presence-rejected and has weight slightly below rejected local_pi_min
    // (pi_prune*2), so it gets pruned and the read is mis-called as 456 unless we
    // implement a strong-evidence rejected-top1 override (prune+gate).
    std::unordered_map<std::string, double> classWeights;
    classWeights["123"] = 1.5e-4; // < 2e-4, so pruned when rejected
    classWeights["456"] = 1e-3;

    ChimeraClassify::TaxDict tax;
    tax.str2id["123"] = 123;
    tax.str2id["456"] = 456;

    ChimeraClassify::PresenceDecision pd;
    pd.threshold = 4.6;
    pd.logPosteriors[123] = -5.0; // rejected: lp <= -tau

    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, &pd, nullptr);

    bool ok = (!results[0].taxidCount.empty() &&
               results[0].taxidCount.front().first == "123");
    if (!expect_true("strong_rejected_top1_override_outputs_full_top1", ok, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    ChimeraClassify::DecisionConfig dc;
    dc.posterior_threshold = 0.56;
    dc.min_class_weight = 2e-4; // pi_prune=min(2e-4,1e-4)=1e-4
    dc.allow_fallback_on_reject = true;
    dc.fallback_gap_min = 0.10;

    ChimeraClassify::classifyResult r;
    r.id = "r_reject_override_tier2";
    r.evaluated = 100.0;
    r.best_taxid_hint = "456"; // typical high-div fallback uses hint
    // Tier2 target: full_p1>=0.75 & full_gap>=0.20, but not Tier1.
    r.posteriors.emplace_back("123", 0.77);
    r.posteriors.emplace_back("456", 0.03);

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    // 123 is presence-rejected and has weight between pi_prune(1e-4) and
    // pi_prune*2 (2e-4), so it is pruned by default but should be recovered by
    // the Tier2 override (floor=pi_prune) when the pruned-top1 is very weak in
    // full posterior.
    std::unordered_map<std::string, double> classWeights;
    classWeights["123"] = 1.5e-4; // < 2e-4, pruned when rejected, but >= 1e-4
    classWeights["456"] = 1e-3;

    ChimeraClassify::TaxDict tax;
    tax.str2id["123"] = 123;
    tax.str2id["456"] = 456;

    ChimeraClassify::PresenceDecision pd;
    pd.threshold = 4.6;
    pd.logPosteriors[123] = -5.0; // rejected: lp <= -tau

    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, &pd, nullptr);

    bool ok = (!results[0].taxidCount.empty() &&
               results[0].taxidCount.front().first == "123");
    if (!expect_true("tier2_rejected_top1_override_outputs_full_top1", ok, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
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
