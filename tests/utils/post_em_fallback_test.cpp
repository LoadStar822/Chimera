#include <iostream>
#include <sstream>
#include <string>
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

bool expect_false(const std::string &name, bool got, std::string &message) {
  if (!got) {
    return true;
  }
  message = name + " 期望 false, 实际 true";
  return false;
}

} // namespace

int main() {
  int failures = 0;
  std::vector<std::string> failure_messages;

  {
    std::string message;
    ChimeraClassify::DecisionConfig dc;
    dc.posterior_threshold = 0.95; // force reject by posterior
    dc.min_class_weight = 0.0;     // disable weight gate
    dc.allow_fallback_on_reject = true;
    dc.fallback_gap_min = 0.10;

    ChimeraClassify::classifyResult r;
    r.id = "r1";
    r.evaluated = 100.0;
    r.posteriors.emplace_back("123", 0.8);
    r.posteriors.emplace_back("456", 0.2);

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    std::unordered_map<std::string, double> classWeights;
    ChimeraClassify::TaxDict tax;
    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, nullptr,
                                    nullptr);

    bool ok = (!results[0].taxidCount.empty() &&
               results[0].taxidCount.front().first == "123");
    if (!expect_true("fallback_keeps_top_taxid_on_reject", ok, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    ChimeraClassify::DecisionConfig dc;
    dc.posterior_threshold = 0.95; // reject by posterior
    dc.min_class_weight = 0.0;
    dc.allow_fallback_on_reject = true;
    dc.fallback_gap_min = 0.10;

    ChimeraClassify::classifyResult r;
    r.id = "r_hint_not_in_topk";
    r.evaluated = 100.0;
    r.best_taxid_hint = "999";
    r.posteriors.emplace_back("123", 0.8);
    r.posteriors.emplace_back("456", 0.2);

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    std::unordered_map<std::string, double> classWeights;
    ChimeraClassify::TaxDict tax;
    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, nullptr,
                                    nullptr);

    // New behavior: if we fallback to hint, require hint to appear in POST_TOPK.
    bool ok = (!results[0].taxidCount.empty() &&
               results[0].taxidCount.front().first == "unclassified");
    if (!expect_true("fallback_blocks_hint_not_in_post_topk", ok, message)) {
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
    r.id = "r_hint_unblock_full_topk";
    r.evaluated = 100.0;
    r.best_taxid_hint = "999";
    // Design the pruned posterior (after dropping 999) to be only slightly above
    // dyn_post so the soft weight gate does not relax pi_min enough to pass.
    // (0.4,0.3)->renorm=(0.571,0.429), gap=0.142
    r.posteriors.emplace_back("123", 0.4);
    r.posteriors.emplace_back("999", 0.3);
    r.posteriors.emplace_back("456", 0.3);

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    // Force reject_by_weight: let top survive prune (w>=pi_prune) but fail weight gate
    // (w<min_class_weight). Hint is present in full posterior topK but pruned away
    // (w<pi_prune), so it is not in pruned POST_TOPK.
    std::unordered_map<std::string, double> classWeights;
    classWeights["123"] = 1.5e-4; // >=pi_prune but <2e-4
    classWeights["456"] = 1e-3;
    classWeights["999"] = 0.0; // pruned away from POST_TOPK

    ChimeraClassify::TaxDict tax;
    tax.str2id["123"] = 123;
    tax.str2id["456"] = 456;
    tax.str2id["999"] = 999;

    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, nullptr,
                                    nullptr);

    bool ok = (!results[0].taxidCount.empty() &&
               results[0].taxidCount.front().first == "999");
    if (!expect_true("fallback_allows_hint_in_full_post_topk_when_strong", ok,
                     message)) {
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
    r.id = "r_hint_unblock_audit_fields";
    r.evaluated = 100.0;
    r.best_taxid_hint = "999";
    r.posteriors.emplace_back("123", 0.4);
    r.posteriors.emplace_back("999", 0.3);
    r.posteriors.emplace_back("456", 0.3);

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    std::unordered_map<std::string, double> classWeights;
    classWeights["123"] = 1.5e-4; // >=pi_prune but <2e-4
    classWeights["456"] = 1e-3;
    classWeights["999"] = 0.0; // pruned away from POST_TOPK

    ChimeraClassify::TaxDict tax;
    tax.str2id["123"] = 123;
    tax.str2id["456"] = 456;
    tax.str2id["999"] = 999;

    std::ostringstream captured;
    auto *old_buf = std::cout.rdbuf(captured.rdbuf());
    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, nullptr,
                                    nullptr);
    std::cout.rdbuf(old_buf);

    const std::string out = captured.str();
    bool ok = (out.find("PostEM fallback gate:") != std::string::npos) &&
              (out.find("log_odds(p50/p90)=") != std::string::npos) &&
              (out.find("log_penalty(p50/p90)=") != std::string::npos);
    if (!expect_true("fallback_gate_prints_margin_breakdown_audit", ok,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    ChimeraClassify::DecisionConfig dc;
    dc.posterior_threshold = 0.95; // reject by posterior -> em_post
    dc.min_class_weight = 0.0;
    dc.allow_fallback_on_reject = true;
    dc.fallback_gap_min = 0.10;

    ChimeraClassify::classifyResult r;
    r.id = "r_em_post_gap_at_threshold";
    r.evaluated = 100.0;
    r.best_taxid_hint = "123"; // in topk
    r.posteriors.emplace_back("123", 0.55);
    r.posteriors.emplace_back("456", 0.45); // gap=0.10

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    std::unordered_map<std::string, double> classWeights;
    ChimeraClassify::TaxDict tax;
    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, nullptr,
                                    nullptr);

    // Desired behavior: em_post uses the same fallback gap threshold.
    bool ok = (!results[0].taxidCount.empty() &&
               results[0].taxidCount.front().first == "123");
    if (!expect_true("fallback_allows_em_post_at_gap_threshold", ok, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    ChimeraClassify::DecisionConfig dc;
    dc.posterior_threshold = 0.95; // force reject
    dc.min_class_weight = 0.0;
    dc.allow_fallback_on_reject = true;
    dc.fallback_gap_min = 0.10;

    ChimeraClassify::NcbiTaxdump taxdump;
    const uint32_t max_id = 300;
    taxdump.parent.resize(static_cast<size_t>(max_id) + 1, 0);
    taxdump.is_species.resize(static_cast<size_t>(max_id) + 1, 0);
    taxdump.is_genus.resize(static_cast<size_t>(max_id) + 1, 0);
    // Two genera: 10 and 20
    taxdump.is_genus[10] = 1;
    taxdump.is_genus[20] = 1;
    // Two species (as candidates) under different genera.
    taxdump.parent[100] = 10;
    taxdump.parent[200] = 20;
    taxdump.is_species[100] = 1;
    taxdump.is_species[200] = 1;

    ChimeraClassify::classifyResult r;
    r.id = "r_genus_conflict_gap_small";
    r.evaluated = 100.0;
    r.posteriors.emplace_back("100", 0.59);
    r.posteriors.emplace_back("200", 0.41); // gap=0.18 (<0.20)

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    std::unordered_map<std::string, double> classWeights;
    ChimeraClassify::TaxDict tax;
    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, nullptr,
                                    &taxdump);

    // Desired behavior: for high-div fallback, block cross-genus ambiguous reads.
    bool ok = (!results[0].taxidCount.empty() &&
               results[0].taxidCount.front().first == "unclassified");
    if (!expect_true("fallback_blocks_genus_conflict_when_gap_small", ok,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    ChimeraClassify::DecisionConfig dc;
    dc.posterior_threshold = 0.95;
    dc.min_class_weight = 0.0;
    dc.allow_fallback_on_reject = false;
    dc.fallback_gap_min = 0.10;

    ChimeraClassify::classifyResult r;
    r.id = "r2";
    r.evaluated = 100.0;
    r.posteriors.emplace_back("123", 0.8);
    r.posteriors.emplace_back("456", 0.2);

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    std::unordered_map<std::string, double> classWeights;
    ChimeraClassify::TaxDict tax;
    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, nullptr,
                                    nullptr);

    bool ok = (!results[0].taxidCount.empty() &&
               results[0].taxidCount.front().first == "unclassified");
    if (!expect_true("no_fallback_rejects_to_unclassified", ok, message)) {
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
