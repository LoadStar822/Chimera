#include <iostream>
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

    ChimeraClassify::classifyResult r;
    r.id = "r1";
    r.evaluated = 100.0;
    r.posteriors.emplace_back("123", 0.8);
    r.posteriors.emplace_back("456", 0.2);

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    std::unordered_map<std::string, double> classWeights;
    ChimeraClassify::TaxDict tax;
    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, nullptr);

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
    dc.posterior_threshold = 0.95;
    dc.min_class_weight = 0.0;
    dc.allow_fallback_on_reject = false;

    ChimeraClassify::classifyResult r;
    r.id = "r2";
    r.evaluated = 100.0;
    r.posteriors.emplace_back("123", 0.8);
    r.posteriors.emplace_back("456", 0.2);

    std::vector<ChimeraClassify::classifyResult> results;
    results.push_back(r);

    std::unordered_map<std::string, double> classWeights;
    ChimeraClassify::TaxDict tax;
    ChimeraClassify::postEmDecision(results, dc, classWeights, tax, nullptr);

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
