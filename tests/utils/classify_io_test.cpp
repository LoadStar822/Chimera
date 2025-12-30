#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "classify/ChimeraClassifyCommon.hpp"

namespace {

bool expect_contains(const std::string &name, const std::string &haystack,
                     const std::string &needle, std::string &message) {
  if (haystack.find(needle) != std::string::npos) {
    return true;
  }
  message = name + " 期望包含: " + needle + " 实际: " + haystack;
  return false;
}

} // namespace

int main() {
  int failures = 0;
  std::vector<std::string> failure_messages;

  std::filesystem::path out =
      std::filesystem::temp_directory_path() / "chimera_classify_io_test.tsv";
  std::error_code ec;
  std::filesystem::remove(out, ec);

  ChimeraClassify::ClassifyConfig config;
  config.outputFile = out.string();
  config.dump_post_topk = 4;
  config.verbose = false;

  ChimeraClassify::classifyResult r;
  r.id = "r1";
  r.taxidCount.emplace_back("unclassified", 1.0);
  r.reject_reason = "posterior_weight";
  r.best_taxid_hint = "1";
  r.posteriors.emplace_back("1", 0.7);
  r.posteriors.emplace_back("2", 0.3);

  std::vector<ChimeraClassify::classifyResult> results;
  results.push_back(r);
  ChimeraClassify::saveResult(results, config);

  std::ifstream is(out);
  std::string line;
  std::getline(is, line);

  {
    std::string message;
    bool ok = true;
    ok = ok && expect_contains("unclassified_line_has_reject_reason", line,
                              "REJECT=posterior_weight", message);
    ok = ok &&
         expect_contains("unclassified_line_has_hint", line, "HINT=1", message);
    ok = ok && expect_contains("unclassified_line_has_post_topk", line,
                               "\tPOST_TOPK=1:0.700000,2:0.300000", message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  std::filesystem::remove(out, ec);

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
