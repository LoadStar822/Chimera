#include <cmath>
#include <iostream>
#include <string>
#include <vector>

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
    std::vector<double> counts{970.0, 10.0, 10.0, 10.0};
    auto stats = ChimeraClassify::compute_low_div_stats(counts, 0, 10);
    bool ok = ChimeraClassify::is_low_diversity(stats, 0.97, 32.0);
    if (!expect_true("low_div_head_dominated", ok, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::vector<double> counts(100, 1.0);
    auto stats = ChimeraClassify::compute_low_div_stats(counts, 0, 10);
    bool ok = ChimeraClassify::is_low_diversity(stats, 0.97, 32.0);
    if (!expect_false("low_div_uniform_high", ok, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    ChimeraClassify::ClassifyConfig config;
    config.preEmTopK = 16;
    config.exclusive_gamma = 1.2;
    config.em_conf_power = 2.0;
    config.em_coexist_penalty = 0.6;
    config.dump_post_topk = 256;
    ChimeraClassify::apply_low_div_overrides(config);
    bool ok = (config.preEmTopK == 256 &&
               config.exclusive_gamma == 0.0 &&
               config.em_conf_power == 1.0 &&
               config.em_coexist_penalty == 0.0 &&
               config.dump_post_topk == 0);
    if (!expect_true("low_div_overrides", ok, message)) {
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
