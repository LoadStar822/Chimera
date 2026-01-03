#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "classify/ChimeraClassifyCommon.hpp"

namespace {

bool expect_near(const std::string &name, double got, double expected, double eps,
                 std::string &message) {
  if (std::fabs(got - expected) <= eps) {
    return true;
  }
  message = name + " 期望 " + std::to_string(expected) + ", 实际 " +
            std::to_string(got);
  return false;
}

} // namespace

int main() {
  int failures = 0;
  std::vector<std::string> failure_messages;

  {
    std::string message;
    auto stats = ChimeraClassify::local_contrast_zero_sum_stats(/*M=*/8,
                                                                /*k=*/1,
                                                                /*delta_total=*/1.0);
    if (!expect_near("net_supported_k1", stats.net_supported, 0.875, 1e-12,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    if (!expect_near("net_unsupported_k1", stats.net_unsupported, -0.125, 1e-12,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    if (!expect_near("l1_sum_k1", stats.l1_sum, 1.75, 1e-12, message)) {
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
