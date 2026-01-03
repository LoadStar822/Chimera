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
    double got = ChimeraClassify::tf_saturation_factor(/*count=*/1,
                                                       /*alpha=*/1.0);
    if (!expect_near("tf_sat_c1", got, 1.0, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    double got = ChimeraClassify::tf_saturation_factor(/*count=*/2,
                                                       /*alpha=*/1.0);
    if (!expect_near("tf_sat_c2", got, 0.5, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    double got = ChimeraClassify::tf_saturation_factor(/*count=*/3,
                                                       /*alpha=*/1.0);
    if (!expect_near("tf_sat_c3", got, 1.0 / 3.0, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    // alpha=2 => 1 / (1 + (c-1)^2)
    double got = ChimeraClassify::tf_saturation_factor(/*count=*/3,
                                                       /*alpha=*/2.0);
    if (!expect_near("tf_sat_alpha2", got, 0.2, 1e-12, message)) {
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

