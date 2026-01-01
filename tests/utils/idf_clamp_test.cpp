#include <cmath>
#include <string>
#include <vector>

#include "classify/ChimeraClassifyCommon.hpp"

namespace {

bool expect_near(const std::string &name, double got, double expected,
                 double eps, std::string &message) {
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
    double got = ChimeraClassify::clamp_idf(0.10, /*low_div_active=*/false,
                                           /*idf_max=*/8.0);
    if (!expect_near("idf_highdiv_allows_below_floor", got, 0.10, 1e-12,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    double got = ChimeraClassify::clamp_idf(0.10, /*low_div_active=*/true,
                                           /*idf_max=*/8.0);
    if (!expect_near("idf_lowdiv_keeps_floor_0p5", got, 0.50, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    double got = ChimeraClassify::clamp_idf(100.0, /*low_div_active=*/false,
                                           /*idf_max=*/8.0);
    if (!expect_near("idf_clamps_to_max", got, 8.00, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  if (failures > 0) {
    for (const auto &msg : failure_messages) {
      std::cerr << msg << std::endl;
    }
    return 1;
  }
  return 0;
}

