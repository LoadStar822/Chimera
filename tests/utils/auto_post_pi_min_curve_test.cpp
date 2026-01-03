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
    auto tuned = ChimeraClassify::tune_post_pi_min_by_avg_len(
        /*pi_hi=*/2e-4, /*pi_lo=*/1e-4, /*avgLen=*/250,
        /*low_div_active=*/true,
        /*L0=*/2000, /*L1=*/800);
    if (!expect_near("low_div_no_change", tuned.tuned, 2e-4, 1e-15, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    auto tuned = ChimeraClassify::tune_post_pi_min_by_avg_len(
        /*pi_hi=*/2e-4, /*pi_lo=*/1e-4, /*avgLen=*/3000,
        /*low_div_active=*/false,
        /*L0=*/2000, /*L1=*/800);
    if (!expect_near("long_read_no_change", tuned.tuned, 2e-4, 1e-15, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    auto tuned = ChimeraClassify::tune_post_pi_min_by_avg_len(
        /*pi_hi=*/2e-4, /*pi_lo=*/1e-4, /*avgLen=*/600,
        /*low_div_active=*/false,
        /*L0=*/2000, /*L1=*/800);
    if (!expect_near("short_read_hits_pi_lo", tuned.tuned, 1e-4, 1e-15,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    auto tuned = ChimeraClassify::tune_post_pi_min_by_avg_len(
        /*pi_hi=*/2e-4, /*pi_lo=*/1e-4, /*avgLen=*/1400,
        /*low_div_active=*/false,
        /*L0=*/2000, /*L1=*/800);
    const double expected = std::sqrt(2e-4 * 1e-4);
    if (!expect_near("mid_read_log_interp", tuned.tuned, expected, 1e-15,
                     message)) {
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

