#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
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
    double got = ChimeraClassify::stopword_tail_factor(
        /*df_est=*/100, /*df_ref=*/10, /*eta=*/0.0, /*min_factor=*/0.1);
    if (!expect_near("stopword_tail_eta0_is_1", got, 1.0, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    double got = ChimeraClassify::stopword_tail_factor(
        /*df_est=*/10, /*df_ref=*/10, /*eta=*/1.0, /*min_factor=*/0.1);
    if (!expect_near("stopword_tail_at_ref_is_1", got, 1.0, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    // stop = (df_est+1)/(df_ref+1) = 4/2 = 2, eta=1 => 0.5
    double got = ChimeraClassify::stopword_tail_factor(
        /*df_est=*/3, /*df_ref=*/1, /*eta=*/1.0, /*min_factor=*/0.1);
    if (!expect_near("stopword_tail_basic", got, 0.5, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    // same as above but clamp to min_factor=0.6
    double got = ChimeraClassify::stopword_tail_factor(
        /*df_est=*/3, /*df_ref=*/1, /*eta=*/1.0, /*min_factor=*/0.6);
    if (!expect_near("stopword_tail_clamps_min", got, 0.6, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    double got = ChimeraClassify::stopword_tail_eta_for_read(
        /*low_div_active=*/true, /*df_high_threshold=*/10, /*readLen=*/10000);
    if (!expect_near("stopword_tail_eta_lowdiv_off", got, 0.0, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    double got = ChimeraClassify::stopword_tail_eta_for_read(
        /*low_div_active=*/false,
        /*df_high_threshold=*/std::numeric_limits<uint32_t>::max(),
        /*readLen=*/10000);
    if (!expect_near("stopword_tail_eta_no_threshold_off", got, 0.0, 1e-12,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    double got = ChimeraClassify::stopword_tail_eta_for_read(
        /*low_div_active=*/false, /*df_high_threshold=*/10, /*readLen=*/2048);
    if (!expect_near("stopword_tail_eta_at_start_is_0", got, 0.0, 1e-12,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    // start=2048 span=4096 -> readLen=4096 gives eta=0.5
    double got = ChimeraClassify::stopword_tail_eta_for_read(
        /*low_div_active=*/false, /*df_high_threshold=*/10, /*readLen=*/4096);
    if (!expect_near("stopword_tail_eta_midpoint", got, 0.5, 1e-12, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    // start=2048 span=4096 -> readLen=6144 gives eta=1.0 (clamped)
    double got = ChimeraClassify::stopword_tail_eta_for_read(
        /*low_div_active=*/false, /*df_high_threshold=*/10, /*readLen=*/6144);
    if (!expect_near("stopword_tail_eta_clamps_to_1", got, 1.0, 1e-12,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    uint32_t got = ChimeraClassify::stopword_tail_df_ref(/*df_high_threshold=*/10,
                                                         /*divisor=*/2);
    if (got != 5u) {
      ++failures;
      failure_messages.push_back("stopword_tail_df_ref_half 期望 5, 实际 " +
                                 std::to_string(got));
    }
  }

  {
    std::string message;
    uint32_t got = ChimeraClassify::stopword_tail_df_ref(/*df_high_threshold=*/1,
                                                         /*divisor=*/2);
    if (got != 1u) {
      ++failures;
      failure_messages.push_back("stopword_tail_df_ref_clamps_to_1 期望 1, 实际 " +
                                 std::to_string(got));
    }
  }

  {
    std::string message;
    uint32_t got = ChimeraClassify::stopword_tail_df_ref(
        /*df_high_threshold=*/std::numeric_limits<uint32_t>::max(),
        /*divisor=*/2);
    if (got != std::numeric_limits<uint32_t>::max()) {
      ++failures;
      failure_messages.push_back(
          "stopword_tail_df_ref_max_preserved 期望 max, 实际 " +
          std::to_string(got));
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
