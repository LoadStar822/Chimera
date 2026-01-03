#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "classify/ChimeraClassifyCommon.hpp"

namespace {

bool expect_equal(const std::string &name, std::size_t got, std::size_t expected,
                  std::string &message) {
  if (got == expected) {
    return true;
  }
  message = name + " 期望 " + std::to_string(expected) + ", 实际 " +
            std::to_string(got);
  return false;
}

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
    std::size_t got = ChimeraClassify::effective_df_bins(
        /*deg_effective=*/1, /*df_bins=*/4, /*low_div_active=*/false);
    if (!expect_equal("effective_df_bins_highdiv_demotes_fragmented_unique", got,
                      1, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::size_t got = ChimeraClassify::effective_df_bins(
        /*deg_effective=*/2, /*df_bins=*/4, /*low_div_active=*/false);
    if (!expect_equal("effective_df_bins_keeps_df_when_deg_gt1", got, 4,
                      message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::size_t got = ChimeraClassify::effective_df_bins(
        /*deg_effective=*/1, /*df_bins=*/4, /*low_div_active=*/true);
    if (!expect_equal("effective_df_bins_lowdiv_no_change", got, 4, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    double got = ChimeraClassify::unique_edge_bonus(
        /*base_bonus=*/3.0, /*df_bins=*/1);
    if (!expect_near("unique_edge_bonus_df1_is_full", got, 3.0, 1e-12,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    // 1 + (3-1) / (1 + log2(2)) == 2.0
    double got = ChimeraClassify::unique_edge_bonus(
        /*base_bonus=*/3.0, /*df_bins=*/2);
    if (!expect_near("unique_edge_bonus_df2_scaled", got, 2.0, 1e-12,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    // 1 + 2 / (1 + log2(4)) == 1 + 2/3
    double got = ChimeraClassify::unique_edge_bonus(
        /*base_bonus=*/3.0, /*df_bins=*/4);
    if (!expect_near("unique_edge_bonus_df4_scaled", got, 1.0 + 2.0 / 3.0,
                     1e-12, message)) {
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

