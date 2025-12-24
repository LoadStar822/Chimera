#include <cstddef>
#include <cstdint>
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

bool expect_vector_equal(const std::string &name,
                         const std::vector<uint64_t> &got,
                         const std::vector<uint64_t> &expected,
                         std::string &message) {
  if (got == expected) {
    return true;
  }
  message = name + " 结果不一致";
  return false;
}

} // namespace

int main() {
  int failures = 0;
  std::vector<std::string> failure_messages;

  {
    std::string message;
    bool got = ChimeraClassify::allow_unique_edge(1, 1, true, 10, 1);
    if (!expect_false("unique_edge_global_blocks_high_df", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    bool got = ChimeraClassify::allow_unique_edge(1, 1, true, 1, 1);
    if (!expect_true("unique_edge_allows_low_df", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    bool got = ChimeraClassify::allow_unique_edge(1, 1, false, 42, 1);
    if (!expect_true("unique_edge_fallback_without_freq", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    bool got = ChimeraClassify::allow_low_df_boost(2, true, 10, 1);
    if (!expect_false("low_df_boost_blocks_high_df", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    bool got = ChimeraClassify::allow_low_df_boost(2, true, 1, 1);
    if (!expect_true("low_df_boost_allows_low_df", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    bool got = ChimeraClassify::allow_low_df_boost(2, false, 10, 1);
    if (!expect_true("low_df_boost_fallback_without_freq", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    bool got = ChimeraClassify::allow_low_df_boost(3, true, 1, 1);
    if (!expect_false("low_df_boost_off_when_df_bins_large", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::vector<std::pair<uint32_t, uint64_t>> scored = {
        {5, 100}, {1, 200}, {2, 300}, {1, 200}};
    auto got = ChimeraClassify::select_rare_route_values(scored, 2);
    std::vector<uint64_t> expected{200, 300};
    if (!expect_vector_equal("select_rare_route_values_basic", got, expected,
                             message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::vector<std::pair<uint32_t, uint32_t>> ranked = {
        {1, 50}, {2, 25}, {3, 25}};
    double head = ChimeraClassify::compute_head_mass(ranked, 2, 100);
    bool ok = (head > 0.74 && head < 0.76);
    if (!expect_true("compute_head_mass_top2", ok, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    size_t cap = ChimeraClassify::compute_candidate_cap(256, 512, 0.3, 0.5, 600);
    bool ok = (cap == 512);
    if (!expect_true("compute_candidate_cap_expand", ok, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    size_t cap = ChimeraClassify::compute_candidate_cap(256, 512, 0.8, 0.5, 600);
    bool ok = (cap == 256);
    if (!expect_true("compute_candidate_cap_no_expand", ok, message)) {
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
