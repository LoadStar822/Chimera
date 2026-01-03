#include <cmath>
#include <cstddef>
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

  // IDF 原始公式 sanity（log2((totalBins+1)/(df+1))）。
  {
    std::string message;
    double got = ChimeraClassify::idf_raw_from_df_bins(/*totalBins=*/100.0,
                                                       /*df_bins=*/10);
    double expected = std::log2(101.0 / 11.0);
    if (!expect_near("idf_raw_from_df_bins_basic", got, expected, 1e-12,
                     message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  // df_for_idf：high-div 下短读允许 df_eff 参与 IDF；长读回退 df_bins。
  {
    const std::size_t df_bins = 8;
    const std::size_t df_eff = 1;
    const std::size_t got = ChimeraClassify::df_for_idf(
        df_bins, df_eff, /*low_div_active=*/false, /*readLen=*/512,
        /*max_len=*/4096);
    if (got != df_eff) {
      ++failures;
      failure_messages.push_back("df_for_idf_short 期望 " +
                                 std::to_string(df_eff) + ", 实际 " +
                                 std::to_string(got));
    }
  }

  {
    const std::size_t df_bins = 8;
    const std::size_t df_eff = 1;
    const std::size_t got = ChimeraClassify::df_for_idf(
        df_bins, df_eff, /*low_div_active=*/false, /*readLen=*/8192,
        /*max_len=*/4096);
    if (got != df_bins) {
      ++failures;
      failure_messages.push_back("df_for_idf_long 期望 " +
                                 std::to_string(df_bins) + ", 实际 " +
                                 std::to_string(got));
    }
  }

  // df_for_idf：low-div 分支强制 df_bins（不引入新变量）。
  {
    const std::size_t df_bins = 8;
    const std::size_t df_eff = 1;
    const std::size_t got = ChimeraClassify::df_for_idf(
        df_bins, df_eff, /*low_div_active=*/true, /*readLen=*/512,
        /*max_len=*/4096);
    if (got != df_bins) {
      ++failures;
      failure_messages.push_back("df_for_idf_lowdiv 期望 " +
                                 std::to_string(df_bins) + ", 实际 " +
                                 std::to_string(got));
    }
  }

  // localUniqueEdge 必须锁回 (deg==1 && df_bins==1)，不能用 df_eff “假装唯一”。
  {
    std::string message;
    bool got = ChimeraClassify::is_local_unique_edge(/*deg_effective=*/1,
                                                     /*df_bins=*/1);
    if (!expect_true("local_unique_edge_true", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    bool got = ChimeraClassify::is_local_unique_edge(/*deg_effective=*/1,
                                                     /*df_bins=*/2);
    if (!expect_false("local_unique_edge_false_df_gt1", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    bool got = ChimeraClassify::is_local_unique_edge(/*deg_effective=*/2,
                                                     /*df_bins=*/1);
    if (!expect_false("local_unique_edge_false_deg_gt1", got, message)) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  // 分片 species (deg=1, df_bins>1) 允许 df_eff=1 参与 bonus 判定，但
  // localUniqueEdge 仍为 false（df_bins 语义不变）。
  {
    std::string message;
    const std::size_t df_bins = 8;
    const std::size_t df_eff = ChimeraClassify::effective_df_bins(
        /*deg_effective=*/1, /*df_bins=*/df_bins, /*low_div_active=*/false);
    if (df_eff != 1) {
      ++failures;
      failure_messages.push_back("effective_df_bins_fragmented 期望 1, 实际 " +
                                 std::to_string(df_eff));
    } else {
      const bool local = ChimeraClassify::is_local_unique_edge(
          /*deg_effective=*/1, /*df_bins=*/df_bins);
      if (!expect_false("local_unique_edge_stays_false_for_fragmented", local,
                        message)) {
        ++failures;
        failure_messages.push_back(std::move(message));
      }
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
