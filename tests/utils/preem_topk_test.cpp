#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "classify/preem_topk.hpp"

namespace {

bool expect_true(const std::string &name, bool got, std::string &message) {
  if (got) {
    return true;
  }
  message = name + " 期望 true, 实际 false";
  return false;
}

} // namespace

int main() {
  int failures = 0;
  std::vector<std::string> failure_messages;

  {
    std::string message;
    std::vector<std::pair<std::string, double>> items = {
        {"a", 2}, {"b", 6}, {"c", 5}, {"d", 4}, {"e", 3}, {"f", 1}};

    ChimeraClassify::normalize_preem_topk(items, 5);

    bool ok = (items.size() == 5);
    if (!expect_true("truncate_to_k", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }

    double max_val = -1.0;
    for (const auto &kv : items) {
      max_val = std::max(max_val, kv.second);
    }
    ok = (!items.empty() && items.front().second == max_val);
    if (!expect_true("front_is_max_after_normalize", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }

    ok = std::is_sorted(items.begin(), items.end(),
                        [](const auto &x, const auto &y) {
                          if (x.second != y.second) {
                            return x.second > y.second;
                          }
                          return x.first < y.first;
                        });
    if (!expect_true("sorted_desc_then_taxid", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> items = {{"b", 1.0},
                                                         {"a", 1.0},
                                                         {"c", 2.0}};
    ChimeraClassify::normalize_preem_topk(items, 16);
    bool ok = (items.size() == 3);
    if (!expect_true("no_expand_when_k_larger", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
    ok = (items.size() >= 3 && items[0].first == "c" && items[1].first == "a" &&
          items[2].first == "b");
    if (!expect_true("stable_tie_break_by_taxid", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> items = {{"x", 1.0}};
    std::vector<std::pair<std::string, double>> ranked = {
        {"a", 10.0}, {"x", 9.0}, {"unclassified", 8.0}, {"b", 7.0}, {"c", 6.0}};

    ChimeraClassify::pad_preem_candidates(items, ranked, 3);
    ChimeraClassify::normalize_preem_topk(items, 3);

    bool ok = (items.size() == 3);
    if (!expect_true("pad_to_target", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
    ok = (items.size() >= 3 && items[0].first == "a" && items[1].first == "b" &&
          items[2].first == "x");
    if (!expect_true("pad_skips_duplicates_and_unclassified", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> items = {{"x", 9.0}};
    std::vector<std::pair<std::string, double>> ranked = {
        {"a", 10.0}, {"x", 9.0}, {"unclassified", 8.0}, {"b", 7.0}, {"c", 0.0}};

    ChimeraClassify::ensure_preem_floor_candidates(items, ranked, 3);
    ChimeraClassify::normalize_preem_topk(items, 3);

    bool ok = (items.size() == 3);
    if (!expect_true("ensure_floor_to_target", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
    }
    ok = (items.size() >= 3 && items[0].first == "a" && items[1].first == "x" &&
          items[2].first == "b");
    if (!expect_true("ensure_floor_uses_ranked_and_keeps_existing", ok, message)) {
      ++failures;
      failure_messages.push_back(message);
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
