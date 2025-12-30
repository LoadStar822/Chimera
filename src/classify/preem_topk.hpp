#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <unordered_set>
#include <vector>

namespace ChimeraClassify {

inline void normalize_preem_topk(std::vector<std::pair<std::string, double>> &items,
                                 size_t k) {
  if (items.empty() || k == 0) {
    return;
  }

  auto cmp = [](const auto &a, const auto &b) {
    if (a.second != b.second) {
      return a.second > b.second;
    }
    return a.first < b.first;
  };

  if (items.size() > k) {
    std::nth_element(items.begin(), items.begin() + static_cast<std::ptrdiff_t>(k),
                     items.end(), cmp);
    items.resize(k);
  }

  std::sort(items.begin(), items.end(), cmp);
}

inline void pad_preem_candidates(
    std::vector<std::pair<std::string, double>> &items,
    const std::vector<std::pair<std::string, double>> &ranked, size_t target) {
  if (target == 0 || items.size() >= target || ranked.empty()) {
    return;
  }

  std::unordered_set<std::string> seen;
  seen.reserve(items.size() + 8);
  for (const auto &kv : items) {
    seen.insert(kv.first);
  }

  for (const auto &kv : ranked) {
    if (items.size() >= target) {
      break;
    }
    if (kv.first == "unclassified") {
      continue;
    }
    if (!(kv.second > 0.0)) {
      continue;
    }
    if (!seen.insert(kv.first).second) {
      continue;
    }
    items.push_back(kv);
  }
}

} // namespace ChimeraClassify
