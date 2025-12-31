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

inline void ensure_preem_floor_candidates(
    std::vector<std::pair<std::string, double>> &items,
    const std::vector<std::pair<std::string, double>> &ranked,
    size_t floor_k) {
  pad_preem_candidates(items, ranked, floor_k);
}

inline bool should_apply_preem_floor(double top1_score, double top2_score,
                                     double min_ratio) {
  if (!(top1_score > 0.0)) {
    return false;
  }
  if (!(top2_score > 0.0)) {
    return false;
  }
  if (!(min_ratio > 0.0)) {
    return false;
  }
  return (top2_score / top1_score) >= min_ratio;
}

enum class KeepaliveResult {
  kSkippedEmpty,
  kSkippedInvalidProtected,
  kSkippedAlreadyPresent,
  kSkippedNoNonUnclassified,
  kSkippedNoTail,
  kBlockedLowAbs,
  kBlockedLowRatio,
  kBlockedLowGain,
  kApplied,
};

inline KeepaliveResult preem_keepalive_replace_tail(
    std::vector<std::pair<std::string, double>> &items,
    const std::pair<std::string, double> &protected_item, double min_ratio,
    double replace_ratio, double abs_min) {
  if (items.empty()) {
    return KeepaliveResult::kSkippedEmpty;
  }
  const std::string &protected_taxid = protected_item.first;
  const double protected_score = protected_item.second;
  if (protected_taxid.empty() || protected_taxid == "unclassified" ||
      !(protected_score > 0.0)) {
    return KeepaliveResult::kSkippedInvalidProtected;
  }

  for (const auto &kv : items) {
    if (kv.first == protected_taxid) {
      return KeepaliveResult::kSkippedAlreadyPresent;
    }
  }

  double top1_score = 0.0;
  size_t top1_idx = items.size();
  bool has_top1 = false;
  for (size_t i = 0; i < items.size(); ++i) {
    const auto &kv = items[i];
    if (kv.first == "unclassified") {
      continue;
    }
    if (!(kv.second > 0.0)) {
      continue;
    }
    top1_score = kv.second;
    top1_idx = i;
    has_top1 = true;
    break;
  }
  if (!has_top1) {
    return KeepaliveResult::kSkippedNoNonUnclassified;
  }

  size_t tail_idx = items.size();
  for (size_t i = items.size(); i-- > 0;) {
    if (items[i].first == "unclassified") {
      continue;
    }
    tail_idx = i;
    break;
  }
  if (tail_idx >= items.size()) {
    return KeepaliveResult::kSkippedNoNonUnclassified;
  }
  if (tail_idx <= top1_idx) {
    return KeepaliveResult::kSkippedNoTail;
  }
  const double tail_score = items[tail_idx].second;

  if (abs_min > 0.0 && protected_score < abs_min) {
    return KeepaliveResult::kBlockedLowAbs;
  }
  if (min_ratio > 0.0 && protected_score < (top1_score * min_ratio)) {
    return KeepaliveResult::kBlockedLowRatio;
  }
  if (replace_ratio > 0.0 && protected_score < (tail_score * replace_ratio)) {
    return KeepaliveResult::kBlockedLowGain;
  }

  items[tail_idx] = protected_item;

  auto cmp = [](const auto &a, const auto &b) {
    if (a.second != b.second) {
      return a.second > b.second;
    }
    return a.first < b.first;
  };
  std::sort(items.begin(), items.end(), cmp);
  return KeepaliveResult::kApplied;
}

} // namespace ChimeraClassify
