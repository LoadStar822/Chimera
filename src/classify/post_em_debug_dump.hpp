// Post-EM debug dump helpers (header-only).
#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace ChimeraClassify {

struct TopkLookup {
  bool found{false};
  std::size_t rank{0};
  double prob{0.0};
};

inline TopkLookup lookup_taxid_in_topk(
    const std::vector<std::pair<std::string, double>> &posterior_sorted,
    const std::string &taxid, std::size_t topk) {
  TopkLookup out;
  if (taxid.empty() || topk == 0 || posterior_sorted.empty()) {
    return out;
  }
  const std::size_t k = std::min<std::size_t>(topk, posterior_sorted.size());
  for (std::size_t i = 0; i < k; ++i) {
    if (posterior_sorted[i].first == taxid) {
      out.found = true;
      out.rank = i;
      out.prob = posterior_sorted[i].second;
      return out;
    }
  }
  return out;
}

} // namespace ChimeraClassify

