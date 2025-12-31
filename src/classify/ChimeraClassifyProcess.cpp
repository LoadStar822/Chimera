#include "ChimeraClassifyCommon.hpp"
#include "preem_topk.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <mutex>
#include <numeric>
#include <thread>
#include <utility>

namespace {

std::ostream *preem_dump_stream() {
  static std::once_flag once;
  static std::mutex mu;
  static std::ofstream ofs;
  static bool enabled = false;
  const char *path = std::getenv("CHIMERA_DUMP_PREEM");
  if (!path || !*path) {
    return nullptr;
  }
  std::string path_copy = path;
  std::call_once(once, [&]() {
    ofs.open(path_copy);
    enabled = ofs.good();
  });
  if (!enabled) {
    return nullptr;
  }
  return &ofs;
}

void dump_preem_line(const std::string &id,
                     const std::vector<std::pair<std::string, double>> &items,
                     bool high_conf_pre) {
  auto *os = preem_dump_stream();
  if (!os) {
    return;
  }
  static std::mutex mu;
  std::lock_guard<std::mutex> lock(mu);
  (*os) << id;
  (*os) << '\t' << "highconf=" << (high_conf_pre ? 1 : 0);
  for (const auto &kv : items) {
    (*os) << '\t' << kv.first << ':' << kv.second;
  }
  (*os) << '\n';
}

std::ostream *preem_stats_stream() {
  static std::once_flag once;
  static std::mutex mu;
  static std::ofstream ofs;
  static bool enabled = false;
  const char *path = std::getenv("CHIMERA_DUMP_PREEM_STATS");
  if (!path || !*path) {
    return nullptr;
  }
  std::string path_copy = path;
  std::call_once(once, [&]() {
    ofs.open(path_copy);
    enabled = ofs.good();
  });
  if (!enabled) {
    return nullptr;
  }
  return &ofs;
}

void dump_preem_stats(const std::string &id, bool high_conf_pre, bool use_em,
                      double eff_eval, double best, double second,
                      double best_ratio, double gap, size_t uniqueCount,
                      double uniqueRatio, size_t thr_beta, size_t thr_eval,
                      size_t thr_min_eval, size_t thr_final) {
  auto *os = preem_stats_stream();
  if (!os) {
    return;
  }
  static std::mutex mu;
  std::lock_guard<std::mutex> lock(mu);
  (*os) << id
        << "\thighconf=" << (high_conf_pre ? 1 : 0)
        << "\tuse_em=" << (use_em ? 1 : 0)
        << "\teff_eval=" << eff_eval
        << "\tbest=" << best
        << "\tsecond=" << second
        << "\tbest_ratio=" << best_ratio
        << "\tgap=" << gap
        << "\tuniq_cnt=" << uniqueCount
        << "\tuniq_ratio=" << uniqueRatio
        << "\tthr_beta=" << thr_beta
        << "\tthr_eval=" << thr_eval
        << "\tthr_min_eval=" << thr_min_eval
        << "\tthr_final=" << thr_final
        << '\n';
}

std::ostream *tidscore_dump_stream() {
  static std::once_flag once;
  static std::mutex mu;
  static std::ofstream ofs;
  static bool enabled = false;
  const char *path = std::getenv("CHIMERA_DUMP_TIDSCORE");
  if (!path || !*path) {
    return nullptr;
  }
  std::string path_copy = path;
  std::call_once(once, [&]() {
    ofs.open(path_copy);
    enabled = ofs.good();
  });
  if (!enabled) {
    return nullptr;
  }
  return &ofs;
}

size_t tidscore_dump_topn() {
  const char *val = std::getenv("CHIMERA_DUMP_TIDSCORE_TOPN");
  if (!val || !*val) {
    return 128;
  }
  try {
    long long v = std::stoll(val);
    if (v <= 0) {
      return 0;
    }
    return static_cast<size_t>(v);
  } catch (...) {
    return 128;
  }
}

void dump_tidscore_line(const std::string &id,
                        const ChimeraClassify::TaxDict &tax,
                        const robin_hood::unordered_flat_map<uint32_t, double> &tidScore) {
  auto *os = tidscore_dump_stream();
  if (!os || tidScore.empty()) {
    return;
  }
  const size_t topn = tidscore_dump_topn();
  if (topn == 0) {
    return;
  }
  std::vector<std::pair<uint32_t, double>> items;
  items.reserve(tidScore.size());
  for (const auto &kv : tidScore) {
    items.emplace_back(kv.first, kv.second);
  }
  if (items.size() > topn) {
    std::nth_element(items.begin(), items.begin() + topn, items.end(),
                     [](const auto &a, const auto &b) { return a.second > b.second; });
    items.resize(topn);
  }
  std::sort(items.begin(), items.end(),
            [](const auto &a, const auto &b) { return a.second > b.second; });

  static std::mutex mu;
  std::lock_guard<std::mutex> lock(mu);
  (*os) << id;
  for (const auto &kv : items) {
    uint32_t tid_id = kv.first;
    if (tid_id >= tax.id2str.size()) {
      continue;
    }
    (*os) << '\t' << tax.id2str[tid_id] << ':' << kv.second;
  }
  (*os) << '\n';
}

} // namespace

namespace ChimeraClassify {

void GroupHeat::ensure(size_t bins) {
  if (score.size() < bins) {
    score.resize(bins, 0);
  }
}

void GroupHeat::decay_if_needed() {
  if (decay_period == 0) {
    return;
  }
  ++counter;
  if (counter >= decay_period) {
    counter = 0;
    for (auto &v : score) {
      v -= (v >> decay_shift);
    }
  }
}

void GroupHeat::boost(uint32_t bin, uint32_t delta) {
  if (bin >= score.size()) {
    score.resize(static_cast<size_t>(bin) + 1, 0);
  }
  uint64_t next =
      static_cast<uint64_t>(score[bin]) + static_cast<uint64_t>(delta);
  score[bin] = static_cast<uint32_t>(
      std::min<uint64_t>(next, std::numeric_limits<uint32_t>::max()));
}

static inline uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  x ^= (x >> 31);
  return x;
}

void processSequence(
    const std::vector<uint64_t> &hashs1, size_t readLen,
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    ClassifyConfig &config, const WeightingContext &weightCtx, GroupHeat &heat,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const std::string &id,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    PresenceAccumulator *presenceAcc, uint64_t decoySeed) {
  size_t hashNum = hashs1.size();
  auto xor_reduce = [](const std::vector<uint64_t> &vals) {
    uint64_t acc = 0;
    for (uint64_t v : vals) {
      acc ^= v;
    }
    return acc;
  };

  const size_t binNumAll = indexToTaxid.size();
  heat.ensure(binNumAll);
  const bool presenceEnabled = (presenceAcc != nullptr);
  const size_t presenceDecoyReps =
      presenceEnabled ? presenceAcc->decoyReps : 0;
  const uint32_t presenceSketchBits =
      presenceEnabled ? presenceAcc->sketchBits : 0;
  const bool presenceSketchPow2 =
      (presenceSketchBits > 0) &&
      ((presenceSketchBits & (presenceSketchBits - 1u)) == 0);
  const uint32_t presenceSketchMask =
      presenceSketchPow2 ? (presenceSketchBits - 1u) : 0u;
  auto sketch_index = [&](uint64_t value) -> uint32_t {
    if (presenceSketchBits == 0) {
      return std::numeric_limits<uint32_t>::max();
    }
    uint64_t mix = splitmix64(value);
    if (presenceSketchPow2) {
      return static_cast<uint32_t>(mix & presenceSketchMask);
    }
    return static_cast<uint32_t>(mix % presenceSketchBits);
  };
  const double exclusiveGamma = std::max(0.0, config.exclusive_gamma);
  constexpr double kUniqueEdgeBonus = 3.0;
  const bool has_sample_weight = weightCtx.has_sample_weights();
  double sample_weight = 1.0;
  if (has_sample_weight && weightCtx.sampleWeights) {
    auto it_w = weightCtx.sampleWeights->find(id);
    if (it_w != weightCtx.sampleWeights->end()) {
      sample_weight = it_w->second;
    }
  }
  if (!(sample_weight > 0.0)) {
    sample_weight = 1.0;
  }
  // Presence model is sensitive to extremely large sample weights (read counts).
  // Use a tempered weight so presence stays a useful gate instead of saturating.
  double presence_weight = sample_weight;
  if (has_sample_weight) {
    presence_weight = std::sqrt(sample_weight);
    if (!(presence_weight > 0.0)) {
      presence_weight = 1.0;
    }
  }

  auto note_reject = [&](const std::string &reason,
                         const std::string &taxid_hint) {
    fileInfo.rejectReasons[reason] += 1;
    if (!taxid_hint.empty()) {
      fileInfo.rejectByTaxid[taxid_hint][reason] += 1;
    }
  };

  // Sample a subset of hashes for coarse candidate selection.
  size_t targetSample = hashNum / 2;
  targetSample = std::clamp<size_t>(targetSample, config.hash_sample_min, config.hash_sample_max);
  const size_t sampleBudget = std::min<size_t>(targetSample, hashs1.size());
  std::vector<uint64_t> sampleVals;
  if (sampleBudget > 0) {
    sampleVals.reserve(sampleBudget);
    if (hashs1.size() <= sampleBudget) {
      sampleVals = hashs1;
    } else if (weightCtx.freqSketch) {
      // Mix: prefer rarer hashes for a sharper candidate set, but keep some
      // uniform coverage to avoid missing genus-level candidates.
      const size_t rareQuota = std::max<size_t>(1, sampleBudget / 2);
      const size_t over =
          std::min<size_t>(hashs1.size(), std::max<size_t>(sampleBudget * 8, sampleBudget));
      std::vector<std::pair<uint32_t, uint64_t>> scored;
      scored.reserve(over);
      size_t step = std::max<size_t>(1, hashs1.size() / over);
      for (size_t i = 0; i < hashs1.size() && scored.size() < over; i += step) {
        uint64_t v = hashs1[i];
        scored.emplace_back(weightCtx.freqSketch->estimate(v), v);
      }
      if (!hashs1.empty() && scored.size() < over) {
        uint64_t v = hashs1.back();
        scored.emplace_back(weightCtx.freqSketch->estimate(v), v);
      }
      std::sort(scored.begin(), scored.end(),
                [](const auto &a, const auto &b) {
                  if (a.first != b.first) {
                    return a.first < b.first;
                  }
                  return a.second < b.second;
                });
      for (const auto &kv : scored) {
        sampleVals.push_back(kv.second);
        if (sampleVals.size() >= rareQuota) {
          break;
        }
      }
      // Top up with a coarse uniform stride.
      size_t step2 = std::max<size_t>(1, hashs1.size() / sampleBudget);
      for (size_t i = 0; i < hashs1.size() && sampleVals.size() < sampleBudget;
           i += step2) {
        sampleVals.push_back(hashs1[i]);
      }
      if (sampleVals.size() < sampleBudget && !hashs1.empty()) {
        sampleVals.push_back(hashs1.back());
      }
      std::sort(sampleVals.begin(), sampleVals.end());
      sampleVals.erase(std::unique(sampleVals.begin(), sampleVals.end()),
                       sampleVals.end());
      if (sampleVals.size() > sampleBudget) {
        sampleVals.resize(sampleBudget);
      }
    } else {
      size_t step = std::max<size_t>(1, hashs1.size() / sampleBudget);
      for (size_t i = 0; i < hashs1.size() && sampleVals.size() < sampleBudget;
           i += step) {
        sampleVals.push_back(hashs1[i]);
      }
      if (sampleVals.size() < sampleBudget && !hashs1.empty()) {
        sampleVals.push_back(hashs1.back());
      }
    }
  }

  std::vector<uint64_t> routeVals = sampleVals;
  bool used_rare_route = false;
  if (!sampleVals.empty() && weightCtx.freqSketch) {
    std::vector<std::pair<uint32_t, uint64_t>> scored;
    scored.reserve(sampleVals.size());
    for (uint64_t v : sampleVals) {
      scored.emplace_back(weightCtx.freqSketch->estimate(v), v);
    }
    const size_t routeBudget =
        std::max<size_t>(1, sampleVals.size() / 2);
    routeVals = select_rare_route_values(scored, routeBudget);
    if (routeVals.empty()) {
      routeVals = sampleVals;
    }
    used_rare_route = (routeVals.size() < sampleVals.size());
  }
  if (sampleBudget > 0) {
    fileInfo.preem_route_reads += 1;
    if (used_rare_route) {
      fileInfo.preem_route_rare_reads += 1;
    }
  }

  std::vector<std::vector<uint32_t>> sampleCount;
  std::vector<std::pair<uint32_t, uint16_t>> touchedS;
  touchedS.reserve(64);
  robin_hood::unordered_flat_map<uint32_t, uint32_t> sampleBinScore;
  uint64_t coarseTotal = 0;
  uint64_t sampleCovered = 0;
  std::vector<uint64_t> deferredEval;
  if (hashs1.size() > 64) {
    deferredEval.reserve(hashs1.size() - 64);
  }

  if (sampleBudget > 0) {
    imcf.bulkCount_sparse(sampleVals, sampleCount, &touchedS);
    sampleBinScore.reserve(touchedS.size());
    for (auto [bi, sp] : touchedS) {
      uint32_t contrib = sampleCount[bi][sp];
      if (contrib == 0) {
        continue;
      }
      coarseTotal += contrib;
      sampleBinScore[bi] += contrib;
    }
  }

  const size_t baseCandidateCap =
      std::min<size_t>(binNumAll, static_cast<size_t>(256));
  const size_t maxCandidateCap =
      std::min<size_t>(binNumAll, static_cast<size_t>(512));
  size_t candidateCap = baseCandidateCap;
  constexpr double kHeadMassThresh = 0.5;
  std::vector<uint32_t> topBins;
  bool fallback_full = (coarseTotal == 0);

  robin_hood::unordered_flat_set<uint32_t> candidateSet;
  candidateSet.reserve(256);
  robin_hood::unordered_flat_set<uint32_t> lowDegPreserve;
  lowDegPreserve.reserve(64);

  if (!routeVals.empty()) {
    std::vector<uint32_t> routed;
    routed.reserve(16);
    for (auto v : routeVals) {
      routed.clear();
      imcf.route(v, routed);
      if (routed.empty()) {
        continue;
      }
      std::sort(routed.begin(), routed.end());
      routed.erase(std::unique(routed.begin(), routed.end()), routed.end());
      if (routed.size() <= 2) {
        for (uint32_t b : routed) {
          if (b < binNumAll) {
            candidateSet.insert(b);
            lowDegPreserve.insert(b);
          }
        }
      }
    }
  }

  constexpr double coverageTarget = 0.92;
  if (coarseTotal > 0 && !sampleBinScore.empty()) {
    std::vector<std::pair<uint32_t, uint32_t>> ranked;
    ranked.reserve(sampleBinScore.size());
    for (const auto &kv : sampleBinScore) {
      ranked.emplace_back(kv.first, kv.second);
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    const double head_mass = compute_head_mass(ranked, 10, coarseTotal);
    candidateCap = compute_candidate_cap(baseCandidateCap, maxCandidateCap,
                                         head_mass, kHeadMassThresh,
                                         ranked.size());
    fileInfo.preem_cap_checks += 1;
    if (candidateCap > baseCandidateCap) {
      fileInfo.preem_cap_expanded += 1;
    }

    uint64_t goal = static_cast<uint64_t>(
        std::ceil(static_cast<double>(coarseTotal) * coverageTarget));
    uint64_t covered = 0;
    for (const auto &[bin, score] : ranked) {
      if (bin >= binNumAll) {
        continue;
      }
      candidateSet.insert(bin);
      covered += score;
      if (goal > 0 && covered >= goal) {
        break;
      }
    }
    sampleCovered = covered;
    fallback_full = false;
  }

  if (!fallback_full) {
    for (const auto &[bin, _] : sampleBinScore) {
      if (bin < binNumAll) {
        candidateSet.insert(bin);
      }
    }
  }

  if (!candidateSet.empty()) {
    std::vector<std::pair<uint32_t, uint64_t>> weighted;
    weighted.reserve(candidateSet.size());
    for (uint32_t bin : candidateSet) {
      uint64_t weight = 0;
      if (auto it = sampleBinScore.find(bin); it != sampleBinScore.end()) {
        weight += it->second * 4ull;
      }
      if (bin < heat.score.size()) {
        weight += heat.score[bin];
      }
      if (lowDegPreserve.find(bin) != lowDegPreserve.end()) {
        weight += (1ull << 32);
      }
      weighted.emplace_back(bin, weight);
    }
    std::sort(weighted.begin(), weighted.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    topBins.reserve(std::min(candidateCap, weighted.size()));
    for (const auto &[bin, _] : weighted) {
      topBins.push_back(bin);
      if (topBins.size() >= candidateCap) {
        break;
      }
    }
  }

  if (topBins.empty()) {
    fallback_full = true;
    topBins.resize(binNumAll);
    std::iota(topBins.begin(), topBins.end(), 0u);
  }

  if (!fallback_full) {
    std::sort(topBins.begin(), topBins.end());
    topBins.erase(std::unique(topBins.begin(), topBins.end()), topBins.end());
  }

  bool candidateEmpty = fallback_full || topBins.empty();
  if (!fallback_full) {
    for (auto bin : topBins) {
      uint32_t delta = 1;
      if (auto it = sampleBinScore.find(bin); it != sampleBinScore.end()) {
        delta = std::max<uint32_t>(delta, it->second);
      }
      heat.boost(bin, delta);
    }
  }
  heat.decay_if_needed();

  if (!fallback_full && topBins.size() != binNumAll) {
    std::sort(topBins.begin(), topBins.end());
    topBins.erase(std::unique(topBins.begin(), topBins.end()), topBins.end());
  }

  robin_hood::unordered_flat_map<uint32_t, double> tidScore;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> uniqueHits;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> consistencyHits;
  tidScore.reserve(128);
  uniqueHits.reserve(128);
  consistencyHits.reserve(128);

  robin_hood::unordered_flat_map<uint32_t, uint32_t> binHitCount;
  binHitCount.reserve(128);

  std::vector<uint32_t> minimizerTids;
  minimizerTids.reserve(16);
  std::vector<uint32_t> minimizerBins;
  minimizerBins.reserve(16);

  double eff_eval = 0.0;
  size_t n_eval = 0;

  const std::vector<uint32_t> *activeSubset =
      (fallback_full || topBins.size() == binNumAll) ? nullptr : &topBins;

  auto bucket_degree = [](size_t d) -> size_t {
    if (d <= 1) {
      return 0;
    }
    if (d <= 3) {
      return 1;
    }
    if (d <= 7) {
      return 2;
    }
    if (d <= 15) {
      return 3;
    }
    if (d <= 31) {
      return 4;
    }
    return 5;
  };

  auto evaluate_minimizer = [&](uint64_t value,
                                const std::vector<uint32_t> *subset) -> double {
    minimizerTids.clear();
    minimizerBins.clear();
    uint32_t last_hit_bin = std::numeric_limits<uint32_t>::max();
    auto emit = [&](uint32_t bin, uint16_t sp) {
      if (bin >= tax.idx2id.size()) {
        return;
      }
      const auto &speciesVec = tax.idx2id[bin];
      if (sp >= speciesVec.size()) {
        return;
      }
      uint32_t tid = speciesVec[sp];
      if (tid >= tax.id2str.size()) {
        return;
      }
      if (weightCtx.tid2speciesRep && tid < weightCtx.tid2speciesRep->size()) {
        tid = (*weightCtx.tid2speciesRep)[tid];
      }
      minimizerTids.push_back(tid);
      if (bin != last_hit_bin) {
        minimizerBins.push_back(bin);
        last_hit_bin = bin;
      }
      ++binHitCount[bin];
    };

    if (!subset) {
      imcf.bulkContain_events(value, emit);
    } else {
      imcf.bulkContain_events_subset(value, *subset, emit);
    }

    if (minimizerTids.empty()) {
      return 0.0;
    }

    std::sort(minimizerTids.begin(), minimizerTids.end());
    minimizerTids.erase(
        std::unique(minimizerTids.begin(), minimizerTids.end()),
        minimizerTids.end());

    const size_t deg_subset = minimizerTids.size();
    if (deg_subset == 0) {
      return 0.0;
    }
    size_t deg_effective = deg_subset;
    if (config.deg_by_species && weightCtx.tid2speciesGroup &&
        deg_subset > 1) {
      std::vector<uint32_t> groups;
      groups.reserve(deg_subset);
      const auto &map = *weightCtx.tid2speciesGroup;
      for (uint32_t tid : minimizerTids) {
        if (tid < map.size()) {
          groups.push_back(map[tid]);
        } else {
          groups.push_back(0x80000000u | tid);
        }
      }
      std::sort(groups.begin(), groups.end());
      groups.erase(std::unique(groups.begin(), groups.end()), groups.end());
      if (!groups.empty()) {
        deg_effective = groups.size();
      }
    }

    size_t df_bins = 0;
    if (!minimizerBins.empty()) {
      std::sort(minimizerBins.begin(), minimizerBins.end());
      minimizerBins.erase(std::unique(minimizerBins.begin(), minimizerBins.end()),
                          minimizerBins.end());
      df_bins = minimizerBins.size();
    }
    if (df_bins == 0) {
      df_bins = deg_subset;
    }

    double exclusivityWeight = 1.0;
    if (exclusiveGamma > 0.0 && deg_effective > 0) {
      exclusivityWeight =
          std::pow(static_cast<double>(deg_effective), -exclusiveGamma);
    }

    double totalBins = subset ? static_cast<double>(subset->size())
                              : static_cast<double>(binNumAll);
    if (totalBins <= 0.0) {
      totalBins = static_cast<double>(binNumAll);
    }

    double idf = std::log2((totalBins + 1.0) /
                           (static_cast<double>(df_bins) + 1.0));
    idf = std::clamp(idf, 0.5, config.idf_max);

    double freqFactor = 1.0;
    const bool has_freq = weightCtx.enabled();
    uint32_t df_est = 0;
    if (has_freq) {
      df_est = weightCtx.freqSketch->estimate(value);
      double vocab = static_cast<double>(
          std::max<uint64_t>(1, weightCtx.freqStats.nonzero_counters));
      double bg_idf =
          std::log((vocab + 1.0) / (static_cast<double>(df_est) + 1.0));
      bg_idf = bg_idf / std::log(2.0);
      bg_idf = std::clamp(bg_idf, 0.25, 6.0);
      if (weightCtx.freqStats.df_high_threshold !=
              std::numeric_limits<uint32_t>::max() &&
          df_est >= weightCtx.freqStats.df_high_threshold) {
        double damp = std::log2(
            static_cast<double>(weightCtx.freqStats.df_high_threshold + 1.0));
        double denom = std::log2(static_cast<double>(df_est + 2.0));
        if (denom > 0.0) {
          damp = std::clamp(damp / denom, 0.05, 0.5);
          freqFactor *= damp;
        } else {
          freqFactor *= 0.25;
        }
      }
      freqFactor *= bg_idf;
      freqFactor = std::clamp(freqFactor, 0.05, 8.0);
    }

    const bool uniqueEdge = allow_unique_edge(
        deg_effective, df_bins, has_freq, weightCtx.freq_trusted, df_est,
        imcfConfig.presenceUniqueDeg);
    const bool low_df_boost = allow_low_df_boost(
        df_bins, has_freq, weightCtx.freq_trusted, df_est,
        imcfConfig.presenceUniqueDeg);
    const bool localUniqueEdge = (deg_effective == 1 && df_bins == 1);
    double denom = std::log2(2.0 + static_cast<double>(deg_effective));
    double weight = denom > 0.0 ? 1.0 / denom : 1.0;
    double contrib = idf * weight * freqFactor * exclusivityWeight;
    if (uniqueEdge) {
      contrib *= kUniqueEdgeBonus;
    } else if (low_df_boost) {
      contrib *= 1.5;
    }
    for (uint32_t tid : minimizerTids) {
      tidScore[tid] += contrib;
      ++consistencyHits[tid];
    }
    if (presenceEnabled && !minimizerTids.empty()) {
      const double hit_weight = presence_weight;
      const double score_weight = exclusivityWeight * hit_weight;
      uint32_t unique_bucket = std::numeric_limits<uint32_t>::max();
      uint32_t breadth_bucket = std::numeric_limits<uint32_t>::max();
      if (localUniqueEdge && presenceSketchBits > 0) {
        unique_bucket = sketch_index(value ^ 0x9e3779b97f4a7c15ULL);
        if (!minimizerBins.empty()) {
          breadth_bucket = sketch_index(
              static_cast<uint64_t>(minimizerBins.front()) ^
              0xd1b54a32d192ed03ULL);
        }
      }
      for (uint32_t tid : minimizerTids) {
        presenceAcc->add_target(tid, hit_weight, score_weight, uniqueEdge,
                                localUniqueEdge, unique_bucket,
                                breadth_bucket);
      }
      if (presenceDecoyReps > 0 && deg_subset > 0) {
        for (size_t rep = 0; rep < presenceDecoyReps; ++rep) {
          uint64_t mix = splitmix64(
              value ^ (static_cast<uint64_t>(rep + 1) << 17) ^ decoySeed);
          size_t pick = static_cast<size_t>(mix % deg_subset);
          uint32_t decoyTid = minimizerTids[pick];
          presenceAcc->add_decoy(rep, decoyTid, score_weight);
        }
      }
    }
    if (uniqueEdge && !minimizerTids.empty()) {
      ++uniqueHits[minimizerTids.front()];
    }

    return contrib;
  };

  auto recompute_subset_state = [&]() {
    activeSubset = (fallback_full || topBins.size() == binNumAll)
                       ? nullptr
                       : &topBins;
  };

  robin_hood::unordered_flat_set<uint32_t> topBinSet;
  if (!fallback_full) {
    topBinSet.reserve(topBins.size());
    for (auto bin : topBins) {
      topBinSet.insert(bin);
    }
  }

  auto compute_top = [&](uint32_t &bestTid, uint32_t &secondTid, double &best,
                         double &second) {
    best = 0.0;
    second = 0.0;
    bestTid = std::numeric_limits<uint32_t>::max();
    secondTid = std::numeric_limits<uint32_t>::max();
    for (const auto &kv : tidScore) {
      double c = kv.second;
      if (c > best) {
        second = best;
        secondTid = bestTid;
        best = c;
        bestTid = kv.first;
      } else if (c > second) {
        second = c;
        secondTid = kv.first;
      }
    }
  };

  auto compute_ratio = [](double best_score, double second_score) {
    if (second_score > 0.0) {
      return best_score /
             std::max(second_score, std::numeric_limits<double>::min());
    }
    return best_score <= 0.0 ? 0.0
                             : std::numeric_limits<double>::infinity();
  };

  struct EvidenceStats {
    uint32_t bestTid = std::numeric_limits<uint32_t>::max();
    uint32_t secondTid = std::numeric_limits<uint32_t>::max();
    double best = 0.0;
    double second = 0.0;
    double ratio = 0.0;
    double gap = 0.0;
    size_t uniqueCount = 0;
    double uniqueRatio = 0.0;
    size_t consistency = 0;
  };

  auto collect_stats = [&]() -> EvidenceStats {
    EvidenceStats stats;
    compute_top(stats.bestTid, stats.secondTid, stats.best, stats.second);
    stats.ratio = compute_ratio(stats.best, stats.second);
    stats.gap = stats.best - stats.second;
    if (stats.bestTid != std::numeric_limits<uint32_t>::max()) {
      if (auto it = uniqueHits.find(stats.bestTid); it != uniqueHits.end()) {
        stats.uniqueCount = it->second;
      }
      if (auto it = consistencyHits.find(stats.bestTid);
          it != consistencyHits.end()) {
        stats.consistency = it->second;
      }
    }
    double denom = eff_eval > 0.0 ? eff_eval : 1.0;
    stats.uniqueRatio = static_cast<double>(stats.uniqueCount) / denom;
    return stats;
  };

  auto consistency_ratio = [&](uint32_t tid) -> double {
    if (tid == std::numeric_limits<uint32_t>::max() || n_eval == 0) {
      return 0.0;
    }
    auto it = consistencyHits.find(tid);
    if (it == consistencyHits.end()) {
      return 0.0;
    }
    return static_cast<double>(it->second) /
           static_cast<double>(std::max<size_t>(1, n_eval));
  };

  size_t n0 = std::min<size_t>(64, hashs1.size());
  for (size_t i = 0; i < n0; ++i) {
    eff_eval += evaluate_minimizer(hashs1[i], activeSubset);
  }
  n_eval = n0;

  auto meets_quick = [&](const EvidenceStats &s) {
    if (s.bestTid == std::numeric_limits<uint32_t>::max()) {
      return false;
    }
    double thr_conf_local = std::ceil(config.shotThreshold * eff_eval);
    double gap_need_local = std::max(0.5, eff_eval / 24.0);
    bool strong = (s.best >= thr_conf_local);
    bool unique_ok = (s.uniqueCount >= 3) ||
                     (s.uniqueCount >= 2 && s.uniqueRatio >= 0.12);
    bool stable = (s.gap >= gap_need_local) && (s.ratio >= 1.35);
    size_t bestRoundedLocal = static_cast<size_t>(
        std::max<double>(0.0, std::llround(s.best)));
    size_t secondRoundedLocal = static_cast<size_t>(
        std::max<double>(0.0, std::llround(s.second)));
    size_t needLocal = static_cast<size_t>(std::ceil(thr_conf_local));
    auto dc = decide_high_conf(bestRoundedLocal, secondRoundedLocal, eff_eval);
    bool margin_ok = (bestRoundedLocal >= needLocal) && dc.accept;
    return strong && unique_ok && stable && margin_ok;
  };

  EvidenceStats stats = collect_stats();
  bool highConfPre = meets_quick(stats);

  if (!highConfPre && hashs1.size() > n0) {
    size_t mask = (static_cast<size_t>(1) << 3) - 1;
    for (size_t i = n0; i < hashs1.size(); ++i) {
      if ((hashs1[i] & mask) != 0) {
        deferredEval.push_back(hashs1[i]);
        continue;
      }
      eff_eval += evaluate_minimizer(hashs1[i], activeSubset);
      ++n_eval;
    }
    stats = collect_stats();
    highConfPre = meets_quick(stats);

    if (!highConfPre && !deferredEval.empty()) {
      for (auto value : deferredEval) {
        eff_eval += evaluate_minimizer(value, activeSubset);
        ++n_eval;
      }
      deferredEval.clear();
      stats = collect_stats();
      highConfPre = meets_quick(stats);
    }
  }

  double thr_conf = std::max(1.0, std::ceil(config.shotThreshold * eff_eval));
  double gap_need = std::max(0.5, eff_eval / 24.0);

  uint32_t bestTid = stats.bestTid;
  double best = stats.best;
  double second = stats.second;
  double best_ratio = stats.ratio;
  double gap = stats.gap;
  size_t uniqueCount = stats.uniqueCount;
  double uniqueRatio = stats.uniqueRatio;
  std::string bestTaxidStr;
  if (bestTid != std::numeric_limits<uint32_t>::max() &&
      bestTid < tax.id2str.size()) {
    bestTaxidStr = tax.id2str[bestTid];
  }

  bool expanded = false;
  if (!fallback_full && bestTid < tax.tid2bin.size()) {
    const auto &shards = tax.tid2bin[bestTid];
    for (uint32_t bin : shards) {
      if (topBinSet.find(bin) == topBinSet.end()) {
        topBins.push_back(bin);
        topBinSet.insert(bin);
        expanded = true;
      }
    }
  }

  if (expanded) {
    std::sort(topBins.begin(), topBins.end());
    topBins.erase(std::unique(topBins.begin(), topBins.end()), topBins.end());
    fallback_full = (topBins.size() == binNumAll);
    recompute_subset_state();

    tidScore.clear();
    binHitCount.clear();
    uniqueHits.clear();
    consistencyHits.clear();
    eff_eval = 0.0;
    n_eval = 0;
    for (auto value : hashs1) {
      eff_eval += evaluate_minimizer(value, activeSubset);
      ++n_eval;
    }
    stats = collect_stats();
    highConfPre = meets_quick(stats);
    thr_conf = std::ceil(config.shotThreshold * eff_eval);
    gap_need = std::max(0.5, eff_eval / 24.0);
    bestTid = stats.bestTid;
    best = stats.best;
    second = stats.second;
    best_ratio = stats.ratio;
    gap = stats.gap;
    uniqueCount = stats.uniqueCount;
    uniqueRatio = stats.uniqueRatio;
  }

  // Optional: dump raw tidScore (pre-EM candidates before thresholds).
  dump_tidscore_line(id, tax, tidScore);

  size_t bestRounded =
      static_cast<size_t>(std::max<double>(0.0, std::llround(best)));
  size_t secondRounded =
      static_cast<size_t>(std::max<double>(0.0, std::llround(second)));
  size_t thrConfNeed = static_cast<size_t>(std::ceil(thr_conf));
  auto dc = decide_high_conf(bestRounded, secondRounded, eff_eval);
  bool marginAccept = (bestRounded >= thrConfNeed) && dc.accept;
  highConfPre = highConfPre && marginAccept;

  // --- EM 模式下的模糊保护 ---
  // 当使用 EM 时，如果 Top1/Top2 的分值比例不足 3 倍，视为不确定，让其进入 EM。
  // 过低的比例容易被近缘物种或 decoy 抢占，需延后由 EM 全局推断。
  if (config.em && best_ratio < 3.0) {
    highConfPre = false;
  }
  // --------------------------------

  classifyResult result;
  result.evaluated = eff_eval;
  result.id = id;
  if (weightCtx.has_sample_weights()) {
    auto it_w = weightCtx.sampleWeights->find(id);
    result.sample_weight = (it_w != weightCtx.sampleWeights->end())
                               ? it_w->second
                               : 1.0;
  }
  result.best_taxid_hint = bestTaxidStr;
  const double effCap = std::max(1.0, eff_eval);
  std::pair<std::string, double> maxCount;
  bool maxCountValid = false;

  if (!maxCountValid && bestTid < tax.id2str.size()) {
    double bestEvidence = std::clamp(best, 0.0, effCap);
    if (bestEvidence > 0.0) {
      maxCount = std::make_pair(tax.id2str[bestTid], bestEvidence);
      maxCountValid = true;
    }
  }

  bool use_em = config.em && !highConfPre;

  double maxEvidence = std::min(best, eff_eval);
  double beta = config.firstFilterBeta;
  bool beta_user = config.firstFilterBeta_user;
  if (beta <= 0.0) {
    beta = 0.8;
    beta_user = false;
  }
  if (config.em && !beta_user) {
    beta = 0.45; // EM 模式放宽初筛门槛，确保真实物种不被挡在外
  }
  beta = std::clamp(beta, 0.0, 1.0);
  size_t thr_beta =
      static_cast<size_t>(std::floor(beta * std::max(0.0, maxEvidence)));
  size_t thr_eval = static_cast<size_t>(std::ceil(
      config.shotThreshold *
      (config.adaptive_shot ? eff_eval
                            : static_cast<double>(hashNum))));
  if (thr_eval == 0) {
    thr_eval = 1;
  }
  if (use_em) {
    double base =
        config.adaptive_shot ? eff_eval : static_cast<double>(hashNum);
    double softened_ratio = std::min(config.shotThreshold, 0.45);
    size_t em_eval =
        static_cast<size_t>(std::ceil(base * softened_ratio));
    if (em_eval == 0 && base > 0.0) {
      em_eval = 1;
    }
    thr_eval = std::min(thr_eval, em_eval);
  }
  size_t thr_min_eval = 0;
  if (n_eval > 0) {
    double factor = use_em ? 0.15 : 0.30;
    thr_min_eval =
        static_cast<size_t>(std::ceil(factor * static_cast<double>(n_eval)));
  }
  if (use_em) {
    thr_min_eval = std::max<size_t>(thr_min_eval, 4);
  } else {
    thr_min_eval = std::max<size_t>(thr_min_eval, 8);
  }

  double TOT = 0.0;
  for (const auto &kv : tidScore) {
    TOT += kv.second;
  }

  size_t M = 0;
  if (!topBins.empty()) {
    for (auto b : topBins) {
      if (b < tax.idx2id.size()) {
        M += tax.idx2id[b].size();
      }
    }
  } else {
    for (size_t b = 0; b < tax.idx2id.size(); ++b) {
      M += tax.idx2id[b].size();
    }
  }
  M = std::max<size_t>(M, 1);

  size_t thr_final = 0;
  if (use_em) {
    size_t thr_beta_eval = std::min(thr_beta, thr_eval);
    size_t min_hits = 1;
    thr_final = std::max({thr_beta_eval, min_hits, thr_min_eval});
  } else {
    thr_final = std::max({thr_eval, thr_beta, thr_min_eval});
  }

  size_t baseTopK =
      config.preEmTopK > 0 ? static_cast<size_t>(config.preEmTopK) : 32;

  const size_t thr_final_raw = thr_final;
  size_t thr_final_used = thr_final_raw;
  bool preem_beta_relax_applied = false;
  double preem_beta_local = beta;
  size_t preem_n_strict = 0;
  size_t preem_n_used = 0;

  if (highConfPre && bestTid < tax.id2str.size()) {
    double bestEvidence = std::clamp(best, 0.0, effCap);
    const std::string &taxid = tax.id2str[bestTid];
    result.taxidCount.emplace_back(taxid, bestEvidence);
    maxCount = std::make_pair(taxid, bestEvidence);
    maxCountValid = true;

  } else {
    if (use_em && config.low_div_active && !tidScore.empty()) {
      const size_t floorK = std::min<size_t>(128, baseTopK);
      const size_t poolK = std::max(baseTopK, floorK);
      std::vector<std::pair<uint32_t, double>> ranked;
      ranked.reserve(tidScore.size());
      for (const auto &kv : tidScore) {
        ranked.emplace_back(kv.first, kv.second);
      }
      size_t want = std::min(poolK, ranked.size());
      if (want > 0) {
        std::partial_sort(
            ranked.begin(), ranked.begin() + want, ranked.end(),
            [](const auto &a, const auto &b) { return a.second > b.second; });
        ranked.resize(want);
      }

      size_t thr_min_eval_lowdiv = thr_min_eval;
      if (thr_min_eval > 0) {
        thr_min_eval_lowdiv = std::max<size_t>(
            1, static_cast<size_t>(std::floor(
                   0.3 * static_cast<double>(thr_min_eval))));
      }
      size_t thr_base = std::max(thr_eval, thr_min_eval_lowdiv);

      std::vector<std::pair<std::string, double>> floor;
      std::vector<std::pair<std::string, double>> tail;
      floor.reserve(floorK);
      tail.reserve(ranked.size());
      const size_t floorCount = std::min(floorK, ranked.size());
      for (size_t i = 0; i < ranked.size(); ++i) {
        const uint32_t tid_id = ranked[i].first;
        if (tid_id >= tax.id2str.size()) {
          continue;
        }
        const double countVal =
            std::clamp(ranked[i].second, 0.0, effCap);
        if (countVal <= 0.0) {
          continue;
        }
        const std::string &taxid = tax.id2str[tid_id];
        if (i < floorCount) {
          floor.emplace_back(taxid, countVal);
        } else if (countVal >= static_cast<double>(thr_base)) {
          tail.emplace_back(taxid, countVal);
        }
      }

      size_t thr_beta_eval = std::min(thr_beta, thr_eval);
      if (floor.size() + tail.size() > baseTopK && thr_beta_eval > 0) {
        const double thr_beta_val = static_cast<double>(thr_beta_eval);
        std::vector<std::pair<std::string, double>> tail_kept;
        tail_kept.reserve(tail.size());
        for (const auto &kv : tail) {
          if (kv.second >= thr_beta_val) {
            tail_kept.push_back(kv);
          }
        }
        tail.swap(tail_kept);
      }

      result.taxidCount.reserve(floor.size() + tail.size());
      result.taxidCount.insert(result.taxidCount.end(), floor.begin(),
                               floor.end());
      result.taxidCount.insert(result.taxidCount.end(), tail.begin(),
                               tail.end());
      if (baseTopK > 0 && result.taxidCount.size() > baseTopK) {
        std::nth_element(
            result.taxidCount.begin(),
            result.taxidCount.begin() + baseTopK, result.taxidCount.end(),
            [](const auto &a, const auto &b) { return a.second > b.second; });
        result.taxidCount.resize(baseTopK);
      }
    } else {
      for (const auto &[tid_id, rawScore] : tidScore) {
        double countVal = std::clamp(rawScore, 0.0, effCap);
        if (countVal >= static_cast<double>(thr_final_raw)) {
          const std::string &taxid = tax.id2str[tid_id];
          result.taxidCount.emplace_back(taxid, countVal);
          if (tid_id == bestTid) {
            maxCount = std::make_pair(taxid, countVal);
            maxCountValid = true;
          }
        }
      }

      // High-div / EM only: read-level pre-EM gate relaxation (fixed budget).
      // Motivation: wrong reads often contain truth branch in tidScore(topN) but
      // are pruned by thr_beta_eval/thr_final, leaving no chance for EM/post.
      // Guard rails: do NOT pad candidates; do NOT expand dynamicTopK due to
      // size>baseTopK if this relaxation is applied (see binOverflow gate).
      preem_n_strict = result.taxidCount.size();
      preem_n_used = preem_n_strict;

      constexpr size_t kPreemBetaRelaxDelta = 8;
      constexpr double kPreemBetaRelaxEffEvalMin = 48.0;
      if (use_em && !config.low_div_active && !beta_user &&
          eff_eval >= kPreemBetaRelaxEffEvalMin && baseTopK > 0) {
        fileInfo.preem_beta_relax_checks += 1;
        auto decision = decide_preem_beta_relax(
            /*use_em=*/use_em, /*low_div_active=*/config.low_div_active,
            /*beta_user=*/beta_user, /*base_beta=*/beta,
            /*maxEvidence=*/maxEvidence, /*eff_eval=*/eff_eval,
            /*best_ratio=*/best_ratio, /*unique_ratio=*/uniqueRatio,
            /*base_topk=*/baseTopK, /*n_strict=*/preem_n_strict,
            /*thr_beta=*/thr_beta, /*thr_eval=*/thr_eval,
            /*thr_min_eval=*/thr_min_eval, /*thr_final_raw=*/thr_final_raw,
            /*delta=*/kPreemBetaRelaxDelta, /*eff_eval_min=*/kPreemBetaRelaxEffEvalMin);
        if (decision.applied) {
          preem_beta_relax_applied = true;
          preem_beta_local = decision.beta_local;
          thr_final_used = decision.thr_final_used;
          fileInfo.preem_beta_relax_applied += 1;
          fileInfo.preem_beta_relax_thr_drop_sum +=
              (thr_final_raw - thr_final_used);

          // Rebuild candidate list using relaxed threshold.
          result.taxidCount.clear();
          for (const auto &[tid_id, rawScore] : tidScore) {
            double countVal = std::clamp(rawScore, 0.0, effCap);
            if (countVal >= static_cast<double>(thr_final_used)) {
              const std::string &taxid = tax.id2str[tid_id];
              result.taxidCount.emplace_back(taxid, countVal);
              if (tid_id == bestTid) {
                maxCount = std::make_pair(taxid, countVal);
                maxCountValid = true;
              }
            }
          }
          preem_n_used = result.taxidCount.size();
        }
      }
    }
  }

  if (use_em && !highConfPre && result.taxidCount.size() == 1 &&
      !tidScore.empty()) {
    const bool unique_ok = (uniqueCount >= 3) ||
                           (uniqueCount >= 2 && uniqueRatio >= 0.12);
    if (!unique_ok) {
      constexpr size_t kMinCand = 4;
      if (result.taxidCount.size() < kMinCand) {
        robin_hood::unordered_flat_set<uint32_t> seen;
        seen.reserve(kMinCand + 2);
        for (const auto &kv : result.taxidCount) {
          auto it = tax.str2id.find(kv.first);
          if (it != tax.str2id.end()) {
            seen.insert(it->second);
          }
        }
        std::vector<std::pair<uint32_t, double>> ranked;
        ranked.reserve(tidScore.size());
        for (const auto &kv : tidScore) {
          ranked.emplace_back(kv.first, kv.second);
        }
        const size_t want =
            std::min<size_t>(ranked.size(), kMinCand * 4);
        if (want > 0) {
          std::partial_sort(
              ranked.begin(), ranked.begin() + want, ranked.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
          for (size_t i = 0;
               i < want && result.taxidCount.size() < kMinCand; ++i) {
            const uint32_t tid_id = ranked[i].first;
            if (seen.find(tid_id) != seen.end()) {
              continue;
            }
            if (tid_id >= tax.id2str.size()) {
              continue;
            }
            const std::string &taxid = tax.id2str[tid_id];
            const double countVal =
                std::clamp(ranked[i].second, 0.0, effCap);
            if (countVal <= 0.0) {
              continue;
            }
            result.taxidCount.emplace_back(taxid, countVal);
            seen.insert(tid_id);
          }
        }
      }
    }
  }

  if (result.taxidCount.empty() && use_em && maxCountValid &&
      maxCount.second > 0.0) {
    result.taxidCount.emplace_back(maxCount);
  }

  size_t dynamicTopK = baseTopK;
  double tieGapNeed = std::max(2.0, eff_eval / 24.0);
  bool nearTie = (gap < tieGapNeed) || (best_ratio < 1.30);
  bool binOverflow = preem_bin_overflow(
      /*fallback_full=*/fallback_full, /*topBins_size=*/topBins.size(),
      /*taxidCount_size=*/result.taxidCount.size(), /*baseTopK=*/baseTopK,
      /*suppress_size_overflow=*/preem_beta_relax_applied);
  if (preem_beta_relax_applied && result.taxidCount.size() > baseTopK) {
    fileInfo.preem_beta_relax_suppressed_overflow += 1;
  }
  if (nearTie || binOverflow) {
    dynamicTopK = std::max<size_t>(96, dynamicTopK);
  } else if (best_ratio >= 2.5 &&
             best >= thr_conf + std::max(1.0, eff_eval / 16.0) &&
             dynamicTopK > 8 &&
             // Respect explicit user requests for larger pre-EM candidate lists.
             // Otherwise, large DBs with many strain/assembly taxids can saturate
             // the top-K with same-species variants, preventing sister species
             // from entering EM/posterior lists.
             config.preEmTopK <= 16) {
    dynamicTopK = 16;
  }

  // NCBI-only (optional): collapse strain/subspecies taxids to ONE representative
  // per species before pre-EM topK truncation, so TopK is not saturated by
  // same-species variants and sister species can enter EM/posterior lists.
  //
  // Important: keep the representative taxid string (not the species taxid),
  // so downstream presence/meta that relies on DB taxids stays consistent.
  if (use_em && config.collapse_strain_candidates && weightCtx.ncbiTaxdump &&
      weightCtx.ncbiTaxdump->enabled() && !result.taxidCount.empty()) {
    auto parse_u32 = [](const std::string &s, uint32_t &out) -> bool {
      if (s.empty()) {
        return false;
      }
      for (unsigned char c : s) {
        if (!std::isdigit(c)) {
          return false;
        }
      }
      try {
        unsigned long v = std::stoul(s);
        if (v > std::numeric_limits<uint32_t>::max()) {
          return false;
        }
        out = static_cast<uint32_t>(v);
        return true;
      } catch (...) {
        return false;
      }
    };

    robin_hood::unordered_flat_map<uint32_t, std::pair<std::string, double>>
        bestBySpecies;
    bestBySpecies.reserve(result.taxidCount.size());
    std::vector<std::pair<std::string, double>> passthrough;
    passthrough.reserve(4);
    double unclassifiedScore = 0.0;
    bool hasUnclassified = false;

    for (const auto &kv : result.taxidCount) {
      const std::string &taxid = kv.first;
      const double score = kv.second;
      if (taxid == "unclassified") {
        hasUnclassified = true;
        unclassifiedScore += score;
        continue;
      }
      uint32_t tid = 0;
      if (!parse_u32(taxid, tid)) {
        passthrough.push_back(kv);
        continue;
      }
      const uint32_t sid = weightCtx.ncbiTaxdump->to_species(tid);
      auto it = bestBySpecies.find(sid);
      if (it == bestBySpecies.end() || score > it->second.second) {
        bestBySpecies[sid] = kv; // keep representative taxid with max score
      }
    }

    std::vector<std::pair<std::string, double>> collapsed;
    collapsed.reserve(bestBySpecies.size() + passthrough.size() +
                      (hasUnclassified ? 1 : 0));
    for (const auto &kv : bestBySpecies) {
      collapsed.push_back(kv.second);
    }
    for (const auto &kv : passthrough) {
      collapsed.push_back(kv);
    }
    if (hasUnclassified) {
      collapsed.emplace_back("unclassified", unclassifiedScore);
    }
    std::sort(collapsed.begin(), collapsed.end(),
              [](const auto &a, const auto &b) {
                if (a.second != b.second) {
                  return a.second > b.second;
                }
                return a.first < b.first;
              });
    result.taxidCount.swap(collapsed);
  }

  // High-div / EM: if the collapsed (species-level) candidate list is still a
  // singleton, keep a tiny floor of additional *species* candidates alive using
  // the same fixed per-read budget (no pad-to-32, no global expansion).
  if (use_em && !config.low_div_active && weightCtx.ncbiTaxdump &&
      weightCtx.ncbiTaxdump->enabled() && config.collapse_strain_candidates &&
      !tidScore.empty()) {
    constexpr size_t kFloorTrigger = 1;

    size_t nonUnclassified = 0;
    bool hasUnclassified = false;
    double topScore = 0.0;
    for (const auto &kv : result.taxidCount) {
      if (kv.first == "unclassified") {
        hasUnclassified = true;
      } else {
        nonUnclassified += 1;
        topScore = std::max(topScore, kv.second);
      }
    }

    const size_t floorTargetNonUncl =
        std::min<size_t>(static_cast<size_t>(config.preem_floor_target),
                         dynamicTopK);
    const size_t floorTargetTotal = floorTargetNonUncl + (hasUnclassified ? 1 : 0);

    if (floorTargetTotal > 0 && nonUnclassified == kFloorTrigger &&
        result.taxidCount.size() < floorTargetTotal && topScore > 0.0) {
      fileInfo.preem_floor_checks += 1;
      const size_t before = result.taxidCount.size();

      robin_hood::unordered_flat_set<uint32_t> seenSpecies;
      seenSpecies.reserve(result.taxidCount.size() + 4);
      for (const auto &kv : result.taxidCount) {
        const std::string &taxid = kv.first;
        if (taxid == "unclassified") {
          continue;
        }
        uint32_t tid = 0;
        bool ok = !taxid.empty();
        for (unsigned char c : taxid) {
          if (!std::isdigit(c)) {
            ok = false;
            break;
          }
          tid = tid * 10 + static_cast<uint32_t>(c - '0');
        }
        if (!ok) {
          continue;
        }
        seenSpecies.insert(weightCtx.ncbiTaxdump->to_species(tid));
      }

      robin_hood::unordered_flat_map<uint32_t, std::pair<std::string, double>>
          bestBySpecies;
      bestBySpecies.reserve(64);
      double bestOtherScore = 0.0;
      for (const auto &kv : tidScore) {
        const uint32_t tid_id = kv.first;
        if (tid_id >= tax.id2str.size()) {
          continue;
        }
        const uint32_t sid = weightCtx.ncbiTaxdump->to_species(tid_id);
        if (seenSpecies.find(sid) != seenSpecies.end()) {
          continue;
        }
        const double countVal = std::clamp(kv.second, 0.0, effCap);
        if (!(countVal > 0.0)) {
          continue;
        }
        bestOtherScore = std::max(bestOtherScore, countVal);
        auto it = bestBySpecies.find(sid);
        if (it == bestBySpecies.end() || countVal > it->second.second) {
          bestBySpecies[sid] = {tax.id2str[tid_id], countVal};
        }
      }

      if (!should_apply_preem_floor(topScore, bestOtherScore,
                                   config.preem_floor_min_ratio)) {
        fileInfo.preem_floor_skipped_dominant += 1;
      } else {
        const double addMinScore = topScore * config.preem_floor_add_min_ratio;

        std::vector<std::pair<std::string, double>> ranked;
        ranked.reserve(bestBySpecies.size());
        for (const auto &kv : bestBySpecies) {
          const auto &cand = kv.second;
          if (!(cand.second > 0.0)) {
            continue;
          }
          if (cand.second < addMinScore) {
            fileInfo.preem_floor_filtered_weak += 1;
            continue;
          }
          ranked.push_back(cand);
        }

        const size_t want =
            std::min<size_t>(ranked.size(), floorTargetNonUncl * 8);
        if (want > 0 && ranked.size() > want) {
          std::nth_element(
              ranked.begin(), ranked.begin() + static_cast<std::ptrdiff_t>(want),
              ranked.end(), [](const auto &a, const auto &b) {
                if (a.second != b.second) {
                  return a.second > b.second;
                }
                return a.first < b.first;
              });
          ranked.resize(want);
        }
        std::sort(ranked.begin(), ranked.end(),
                  [](const auto &a, const auto &b) {
                    if (a.second != b.second) {
                      return a.second > b.second;
                    }
                    return a.first < b.first;
                  });

        ensure_preem_floor_candidates(result.taxidCount, ranked, floorTargetTotal);
        if (result.taxidCount.size() > before) {
          fileInfo.preem_floor_applied += 1;
          fileInfo.preem_floor_added += (result.taxidCount.size() - before);
        }
      }
    }
  }

  if (!result.taxidCount.empty() && use_em) {
    normalize_preem_topk(result.taxidCount, dynamicTopK);
  }

  dump_preem_stats(id, highConfPre, use_em, eff_eval, best, second, best_ratio,
                   gap, uniqueCount, uniqueRatio, thr_beta, thr_eval,
                   thr_min_eval, thr_final_used);

  // High-div / EM only: fixed-budget keepalive of a strong per-read hint
  // candidate. Do NOT expand K; only replace the weakest tail candidate.
  if (use_em && !config.low_div_active && !result.taxidCount.empty() &&
      config.preem_keepalive_min_ratio > 0.0 &&
      config.preem_keepalive_replace_ratio > 0.0 &&
      config.preem_keepalive_abs_min > 0.0 &&
      !result.best_taxid_hint.empty() &&
      result.best_taxid_hint != "unclassified") {
    bool hint_in_topk = false;
    for (const auto &kv : result.taxidCount) {
      if (kv.first == result.best_taxid_hint) {
        hint_in_topk = true;
        break;
      }
    }
    if (!hint_in_topk) {
      uint32_t hint_id = std::numeric_limits<uint32_t>::max();
      if (auto it = tax.str2id.find(result.best_taxid_hint);
          it != tax.str2id.end()) {
        hint_id = it->second;
      }
      double hint_score = 0.0;
      if (hint_id != std::numeric_limits<uint32_t>::max()) {
        if (auto it = tidScore.find(hint_id); it != tidScore.end()) {
          hint_score = std::clamp(it->second, 0.0, effCap);
        }
      }

      fileInfo.preem_keepalive_attempt += 1;
      auto res = preem_keepalive_replace_tail(
          result.taxidCount, {result.best_taxid_hint, hint_score},
          config.preem_keepalive_min_ratio, config.preem_keepalive_replace_ratio,
          config.preem_keepalive_abs_min);
      if (res == KeepaliveResult::kApplied) {
        fileInfo.preem_keepalive_applied += 1;
      } else if (res == KeepaliveResult::kBlockedLowRatio) {
        fileInfo.preem_keepalive_blocked_low_ratio += 1;
      } else if (res == KeepaliveResult::kBlockedLowGain) {
        fileInfo.preem_keepalive_blocked_low_gain += 1;
      } else if (res == KeepaliveResult::kBlockedLowAbs) {
        fileInfo.preem_keepalive_blocked_low_abs += 1;
      }
    }
  }

  // Optional: dump pre-EM candidate list for P1/P2 analysis.
  dump_preem_line(result.id, result.taxidCount, highConfPre);

  if (!result.taxidCount.empty()) {
    if (presenceEnabled && presenceAcc) {
      const auto &top = result.taxidCount.front();
      if (top.first != "unclassified") {
        auto itid = tax.str2id.find(top.first);
        if (itid != tax.str2id.end()) {
          bool unique_ok = (uniqueCount >= 3) ||
                           (uniqueCount >= 2 && uniqueRatio >= 0.12);
          presenceAcc->add_read_support(itid->second, unique_ok);
        }
      }
    }
    for (const auto &entry : result.taxidCount) {
      fileInfo.taxidTotalMatches[entry.first] += 1;
    }

    bool isUniqueMapping = (result.taxidCount.size() == 1);
    if (isUniqueMapping) {
      fileInfo.uniqueTaxids.insert(result.taxidCount.front().first);
      const std::string &taxid = result.taxidCount.front().first;
      fileInfo.taxidUniqueMatches[taxid] += 1;
    }

    fileInfo.classifiedNum++;
    if (result.taxidCount.size() == 1) {
      if (!maxCountValid) {
        maxCount = result.taxidCount.front();
        maxCountValid = true;
      }
      result.taxidCount.clear();
      result.taxidCount.emplace_back(maxCount);
    }
  } else {
    fileInfo.unclassifiedNum++;
    result.taxidCount.emplace_back("unclassified", 1.0);
    if (result.reject_reason.empty()) {
      if (!maxCountValid && tidScore.empty()) {
        result.reject_reason = "no_candidate";
      } else {
        result.reject_reason = "low_eval";
      }
    }
    note_reject(result.reject_reason,
                result.best_taxid_hint.empty() ? bestTaxidStr
                                               : result.best_taxid_hint);
  }

  classifyResults.emplace_back(std::move(result));
}

void processBatch(
    batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    ClassifyConfig &config, chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<classifyResult> &classifyResults,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    FileInfo &fileInfo, GroupHeat &heat, const WeightingContext &weightCtx,
    PresenceAccumulator *presenceAcc, uint64_t decoySeed) {
  std::vector<uint64_t> hashs1;
  hashs1.reserve(2048);
  if (!batch.seqs2.empty()) {
    for (size_t i = 0; i < batch.ids.size(); ++i) {
      hashs1.clear();
      size_t len1 = (i < batch.seqs.size()) ? batch.seqs[i].size() : 0;
      size_t len2 = (i < batch.seqs2.size()) ? batch.seqs2[i].size() : 0;
      size_t readLen = len1 + len2;
      if (readLen > 0) {
        if (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength ||
            readLen < fileInfo.minLen)
          fileInfo.minLen = readLen;
        if (readLen > fileInfo.maxLen)
          fileInfo.maxLen = readLen;
        fileInfo.bpLength += readLen;
      }
      if (i < batch.seqs.size() && batch.seqs[i].size() >= feature_min_len) {
        chimera::feature::compute_hashes_append(batch.seqs[i], feature_params,
                                                hashs1);
      }
      if (i < batch.seqs2.size() && batch.seqs2[i].size() >= feature_min_len) {
        chimera::feature::compute_hashes_append(batch.seqs2[i], feature_params,
                                                hashs1);
      }
      if (hashs1.size() > 2048) {
        std::sort(hashs1.begin(), hashs1.end());
        hashs1.erase(std::unique(hashs1.begin(), hashs1.end()), hashs1.end());
      }
      processSequence(hashs1, readLen, imcfConfig, indexToTaxid, tax, config,
                      weightCtx, heat, imcf, batch.ids[i], classifyResults,
                      fileInfo, presenceAcc, decoySeed);
    }
  } else {
    for (size_t i = 0; i < batch.seqs.size(); i++) {
      hashs1.clear();
      size_t readLen = batch.seqs[i].size();
      if (readLen > 0) {
        if (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength ||
            readLen < fileInfo.minLen)
          fileInfo.minLen = readLen;
        if (readLen > fileInfo.maxLen)
          fileInfo.maxLen = readLen;
        fileInfo.bpLength += readLen;
      }
      if (batch.seqs[i].size() >= feature_min_len) {
        chimera::feature::compute_hashes_append(batch.seqs[i], feature_params,
                                                hashs1);
      }
      if (hashs1.size() > 2048) {
        std::sort(hashs1.begin(), hashs1.end());
        hashs1.erase(std::unique(hashs1.begin(), hashs1.end()), hashs1.end());
      }
      processSequence(hashs1, readLen, imcfConfig, indexToTaxid, tax, config,
                      weightCtx, heat, imcf, batch.ids[i], classifyResults,
                      fileInfo, presenceAcc, decoySeed);
    }
  }
}

void classify_streaming(
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
    ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    std::atomic<bool> &producer_done,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, PresenceSummary *presenceSummary,
    uint64_t decoySeed) {

#pragma omp parallel
  {
#ifdef _OPENMP
    const int thread_id = omp_get_thread_num();
#else
    const int thread_id = 0;
#endif
    moodycamel::ConcurrentQueue<batchReads> &readQueue =
        readQueues[static_cast<size_t>(
            std::clamp<int>(thread_id, 0,
                            static_cast<int>(readQueues.size() - 1)))];

    batchReads batch;
    std::vector<classifyResult> localClassifyResults;
    FileInfo localFileInfo;
    localFileInfo.minLen = kInvalidLength;
    localFileInfo.maxLen = 0;
    localFileInfo.bpLength = 0;
    GroupHeat heat;
    heat.ensure(indexToTaxid.size());
    PresenceAccumulator presenceLocal(
        presenceSummary ? presenceSummary->decoyReps : 0,
        presenceSummary ? presenceSummary->sketchBits : 0);
    PresenceAccumulator *presencePtr =
        presenceSummary ? &presenceLocal : nullptr;

    for (;;) {
      if (readQueue.try_dequeue(batch)) {
        processBatch(batch, imcfConfig, indexToTaxid, tax, config, imcf,
                     localClassifyResults, feature_params, feature_min_len,
                     localFileInfo, heat, weightCtx, presencePtr, decoySeed);
        continue;
      }
      if (producer_done.load(std::memory_order_acquire)) {
        break;
      }
      std::this_thread::yield();
    }

#pragma omp critical
    {
      classifyResults.insert(classifyResults.end(),
                             localClassifyResults.begin(),
                             localClassifyResults.end());

      fileInfo.classifiedNum += localFileInfo.classifiedNum;
      fileInfo.unclassifiedNum += localFileInfo.unclassifiedNum;

      fileInfo.uniqueTaxids.insert(localFileInfo.uniqueTaxids.begin(),
                                   localFileInfo.uniqueTaxids.end());

      for (const auto &[taxid, count] : localFileInfo.taxidTotalMatches) {
        fileInfo.taxidTotalMatches[taxid] += count;
      }

      for (const auto &[taxid, count] : localFileInfo.taxidUniqueMatches) {
        fileInfo.taxidUniqueMatches[taxid] += count;
      }
      size_t localMin =
          (localFileInfo.minLen == kInvalidLength) ? 0 : localFileInfo.minLen;
      if (localMin > 0 &&
          (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength ||
           localMin < fileInfo.minLen)) {
        fileInfo.minLen = localMin;
      }
      if (localFileInfo.maxLen > fileInfo.maxLen) {
        fileInfo.maxLen = localFileInfo.maxLen;
      }
      fileInfo.bpLength += localFileInfo.bpLength;

      fileInfo.preem_route_reads += localFileInfo.preem_route_reads;
      fileInfo.preem_route_rare_reads += localFileInfo.preem_route_rare_reads;
      fileInfo.preem_cap_checks += localFileInfo.preem_cap_checks;
      fileInfo.preem_cap_expanded += localFileInfo.preem_cap_expanded;
      fileInfo.preem_floor_checks += localFileInfo.preem_floor_checks;
      fileInfo.preem_floor_applied += localFileInfo.preem_floor_applied;
      fileInfo.preem_floor_added += localFileInfo.preem_floor_added;
      fileInfo.preem_floor_skipped_dominant +=
          localFileInfo.preem_floor_skipped_dominant;
      fileInfo.preem_floor_filtered_weak += localFileInfo.preem_floor_filtered_weak;
      fileInfo.preem_keepalive_attempt += localFileInfo.preem_keepalive_attempt;
      fileInfo.preem_keepalive_applied += localFileInfo.preem_keepalive_applied;
      fileInfo.preem_keepalive_blocked_low_ratio +=
          localFileInfo.preem_keepalive_blocked_low_ratio;
      fileInfo.preem_keepalive_blocked_low_gain +=
          localFileInfo.preem_keepalive_blocked_low_gain;
      fileInfo.preem_keepalive_blocked_low_abs +=
          localFileInfo.preem_keepalive_blocked_low_abs;
      if (presenceSummary) {
        presenceSummary->merge(presenceLocal);
      }
    }
  }
}

void classify(
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
    ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx) {

#pragma omp parallel
  {
#ifdef _OPENMP
    const int thread_id = omp_get_thread_num();
#else
    const int thread_id = 0;
#endif
    moodycamel::ConcurrentQueue<batchReads> &readQueue =
        readQueues[static_cast<size_t>(
            std::clamp<int>(thread_id, 0,
                            static_cast<int>(readQueues.size() - 1)))];

    batchReads batch;
    std::vector<classifyResult> localClassifyResults;
    FileInfo localFileInfo;
    localFileInfo.minLen = kInvalidLength;
    localFileInfo.maxLen = 0;
    localFileInfo.bpLength = 0;
    GroupHeat heat;
    heat.ensure(indexToTaxid.size());
    while (readQueue.try_dequeue(batch)) {
      processBatch(batch, imcfConfig, indexToTaxid, tax, config, imcf,
                   localClassifyResults, feature_params, feature_min_len,
                   localFileInfo, heat, weightCtx, nullptr, 0);
    }
#pragma omp critical
    {
      classifyResults.insert(classifyResults.end(),
                             localClassifyResults.begin(),
                             localClassifyResults.end());

      fileInfo.classifiedNum += localFileInfo.classifiedNum;
      fileInfo.unclassifiedNum += localFileInfo.unclassifiedNum;

      fileInfo.uniqueTaxids.insert(localFileInfo.uniqueTaxids.begin(),
                                   localFileInfo.uniqueTaxids.end());

      for (const auto &[taxid, count] : localFileInfo.taxidTotalMatches) {
        fileInfo.taxidTotalMatches[taxid] += count;
      }

      for (const auto &[taxid, count] : localFileInfo.taxidUniqueMatches) {
        fileInfo.taxidUniqueMatches[taxid] += count;
      }
      size_t localMin =
          (localFileInfo.minLen == kInvalidLength) ? 0 : localFileInfo.minLen;
      if (localMin > 0 &&
          (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength ||
           localMin < fileInfo.minLen)) {
        fileInfo.minLen = localMin;
      }
      if (localFileInfo.maxLen > fileInfo.maxLen) {
        fileInfo.maxLen = localFileInfo.maxLen;
      }
      fileInfo.bpLength += localFileInfo.bpLength;

      fileInfo.preem_route_reads += localFileInfo.preem_route_reads;
      fileInfo.preem_route_rare_reads += localFileInfo.preem_route_rare_reads;
      fileInfo.preem_cap_checks += localFileInfo.preem_cap_checks;
      fileInfo.preem_cap_expanded += localFileInfo.preem_cap_expanded;
      fileInfo.preem_floor_checks += localFileInfo.preem_floor_checks;
      fileInfo.preem_floor_applied += localFileInfo.preem_floor_applied;
      fileInfo.preem_floor_added += localFileInfo.preem_floor_added;
      fileInfo.preem_floor_skipped_dominant +=
          localFileInfo.preem_floor_skipped_dominant;
      fileInfo.preem_floor_filtered_weak += localFileInfo.preem_floor_filtered_weak;
      fileInfo.preem_keepalive_attempt += localFileInfo.preem_keepalive_attempt;
      fileInfo.preem_keepalive_applied += localFileInfo.preem_keepalive_applied;
      fileInfo.preem_keepalive_blocked_low_ratio +=
          localFileInfo.preem_keepalive_blocked_low_ratio;
      fileInfo.preem_keepalive_blocked_low_gain +=
          localFileInfo.preem_keepalive_blocked_low_gain;
      fileInfo.preem_keepalive_blocked_low_abs +=
          localFileInfo.preem_keepalive_blocked_low_abs;
      fileInfo.preem_beta_relax_checks += localFileInfo.preem_beta_relax_checks;
      fileInfo.preem_beta_relax_applied += localFileInfo.preem_beta_relax_applied;
      fileInfo.preem_beta_relax_thr_drop_sum +=
          localFileInfo.preem_beta_relax_thr_drop_sum;
      fileInfo.preem_beta_relax_suppressed_overflow +=
          localFileInfo.preem_beta_relax_suppressed_overflow;
    }
  }
}

} // namespace ChimeraClassify
