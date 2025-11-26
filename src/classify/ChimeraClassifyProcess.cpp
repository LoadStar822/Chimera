#include "ChimeraClassifyCommon.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <thread>
#include <utility>

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
  const double exclusiveGamma = std::max(0.0, config.exclusive_gamma);

  auto note_reject = [&](const std::string &reason,
                         const std::string &taxid_hint) {
    fileInfo.rejectReasons[reason] += 1;
    if (!taxid_hint.empty()) {
      fileInfo.rejectByTaxid[taxid_hint][reason] += 1;
    }
  };

  size_t targetSample = hashNum / 4;
  targetSample = std::clamp<size_t>(targetSample, 16, 96);
  const size_t sampleBudget = std::min<size_t>(targetSample, hashs1.size());
  std::vector<uint64_t> sampleVals;
  if (sampleBudget > 0) {
    sampleVals.reserve(sampleBudget);
    size_t step = std::max<size_t>(1, hashs1.size() / sampleBudget);
    for (size_t i = 0; i < hashs1.size() && sampleVals.size() < sampleBudget;
         i += step) {
      sampleVals.push_back(hashs1[i]);
    }
    if (sampleVals.size() < sampleBudget && !hashs1.empty()) {
      sampleVals.push_back(hashs1.back());
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

  const size_t candidateCap =
      std::min<size_t>(binNumAll, static_cast<size_t>(256));
  std::vector<uint32_t> topBins;
  bool fallback_full = (coarseTotal == 0);

  robin_hood::unordered_flat_set<uint32_t> candidateSet;
  candidateSet.reserve(256);
  robin_hood::unordered_flat_set<uint32_t> lowDegPreserve;
  lowDegPreserve.reserve(64);

  if (!sampleVals.empty()) {
    std::vector<uint32_t> routed;
    routed.reserve(16);
    for (auto v : sampleVals) {
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
  std::vector<uint32_t> routedBinsBuf;
  routedBinsBuf.reserve(16);

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
      minimizerTids.push_back(tid);
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

    size_t deg = minimizerTids.size();
    if (deg == 0) {
      return 0.0;
    }
    double exclusivityWeight = 1.0;
    if (exclusiveGamma > 0.0 && deg > 0) {
      exclusivityWeight = std::pow(static_cast<double>(deg), -exclusiveGamma);
    }

    double totalBins = subset ? static_cast<double>(subset->size())
                              : static_cast<double>(binNumAll);
    if (totalBins <= 0.0) {
      totalBins = static_cast<double>(binNumAll);
    }

    routedBinsBuf.clear();
    imcf.route(value, routedBinsBuf);
    if (!routedBinsBuf.empty()) {
      std::sort(routedBinsBuf.begin(), routedBinsBuf.end());
      routedBinsBuf.erase(
          std::unique(routedBinsBuf.begin(), routedBinsBuf.end()),
          routedBinsBuf.end());
    }

    size_t df = routedBinsBuf.empty() ? deg : routedBinsBuf.size();
    double idf = std::log2((totalBins + 1.0) /
                           (static_cast<double>(df) + 1.0));
    idf = std::clamp(idf, 0.5, 5.0);

    double freqFactor = 1.0;
    if (weightCtx.enabled()) {
      uint32_t df_est = weightCtx.freqSketch->estimate(value);
      double vocab = static_cast<double>(
          std::max<uint64_t>(1, weightCtx.freqStats.nonzero_counters));
      double bg_idf =
          std::log((vocab + 1.0) / (static_cast<double>(df_est) + 1.0));
      bg_idf = bg_idf / std::log(2.0);
      bg_idf = std::clamp(bg_idf, 0.25, 6.0);
      if (weightCtx.freqStats.df_high_threshold !=
              std::numeric_limits<uint32_t>::max() &&
          df_est >= weightCtx.freqStats.df_high_threshold) {
        double damp =
            std::log2(static_cast<double>(weightCtx.freqStats.df_high_threshold + 1.0));
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

    double denom = std::log2(2.0 + static_cast<double>(deg));
    double weight = denom > 0.0 ? 1.0 / denom : 1.0;
    double contrib = idf * weight * freqFactor * exclusivityWeight;
    for (uint32_t tid : minimizerTids) {
      tidScore[tid] += contrib;
      ++consistencyHits[tid];
    }
    if (presenceEnabled && !minimizerTids.empty()) {
      bool uniqueEdge = (deg == 1);
      for (uint32_t tid : minimizerTids) {
        presenceAcc->add_target(tid, exclusivityWeight, uniqueEdge);
      }
      if (presenceDecoyReps > 0 && deg > 0) {
        for (size_t rep = 0; rep < presenceDecoyReps; ++rep) {
          uint64_t mix = splitmix64(
              value ^ (static_cast<uint64_t>(rep + 1) << 17) ^ decoySeed);
          size_t pick = static_cast<size_t>(mix % deg);
          uint32_t decoyTid = minimizerTids[pick];
          presenceAcc->add_decoy(rep, decoyTid, exclusivityWeight);
        }
      }
    }
    if (deg == 1 && !minimizerTids.empty()) {
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

  size_t bestRounded =
      static_cast<size_t>(std::max<double>(0.0, std::llround(best)));
  size_t secondRounded =
      static_cast<size_t>(std::max<double>(0.0, std::llround(second)));
  size_t thrConfNeed = static_cast<size_t>(std::ceil(thr_conf));
  auto dc = decide_high_conf(bestRounded, secondRounded, eff_eval);
  bool marginAccept = (bestRounded >= thrConfNeed) && dc.accept;
  highConfPre = highConfPre && marginAccept;

  classifyResult result;
  result.evaluated = eff_eval;
  result.id = id;
  result.best_taxid_hint = bestTaxidStr;
  std::pair<std::string, std::size_t> maxCount;
  bool maxCountValid = false;

  size_t effRounded =
      static_cast<size_t>(std::max<double>(1.0, std::llround(eff_eval)));
  if (!maxCountValid && bestRounded > 0 &&
      bestTid < tax.id2str.size()) {
    maxCount = std::make_pair(
        tax.id2str[bestTid], std::min(bestRounded, effRounded));
    maxCountValid = true;
  }

  bool use_em = config.em && !highConfPre;

  double maxEvidence = std::min(best, eff_eval);
  double beta = config.firstFilterBeta;
  if (beta <= 0.0) {
    beta = 0.8;
  }
  beta = std::clamp(beta, 0.0, 1.0);
  if (use_em && beta < 0.50) {
    beta = 0.50;
  }
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

  if (highConfPre && bestTid < tax.id2str.size()) {
    size_t bestRoundedDirect =
        static_cast<size_t>(std::max<double>(1.0, std::llround(best)));
    const std::string &taxid = tax.id2str[bestTid];
    result.taxidCount.emplace_back(taxid, bestRoundedDirect);
    maxCount = std::make_pair(taxid, bestRoundedDirect);
    maxCountValid = true;

  } else {
    for (const auto &[tid_id, rawScore] : tidScore) {
      size_t rounded =
          static_cast<size_t>(std::max<double>(0.0, std::llround(rawScore)));
      size_t countVal = std::min(rounded, effRounded);
      if (countVal >= thr_final) {
        const std::string &taxid = tax.id2str[tid_id];
        result.taxidCount.emplace_back(taxid, countVal);
        if (countVal == bestRounded) {
          maxCount = std::make_pair(taxid, countVal);
          maxCountValid = true;
        }
      }
    }
  }

  if (result.taxidCount.empty() && use_em && maxCountValid &&
      maxCount.second > 0) {
    result.taxidCount.emplace_back(maxCount);
  }

  size_t baseTopK =
      config.preEmTopK > 0 ? static_cast<size_t>(config.preEmTopK) : 32;
  size_t dynamicTopK = baseTopK;
  double tieGapNeed = std::max(2.0, eff_eval / 24.0);
  bool nearTie = (gap < tieGapNeed) || (best_ratio < 1.30);
  bool binOverflow = (!fallback_full && topBins.size() > 40) ||
                     (result.taxidCount.size() > baseTopK);
  if (nearTie || binOverflow) {
    dynamicTopK = std::max<size_t>(96, dynamicTopK);
  } else if (best_ratio >= 2.5 &&
             best >= thr_conf + std::max(1.0, eff_eval / 16.0) &&
             dynamicTopK > 8) {
    dynamicTopK = 16;
  }

  if (!result.taxidCount.empty() && use_em) {
    size_t K = dynamicTopK;
    if (K > 0 && result.taxidCount.size() > K) {
      std::nth_element(
          result.taxidCount.begin(), result.taxidCount.begin() + K,
          result.taxidCount.end(),
          [](const auto &a, const auto &b) { return a.second > b.second; });
      result.taxidCount.resize(K);
    }
  }

  if (!result.taxidCount.empty()) {
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
      result.taxidCount.clear();
      result.taxidCount.emplace_back(maxCount);
    }
  } else {
    fileInfo.unclassifiedNum++;
    result.taxidCount.emplace_back("unclassified", 1);
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
    moodycamel::ConcurrentQueue<batchReads> &readQueue, ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    std::atomic<bool> &producer_done,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, PresenceSummary *presenceSummary,
    uint64_t decoySeed) {

#pragma omp parallel
  {
    batchReads batch;
    std::vector<classifyResult> localClassifyResults;
    FileInfo localFileInfo;
    localFileInfo.minLen = kInvalidLength;
    localFileInfo.maxLen = 0;
    localFileInfo.bpLength = 0;
    GroupHeat heat;
    heat.ensure(indexToTaxid.size());
    PresenceAccumulator presenceLocal(
        presenceSummary ? presenceSummary->decoyReps : 0);
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
      if (presenceSummary) {
        presenceSummary->merge(presenceLocal);
      }
    }
  }
}

void classify(
    ChimeraBuild::IMCFConfig &imcfConfig,
    moodycamel::ConcurrentQueue<batchReads> &readQueue, ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx) {

#pragma omp parallel
  {
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
    }
  }
}

} // namespace ChimeraClassify
