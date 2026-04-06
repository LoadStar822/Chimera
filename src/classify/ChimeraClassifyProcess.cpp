#include "ChimeraClassifyCommon.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <mutex>
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
    ChimeraBuild::IMCFConfig &imcfConfig, const TaxDict &tax,
    ClassifyConfig &config, const WeightingContext &weightCtx,
    const AutoClassifyPolicy &autoPolicy, GroupHeat &heat,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const std::string &id,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    PresenceAccumulator *presenceAcc, ProcessScratch &scratch) {
  size_t hashNum = hashs1.size();
  auto xor_reduce = [](const std::vector<uint64_t> &vals) {
    uint64_t acc = 0;
    for (uint64_t v : vals) {
      acc ^= v;
    }
    return acc;
  };
  const size_t binNumAll = tax.idx2id.size();
  heat.ensure(binNumAll);
  const double tail_risk_s = clamp01(config.tail_risk_s);
  const double tail_risk_hi = smoothstep(tail_risk_s, 0.60, 0.95);
  const bool presenceEnabled = (presenceAcc != nullptr);
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
  constexpr double kUniqueEdgeBonus = 3.0;
  double presence_weight = 1.0;

  const bool taxpool_enabled = autoPolicy.candidate.enable_taxpool;
  const char *dump_preem_env = std::getenv("CHIMERA_DUMP_PREEM_TSV");
  const bool dump_preem_enabled = (dump_preem_env && dump_preem_env[0] != '\0');
  const std::string dump_preem_path =
      dump_preem_enabled ? std::string(dump_preem_env) : std::string();
  auto map_species_rep = [&](uint32_t tid_id) -> uint32_t {
    if (weightCtx.tid2speciesRep && tid_id < weightCtx.tid2speciesRep->size()) {
      return (*weightCtx.tid2speciesRep)[tid_id];
    }
    return tid_id;
  };
  auto map_genus = [&](uint32_t tid_id) -> uint32_t {
    if (weightCtx.tid2genus && tid_id < weightCtx.tid2genus->size()) {
      return (*weightCtx.tid2genus)[tid_id];
    }
    return 0u;
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
      for (size_t i = 0; i < hashs1.size() && scored.size() < over;
           i += step) {
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
    std::vector<std::pair<uint32_t, uint64_t>> routeScored;
    routeScored.reserve(sampleVals.size());
    for (uint64_t v : sampleVals) {
      routeScored.emplace_back(weightCtx.freqSketch->estimate(v), v);
    }
    const double route_ratio = lerp(0.50, 0.50, tail_risk_hi);
    const size_t routeBudget =
        std::max<size_t>(1, static_cast<size_t>(
                                std::llround(route_ratio * sampleVals.size())));
    routeVals = select_rare_route_values(routeScored, routeBudget);
    if (routeVals.empty()) {
      routeVals = sampleVals;
    }
    used_rare_route = (routeVals.size() < sampleVals.size());
  }

  auto &sampleCount = scratch.sampleCount;
  sampleCount.clear();
  auto &touchedS = scratch.touchedS;
  touchedS.clear();
  touchedS.reserve(64);
  robin_hood::unordered_flat_map<uint32_t, uint32_t> sampleBinScore;
  robin_hood::unordered_flat_map<uint32_t, uint64_t> repCoarseScore;
  robin_hood::unordered_flat_map<uint32_t, uint64_t> genusCoarseScore;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> repRareSupport;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> genusRareSupport;
  auto &rankedBins = scratch.rankedBins;
  rankedBins.clear();
  uint64_t coarseTotal = 0;
  uint64_t sampleCovered = 0;
  auto &deferredEval = scratch.deferredEval;
  deferredEval.clear();
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
      if (bi < tax.idx2id.size() && sp < tax.idx2id[bi].size()) {
        const uint32_t tid_id = tax.idx2id[bi][sp];
        const uint32_t rep_id = map_species_rep(tid_id);
        repCoarseScore[rep_id] += static_cast<uint64_t>(contrib);
        const uint32_t genus_id = map_genus(rep_id);
        if (genus_id != 0u) {
          genusCoarseScore[genus_id] += static_cast<uint64_t>(contrib);
        }
      }
    }
  }

  const size_t baseCandidateCap = std::min<size_t>(
      binNumAll,
      static_cast<size_t>(std::llround(lerp(256.0, 320.0, tail_risk_hi))));
  const size_t maxCandidateCap = std::min<size_t>(
      binNumAll,
      static_cast<size_t>(std::llround(lerp(512.0, 640.0, tail_risk_hi))));
  size_t candidateCap = baseCandidateCap;
  constexpr double kHeadMassThresh = 0.5;
  auto &topBins = scratch.topBins;
  topBins.clear();
  bool fallback_full = (coarseTotal == 0);

  robin_hood::unordered_flat_set<uint32_t> candidateSet;
  candidateSet.reserve(256);
  robin_hood::unordered_flat_set<uint32_t> lowDegPreserve;
  lowDegPreserve.reserve(64);

  uint32_t rare_freq_threshold = 2u;
  if (weightCtx.freqSketch) {
    const uint32_t q = static_cast<uint32_t>(
        std::llround(std::max(1.0, weightCtx.freqQuantile * 0.25)));
    rare_freq_threshold = std::clamp<uint32_t>(q, 2u, 32u);
  }
  if (!routeVals.empty()) {
    auto &routed = scratch.routed;
    auto &rareRepHits = scratch.rareRepHits;
    auto &rareGenusHits = scratch.rareGenusHits;
    routed.reserve(16);
    rareRepHits.reserve(16);
    rareGenusHits.reserve(16);
    for (auto v : routeVals) {
      routed.clear();
      imcf.route(v, routed);
      if (routed.empty()) {
        continue;
      }
      const bool low_fanout = (routed.size() <= 2);
      const bool low_freq =
          (weightCtx.freqSketch &&
           weightCtx.freqSketch->estimate(v) <= rare_freq_threshold);
      if (low_fanout) {
        for (uint32_t b : routed) {
          if (b < binNumAll) {
            candidateSet.insert(b);
            lowDegPreserve.insert(b);
          }
        }
      }
      if (low_fanout || low_freq) {
        rareRepHits.clear();
        rareGenusHits.clear();
        imcf.bulkContain_events_subset(
            v, routed, [&](uint32_t bin, uint16_t sp) {
              if (bin >= tax.idx2id.size()) {
                return;
              }
              const auto &speciesVec = tax.idx2id[bin];
              if (sp >= speciesVec.size()) {
                return;
              }
              const uint32_t tid_id = speciesVec[sp];
              const uint32_t rep_id = map_species_rep(tid_id);
              rareRepHits.push_back(rep_id);
              const uint32_t genus_id = map_genus(rep_id);
              if (genus_id != 0u) {
                rareGenusHits.push_back(genus_id);
              }
            });
        if (!rareRepHits.empty()) {
          std::sort(rareRepHits.begin(), rareRepHits.end());
          rareRepHits.erase(std::unique(rareRepHits.begin(), rareRepHits.end()),
                            rareRepHits.end());
          for (uint32_t rep_id : rareRepHits) {
            repRareSupport[rep_id] += 1u;
          }
        }
        if (!rareGenusHits.empty()) {
          std::sort(rareGenusHits.begin(), rareGenusHits.end());
          rareGenusHits.erase(
              std::unique(rareGenusHits.begin(), rareGenusHits.end()),
              rareGenusHits.end());
          for (uint32_t genus_id : rareGenusHits) {
            genusRareSupport[genus_id] += 1u;
          }
        }
      }
    }
  }

  constexpr double coverageTarget = 0.92;
  if (coarseTotal > 0 && !sampleBinScore.empty()) {
    rankedBins.clear();
    rankedBins.reserve(sampleBinScore.size());
    for (const auto &kv : sampleBinScore) {
      rankedBins.emplace_back(kv.first, kv.second);
    }
    std::sort(rankedBins.begin(), rankedBins.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    const double head_mass = compute_head_mass(rankedBins, 10, coarseTotal);
    candidateCap = compute_candidate_cap(baseCandidateCap, maxCandidateCap,
                                         head_mass, kHeadMassThresh,
                                         rankedBins.size());

    uint64_t goal = static_cast<uint64_t>(
        std::ceil(static_cast<double>(coarseTotal) * coverageTarget));
    uint64_t covered = 0;
    for (const auto &[bin, score] : rankedBins) {
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
    auto &weighted = scratch.weighted;
    weighted.clear();
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

  bool taxpool_mode = false;
  const bool short_like_read = (readLen <= 600);
  const double taxpool_tail_gate = short_like_read ? 0.20 : 0.35;
  if (taxpool_enabled && !fallback_full && tail_risk_hi > taxpool_tail_gate &&
      coarseTotal > 0 && !repCoarseScore.empty() && !genusCoarseScore.empty() &&
      weightCtx.tid2genus) {
    auto &repRanked = scratch.repRanked;
    auto &genusRanked = scratch.genusRanked;
    repRanked.clear();
    genusRanked.clear();
    repRanked.reserve(repCoarseScore.size());
    genusRanked.reserve(genusCoarseScore.size());
    uint64_t repTotal = 0;
    for (const auto &kv : repCoarseScore) {
      if (kv.second == 0) {
        continue;
      }
      repRanked.emplace_back(kv.first, kv.second);
      repTotal += kv.second;
    }
    for (const auto &kv : genusCoarseScore) {
      if (kv.second == 0 || kv.first == 0u) {
        continue;
      }
      genusRanked.emplace_back(kv.first, kv.second);
    }
    std::sort(repRanked.begin(), repRanked.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    std::sort(genusRanked.begin(), genusRanked.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    if (!repRanked.empty() && !genusRanked.empty() && repTotal > 0) {
      const uint32_t dominant_genus = genusRanked.front().first;
      const double dominant_genus_share =
          static_cast<double>(genusRanked.front().second) /
          static_cast<double>(coarseTotal);
      const size_t repTopN = std::min<size_t>(10, repRanked.size());
      uint64_t repHeadSum = 0;
      for (size_t i = 0; i < repTopN; ++i) {
        repHeadSum += repRanked[i].second;
      }
      const double head_mass_rep =
          static_cast<double>(repHeadSum) / static_cast<double>(repTotal);
      const uint64_t top1_rep_score = repRanked.front().second;
      size_t dominant_crowded = 0;
      for (const auto &[rep_id, score] : repRanked) {
        if (map_genus(rep_id) != dominant_genus) {
          continue;
        }
        if (score * 10ull >= top1_rep_score * 4ull) {
          ++dominant_crowded;
        }
      }
      const bool short_like = short_like_read;
      taxpool_mode =
          (dominant_genus_share >= 0.72) ||
          (short_like && dominant_genus_share >= 0.60) ||
          (head_mass_rep >= 0.55 && dominant_genus_share >= 0.65) ||
          (dominant_crowded >= 3);
      if (taxpool_mode) {
        const size_t rep_pool_target = short_like ? 16u : 6u;
        const size_t dominant_limit = short_like ? 5u : 3u;
        const size_t non_dominant_limit = short_like ? 5u : 2u;
        auto &repPool = scratch.repPool;
        repPool.clear();
        repPool.reserve(rep_pool_target);
        robin_hood::unordered_flat_set<uint32_t> repSeen;
        repSeen.reserve(rep_pool_target * 2 + 4);
        auto add_rep = [&](uint32_t rep_id) -> bool {
          if (repSeen.find(rep_id) != repSeen.end()) {
            return false;
          }
          repSeen.insert(rep_id);
          repPool.push_back(rep_id);
          return true;
        };

        add_rep(repRanked.front().first);
        size_t dominant_added = 1;
        auto &dominantCandidates = scratch.dominantCandidates;
        dominantCandidates.clear();
        dominantCandidates.reserve(repRanked.size());
        for (const auto &[rep_id, _] : repRanked) {
          if (map_genus(rep_id) != dominant_genus) {
            continue;
          }
          if (rep_id == repRanked.front().first) {
            continue;
          }
          dominantCandidates.push_back(rep_id);
        }
        std::sort(dominantCandidates.begin(), dominantCandidates.end(),
                  [&](uint32_t a, uint32_t b) {
                    const uint32_t ar = repRareSupport.count(a) ? repRareSupport[a] : 0u;
                    const uint32_t br = repRareSupport.count(b) ? repRareSupport[b] : 0u;
                    if (ar != br) {
                      return ar > br;
                    }
                    const uint64_t as = repCoarseScore.count(a) ? repCoarseScore[a] : 0ull;
                    const uint64_t bs = repCoarseScore.count(b) ? repCoarseScore[b] : 0ull;
                    if (as != bs) {
                      return as > bs;
                    }
                    return a < b;
                  });
        for (uint32_t rep_id : dominantCandidates) {
          if (repPool.size() >= rep_pool_target || dominant_added >= dominant_limit) {
            break;
          }
          if (add_rep(rep_id)) {
            ++dominant_added;
          }
        }

        size_t non_dominant_added = 0;
        auto &nonDominantGenusPicks = scratch.nonDominantGenusPicks;
        nonDominantGenusPicks.clear();
        nonDominantGenusPicks.reserve(genusCoarseScore.size());
        for (const auto &[genus_id, coarse_score] : genusCoarseScore) {
          if (genus_id == 0u || genus_id == dominant_genus || coarse_score == 0u) {
            continue;
          }
          uint32_t best_rep = 0u;
          uint32_t best_rep_rare = 0u;
          uint64_t best_rep_coarse = 0u;
          bool have_best_rep = false;
          for (const auto &[rep_id, rep_score] : repCoarseScore) {
            if (rep_score == 0u || map_genus(rep_id) != genus_id) {
              continue;
            }
            const uint32_t rep_rare =
                repRareSupport.count(rep_id) ? repRareSupport[rep_id] : 0u;
            if (!have_best_rep || rep_rare > best_rep_rare ||
                (rep_rare == best_rep_rare && rep_score > best_rep_coarse) ||
                (rep_rare == best_rep_rare && rep_score == best_rep_coarse &&
                 rep_id < best_rep)) {
              have_best_rep = true;
              best_rep = rep_id;
              best_rep_rare = rep_rare;
              best_rep_coarse = rep_score;
            }
          }
          if (!have_best_rep) {
            continue;
          }
          const uint32_t g_rare =
              genusRareSupport.count(genus_id) ? genusRareSupport[genus_id] : 0u;
          nonDominantGenusPicks.push_back(
              {genus_id, best_rep, g_rare, coarse_score});
        }
        std::sort(nonDominantGenusPicks.begin(), nonDominantGenusPicks.end(),
                  [](const TaxpoolGenusPick &a, const TaxpoolGenusPick &b) {
                    if (a.rare_support != b.rare_support) {
                      return a.rare_support > b.rare_support;
                    }
                    if (a.coarse_score != b.coarse_score) {
                      return a.coarse_score > b.coarse_score;
                    }
                    if (a.genus_id != b.genus_id) {
                      return a.genus_id < b.genus_id;
                    }
                    return a.rep_id < b.rep_id;
                  });
        for (const auto &pick : nonDominantGenusPicks) {
          if (repPool.size() >= rep_pool_target ||
              non_dominant_added >= non_dominant_limit) {
            break;
          }
          if (add_rep(pick.rep_id)) {
            ++non_dominant_added;
          }
        }
        for (const auto &[rep_id, _] : repRanked) {
          if (repPool.size() >= rep_pool_target) {
            break;
          }
          add_rep(rep_id);
        }

        robin_hood::unordered_flat_set<uint32_t> taxpoolBins;
        taxpoolBins.reserve(512);
        for (uint32_t rep_id : repPool) {
          if (rep_id >= tax.tid2bin.size()) {
            continue;
          }
          for (uint32_t bin : tax.tid2bin[rep_id]) {
            if (bin < binNumAll) {
              taxpoolBins.insert(bin);
            }
          }
        }
        for (uint32_t bin : lowDegPreserve) {
          if (bin < binNumAll) {
            taxpoolBins.insert(bin);
          }
        }
        constexpr size_t kLegacyGlobalTopBins = 16;
        for (size_t i = 0; i < rankedBins.size() && i < kLegacyGlobalTopBins; ++i) {
          if (rankedBins[i].first < binNumAll) {
            taxpoolBins.insert(rankedBins[i].first);
          }
        }
        if (!taxpoolBins.empty()) {
          topBins.assign(taxpoolBins.begin(), taxpoolBins.end());
          std::sort(topBins.begin(), topBins.end());
          topBins.erase(std::unique(topBins.begin(), topBins.end()),
                        topBins.end());
          fallback_full = (topBins.size() == binNumAll);
        }
      }
    }
  }
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

  robin_hood::unordered_flat_map<uint32_t, double> tidScore;
  robin_hood::unordered_flat_map<uint32_t, double> uniqueHits;
  tidScore.reserve(128);
  uniqueHits.reserve(128);

  auto &minimizerTids = scratch.minimizerTids;
  minimizerTids.clear();
  minimizerTids.reserve(64);
  auto &minimizerBins = scratch.minimizerBins;
  minimizerBins.clear();
  minimizerBins.reserve(64);

  double eff_eval = 0.0;
  size_t n_eval = 0;

  const std::vector<uint32_t> *activeSubset =
      (fallback_full || topBins.size() == binNumAll) ? nullptr : &topBins;
  static thread_local std::vector<uint32_t> topBinMarks;
  static thread_local uint32_t topBinMarkEpoch = 0u;
  bool topBinMarksActive = false;
  auto rebuild_top_bin_marks = [&]() {
    topBinMarksActive = false;
    if (fallback_full || topBins.size() == binNumAll) {
      return;
    }
    if (topBinMarks.size() != binNumAll) {
      topBinMarks.assign(binNumAll, 0u);
      topBinMarkEpoch = 1u;
    } else {
      ++topBinMarkEpoch;
      if (topBinMarkEpoch == 0u) {
        std::fill(topBinMarks.begin(), topBinMarks.end(), 0u);
        topBinMarkEpoch = 1u;
      }
    }
    for (uint32_t bin : topBins) {
      if (bin < topBinMarks.size()) {
        topBinMarks[bin] = topBinMarkEpoch;
      }
    }
    topBinMarksActive = true;
  };
  rebuild_top_bin_marks();

  double idf_power = 1.0 + (0.25 * tail_risk_hi);

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
    };

    if (!subset) {
      imcf.bulkContain_events(value, emit);
    } else if (topBinMarksActive && subset == &topBins) {
      imcf.bulkContain_events_subset_marked(value, topBinMarks, topBinMarkEpoch,
                                            emit);
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
    const size_t df_eff = effective_df_bins(deg_effective, df_bins);

    double totalBins = subset ? static_cast<double>(subset->size())
                              : static_cast<double>(binNumAll);
    if (totalBins <= 0.0) {
      totalBins = static_cast<double>(binNumAll);
    }

    // Hybrid IDF DF selection:
    // - short reads/contigs: allow df_eff (species-fragmentation correction)
    // - long reads: keep df_bins to avoid FP sensitivity under subset/topBins
    constexpr size_t kDfEffIdfMaxLen = 4096;
    const double df_idf = df_for_idf(df_bins, df_eff, readLen, tail_risk_hi,
                                     kDfEffIdfMaxLen);
    double idf_raw = idf_raw_from_df_bins(totalBins, df_idf);
    double idf = clamp_idf(idf_raw, tail_risk_hi, config.idf_max, idf_power);

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
      freqFactor *= bg_idf;
      freqFactor = std::clamp(freqFactor, 0.05, 8.0);
    }

    const bool uniqueEdge = allow_unique_edge(
        deg_effective, df_bins, has_freq, weightCtx.freq_trusted, df_est,
        imcfConfig.presenceUniqueDeg);
    const bool fragmentedUnique =
        (!uniqueEdge && deg_effective == 1 && df_eff == 1 && df_bins > 1);
    const double unique_strength =
        uniqueEdge ? 1.0 : (fragmentedUnique ? tail_risk_hi : 0.0);
    const bool localUniqueEdge = is_local_unique_edge(deg_effective, df_bins);
    double denom = std::log2(2.0 + static_cast<double>(deg_effective));
    double weight = denom > 0.0 ? 1.0 / denom : 1.0;
    double base = weight * freqFactor;
    double bonus = 1.0;
    if (unique_strength > 0.0) {
      const double target_bonus =
          uniqueEdge ? kUniqueEdgeBonus : unique_edge_bonus(kUniqueEdgeBonus, df_bins);
      bonus = 1.0 + unique_strength * (target_bonus - 1.0);
    }
    double contrib = idf * base * bonus;

    for (uint32_t tid : minimizerTids) {
      tidScore[tid] += contrib;
    }
    if (presenceEnabled && !minimizerTids.empty()) {
      const double hit_weight = presence_weight;
      const double score_weight = hit_weight;
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
        presenceAcc->add_target(tid, hit_weight, score_weight, unique_strength,
                                localUniqueEdge, unique_bucket,
                                breadth_bucket);
      }
    }
    if (unique_strength > 0.0 && !minimizerTids.empty()) {
      uniqueHits[minimizerTids.front()] += unique_strength;
    }

    return contrib;
  };

  auto recompute_subset_state = [&]() {
    activeSubset = (fallback_full || topBins.size() == binNumAll)
                       ? nullptr
                       : &topBins;
    rebuild_top_bin_marks();
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
    double uniqueCount = 0.0;
    double uniqueRatio = 0.0;
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
    }
    double denom = eff_eval > 0.0 ? eff_eval : 1.0;
    stats.uniqueRatio = stats.uniqueCount / denom;
    return stats;
  };

  size_t n0 = std::min<size_t>(64, hashs1.size());
  for (size_t i = 0; i < n0; ++i) {
    const double c = evaluate_minimizer(hashs1[i], activeSubset);
    eff_eval += c;
  }
  n_eval = n0;

  auto meets_quick = [&](const EvidenceStats &s) {
    if (s.bestTid == std::numeric_limits<uint32_t>::max()) {
      return false;
    }
    double thr_conf_local = std::ceil(config.shotThreshold * eff_eval);
    double gap_need_local = std::max(0.5, eff_eval / 24.0);
    bool strong = (s.best >= thr_conf_local);
    bool unique_ok = (s.uniqueCount >= 3.0) ||
                     (s.uniqueCount >= 2.0 && s.uniqueRatio >= 0.12);
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
        deferredEval.push_back(i);
        continue;
      }
      const double c = evaluate_minimizer(hashs1[i], activeSubset);
      eff_eval += c;
      ++n_eval;
    }
    stats = collect_stats();
    highConfPre = meets_quick(stats);

    if (!highConfPre && !deferredEval.empty()) {
      for (size_t idx : deferredEval) {
        if (idx >= hashs1.size()) {
          continue;
        }
        const double c = evaluate_minimizer(hashs1[idx], activeSubset);
        eff_eval += c;
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
  uint32_t secondTid = stats.secondTid;
  double best = stats.best;
  double second = stats.second;
  double best_ratio = stats.ratio;
  double gap = stats.gap;
  double uniqueCount = stats.uniqueCount;
  double uniqueRatio = stats.uniqueRatio;
  std::string bestTaxidStr;
  if (bestTid != std::numeric_limits<uint32_t>::max() &&
      bestTid < tax.id2str.size()) {
    bestTaxidStr = tax.id2str[bestTid];
  }

  bool expanded = false;
  if (!taxpool_mode && !fallback_full && bestTid < tax.tid2bin.size()) {
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
    uniqueHits.clear();
    eff_eval = 0.0;
    n_eval = 0;
    for (size_t i = 0; i < hashs1.size(); ++i) {
      const double c = evaluate_minimizer(hashs1[i], activeSubset);
      eff_eval += c;
      ++n_eval;
    }
    stats = collect_stats();
    highConfPre = meets_quick(stats);
    thr_conf = std::ceil(config.shotThreshold * eff_eval);
    gap_need = std::max(0.5, eff_eval / 24.0);
    bestTid = stats.bestTid;
    secondTid = stats.secondTid;
    best = stats.best;
    second = stats.second;
    best_ratio = stats.ratio;
    gap = stats.gap;
    uniqueCount = stats.uniqueCount;
    uniqueRatio = stats.uniqueRatio;
  }


  // Optional: dump raw tidScore (pre-EM candidates before thresholds).

  size_t bestRounded =
      static_cast<size_t>(std::max<double>(0.0, std::llround(best)));
  size_t secondRounded =
      static_cast<size_t>(std::max<double>(0.0, std::llround(second)));
  size_t thrConfNeed = static_cast<size_t>(std::ceil(thr_conf));
  auto dc = decide_high_conf(bestRounded, secondRounded, eff_eval);
  bool marginAccept = (bestRounded >= thrConfNeed) && dc.accept;
  highConfPre = highConfPre && marginAccept;

  // --- Ambiguity guard before EM ---
  // When EM is enabled, send the read into EM unless Top1 is at least 3x Top2.
  // Lower ratios are easily hijacked by near-neighbor species and are better
  // resolved by the global EM pass.
  if (best_ratio < 3.0) {
    highConfPre = false;
  }
  // --------------------------------

  classifyResult result;
  result.evaluated = eff_eval;
  result.id = id;
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

  bool use_em = !highConfPre;

  double maxEvidence = std::min(best, eff_eval);
  double beta = config.firstFilterBeta;
  bool beta_user = config.firstFilterBeta_user;
  if (!(beta > 0.0)) {
    beta = 0.8;
  }
  if (!beta_user) {
    const double beta_em = lerp(0.50, 0.45, tail_risk_hi);
    const double beta_no_em = lerp(0.50, 0.80, tail_risk_hi);
    beta = use_em ? beta_em : beta_no_em;
  }
  beta = std::clamp(beta, 0.0, 1.0);
  size_t thr_beta =
      static_cast<size_t>(std::floor(beta * std::max(0.0, maxEvidence)));
  size_t thr_eval = static_cast<size_t>(
      std::ceil(config.shotThreshold * eff_eval));
  if (thr_eval == 0) {
    thr_eval = 1;
  }
  if (use_em) {
    double base = eff_eval;
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

  size_t thr_final = 0;
  if (use_em) {
    size_t thr_beta_eval = std::min(thr_beta, thr_eval);
    size_t min_hits = 1;
    thr_final = std::max({thr_beta_eval, min_hits, thr_min_eval});
  } else {
    thr_final = std::max({thr_eval, thr_beta, thr_min_eval});
  }

  const size_t thr_final_raw = thr_final;
  if (highConfPre && bestTid < tax.id2str.size()) {
    double bestEvidence = std::clamp(best, 0.0, effCap);
    const std::string &taxid = tax.id2str[bestTid];
    result.taxidCount.emplace_back(taxid, bestEvidence);
    maxCount = std::make_pair(taxid, bestEvidence);
    maxCountValid = true;

  } else {
    if (use_em && !tidScore.empty()) {
      auto &ranked = scratch.rankedTidScores;
      ranked.clear();
      ranked.reserve(tidScore.size());
      for (const auto &kv : tidScore) {
        ranked.emplace_back(kv.first, kv.second);
      }
      std::sort(ranked.begin(), ranked.end(),
                [](const auto &a, const auto &b) { return a.second > b.second; });

      size_t thr_min_eval_low = thr_min_eval;
      if (thr_min_eval > 0) {
        thr_min_eval_low = std::max<size_t>(
            1, static_cast<size_t>(
                   std::floor(0.3 * static_cast<double>(thr_min_eval))));
      }
      const double thr_low =
          static_cast<double>(std::max(thr_eval, thr_min_eval_low));
      const double thr_high = static_cast<double>(thr_final_raw);
      const double thr_base = lerp(thr_low, thr_high, tail_risk_hi);

      const size_t floorK_max = std::min<size_t>(128, ranked.size());
      const size_t floorCount = std::min<size_t>(
          floorK_max,
          static_cast<size_t>(std::llround((1.0 - tail_risk_hi) *
                                           static_cast<double>(floorK_max))));

      const uint32_t best_genus_id = map_genus(bestTid);

      if (dump_preem_enabled) {
        static std::mutex preem_dump_mu;
        static bool preem_dump_header_written = false;
        std::lock_guard<std::mutex> lock(preem_dump_mu);
        std::ofstream dump_os(dump_preem_path,
                              std::ios::out | std::ios::app);
        if (dump_os.is_open()) {
          if (!preem_dump_header_written) {
            dump_os
                << "read_id\ttid\trank\tscore\tthr_base\tfloor_count\tbest_tid\tbest_genus\tcand_genus\tbaseline_kept\n";
            preem_dump_header_written = true;
          }
          for (size_t i = 0; i < ranked.size(); ++i) {
            const uint32_t tid_id = ranked[i].first;
            if (tid_id >= tax.id2str.size()) {
              continue;
            }
            const double countVal = std::clamp(ranked[i].second, 0.0, effCap);
            if (countVal <= 0.0) {
              continue;
            }
            const bool baseline_kept = (i < floorCount || countVal >= thr_base);
            const uint32_t cand_genus_id = map_genus(tid_id);
            dump_os << id << '\t' << tax.id2str[tid_id] << '\t' << i << '\t'
                    << countVal << '\t' << thr_base << '\t' << floorCount << '\t'
                    << bestTaxidStr << '\t' << best_genus_id << '\t'
                    << cand_genus_id << '\t' << (baseline_kept ? 1 : 0)
                    << '\n';
          }
        }
      }

      result.taxidCount.reserve(ranked.size());
      for (size_t i = 0; i < ranked.size(); ++i) {
        const uint32_t tid_id = ranked[i].first;
        if (tid_id >= tax.id2str.size()) {
          continue;
        }
        const double countVal = std::clamp(ranked[i].second, 0.0, effCap);
        if (countVal <= 0.0) {
          continue;
        }
        const bool keep = (i < floorCount || countVal >= thr_base);
        if (keep) {
          const std::string &taxid = tax.id2str[tid_id];
          result.taxidCount.emplace_back(taxid, countVal);
          if (tid_id == bestTid) {
            maxCount = std::make_pair(taxid, countVal);
            maxCountValid = true;
          }
        }
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

    }
  }

  if (use_em && !highConfPre && result.taxidCount.size() == 1 &&
      !tidScore.empty()) {
    const bool unique_ok = (uniqueCount >= 3.0) ||
                           (uniqueCount >= 2.0 && uniqueRatio >= 0.12);
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
        auto &ranked = scratch.rankedTidScores;
        ranked.clear();
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

  if (!result.taxidCount.empty()) {
    if (presenceEnabled && presenceAcc) {
      const auto &top = result.taxidCount.front();
      if (top.first != "unclassified") {
        auto itid = tax.str2id.find(top.first);
        if (itid != tax.str2id.end()) {
          bool unique_ok = (uniqueCount >= 3.0) ||
                           (uniqueCount >= 2.0 && uniqueRatio >= 0.12);
          presenceAcc->add_read_support(itid->second, unique_ok);
        }
      }
    }
    bool isUniqueMapping = (result.taxidCount.size() == 1);
    if (isUniqueMapping) {
      fileInfo.uniqueTaxids.insert(result.taxidCount.front().first);
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
  }

  classifyResults.emplace_back(std::move(result));
}

void processBatch(
    batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
    const TaxDict &tax, ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<classifyResult> &classifyResults,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    FileInfo &fileInfo, GroupHeat &heat, const WeightingContext &weightCtx,
    const AutoClassifyPolicy &autoPolicy, PresenceAccumulator *presenceAcc,
    ProcessScratch &scratch) {
  auto &hashs1 = scratch.hashs1;
  hashs1.clear();
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
      processSequence(hashs1, readLen, imcfConfig, tax, config, weightCtx,
                      autoPolicy, heat, imcf, batch.ids[i], classifyResults,
                      fileInfo, presenceAcc, scratch);
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
      processSequence(hashs1, readLen, imcfConfig, tax, config, weightCtx,
                      autoPolicy, heat, imcf, batch.ids[i], classifyResults,
                      fileInfo, presenceAcc, scratch);
    }
  }
}

void classify_streaming(
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
    ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    std::atomic<bool> &producer_done,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, const AutoClassifyPolicy &autoPolicy,
    PresenceSummary *presenceSummary) {

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
    localFileInfo.avgLen = fileInfo.avgLen;
    localFileInfo.minLen = kInvalidLength;
    localFileInfo.maxLen = 0;
    localFileInfo.bpLength = 0;
    GroupHeat heat;
    heat.ensure(tax.idx2id.size());
    ProcessScratch scratch;
    PresenceAccumulator presenceLocal(
        presenceSummary ? presenceSummary->sketchBits : 0);
    PresenceAccumulator *presencePtr =
        presenceSummary ? &presenceLocal : nullptr;

    for (;;) {
      if (readQueue.try_dequeue(batch)) {
        processBatch(batch, imcfConfig, tax, config, imcf,
                     localClassifyResults, feature_params, feature_min_len,
                     localFileInfo, heat, weightCtx, autoPolicy, presencePtr,
                     scratch);
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
    std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
    ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, const AutoClassifyPolicy &autoPolicy) {

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
    localFileInfo.avgLen = fileInfo.avgLen;
    localFileInfo.minLen = kInvalidLength;
    localFileInfo.maxLen = 0;
    localFileInfo.bpLength = 0;
    GroupHeat heat;
    heat.ensure(tax.idx2id.size());
    ProcessScratch scratch;
    while (readQueue.try_dequeue(batch)) {
      processBatch(batch, imcfConfig, tax, config, imcf, localClassifyResults,
                   feature_params, feature_min_len, localFileInfo, heat,
                   weightCtx, autoPolicy, nullptr, scratch);
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
