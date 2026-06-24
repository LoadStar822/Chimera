#include "ChimeraClassifyCommon.hpp"

#include <utils/Parse.hpp>
#include <utils/NativeBoundedIndex.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

namespace ChimeraClassify {

static uint32_t species_rep_for_tid(const WeightingContext &weightCtx,
                                    uint32_t tid_id);
static uint32_t genus_for_tid(const WeightingContext &weightCtx,
                              uint32_t tid_id);

namespace {

static uint32_t profile_response_taxid_from_candidates(
    const std::vector<SpoolCandidate> &candidates, const TaxDict &tax,
    const WeightingContext &weightCtx) {
  if (candidates.empty()) {
    return 0;
  }
  const uint32_t tid_id = candidates.front().tid;
  if (tid_id == kSpoolUnclassifiedTid || tid_id >= tax.id2str.size()) {
    return 0;
  }
  uint32_t taxid = 0;
  if (!chimera::utils::try_parse_u32(tax.id2str[tid_id], taxid) ||
      taxid == 0) {
    return 0;
  }
  if (weightCtx.ncbiTaxdump != nullptr && weightCtx.ncbiTaxdump->enabled()) {
    taxid = weightCtx.ncbiTaxdump->to_species(taxid);
  }
  return taxid;
}

} // namespace

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

static uint64_t xor_reduce(const std::vector<uint64_t> &vals) {
  uint64_t acc = 0;
  for (uint64_t v : vals) {
    acc ^= v;
  }
  return acc;
}

static void ensure_epoch_vector(std::vector<uint32_t> &epochs, size_t size) {
  if (epochs.size() < size) {
    epochs.resize(size, 0);
  }
}

static void ensure_double_vector(std::vector<double> &values, size_t size) {
  if (values.size() < size) {
    values.resize(size, 0.0);
  }
}

static uint32_t next_epoch(std::vector<uint32_t> &epochs, uint32_t &epoch) {
  ++epoch;
  if (epoch == 0u) {
    std::fill(epochs.begin(), epochs.end(), 0u);
    epoch = 1u;
  }
  return epoch;
}

static uint32_t species_rep_for_tid(const WeightingContext &weightCtx,
                                    uint32_t tid_id) {
  if (weightCtx.tid2speciesRep && tid_id < weightCtx.tid2speciesRep->size()) {
    return (*weightCtx.tid2speciesRep)[tid_id];
  }
  return tid_id;
}

static uint32_t genus_for_tid(const WeightingContext &weightCtx,
                              uint32_t tid_id) {
  if (weightCtx.tid2genus && tid_id < weightCtx.tid2genus->size()) {
    return (*weightCtx.tid2genus)[tid_id];
  }
  return 0u;
}

static bool expand_adaptive_candidate_surface(
    size_t readLen, double dispersion_hi, uint64_t coarseTotal,
    robin_hood::unordered_flat_map<uint32_t, uint64_t> &repCoarseScore,
    robin_hood::unordered_flat_map<uint32_t, uint64_t> &genusCoarseScore,
    robin_hood::unordered_flat_map<uint32_t, uint32_t> &repRareSupport,
    robin_hood::unordered_flat_map<uint32_t, uint32_t> &genusRareSupport,
    const robin_hood::unordered_flat_set<uint32_t> &lowDegPreserve,
    const std::vector<std::pair<uint32_t, uint32_t>> &rankedBins,
    const TaxDict &tax, const WeightingContext &weightCtx, size_t binNumAll,
    bool &full_surface_mode, std::vector<uint32_t> &topBins,
    ProcessScratch &scratch) {
  // Adaptive candidate surface expansion keeps genus-structured evidence in
  // the candidate set before readout. It is not a readout-level genus gate.
  bool adaptive_surface_mode = false;
  const bool short_like_read = (readLen <= 600);
  const double adaptive_surface_dispersion_gate =
      short_like_read ? 0.20 : 0.35;
  if (!full_surface_mode && dispersion_hi > adaptive_surface_dispersion_gate &&
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
        if (genus_for_tid(weightCtx, rep_id) != dominant_genus) {
          continue;
        }
        if (score * 10ull >= top1_rep_score * 4ull) {
          ++dominant_crowded;
        }
      }
      const bool short_like = short_like_read;
      adaptive_surface_mode =
          (dominant_genus_share >= 0.72) ||
          (short_like && dominant_genus_share >= 0.60) ||
          (head_mass_rep >= 0.55 && dominant_genus_share >= 0.65) ||
          (dominant_crowded >= 3);
      if (adaptive_surface_mode) {
        const size_t rep_surface_target = short_like ? 16u : 6u;
        const size_t dominant_limit = short_like ? 5u : 3u;
        const size_t non_dominant_limit = short_like ? 5u : 2u;
        auto &repPool = scratch.repPool;
        repPool.clear();
        repPool.reserve(rep_surface_target);
        robin_hood::unordered_flat_set<uint32_t> repSeen;
        repSeen.reserve(rep_surface_target * 2 + 4);
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
          if (genus_for_tid(weightCtx, rep_id) != dominant_genus) {
            continue;
          }
          if (rep_id == repRanked.front().first) {
            continue;
          }
          dominantCandidates.push_back(rep_id);
        }
        std::sort(dominantCandidates.begin(), dominantCandidates.end(),
                  [&](uint32_t a, uint32_t b) {
                    const uint32_t ar =
                        repRareSupport.count(a) ? repRareSupport[a] : 0u;
                    const uint32_t br =
                        repRareSupport.count(b) ? repRareSupport[b] : 0u;
                    if (ar != br) {
                      return ar > br;
                    }
                    const uint64_t as =
                        repCoarseScore.count(a) ? repCoarseScore[a] : 0ull;
                    const uint64_t bs =
                        repCoarseScore.count(b) ? repCoarseScore[b] : 0ull;
                    if (as != bs) {
                      return as > bs;
                    }
                    return a < b;
                  });
        for (uint32_t rep_id : dominantCandidates) {
          if (repPool.size() >= rep_surface_target ||
              dominant_added >= dominant_limit) {
            break;
          }
          if (add_rep(rep_id)) {
            ++dominant_added;
          }
        }

        size_t non_dominant_added = 0;
        auto &secondaryGenusPicks = scratch.secondaryGenusPicks;
        secondaryGenusPicks.clear();
        secondaryGenusPicks.reserve(genusCoarseScore.size());
        for (const auto &[genus_id, coarse_score] : genusCoarseScore) {
          if (genus_id == 0u || genus_id == dominant_genus ||
              coarse_score == 0u) {
            continue;
          }
          uint32_t best_rep = 0u;
          uint32_t best_rep_rare = 0u;
          uint64_t best_rep_coarse = 0u;
          bool have_best_rep = false;
          for (const auto &[rep_id, rep_score] : repCoarseScore) {
            if (rep_score == 0u ||
                genus_for_tid(weightCtx, rep_id) != genus_id) {
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
              genusRareSupport.count(genus_id) ? genusRareSupport[genus_id]
                                               : 0u;
          secondaryGenusPicks.push_back(
              {genus_id, best_rep, g_rare, coarse_score});
        }
        std::sort(secondaryGenusPicks.begin(), secondaryGenusPicks.end(),
                  [](const AdaptiveSurfaceGenusPick &a,
                     const AdaptiveSurfaceGenusPick &b) {
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
        for (const auto &pick : secondaryGenusPicks) {
          if (repPool.size() >= rep_surface_target ||
              non_dominant_added >= non_dominant_limit) {
            break;
          }
          if (add_rep(pick.rep_id)) {
            ++non_dominant_added;
          }
        }
        for (const auto &[rep_id, _] : repRanked) {
          if (repPool.size() >= rep_surface_target) {
            break;
          }
          add_rep(rep_id);
        }

        robin_hood::unordered_flat_set<uint32_t> candidateSurfaceBins;
        candidateSurfaceBins.reserve(512);
        for (uint32_t rep_id : repPool) {
          if (rep_id >= tax.tid2bin.size()) {
            continue;
          }
          for (uint32_t bin : tax.tid2bin[rep_id]) {
            if (bin < binNumAll) {
              candidateSurfaceBins.insert(bin);
            }
          }
        }
        for (uint32_t bin : lowDegPreserve) {
          if (bin < binNumAll) {
            candidateSurfaceBins.insert(bin);
          }
        }
        constexpr size_t kGlobalEvidenceAnchorBins = 16;
        for (size_t i = 0;
             i < rankedBins.size() && i < kGlobalEvidenceAnchorBins; ++i) {
          if (rankedBins[i].first < binNumAll) {
            candidateSurfaceBins.insert(rankedBins[i].first);
          }
        }
        if (!candidateSurfaceBins.empty()) {
          topBins.assign(candidateSurfaceBins.begin(),
                         candidateSurfaceBins.end());
          std::sort(topBins.begin(), topBins.end());
          topBins.erase(std::unique(topBins.begin(), topBins.end()),
                        topBins.end());
          full_surface_mode = (topBins.size() == binNumAll);
        }
      }
    }
  }
  return adaptive_surface_mode;
}

static std::vector<SpoolCandidate>
capture_ranked_score_candidates(const ProcessScratch &scratch,
                                size_t taxid_count) {
  std::vector<SpoolCandidate> ranked;
  ranked.reserve(scratch.activeTidScores.size());
  for (uint32_t tid : scratch.activeTidScores) {
    if (tid < taxid_count) {
      const double score = scratch.tidScoreDense[tid];
      if (score > 0.0) {
        ranked.push_back(SpoolCandidate{tid, score, score});
      }
    }
  }
  std::sort(ranked.begin(), ranked.end(), [](const auto &a, const auto &b) {
    if (a.score != b.score) {
      return a.score > b.score;
    }
    return a.tid < b.tid;
  });
  return ranked;
}

static std::vector<SpoolCandidate> capture_abundance_candidates_from_scores(
    const ProcessScratch &scratch, size_t taxid_count, double eff_eval,
    size_t n_eval) {
  std::vector<SpoolCandidate> out;
  auto ranked = capture_ranked_score_candidates(scratch, taxid_count);
  if (ranked.empty()) {
    return out;
  }

  const double evidenceN =
      std::max({2.0, eff_eval, static_cast<double>(n_eval)});
  const double candidateCost =
      (0.5 * std::log(evidenceN)) +
      std::log(std::max<double>(1.0, static_cast<double>(ranked.size())));
  out.reserve(ranked.size());
  out.push_back(ranked.front());
  for (size_t i = 1; i < ranked.size(); ++i) {
    const auto &cand = ranked[i];
    if (cand.tid >= taxid_count || !(cand.score > 0.0)) {
      continue;
    }
    if (cand.score < candidateCost) {
      break;
    }
    out.push_back(cand);
  }
  return out;
}

static std::vector<SpoolCandidate>
capture_sample_mixture_candidates_from_base_scores(const ProcessScratch &scratch,
                                                   size_t taxid_count) {
  std::vector<SpoolCandidate> out;
  out.reserve(scratch.activeTidScores.size());
  for (uint32_t tid : scratch.activeTidScores) {
    if (tid >= taxid_count) {
      continue;
    }
    const double baseScore = scratch.tidBaseScoreDense[tid];
    if (!(baseScore > 0.0)) {
      continue;
    }
    out.push_back(SpoolCandidate{tid, baseScore, baseScore});
  }
  if (out.empty()) {
    return out;
  }
  std::sort(out.begin(), out.end(), [](const auto &a, const auto &b) {
    if (a.score != b.score) {
      return a.score > b.score;
    }
    return a.tid < b.tid;
  });
  double total = 0.0;
  double square = 0.0;
  for (const auto &cand : out) {
    total += cand.score;
    square += cand.score * cand.score;
  }
  if (total > 0.0 && square > 0.0) {
    const double eff = (total * total) / square;
    const double keepD = std::ceil(eff * eff);
    if (std::isfinite(keepD) && keepD >= 1.0 &&
        keepD < static_cast<double>(out.size())) {
      out.resize(static_cast<size_t>(keepD));
    }
  }
  return out;
}

static void prepare_hash_sample_and_route(
    const std::vector<uint64_t> &hashs, const ClassifyConfig &config,
    const WeightingContext &weightCtx, ProcessScratch &scratch) {
  size_t targetSample = hashs.size() / 2;
  targetSample = std::clamp<size_t>(
      targetSample, config.hash_sample_min, config.hash_sample_max);
  const size_t sampleBudget = std::min<size_t>(targetSample, hashs.size());

  auto &sampleVals = scratch.sampleVals;
  sampleVals.clear();
  if (sampleBudget > 0) {
    sampleVals.reserve(sampleBudget);
    if (hashs.size() <= sampleBudget) {
      sampleVals = hashs;
    } else if (weightCtx.freqSketch) {
      // Mix rare-biased hashes with uniform coverage so coarse routing stays
      // sharp without losing genus-level candidates.
      const size_t rareQuota = std::max<size_t>(1, sampleBudget / 2);
      const size_t over = std::min<size_t>(
          hashs.size(), std::max<size_t>(sampleBudget * 8, sampleBudget));
      auto &scored = scratch.sampleScored;
      scored.clear();
      scored.reserve(over);
      const size_t step = std::max<size_t>(1, hashs.size() / over);
      for (size_t i = 0; i < hashs.size() && scored.size() < over; i += step) {
        const uint64_t v = hashs[i];
        scored.emplace_back(weightCtx.freqSketch->estimate(v), v);
      }
      if (!hashs.empty() && scored.size() < over) {
        const uint64_t v = hashs.back();
        scored.emplace_back(weightCtx.freqSketch->estimate(v), v);
      }
      std::sort(scored.begin(), scored.end(), [](const auto &a,
                                                 const auto &b) {
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

      const size_t step2 = std::max<size_t>(1, hashs.size() / sampleBudget);
      for (size_t i = 0; i < hashs.size() && sampleVals.size() < sampleBudget;
           i += step2) {
        sampleVals.push_back(hashs[i]);
      }
      if (sampleVals.size() < sampleBudget && !hashs.empty()) {
        sampleVals.push_back(hashs.back());
      }
      std::sort(sampleVals.begin(), sampleVals.end());
      sampleVals.erase(std::unique(sampleVals.begin(), sampleVals.end()),
                       sampleVals.end());
      if (sampleVals.size() > sampleBudget) {
        sampleVals.resize(sampleBudget);
      }
    } else {
      const size_t step = std::max<size_t>(1, hashs.size() / sampleBudget);
      for (size_t i = 0; i < hashs.size() && sampleVals.size() < sampleBudget;
           i += step) {
        sampleVals.push_back(hashs[i]);
      }
      if (sampleVals.size() < sampleBudget && !hashs.empty()) {
        sampleVals.push_back(hashs.back());
      }
    }
  }

  auto &routeVals = scratch.routeVals;
  routeVals = sampleVals;
  if (!sampleVals.empty() && weightCtx.freqSketch) {
    auto &routeScored = scratch.routeScored;
    routeScored.clear();
    routeScored.reserve(sampleVals.size());
    for (uint64_t v : sampleVals) {
      routeScored.emplace_back(weightCtx.freqSketch->estimate(v), v);
    }
    constexpr double route_ratio = 0.50;
    const size_t routeBudget = std::max<size_t>(
        1, static_cast<size_t>(std::llround(route_ratio * sampleVals.size())));
    routeVals = select_rare_route_values(routeScored, routeBudget);
    if (routeVals.empty()) {
      routeVals = sampleVals;
    }
  }

}

static uint64_t collect_coarse_evidence(
    const std::vector<uint64_t> &sampleVals,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const TaxDict &tax,
    const WeightingContext &weightCtx, ProcessScratch &scratch,
    robin_hood::unordered_flat_map<uint32_t, uint32_t> &sampleBinScore,
    robin_hood::unordered_flat_map<uint32_t, uint64_t> &repCoarseScore,
    robin_hood::unordered_flat_map<uint32_t, uint64_t> &genusCoarseScore) {
  using PairCode = uint64_t;
  auto &coarsePairCodes = scratch.coarsePairCodes;
  coarsePairCodes.clear();
  coarsePairCodes.reserve(64);
  auto &coarsePairCounts = scratch.coarsePairCounts;
  coarsePairCounts.clear();
  auto &coarseTouchedPairs = scratch.coarseTouchedPairs;
  coarseTouchedPairs.clear();

  if (sampleVals.empty()) {
    return 0;
  }

  for (uint64_t value : sampleVals) {
    coarsePairCodes.clear();
    imcf.bulkContain_events(value, [&](uint32_t bin, uint16_t sp) {
      coarsePairCodes.push_back((PairCode(bin) << 16) | PairCode(sp));
    });
    if (coarsePairCodes.empty()) {
      continue;
    }
    std::sort(coarsePairCodes.begin(), coarsePairCodes.end());
    coarsePairCodes.erase(
        std::unique(coarsePairCodes.begin(), coarsePairCodes.end()),
        coarsePairCodes.end());
    for (PairCode code : coarsePairCodes) {
      auto [it, inserted] = coarsePairCounts.try_emplace(code, 0u);
      if (inserted) {
        coarseTouchedPairs.push_back(code);
      }
      if (it->second < std::numeric_limits<uint32_t>::max()) {
        ++it->second;
      }
    }
  }

  uint64_t coarseTotal = 0;
  sampleBinScore.reserve(coarseTouchedPairs.size());
  for (PairCode code : coarseTouchedPairs) {
    auto countIt = coarsePairCounts.find(code);
    if (countIt == coarsePairCounts.end() || countIt->second == 0u) {
      continue;
    }
    const uint32_t bi = static_cast<uint32_t>(code >> 16);
    const uint16_t sp = static_cast<uint16_t>(code & 0xFFFFu);
    const uint32_t contrib = countIt->second;
    coarseTotal += contrib;
    sampleBinScore[bi] += contrib;
    if (bi < tax.idx2id.size() && sp < tax.idx2id[bi].size()) {
      const uint32_t tid_id = tax.idx2id[bi][sp];
      const uint32_t rep_id = species_rep_for_tid(weightCtx, tid_id);
      repCoarseScore[rep_id] += static_cast<uint64_t>(contrib);
      const uint32_t genus_id = genus_for_tid(weightCtx, rep_id);
      if (genus_id != 0u) {
        genusCoarseScore[genus_id] += static_cast<uint64_t>(contrib);
      }
    }
  }
  return coarseTotal;
}

static bool build_initial_candidate_surface(
    const std::vector<uint64_t> &routeVals,
    const robin_hood::unordered_flat_map<uint32_t, uint32_t> &sampleBinScore,
    uint64_t coarseTotal, double dispersion_hi, size_t binNumAll,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const TaxDict &tax,
    const WeightingContext &weightCtx, const GroupHeat &heat,
    ProcessScratch &scratch,
    robin_hood::unordered_flat_map<uint32_t, uint32_t> &repRareSupport,
    robin_hood::unordered_flat_map<uint32_t, uint32_t> &genusRareSupport,
    robin_hood::unordered_flat_set<uint32_t> &lowDegPreserve) {
  const size_t baseCandidateCap = std::min<size_t>(
      binNumAll,
      static_cast<size_t>(std::llround(lerp(256.0, 320.0, dispersion_hi))));
  const size_t maxCandidateCap = std::min<size_t>(
      binNumAll,
      static_cast<size_t>(std::llround(lerp(512.0, 640.0, dispersion_hi))));
  size_t candidateCap = baseCandidateCap;

  auto &topBins = scratch.topBins;
  topBins.clear();
  bool full_surface_mode = (coarseTotal == 0);

  robin_hood::unordered_flat_set<uint32_t> candidateSet;
  candidateSet.reserve(256);
  lowDegPreserve.clear();
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
              const uint32_t rep_id = species_rep_for_tid(weightCtx, tid_id);
              rareRepHits.push_back(rep_id);
              const uint32_t genus_id = genus_for_tid(weightCtx, rep_id);
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
  auto &rankedBins = scratch.rankedBins;
  if (coarseTotal > 0 && !sampleBinScore.empty()) {
    rankedBins.clear();
    rankedBins.reserve(sampleBinScore.size());
    for (const auto &kv : sampleBinScore) {
      rankedBins.emplace_back(kv.first, kv.second);
    }
    std::sort(rankedBins.begin(), rankedBins.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    const double core_mass =
        compute_simpson_core_mass(rankedBins, coarseTotal);
    candidateCap = compute_candidate_cap(baseCandidateCap, maxCandidateCap,
                                         core_mass, rankedBins.size());

    const uint64_t goal = static_cast<uint64_t>(
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
    full_surface_mode = false;
  }

  if (!full_surface_mode) {
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
    full_surface_mode = true;
    topBins.resize(binNumAll);
    std::iota(topBins.begin(), topBins.end(), 0u);
  }

  if (!full_surface_mode) {
    std::sort(topBins.begin(), topBins.end());
    topBins.erase(std::unique(topBins.begin(), topBins.end()), topBins.end());
  }
  if (full_surface_mode) {
    scratch.baseTopBins.clear();
  } else {
    scratch.baseTopBins = topBins;
  }

  return full_surface_mode;
}

static void dump_preem_candidate_retention(
    const std::string &dump_path, const std::string &read_id,
    const TaxDict &tax, const WeightingContext &weightCtx,
    const std::vector<std::pair<uint32_t, double>> &ranked, double effCap,
    double thr_base, size_t floorCount, const std::string &bestTaxidStr,
    uint32_t best_genus_id) {
  static std::mutex preem_dump_mu;
  static bool preem_dump_header_written = false;
  auto map_genus = [&](uint32_t tid_id) -> uint32_t {
    if (weightCtx.tid2genus && tid_id < weightCtx.tid2genus->size()) {
      return (*weightCtx.tid2genus)[tid_id];
    }
    return 0u;
  };

  std::lock_guard<std::mutex> lock(preem_dump_mu);
  std::ofstream dump_os(dump_path, std::ios::out | std::ios::app);
  if (!dump_os.is_open()) {
    return;
  }
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
    dump_os << read_id << '\t' << tax.id2str[tid_id] << '\t' << i << '\t'
            << countVal << '\t' << thr_base << '\t' << floorCount << '\t'
            << bestTaxidStr << '\t' << best_genus_id << '\t' << cand_genus_id
            << '\t' << (baseline_kept ? 1 : 0) << '\n';
  }
}

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

static double score_ratio(double best_score, double second_score) {
  if (second_score > 0.0) {
    return best_score /
           std::max(second_score, std::numeric_limits<double>::min());
  }
  return best_score <= 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
}

static EvidenceStats collect_evidence_stats(
    const robin_hood::unordered_flat_map<uint32_t, double> &tidScore,
    const robin_hood::unordered_flat_map<uint32_t, double> &uniqueHits,
    double eff_eval) {
  EvidenceStats stats;
  for (const auto &kv : tidScore) {
    const double c = kv.second;
    if (c > stats.best) {
      stats.second = stats.best;
      stats.secondTid = stats.bestTid;
      stats.best = c;
      stats.bestTid = kv.first;
    } else if (c > stats.second) {
      stats.second = c;
      stats.secondTid = kv.first;
    }
  }
  stats.ratio = score_ratio(stats.best, stats.second);
  stats.gap = stats.best - stats.second;
  if (stats.bestTid != std::numeric_limits<uint32_t>::max()) {
    if (auto it = uniqueHits.find(stats.bestTid); it != uniqueHits.end()) {
      stats.uniqueCount = it->second;
    }
  }
  const double denom = eff_eval > 0.0 ? eff_eval : 1.0;
  stats.uniqueRatio = stats.uniqueCount / denom;
  return stats;
}

static bool meets_quick_high_conf(const ClassifyConfig &config,
                                  const EvidenceStats &stats,
                                  double eff_eval) {
  if (stats.bestTid == std::numeric_limits<uint32_t>::max()) {
    return false;
  }
  const double thr_conf_local = std::ceil(config.shotThreshold * eff_eval);
  const double gap_need_local = std::max(0.5, eff_eval / 24.0);
  const bool strong = (stats.best >= thr_conf_local);
  const bool unique_ok = (stats.uniqueCount >= 3.0) ||
                         (stats.uniqueCount >= 2.0 &&
                          stats.uniqueRatio >= 0.12);
  const bool stable =
      (stats.gap >= gap_need_local) && (stats.ratio >= 1.35);
  const size_t bestRoundedLocal = static_cast<size_t>(
      std::max<double>(0.0, std::llround(stats.best)));
  const size_t secondRoundedLocal = static_cast<size_t>(
      std::max<double>(0.0, std::llround(stats.second)));
  const size_t needLocal = static_cast<size_t>(std::ceil(thr_conf_local));
  const auto dc =
      decide_high_conf(bestRoundedLocal, secondRoundedLocal, eff_eval);
  const bool margin_ok = (bestRoundedLocal >= needLocal) && dc.accept;
  return strong && unique_ok && stable && margin_ok;
}

struct ReadScoringState {
  ProcessScratch &scratch;
  const WeightingContext &weightCtx;
  size_t tid_count{0};
  bool adaptive_surface_mode{false};
  uint32_t tid_score_epoch{0};
  uint32_t unique_hits_epoch{0};
  robin_hood::unordered_flat_map<uint32_t, double> tid_score;
  robin_hood::unordered_flat_map<uint32_t, double> unique_hits;

  ReadScoringState(ProcessScratch &scratch_in,
                   const WeightingContext &weightCtx_in, size_t tid_count_in,
                   bool adaptive_surface_mode_in)
      : scratch(scratch_in), weightCtx(weightCtx_in), tid_count(tid_count_in),
        adaptive_surface_mode(adaptive_surface_mode_in) {
    ensure_epoch_vector(scratch.tidScoreEpoch, tid_count);
    ensure_double_vector(scratch.tidScoreDense, tid_count);
    ensure_double_vector(scratch.tidBaseScoreDense, tid_count);
    ensure_double_vector(scratch.tidCompletionScoreDense, tid_count);
    ensure_epoch_vector(scratch.uniqueHitsEpoch, tid_count);
    ensure_double_vector(scratch.uniqueHitsDense, tid_count);
    scratch.activeTidScores.clear();
    scratch.activeUniqueHits.clear();
    tid_score.reserve(128);
    unique_hits.reserve(128);
    tid_score_epoch =
        next_epoch(scratch.tidScoreEpoch, scratch.tidScoreEpochValue);
    unique_hits_epoch =
        next_epoch(scratch.uniqueHitsEpoch, scratch.uniqueHitsEpochValue);
  }

  void add_tid_score(uint32_t tid, double contrib, bool baseSource) {
    if (tid >= tid_count) {
      return;
    }
    if (scratch.tidScoreEpoch[tid] != tid_score_epoch) {
      scratch.tidScoreEpoch[tid] = tid_score_epoch;
      scratch.tidScoreDense[tid] = 0.0;
      scratch.tidBaseScoreDense[tid] = 0.0;
      scratch.tidCompletionScoreDense[tid] = 0.0;
      scratch.activeTidScores.push_back(tid);
    }
    scratch.tidScoreDense[tid] += contrib;
    if (baseSource) {
      scratch.tidBaseScoreDense[tid] += contrib;
    } else {
      scratch.tidCompletionScoreDense[tid] += contrib;
    }
  }

  void add_unique_hit(uint32_t tid, double value) {
    if (tid >= tid_count) {
      return;
    }
    if (scratch.uniqueHitsEpoch[tid] != unique_hits_epoch) {
      scratch.uniqueHitsEpoch[tid] = unique_hits_epoch;
      scratch.uniqueHitsDense[tid] = 0.0;
      scratch.activeUniqueHits.push_back(tid);
    }
    scratch.uniqueHitsDense[tid] += value;
  }

  void reset_dense_scores() {
    tid_score.clear();
    unique_hits.clear();
    scratch.activeTidScores.clear();
    scratch.activeUniqueHits.clear();
    tid_score_epoch =
        next_epoch(scratch.tidScoreEpoch, scratch.tidScoreEpochValue);
    unique_hits_epoch =
        next_epoch(scratch.uniqueHitsEpoch, scratch.uniqueHitsEpochValue);
  }

  void materialize_score_maps() {
    tid_score.clear();
    robin_hood::unordered_flat_map<uint32_t, double> genusBaseSum;
    genusBaseSum.reserve(scratch.activeTidScores.size());
    if (adaptive_surface_mode && weightCtx.tid2genus) {
      for (uint32_t tid : scratch.activeTidScores) {
        const double baseScore = scratch.tidBaseScoreDense[tid];
        if (!(baseScore > 0.0)) {
          continue;
        }
        const uint32_t genus = genus_for_tid(weightCtx, tid);
        if (genus != 0u) {
          genusBaseSum[genus] += baseScore;
        }
      }
    }
    for (uint32_t tid : scratch.activeTidScores) {
      const double totalScore = scratch.tidScoreDense[tid];
      double primaryScore = totalScore;
      if (adaptive_surface_mode && weightCtx.tid2genus) {
        const double baseScore = scratch.tidBaseScoreDense[tid];
        const double completionScore = scratch.tidCompletionScoreDense[tid];
        if (completionScore > 0.0) {
          const uint32_t genus = genus_for_tid(weightCtx, tid);
          const auto genusIt =
              (genus != 0u) ? genusBaseSum.find(genus) : genusBaseSum.end();
          const double support =
              (genusIt != genusBaseSum.end()) ? genusIt->second : 0.0;
          const double completion =
              support > 0.0
                  ? completionScore * support / (support + completionScore)
                  : 0.0;
          primaryScore = baseScore + completion;
        }
      }
      if (primaryScore > 0.0) {
        tid_score[tid] = primaryScore;
      }
    }
    unique_hits.clear();
    for (uint32_t tid : scratch.activeUniqueHits) {
      unique_hits[tid] = scratch.uniqueHitsDense[tid];
    }
  }

  EvidenceStats collect_stats(double eff_eval) {
    materialize_score_maps();
    return collect_evidence_stats(tid_score, unique_hits, eff_eval);
  }
};

struct ScoreStateSnapshot {
  double eff_eval{0.0};
  size_t n_eval{0};
  std::vector<size_t> deferred_eval;
  EvidenceStats stats;
  bool high_conf_pre{false};
  double thr_conf{0.0};
  double gap_need{0.0};
  uint32_t best_tid{std::numeric_limits<uint32_t>::max()};
  uint32_t second_tid{std::numeric_limits<uint32_t>::max()};
  double best{0.0};
  double second{0.0};
  double best_ratio{0.0};
  double gap{0.0};
  double unique_count{0.0};
  double unique_ratio{0.0};
  std::vector<uint32_t> active_tid_scores;
  std::vector<double> active_tid_values;
  std::vector<double> active_tid_base_values;
  std::vector<double> active_tid_completion_values;
  std::vector<uint32_t> active_unique_hits;
  std::vector<double> active_unique_values;
};

static ScoreStateSnapshot snapshot_score_state(
    const ProcessScratch &scratch, double eff_eval, size_t n_eval,
    const std::vector<size_t> &deferredEval, const EvidenceStats &stats,
    bool highConfPre, double thr_conf, double gap_need, uint32_t bestTid,
    uint32_t secondTid, double best, double second, double best_ratio,
    double gap, double uniqueCount, double uniqueRatio) {
  ScoreStateSnapshot snap;
  snap.eff_eval = eff_eval;
  snap.n_eval = n_eval;
  snap.deferred_eval = deferredEval;
  snap.stats = stats;
  snap.high_conf_pre = highConfPre;
  snap.thr_conf = thr_conf;
  snap.gap_need = gap_need;
  snap.best_tid = bestTid;
  snap.second_tid = secondTid;
  snap.best = best;
  snap.second = second;
  snap.best_ratio = best_ratio;
  snap.gap = gap;
  snap.unique_count = uniqueCount;
  snap.unique_ratio = uniqueRatio;
  snap.active_tid_scores = scratch.activeTidScores;
  snap.active_tid_values.reserve(snap.active_tid_scores.size());
  snap.active_tid_base_values.reserve(snap.active_tid_scores.size());
  snap.active_tid_completion_values.reserve(snap.active_tid_scores.size());
  for (uint32_t tid : snap.active_tid_scores) {
    snap.active_tid_values.push_back(scratch.tidScoreDense[tid]);
    snap.active_tid_base_values.push_back(scratch.tidBaseScoreDense[tid]);
    snap.active_tid_completion_values.push_back(
        scratch.tidCompletionScoreDense[tid]);
  }
  snap.active_unique_hits = scratch.activeUniqueHits;
  snap.active_unique_values.reserve(snap.active_unique_hits.size());
  for (uint32_t tid : snap.active_unique_hits) {
    snap.active_unique_values.push_back(scratch.uniqueHitsDense[tid]);
  }
  return snap;
}

static void restore_score_state(
    const ScoreStateSnapshot &snap, ProcessScratch &scratch,
    uint32_t tidScoreEpoch, uint32_t uniqueHitsEpoch, double &eff_eval,
    size_t &n_eval, std::vector<size_t> &deferredEval, EvidenceStats &stats,
    bool &highConfPre, double &thr_conf, double &gap_need, uint32_t &bestTid,
    uint32_t &secondTid, double &best, double &second, double &best_ratio,
    double &gap, double &uniqueCount, double &uniqueRatio) {
  robin_hood::unordered_flat_set<uint32_t> oldTidScores;
  oldTidScores.reserve(snap.active_tid_scores.size());
  for (uint32_t tid : snap.active_tid_scores) {
    oldTidScores.insert(tid);
  }
  for (uint32_t tid : scratch.activeTidScores) {
    if (oldTidScores.find(tid) == oldTidScores.end() &&
        tid < scratch.tidScoreEpoch.size() &&
        scratch.tidScoreEpoch[tid] == tidScoreEpoch) {
      scratch.tidScoreEpoch[tid] = 0u;
      scratch.tidScoreDense[tid] = 0.0;
      scratch.tidBaseScoreDense[tid] = 0.0;
      scratch.tidCompletionScoreDense[tid] = 0.0;
    }
  }
  scratch.activeTidScores = snap.active_tid_scores;
  for (size_t i = 0; i < snap.active_tid_scores.size(); ++i) {
    const uint32_t tid = snap.active_tid_scores[i];
    scratch.tidScoreEpoch[tid] = tidScoreEpoch;
    scratch.tidScoreDense[tid] = snap.active_tid_values[i];
    scratch.tidBaseScoreDense[tid] = snap.active_tid_base_values[i];
    scratch.tidCompletionScoreDense[tid] =
        snap.active_tid_completion_values[i];
  }

  robin_hood::unordered_flat_set<uint32_t> oldUniqueHits;
  oldUniqueHits.reserve(snap.active_unique_hits.size());
  for (uint32_t tid : snap.active_unique_hits) {
    oldUniqueHits.insert(tid);
  }
  for (uint32_t tid : scratch.activeUniqueHits) {
    if (oldUniqueHits.find(tid) == oldUniqueHits.end() &&
        tid < scratch.uniqueHitsEpoch.size() &&
        scratch.uniqueHitsEpoch[tid] == uniqueHitsEpoch) {
      scratch.uniqueHitsEpoch[tid] = 0u;
      scratch.uniqueHitsDense[tid] = 0.0;
    }
  }
  scratch.activeUniqueHits = snap.active_unique_hits;
  for (size_t i = 0; i < snap.active_unique_hits.size(); ++i) {
    const uint32_t tid = snap.active_unique_hits[i];
    scratch.uniqueHitsEpoch[tid] = uniqueHitsEpoch;
    scratch.uniqueHitsDense[tid] = snap.active_unique_values[i];
  }
  eff_eval = snap.eff_eval;
  n_eval = snap.n_eval;
  deferredEval = snap.deferred_eval;
  stats = snap.stats;
  highConfPre = snap.high_conf_pre;
  thr_conf = snap.thr_conf;
  gap_need = snap.gap_need;
  bestTid = snap.best_tid;
  secondTid = snap.second_tid;
  best = snap.best;
  second = snap.second;
  best_ratio = snap.best_ratio;
  gap = snap.gap;
  uniqueCount = snap.unique_count;
  uniqueRatio = snap.unique_ratio;
}

struct ResultCandidateSelection {
  std::vector<SpoolCandidate> candidates;
  uint32_t max_count_tid{kSpoolUnclassifiedTid};
  double max_count_score{0.0};
  double max_count_raw_score{0.0};
  bool max_count_valid{false};
};

struct CandidateThresholds {
  double eff_cap{1.0};
  bool use_em{false};
  size_t thr_eval{1};
  size_t thr_min_eval{0};
  size_t thr_final_raw{0};
  double thr_base{0.0};
};

static CandidateThresholds compute_candidate_thresholds(
    const ClassifyConfig &config, bool highConfPre, double dispersion_hi,
    double eff_eval, size_t n_eval, double best) {
  CandidateThresholds out;
  out.eff_cap = std::max(1.0, eff_eval);
  out.use_em = !highConfPre;

  const double maxEvidence = std::min(best, eff_eval);
  double beta = config.firstFilterBeta;
  if (!(beta > 0.0)) {
    beta = 0.8;
  }
  if (!config.firstFilterBeta_user) {
    const double beta_em = lerp(0.50, 0.45, dispersion_hi);
    const double beta_no_em = lerp(0.50, 0.80, dispersion_hi);
    beta = out.use_em ? beta_em : beta_no_em;
  }
  beta = std::clamp(beta, 0.0, 1.0);
  const size_t thr_beta =
      static_cast<size_t>(std::floor(beta * std::max(0.0, maxEvidence)));

  out.thr_eval =
      static_cast<size_t>(std::ceil(config.shotThreshold * eff_eval));
  if (out.thr_eval == 0) {
    out.thr_eval = 1;
  }
  if (out.use_em) {
    const double softened_ratio = std::min(config.shotThreshold, 0.45);
    size_t em_eval =
        static_cast<size_t>(std::ceil(eff_eval * softened_ratio));
    if (em_eval == 0 && eff_eval > 0.0) {
      em_eval = 1;
    }
    out.thr_eval = std::min(out.thr_eval, em_eval);
  }

  if (n_eval > 0) {
    const double factor = out.use_em ? 0.15 : 0.30;
    out.thr_min_eval =
        static_cast<size_t>(std::ceil(factor * static_cast<double>(n_eval)));
  }
  out.thr_min_eval =
      std::max<size_t>(out.thr_min_eval, out.use_em ? 4 : 8);

  if (out.use_em) {
    const size_t thr_beta_eval = std::min(thr_beta, out.thr_eval);
    out.thr_final_raw =
        std::max({thr_beta_eval, static_cast<size_t>(1), out.thr_min_eval});
  } else {
    out.thr_final_raw =
        std::max({out.thr_eval, thr_beta, out.thr_min_eval});
  }

  size_t thr_min_eval_low = out.thr_min_eval;
  if (out.thr_min_eval > 0) {
    thr_min_eval_low = std::max<size_t>(
        1, static_cast<size_t>(
               std::floor(0.3 * static_cast<double>(out.thr_min_eval))));
  }
  const double thr_low =
      static_cast<double>(std::max(out.thr_eval, thr_min_eval_low));
  const double thr_high = static_cast<double>(out.thr_final_raw);
  out.thr_base = lerp(thr_low, thr_high, dispersion_hi);
  return out;
}

static size_t count_same_genus_competitors(
    const std::vector<std::pair<uint32_t, double>> &ranked,
    const TaxDict &tax, const WeightingContext &weightCtx, uint32_t bestTid,
    uint32_t best_genus_id, size_t scan_limit) {
  const uint32_t best_rep_tid_id = species_rep_for_tid(weightCtx, bestTid);
  size_t same_genus_competitors = 0;
  robin_hood::unordered_flat_set<uint32_t> scanned_rep_ids;
  scanned_rep_ids.reserve(scan_limit);
  for (size_t i = 0; i < scan_limit; ++i) {
    const uint32_t tid_id = ranked[i].first;
    if (tid_id >= tax.id2str.size()) {
      continue;
    }
    const uint32_t rep_id = species_rep_for_tid(weightCtx, tid_id);
    if (!scanned_rep_ids.insert(rep_id).second) {
      continue;
    }
    if (rep_id == best_rep_tid_id) {
      continue;
    }
    if (best_genus_id != 0u &&
        genus_for_tid(weightCtx, tid_id) == best_genus_id) {
      ++same_genus_competitors;
    }
  }
  return same_genus_competitors;
}

static double effective_score_count(
    const std::vector<std::pair<uint32_t, double>> &ranked) {
  double score_total = 0.0;
  double score_square_total = 0.0;
  for (const auto &kv : ranked) {
    const double score = std::max(0.0, kv.second);
    score_total += score;
    score_square_total += score * score;
  }
  if (score_total > 0.0 && score_square_total > 0.0) {
    return (score_total * score_total) / score_square_total;
  }
  return 1.0;
}

struct LocalRescueCandidate {
  uint32_t tid{std::numeric_limits<uint32_t>::max()};
  uint32_t rep{std::numeric_limits<uint32_t>::max()};
  uint32_t priority{2u};
  double score{0.0};
  double raw_score{0.0};
};

static void append_local_rescue_candidates(
    const TaxDict &tax, const WeightingContext &weightCtx,
    const std::vector<std::pair<uint32_t, double>> &ranked, double effCap,
    double thr_base, uint32_t bestTid, uint32_t best_genus_id, double best,
    double uniqueRatio, std::vector<SpoolCandidate> &candidates) {
  const uint32_t best_rep_tid_id = species_rep_for_tid(weightCtx, bestTid);
  robin_hood::unordered_flat_set<uint32_t> kept_rep_ids;
  kept_rep_ids.reserve(candidates.size() + 8);
  for (const auto &cand : candidates) {
    if (cand.tid < tax.id2str.size()) {
      kept_rep_ids.insert(species_rep_for_tid(weightCtx, cand.tid));
    }
  }

  const size_t scan_limit = std::min<size_t>(64, ranked.size());
  const size_t same_genus_competitors = count_same_genus_competitors(
      ranked, tax, weightCtx, bestTid, best_genus_id, scan_limit);
  const double best_score = std::clamp(best, 0.0, effCap);
  const double second_score =
      (ranked.size() > 1) ? std::clamp(ranked[1].second, 0.0, effCap) : 0.0;
  const double gap_term =
      1.0 - std::clamp((best_score - second_score) /
                           std::max(best_score, 1e-6),
                       0.0, 1.0);
  const double genus_term =
      std::clamp(static_cast<double>(same_genus_competitors) / 4.0, 0.0,
                 1.0);
  const double unique_term =
      1.0 - std::clamp(uniqueRatio / 0.25, 0.0, 1.0);
  const double uncertainty =
      std::clamp(0.45 * gap_term + 0.35 * genus_term + 0.20 * unique_term,
                 0.0, 1.0);
  const size_t keep_budget =
      1u + static_cast<size_t>(std::ceil(7.0 * uncertainty));
  const double keep_floor = std::max(0.25 * best_score, 0.50 * thr_base);

  std::vector<LocalRescueCandidate> rescue_candidates;
  rescue_candidates.reserve(scan_limit);
  for (size_t i = 0; i < scan_limit; ++i) {
    const uint32_t tid_id = ranked[i].first;
    if (tid_id >= tax.id2str.size()) {
      continue;
    }
    const double raw_score = ranked[i].second;
    const double score = std::clamp(raw_score, 0.0, effCap);
    if (score < keep_floor) {
      continue;
    }
    const uint32_t rep_id = species_rep_for_tid(weightCtx, tid_id);
    if (kept_rep_ids.find(rep_id) != kept_rep_ids.end()) {
      continue;
    }
    uint32_t priority = 2u;
    if (best_genus_id != 0u &&
        genus_for_tid(weightCtx, tid_id) == best_genus_id) {
      priority = 0u;
    } else if (best_score > 0.0 && score >= 0.50 * best_score) {
      priority = 1u;
    }
    rescue_candidates.push_back(
        LocalRescueCandidate{tid_id, rep_id, priority, score, raw_score});
  }
  std::sort(rescue_candidates.begin(), rescue_candidates.end(),
            [](const LocalRescueCandidate &lhs,
               const LocalRescueCandidate &rhs) {
              if (lhs.priority != rhs.priority) {
                return lhs.priority < rhs.priority;
              }
              if (lhs.score != rhs.score) {
                return lhs.score > rhs.score;
              }
              return lhs.tid < rhs.tid;
            });
  for (const auto &cand : rescue_candidates) {
    if (kept_rep_ids.size() >= keep_budget) {
      break;
    }
    if (!kept_rep_ids.insert(cand.rep).second) {
      continue;
    }
    candidates.push_back(SpoolCandidate{cand.tid, cand.score, cand.raw_score});
  }
}

static ResultCandidateSelection select_result_candidates_from_scores(
    const std::string &id, const TaxDict &tax,
    const WeightingContext &weightCtx, const ClassifyConfig &config,
    ProcessScratch &scratch,
    const robin_hood::unordered_flat_map<uint32_t, double> &tidScore,
    bool highConfPre, double dispersion_hi, double eff_eval, size_t n_eval,
    uint32_t bestTid, double best, double uniqueRatio,
    const std::string &bestTaxidStr, bool dump_preem_enabled,
    const std::string &dump_preem_path) {
  const CandidateThresholds thresholds = compute_candidate_thresholds(
      config, highConfPre, dispersion_hi, eff_eval, n_eval, best);
  const double effCap = thresholds.eff_cap;
  const bool use_em = thresholds.use_em;
  const size_t thr_final_raw = thresholds.thr_final_raw;
  ResultCandidateSelection selection;
  auto add_candidate = [&](uint32_t tid_id, double score, double raw_score) {
    if (tid_id < tax.id2str.size() && score > 0.0) {
      selection.candidates.push_back(SpoolCandidate{tid_id, score, raw_score});
    }
  };

  if (!selection.max_count_valid && bestTid < tax.id2str.size()) {
    const double bestEvidence = std::clamp(best, 0.0, effCap);
    if (bestEvidence > 0.0) {
      selection.max_count_tid = bestTid;
      selection.max_count_score = bestEvidence;
      selection.max_count_raw_score = best;
      selection.max_count_valid = true;
    }
  }

  if (highConfPre && bestTid < tax.id2str.size()) {
    const double bestEvidence = std::clamp(best, 0.0, effCap);
    add_candidate(bestTid, bestEvidence, best);
    selection.max_count_tid = bestTid;
    selection.max_count_score = bestEvidence;
    selection.max_count_raw_score = best;
    selection.max_count_valid = true;
  } else if (use_em && !tidScore.empty()) {
    auto &ranked = scratch.rankedTidScores;
    ranked.clear();
    ranked.reserve(tidScore.size());
    for (const auto &kv : tidScore) {
      ranked.emplace_back(kv.first, kv.second);
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    const uint32_t best_genus_id = genus_for_tid(weightCtx, bestTid);
    size_t same_genus_competitors_for_calibration = 0;
    double effective_candidate_count_for_calibration = 1.0;
    if (config.sample_state_calibration) {
      const size_t scan_limit = std::min<size_t>(64, ranked.size());
      same_genus_competitors_for_calibration =
          count_same_genus_competitors(ranked, tax, weightCtx, bestTid,
                                       best_genus_id, scan_limit);
      effective_candidate_count_for_calibration =
          effective_score_count(ranked);
    }

    const size_t floorK_max = std::min<size_t>(128, ranked.size());
    const size_t dispersionFloor = std::min<size_t>(
        floorK_max,
        static_cast<size_t>(std::llround((1.0 - dispersion_hi) *
                                         static_cast<double>(floorK_max))));
    const size_t calibrationSupportFloor =
        config.sample_state_calibration
            ? sample_state_calibration_support_floor(
                  effective_candidate_count_for_calibration,
                  same_genus_competitors_for_calibration, floorK_max)
            : 0;
    const size_t floorCount =
        std::max(dispersionFloor, calibrationSupportFloor);

    if (dump_preem_enabled) {
      dump_preem_candidate_retention(dump_preem_path, id, tax, weightCtx,
                                     ranked, effCap, thresholds.thr_base,
                                     floorCount, bestTaxidStr, best_genus_id);
    }

    selection.candidates.reserve(ranked.size());
    for (size_t i = 0; i < ranked.size(); ++i) {
      const uint32_t tid_id = ranked[i].first;
      if (tid_id >= tax.id2str.size()) {
        continue;
      }
      const double rawScore = ranked[i].second;
      const double countVal = std::clamp(rawScore, 0.0, effCap);
      if (countVal <= 0.0) {
        continue;
      }
      const bool keep = (i < floorCount || countVal >= thresholds.thr_base);
      if (keep) {
        add_candidate(tid_id, countVal, rawScore);
        if (tid_id == bestTid) {
          selection.max_count_tid = tid_id;
          selection.max_count_score = countVal;
          selection.max_count_raw_score = rawScore;
          selection.max_count_valid = true;
        }
      }
    }

    if (!highConfPre && !ranked.empty()) {
      append_local_rescue_candidates(
          tax, weightCtx, ranked, effCap, thresholds.thr_base, bestTid,
          best_genus_id, best, uniqueRatio, selection.candidates);
    }
  } else {
    for (const auto &[tid_id, rawScore] : tidScore) {
      const double countVal = std::clamp(rawScore, 0.0, effCap);
      if (countVal >= static_cast<double>(thr_final_raw)) {
        add_candidate(tid_id, countVal, rawScore);
        if (tid_id == bestTid) {
          selection.max_count_tid = tid_id;
          selection.max_count_score = countVal;
          selection.max_count_raw_score = rawScore;
          selection.max_count_valid = true;
        }
      }
    }
  }

  if (selection.candidates.empty() && use_em && selection.max_count_valid &&
      selection.max_count_score > 0.0) {
    selection.candidates.push_back(SpoolCandidate{selection.max_count_tid,
                                                  selection.max_count_score,
                                                  selection.max_count_raw_score});
  }
  return selection;
}

static void finalize_read_record(
    const std::string &id, const TaxDict &tax, FileInfo &fileInfo,
    PresenceAccumulator *presenceAcc, double uniqueCount, double uniqueRatio,
    double eff_eval, size_t readLen, uint32_t bestTaxidHintTid,
    const std::string &bestTaxidStr, bool maxCountValid,
    uint32_t maxCountTid, double maxCountScore, double maxCountRawScore,
    bool tidScoreEmpty, std::string rejectReason,
    uint32_t profileResponseTaxid,
    std::vector<SpoolCandidate> resultCandidates,
    std::vector<SpoolCandidate> abundanceCandidates,
    std::vector<SpoolCandidate> sampleMixtureCandidates,
    std::vector<classifyResult> *classifyResults,
    std::vector<CompactClassifyResult> *compactResults) {
  if (!resultCandidates.empty()) {
    if (presenceAcc != nullptr) {
      const auto &top = resultCandidates.front();
      if (top.tid != kSpoolUnclassifiedTid) {
        const bool unique_ok = (uniqueCount >= 3.0) ||
                               (uniqueCount >= 2.0 && uniqueRatio >= 0.12);
        presenceAcc->add_read_support(top.tid, unique_ok);
      }
    }
    const bool isUniqueMapping = (resultCandidates.size() == 1);
    if (isUniqueMapping && resultCandidates.front().tid < tax.id2str.size()) {
      fileInfo.uniqueTaxids.insert(tax.id2str[resultCandidates.front().tid]);
    }

    fileInfo.classifiedNum++;
    if (resultCandidates.size() == 1) {
      if (!maxCountValid) {
        maxCountTid = resultCandidates.front().tid;
        maxCountScore = resultCandidates.front().score;
        maxCountRawScore = resultCandidates.front().raw_score > 0.0
                               ? resultCandidates.front().raw_score
                               : resultCandidates.front().score;
        maxCountValid = true;
      }
      resultCandidates.clear();
      resultCandidates.push_back(
          SpoolCandidate{maxCountTid, maxCountScore, maxCountRawScore});
    }
  } else {
    fileInfo.unclassifiedNum++;
    resultCandidates.push_back(SpoolCandidate{kSpoolUnclassifiedTid, 1.0, 1.0});
    if (rejectReason.empty()) {
      rejectReason = (!maxCountValid && tidScoreEmpty) ? "no_candidate"
                                                       : "low_eval";
    }
  }

  if (compactResults != nullptr) {
    CompactClassifyResult result;
    result.evaluated = eff_eval;
    result.query_length = static_cast<uint32_t>(
        std::min<size_t>(readLen, std::numeric_limits<uint32_t>::max()));
    result.id = id;
    result.best_taxid_hint = bestTaxidHintTid;
    result.profile_response_taxid = profileResponseTaxid;
    result.reject_reason = std::move(rejectReason);
    result.candidates = std::move(resultCandidates);
    result.abundance_candidates = std::move(abundanceCandidates);
    result.sample_mixture_candidates = std::move(sampleMixtureCandidates);
    compactResults->emplace_back(std::move(result));
  } else if (classifyResults != nullptr) {
    classifyResult result;
    result.evaluated = eff_eval;
    result.query_length = static_cast<uint32_t>(
        std::min<size_t>(readLen, std::numeric_limits<uint32_t>::max()));
    result.id = id;
    result.best_taxid_hint = bestTaxidStr;
    result.profile_response_taxid = profileResponseTaxid;
    result.reject_reason = std::move(rejectReason);
    result.taxidCount.reserve(resultCandidates.size());
    for (const auto &candidate : resultCandidates) {
      if (candidate.tid == kSpoolUnclassifiedTid) {
        result.taxidCount.emplace_back("unclassified", candidate.score);
      } else if (candidate.tid < tax.id2str.size()) {
        result.taxidCount.emplace_back(tax.id2str[candidate.tid],
                                       candidate.score);
      }
    }
    result.abundanceCount.reserve(abundanceCandidates.size());
    for (const auto &candidate : abundanceCandidates) {
      if (candidate.tid < tax.id2str.size() && candidate.score > 0.0) {
        result.abundanceCount.emplace_back(tax.id2str[candidate.tid],
                                           candidate.score);
      }
    }
    classifyResults->emplace_back(std::move(result));
  }
}

void processSequence(
    const std::vector<uint64_t> &hashs1, size_t readLen,
    ChimeraBuild::IMCFConfig &imcfConfig, const TaxDict &tax,
    ClassifyConfig &config, const WeightingContext &weightCtx, GroupHeat &heat,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const std::string &id,
    std::vector<classifyResult> *classifyResults,
    std::vector<CompactClassifyResult> *compactResults, FileInfo &fileInfo,
    PresenceAccumulator *presenceAcc, ProcessScratch &scratch) {
  const size_t binNumAll = tax.idx2id.size();
  heat.ensure(binNumAll);
  const double community_dispersion_s =
      clamp01(config.community_dispersion_s);
  const double dispersion_hi = community_dispersion_s;
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

  const char *dump_preem_env = std::getenv("CHIMERA_DUMP_PREEM_TSV");
  const bool dump_preem_enabled = (dump_preem_env && dump_preem_env[0] != '\0');
  const std::string dump_preem_path =
      dump_preem_enabled ? std::string(dump_preem_env) : std::string();
  prepare_hash_sample_and_route(hashs1, config, weightCtx, scratch);
  uint32_t profileResponseTaxid = 0;
  auto &sampleVals = scratch.sampleVals;
  auto &routeVals = scratch.routeVals;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> sampleBinScore;
  robin_hood::unordered_flat_map<uint32_t, uint64_t> repCoarseScore;
  robin_hood::unordered_flat_map<uint32_t, uint64_t> genusCoarseScore;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> repRareSupport;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> genusRareSupport;
  auto &rankedBins = scratch.rankedBins;
  rankedBins.clear();
  uint64_t coarseTotal = 0;
  auto &deferredEval = scratch.deferredEval;
  deferredEval.clear();
  if (hashs1.size() > 64) {
    deferredEval.reserve(hashs1.size() - 64);
  }

  coarseTotal =
      collect_coarse_evidence(sampleVals, imcf, tax, weightCtx, scratch,
                              sampleBinScore, repCoarseScore, genusCoarseScore);
  auto &topBins = scratch.topBins;
  robin_hood::unordered_flat_set<uint32_t> lowDegPreserve;
  bool full_surface_mode = build_initial_candidate_surface(
      routeVals, sampleBinScore, coarseTotal, dispersion_hi, binNumAll, imcf,
      tax, weightCtx, heat, scratch, repRareSupport, genusRareSupport,
      lowDegPreserve);

  const bool adaptive_surface_mode = expand_adaptive_candidate_surface(
      readLen, dispersion_hi, coarseTotal, repCoarseScore, genusCoarseScore,
      repRareSupport, genusRareSupport, lowDegPreserve, rankedBins, tax,
      weightCtx, binNumAll, full_surface_mode, topBins, scratch);
  if (!full_surface_mode) {
    for (auto bin : topBins) {
      uint32_t delta = 1;
      if (auto it = sampleBinScore.find(bin); it != sampleBinScore.end()) {
        delta = std::max<uint32_t>(delta, it->second);
      }
      heat.boost(bin, delta);
    }
  }
  heat.decay_if_needed();

  auto &minimizerTids = scratch.minimizerTids;
  minimizerTids.clear();
  minimizerTids.reserve(64);
  auto &minimizerBins = scratch.minimizerBins;
  minimizerBins.clear();
  minimizerBins.reserve(64);
  const size_t tidCountAll = tax.id2str.size();
  ensure_epoch_vector(scratch.minimizerTidEpoch, tidCountAll);
  ensure_epoch_vector(scratch.minimizerTidBaseEpoch, tidCountAll);
  ensure_epoch_vector(scratch.minimizerTidCompletionEpoch, tidCountAll);
  ensure_epoch_vector(scratch.minimizerBinEpoch, binNumAll);
  ReadScoringState scoring(scratch, weightCtx, tidCountAll,
                           adaptive_surface_mode);

  double eff_eval = 0.0;
  size_t n_eval = 0;

  const std::vector<uint32_t> *activeSubset =
      (full_surface_mode || topBins.size() == binNumAll) ? nullptr : &topBins;
  static thread_local std::vector<uint32_t> baseBinMarks;
  static thread_local uint32_t baseBinMarkEpoch = 0u;
  bool baseBinMarksActive = false;
  auto rebuild_base_bin_marks = [&]() {
    baseBinMarksActive = false;
    if (!adaptive_surface_mode || full_surface_mode ||
        scratch.baseTopBins.empty()) {
      return;
    }
    if (baseBinMarks.size() != binNumAll) {
      baseBinMarks.assign(binNumAll, 0u);
      baseBinMarkEpoch = 1u;
    } else {
      ++baseBinMarkEpoch;
      if (baseBinMarkEpoch == 0u) {
        std::fill(baseBinMarks.begin(), baseBinMarks.end(), 0u);
        baseBinMarkEpoch = 1u;
      }
    }
    for (uint32_t bin : scratch.baseTopBins) {
      if (bin < baseBinMarks.size()) {
        baseBinMarks[bin] = baseBinMarkEpoch;
      }
    }
    baseBinMarksActive = true;
  };
  rebuild_base_bin_marks();
  static thread_local std::vector<uint32_t> topBinMarks;
  static thread_local uint32_t topBinMarkEpoch = 0u;
  bool topBinMarksActive = false;
  auto rebuild_top_bin_marks = [&]() {
    topBinMarksActive = false;
    if (full_surface_mode || topBins.size() == binNumAll) {
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

  auto is_base_evidence_bin = [&](uint32_t bin) {
    if (!adaptive_surface_mode || full_surface_mode) {
      return true;
    }
    if (!baseBinMarksActive || bin >= baseBinMarks.size()) {
      return false;
    }
    return baseBinMarks[bin] == baseBinMarkEpoch;
  };

  double idf_power = 1.0 + (0.25 * dispersion_hi);
  bool emitPresenceEvidence = true;

  auto evaluate_minimizer = [&](uint64_t value,
                                const std::vector<uint32_t> *subset) -> double {
    const uint32_t minimizerEpoch = next_epoch(
        scratch.minimizerTidEpoch, scratch.minimizerTidEpochValue);
    const uint32_t minimizerBinEpoch = next_epoch(
        scratch.minimizerBinEpoch, scratch.minimizerBinEpochValue);
    minimizerTids.clear();
    minimizerBins.clear();
    uint32_t last_hit_bin = std::numeric_limits<uint32_t>::max();
    uint32_t min_hit_bin = std::numeric_limits<uint32_t>::max();
    auto emit = [&](uint32_t bin, uint16_t sp) {
      uint32_t tid = tax.rep_tid_for_bin_slot(bin, sp);
      if (tid == kInvalidTidId) {
        return;
      }
      if (scratch.minimizerTidEpoch[tid] != minimizerEpoch) {
        scratch.minimizerTidEpoch[tid] = minimizerEpoch;
        minimizerTids.push_back(tid);
      }
      if (is_base_evidence_bin(bin)) {
        scratch.minimizerTidBaseEpoch[tid] = minimizerEpoch;
      } else {
        scratch.minimizerTidCompletionEpoch[tid] = minimizerEpoch;
      }
      if (bin != last_hit_bin) {
        last_hit_bin = bin;
        if (scratch.minimizerBinEpoch[bin] != minimizerBinEpoch) {
          scratch.minimizerBinEpoch[bin] = minimizerBinEpoch;
          minimizerBins.push_back(bin);
          min_hit_bin = std::min(min_hit_bin, bin);
        }
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

    const size_t deg_subset = minimizerTids.size();
    if (deg_subset == 0) {
      return 0.0;
    }
    size_t deg_effective = deg_subset;

    size_t df_bins = 0;
    if (!minimizerBins.empty()) {
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
    const double df_idf = df_for_idf(df_bins, df_eff, readLen, dispersion_hi,
                                     kDfEffIdfMaxLen);
    double idf_raw = idf_raw_from_df_bins(totalBins, df_idf);
    double idf = clamp_idf(idf_raw, dispersion_hi, config.idf_max, idf_power);

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
        uniqueEdge ? 1.0 : (fragmentedUnique ? dispersion_hi : 0.0);
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
      const bool hasBase =
          scratch.minimizerTidBaseEpoch[tid] == minimizerEpoch;
      const bool hasCompletion =
          scratch.minimizerTidCompletionEpoch[tid] == minimizerEpoch;
      scoring.add_tid_score(tid, contrib, hasBase || !hasCompletion);
    }
    if (emitPresenceEvidence && presenceEnabled && !minimizerTids.empty()) {
      const double hit_weight = presence_weight;
      const double score_weight = hit_weight;
      uint32_t unique_bucket = std::numeric_limits<uint32_t>::max();
      uint32_t breadth_bucket = std::numeric_limits<uint32_t>::max();
      if (localUniqueEdge && presenceSketchBits > 0) {
        unique_bucket = sketch_index(value ^ 0x9e3779b97f4a7c15ULL);
        if (min_hit_bin != std::numeric_limits<uint32_t>::max()) {
          breadth_bucket = sketch_index(
              static_cast<uint64_t>(min_hit_bin) ^
              0xd1b54a32d192ed03ULL);
        }
      }
      presenceAcc->add_targets(minimizerTids, hit_weight, score_weight,
                               unique_strength, localUniqueEdge, unique_bucket,
                               breadth_bucket);
    }
    if (unique_strength > 0.0 && !minimizerTids.empty()) {
      scoring.add_unique_hit(minimizerTids.front(), unique_strength);
    }

    return contrib;
  };

  auto recompute_subset_state = [&]() {
    activeSubset = (full_surface_mode || topBins.size() == binNumAll)
                       ? nullptr
                       : &topBins;
    rebuild_top_bin_marks();
  };

  auto top_bin_marked = [&](uint32_t bin) {
    if (topBinMarksActive && bin < topBinMarks.size()) {
      return topBinMarks[bin] == topBinMarkEpoch;
    }
    return std::binary_search(topBins.begin(), topBins.end(), bin);
  };
  auto mark_top_bin = [&](uint32_t bin) {
    if (topBinMarksActive && bin < topBinMarks.size()) {
      topBinMarks[bin] = topBinMarkEpoch;
    }
  };

  auto collect_stats = [&]() -> EvidenceStats {
    return scoring.collect_stats(eff_eval);
  };

  size_t n0 = std::min<size_t>(64, hashs1.size());
  for (size_t i = 0; i < n0; ++i) {
    const double c = evaluate_minimizer(hashs1[i], activeSubset);
    eff_eval += c;
  }
  n_eval = n0;

  EvidenceStats stats = collect_stats();
  bool highConfPre = meets_quick_high_conf(config, stats, eff_eval);

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
    highConfPre = meets_quick_high_conf(config, stats, eff_eval);

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
      highConfPre = meets_quick_high_conf(config, stats, eff_eval);
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
  auto apply_stats = [&](const EvidenceStats &s) {
    bestTid = s.bestTid;
    secondTid = s.secondTid;
    best = s.best;
    second = s.second;
    best_ratio = s.ratio;
    gap = s.gap;
    uniqueCount = s.uniqueCount;
    uniqueRatio = s.uniqueRatio;
  };
  auto complete_original_subset_eval = [&]() {
    if (n_eval >= hashs1.size()) {
      return;
    }
    if (!deferredEval.empty()) {
      for (size_t idx : deferredEval) {
        if (idx >= hashs1.size()) {
          continue;
        }
        const double c = evaluate_minimizer(hashs1[idx], activeSubset);
        eff_eval += c;
        ++n_eval;
      }
      deferredEval.clear();
    } else {
      for (size_t i = n_eval; i < hashs1.size(); ++i) {
        const double c = evaluate_minimizer(hashs1[i], activeSubset);
        eff_eval += c;
        ++n_eval;
      }
    }
    stats = collect_stats();
    highConfPre = meets_quick_high_conf(config, stats, eff_eval);
    thr_conf = std::ceil(config.shotThreshold * eff_eval);
    gap_need = std::max(0.5, eff_eval / 24.0);
    apply_stats(stats);
  };
  std::string bestTaxidStr;
  if (bestTid != std::numeric_limits<uint32_t>::max() &&
      bestTid < tax.id2str.size()) {
    bestTaxidStr = tax.id2str[bestTid];
  }
  uint32_t bestTaxidHintTid =
      bestTaxidStr.empty() ? kSpoolUnclassifiedTid : bestTid;

  std::vector<SpoolCandidate> abundanceCandidates;

  if (highConfPre && n_eval < hashs1.size()) {
    const ScoreStateSnapshot primarySnapshot = snapshot_score_state(
        scratch, eff_eval, n_eval, deferredEval, stats, highConfPre, thr_conf,
        gap_need, bestTid, secondTid, best, second, best_ratio, gap,
        uniqueCount, uniqueRatio);
    const bool prevEmitPresenceEvidence = emitPresenceEvidence;
    emitPresenceEvidence = false;
    complete_original_subset_eval();
    abundanceCandidates = capture_abundance_candidates_from_scores(
        scratch, tax.id2str.size(), eff_eval, n_eval);
    emitPresenceEvidence = prevEmitPresenceEvidence;
    restore_score_state(primarySnapshot, scratch, scoring.tid_score_epoch,
                        scoring.unique_hits_epoch, eff_eval, n_eval,
                        deferredEval, stats, highConfPre, thr_conf, gap_need,
                        bestTid, secondTid, best, second, best_ratio, gap,
                        uniqueCount, uniqueRatio);
    scoring.materialize_score_maps();
  }

  bool expanded = false;
  if (!adaptive_surface_mode && !full_surface_mode &&
      bestTid < tax.tid2bin.size()) {
    const auto &shards = tax.tid2bin[bestTid];
    for (uint32_t bin : shards) {
      if (!top_bin_marked(bin)) {
        topBins.push_back(bin);
        mark_top_bin(bin);
        expanded = true;
      }
    }
  }

  if (expanded) {
    std::sort(topBins.begin(), topBins.end());
    topBins.erase(std::unique(topBins.begin(), topBins.end()), topBins.end());
    full_surface_mode = (topBins.size() == binNumAll);
    recompute_subset_state();

    scoring.reset_dense_scores();
    eff_eval = 0.0;
    n_eval = 0;
    for (size_t i = 0; i < hashs1.size(); ++i) {
      const double c = evaluate_minimizer(hashs1[i], activeSubset);
      eff_eval += c;
      ++n_eval;
    }
    stats = collect_stats();
    highConfPre = meets_quick_high_conf(config, stats, eff_eval);
    thr_conf = std::ceil(config.shotThreshold * eff_eval);
    gap_need = std::max(0.5, eff_eval / 24.0);
    apply_stats(stats);
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

  std::string rejectReason;
  auto selection = select_result_candidates_from_scores(
      id, tax, weightCtx, config, scratch, scoring.tid_score, highConfPre,
      dispersion_hi, eff_eval, n_eval, bestTid, best, uniqueRatio,
      bestTaxidStr, dump_preem_enabled, dump_preem_path);
  std::vector<SpoolCandidate> resultCandidates =
      std::move(selection.candidates);
  profileResponseTaxid =
      profile_response_taxid_from_candidates(resultCandidates, tax, weightCtx);
  uint32_t maxCountTid = selection.max_count_tid;
  double maxCountScore = selection.max_count_score;
  double maxCountRawScore = selection.max_count_raw_score;
  bool maxCountValid = selection.max_count_valid;

  if (abundanceCandidates.empty()) {
    abundanceCandidates = capture_abundance_candidates_from_scores(
        scratch, tax.id2str.size(), eff_eval, n_eval);
  }
  std::vector<SpoolCandidate> sampleMixtureCandidates =
      capture_sample_mixture_candidates_from_base_scores(scratch,
                                                         tax.id2str.size());

  finalize_read_record(
      id, tax, fileInfo, presenceAcc, uniqueCount, uniqueRatio, eff_eval,
      readLen,
      bestTaxidHintTid, bestTaxidStr, maxCountValid, maxCountTid,
      maxCountScore, maxCountRawScore, scoring.tid_score.empty(),
      std::move(rejectReason), profileResponseTaxid,
      std::move(resultCandidates), std::move(abundanceCandidates),
      std::move(sampleMixtureCandidates), classifyResults, compactResults);
}

void processBatch(
    batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
    const TaxDict &tax, ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<classifyResult> &classifyResults,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    FileInfo &fileInfo, GroupHeat &heat, const WeightingContext &weightCtx,
    PresenceAccumulator *presenceAcc, ProcessScratch &scratch) {
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
                      heat, imcf, batch.ids[i],
                      &classifyResults,
                      nullptr, fileInfo, presenceAcc, scratch);
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
                      heat, imcf, batch.ids[i],
                      &classifyResults,
                      nullptr, fileInfo, presenceAcc, scratch);
	    }
  }
}

void processBatchCompact(
    batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
    const TaxDict &tax, ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    std::vector<CompactClassifyResult> &classifyResults,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    FileInfo &fileInfo, GroupHeat &heat, const WeightingContext &weightCtx,
    PresenceAccumulator *presenceAcc, ProcessScratch &scratch) {
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
                      heat, imcf, batch.ids[i], nullptr,
                      &classifyResults, fileInfo, presenceAcc, scratch);
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
                      heat, imcf, batch.ids[i], nullptr,
                      &classifyResults, fileInfo, presenceAcc, scratch);
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
    const WeightingContext &weightCtx, PresenceSummary *presenceSummary,
    std::vector<QueueThrottle> *queueThrottles) {

#pragma omp parallel
  {
#ifdef _OPENMP
    const int thread_id = omp_get_thread_num();
#else
    const int thread_id = 0;
#endif
    const size_t queue_index = static_cast<size_t>(
        std::clamp<int>(thread_id, 0, static_cast<int>(readQueues.size() - 1)));
    moodycamel::ConcurrentQueue<batchReads> &readQueue =
        readQueues[queue_index];
    QueueThrottle *queueThrottle =
        (queueThrottles != nullptr && queue_index < queueThrottles->size())
            ? &(*queueThrottles)[queue_index]
            : nullptr;

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
        release_queue_slot(queueThrottle, estimate_batch_bytes(batch));
        processBatch(batch, imcfConfig, tax, config, imcf,
                     localClassifyResults, feature_params, feature_min_len,
                     localFileInfo, heat, weightCtx, presencePtr, scratch);
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
                             std::make_move_iterator(
                                 localClassifyResults.begin()),
                             std::make_move_iterator(
                                 localClassifyResults.end()));
      std::vector<classifyResult>().swap(localClassifyResults);

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

namespace {

void write_spool_results(const std::vector<CompactClassifyResult> &results,
                         std::ostream &os) {
  for (const auto &result : results) {
    write_spool_record(os, result);
  }
}

} // namespace

void classify_streaming_spool(
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
    ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const TaxDict &tax,
    const std::vector<std::string> &spoolPaths, FileInfo &fileInfo,
    std::atomic<bool> &producer_done,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const WeightingContext &weightCtx, PresenceSummary *presenceSummary,
    std::vector<QueueThrottle> *queueThrottles,
    ClassifyProgressCounters *progress) {

#pragma omp parallel
  {
#ifdef _OPENMP
    const int thread_id = omp_get_thread_num();
#else
    const int thread_id = 0;
#endif
    const size_t queue_index = static_cast<size_t>(
        std::clamp<int>(thread_id, 0, static_cast<int>(readQueues.size() - 1)));
    moodycamel::ConcurrentQueue<batchReads> &readQueue =
        readQueues[queue_index];
    QueueThrottle *queueThrottle =
        (queueThrottles != nullptr && queue_index < queueThrottles->size())
            ? &(*queueThrottles)[queue_index]
            : nullptr;

    const size_t spool_index = static_cast<size_t>(std::clamp<int>(
        thread_id, 0, static_cast<int>(spoolPaths.size() - 1)));
    std::vector<char> spoolBuffer(1 << 20);
    std::ofstream spool;
    spool.rdbuf()->pubsetbuf(spoolBuffer.data(),
                             static_cast<std::streamsize>(spoolBuffer.size()));
    spool.open(spoolPaths[spool_index], std::ios::binary);
    if (!spool.is_open()) {
      throw std::runtime_error("Failed to open classify spool: " +
                               spoolPaths[spool_index]);
    }
    write_spool_header(spool);

    batchReads batch;
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
    std::vector<CompactClassifyResult> localClassifyResults;

    for (;;) {
      if (readQueue.try_dequeue(batch)) {
        const size_t batch_size = batch.ids.size();
        release_queue_slot(queueThrottle, estimate_batch_bytes(batch));
        processBatchCompact(batch, imcfConfig, tax, config, imcf,
                            localClassifyResults, feature_params,
                            feature_min_len, localFileInfo, heat, weightCtx,
                            presencePtr, scratch);
        write_spool_results(localClassifyResults, spool);
        localClassifyResults.clear();
        if (progress != nullptr) {
          progress->processed_reads.fetch_add(batch_size,
                                              std::memory_order_relaxed);
        }
        continue;
      }
      if (producer_done.load(std::memory_order_acquire)) {
        break;
      }
      std::this_thread::yield();
    }

#pragma omp critical
    {
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
                   weightCtx, nullptr, scratch);
    }
#pragma omp critical
    {
      classifyResults.insert(classifyResults.end(),
                             std::make_move_iterator(
                                 localClassifyResults.begin()),
                             std::make_move_iterator(
                                 localClassifyResults.end()));
      std::vector<classifyResult>().swap(localClassifyResults);

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
