#include "ChimeraClassifyReadout.hpp"

#include "ChimeraClassifyCommon.hpp"

#include <utils/Parse.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace ChimeraClassify::readout {

bool apply_deferred_abstention_release(classifyResult &result,
                                       const NcbiTaxdump *ncbiTaxdump) {
  if (!ncbiTaxdump || !ncbiTaxdump->enabled() ||
      result.sampleMixturePosteriors.empty() || result.taxidCount.empty() ||
      result.taxidCount.front().first != "unclassified" ||
      result.reject_reason != "posterior_weight") {
    return false;
  }

  std::unordered_map<uint32_t, double> speciesMass;
  std::unordered_map<uint32_t, std::pair<std::string, double>>
      bestMemberBySpecies;
  speciesMass.reserve(result.sampleMixturePosteriors.size());
  bestMemberBySpecies.reserve(result.sampleMixturePosteriors.size());

  for (const auto &[taxidText, posterior] : result.sampleMixturePosteriors) {
    if (!(posterior > 0.0)) {
      continue;
    }
    uint32_t taxid = 0;
    if (!chimera::utils::try_parse_u32(taxidText, taxid)) {
      continue;
    }
    const uint32_t species = ncbiTaxdump->to_species(taxid);
    if (species == 0) {
      continue;
    }
    speciesMass[species] += posterior;
    auto it = bestMemberBySpecies.find(species);
    if (it == bestMemberBySpecies.end() || posterior > it->second.second ||
        (posterior == it->second.second && taxidText < it->second.first)) {
      bestMemberBySpecies[species] = {taxidText, posterior};
    }
  }
  if (speciesMass.empty()) {
    return false;
  }

  uint32_t topSpecies = 0;
  double topSpeciesMass = 0.0;
  for (const auto &[species, mass] : speciesMass) {
    if (mass > topSpeciesMass ||
        (mass == topSpeciesMass && species < topSpecies)) {
      topSpecies = species;
      topSpeciesMass = mass;
    }
  }
  auto memberIt = bestMemberBySpecies.find(topSpecies);
  if (memberIt == bestMemberBySpecies.end() ||
      memberIt->second.first.empty()) {
    return false;
  }

  constexpr double eps = 1e-12;
  const double p = std::clamp(topSpeciesMass, eps, 1.0 - eps);
  const double releaseGain = std::log(p / (1.0 - p));
  const double releaseCost =
      0.5 * std::log1p(std::max(0.0, result.sampleMixtureTopScore));
  if (!(releaseGain > releaseCost)) {
    return false;
  }

  result.taxidCount.clear();
  result.taxidCount.emplace_back(memberIt->second.first, topSpeciesMass);
  result.posteriors = result.sampleMixturePosteriors;
  result.reject_reason.clear();
  return true;
}

bool apply_local_correction(classifyResult &result,
                            const NcbiTaxdump *ncbiTaxdump) {
  if (!ncbiTaxdump || !ncbiTaxdump->enabled() ||
      result.sampleMixturePosteriors.empty() || result.taxidCount.empty()) {
    return false;
  }

  std::unordered_map<uint32_t, double> speciesMass;
  std::unordered_map<uint32_t, double> genusMass;
  std::unordered_map<uint32_t, std::pair<std::string, double>>
      bestMemberBySpecies;
  speciesMass.reserve(result.sampleMixturePosteriors.size());
  genusMass.reserve(result.sampleMixturePosteriors.size());
  bestMemberBySpecies.reserve(result.sampleMixturePosteriors.size());

  for (const auto &[taxidText, posterior] : result.sampleMixturePosteriors) {
    if (!(posterior > 0.0)) {
      continue;
    }
    uint32_t taxid = 0;
    if (!chimera::utils::try_parse_u32(taxidText, taxid)) {
      continue;
    }
    const uint32_t species = ncbiTaxdump->to_species(taxid);
    const uint32_t genus = ncbiTaxdump->to_genus(taxid);
    if (species != 0) {
      speciesMass[species] += posterior;
      auto it = bestMemberBySpecies.find(species);
      if (it == bestMemberBySpecies.end() || posterior > it->second.second ||
          (posterior == it->second.second && taxidText < it->second.first)) {
        bestMemberBySpecies[species] = {taxidText, posterior};
      }
    }
    if (genus != 0) {
      genusMass[genus] += posterior;
    }
  }
  if (genusMass.empty() || speciesMass.empty()) {
    return false;
  }

  uint32_t topGenus = 0;
  double topGenusMass = 0.0;
  for (const auto &[genus, mass] : genusMass) {
    if (mass > topGenusMass || (mass == topGenusMass && genus < topGenus)) {
      topGenus = genus;
      topGenusMass = mass;
    }
  }
  if (topGenus == 0) {
    return false;
  }

  uint32_t topSpecies = 0;
  double topSpeciesMass = 0.0;
  for (const auto &[species, mass] : speciesMass) {
    if (ncbiTaxdump->to_genus(species) != topGenus) {
      continue;
    }
    if (mass > topSpeciesMass ||
        (mass == topSpeciesMass && species < topSpecies)) {
      topSpecies = species;
      topSpeciesMass = mass;
    }
  }
  if (topSpecies == 0 || !(topSpeciesMass > 0.0)) {
    return false;
  }

  uint32_t oldTaxid = 0;
  if (!chimera::utils::try_parse_u32(result.taxidCount.front().first,
                                    oldTaxid)) {
    return false;
  }
  if (oldTaxid == 0) {
    return false;
  }
  const uint32_t oldSpecies = ncbiTaxdump->to_species(oldTaxid);
  if (oldSpecies == topSpecies) {
    return false;
  }

  const auto oldMassIt = speciesMass.find(oldSpecies);
  const double oldSpeciesMass =
      (oldMassIt == speciesMass.end()) ? 0.0 : oldMassIt->second;

  std::vector<double> localScores;
  localScores.reserve(result.sampleMixtureLocalScores.size());
  for (const auto &[_, score] : result.sampleMixtureLocalScores) {
    if (score > 0.0) {
      localScores.push_back(score);
    }
  }
  std::sort(localScores.begin(), localScores.end(), std::greater<double>());
  double localRatio = 1.0;
  if (localScores.size() >= 2 && localScores[1] > 0.0) {
    localRatio = localScores[0] / localScores[1];
  } else if (!localScores.empty()) {
    localRatio = 1.0 + std::max(0.0, localScores[0]);
  }
  if (!std::isfinite(localRatio)) {
    localRatio = 1.0;
  }
  constexpr double eps = 1e-12;
  const double switchGain =
      std::log(std::max(topSpeciesMass, eps) / std::max(oldSpeciesMass, eps));
  const double localCost = std::log(std::max(1.0, localRatio));
  const double speciesModelCost =
      0.5 * std::log(std::max(1.0, static_cast<double>(speciesMass.size())));
  const double switchCost = std::max(localCost, speciesModelCost);
  if (!(switchGain > switchCost)) {
    return false;
  }

  auto memberIt = bestMemberBySpecies.find(topSpecies);
  if (memberIt == bestMemberBySpecies.end() || memberIt->second.first.empty() ||
      memberIt->second.first == result.taxidCount.front().first) {
    return false;
  }

  result.taxidCount.clear();
  result.taxidCount.emplace_back(memberIt->second.first, topSpeciesMass);
  result.posteriors = result.sampleMixturePosteriors;
  return true;
}

bool apply_selective_readout(classifyResult &result,
                             const NcbiTaxdump *ncbiTaxdump) {
  apply_deferred_abstention_release(result, ncbiTaxdump);
  return apply_local_correction(result, ncbiTaxdump);
}

} // namespace ChimeraClassify::readout
