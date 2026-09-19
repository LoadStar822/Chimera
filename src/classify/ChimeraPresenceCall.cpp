#include "ChimeraPresenceCall.hpp"

#include "ChimeraClassifyCommon.hpp"

#include <robin_hood.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace ChimeraClassify::presence_call {

namespace psk = chimera::presence_sketch;

void SampleSketchCollector::absorb(std::vector<uint64_t> &&hashes,
                                   uint64_t sequences, uint64_t bases) {
  std::lock_guard<std::mutex> lock(mutex_);
  sequences_ += sequences;
  bases_ += bases;
  if (!hashes.empty()) {
    parts_.push_back(std::move(hashes));
  }
}

SampleSketch SampleSketchCollector::finalize() {
  std::lock_guard<std::mutex> lock(mutex_);
  SampleSketch sketch;
  sketch.sequences = sequences_;
  sketch.bases = bases_;
  size_t total = 0;
  for (const auto &part : parts_) {
    total += part.size();
  }
  std::vector<uint64_t> all;
  all.reserve(total);
  for (auto &part : parts_) {
    all.insert(all.end(), part.begin(), part.end());
    std::vector<uint64_t>().swap(part);
  }
  parts_.clear();
  sketch.sampled_kmers = all.size();
  std::sort(all.begin(), all.end());
  sketch.keys.reserve(all.size());
  sketch.counts.reserve(all.size());
  for (size_t i = 0; i < all.size();) {
    size_t j = i;
    while (j < all.size() && all[j] == all[i]) {
      ++j;
    }
    sketch.keys.push_back(all[i]);
    sketch.counts.push_back(static_cast<uint32_t>(
        std::min<size_t>(j - i, std::numeric_limits<uint32_t>::max())));
    i = j;
  }
  return sketch;
}

namespace {

// Solve mean_positive = lambda / (1 - exp(-lambda)) for lambda (zero-truncated Poisson).
double ztp_lambda(double mean_positive) {
  if (!(mean_positive > 1.0)) {
    return -1.0;
  }
  double lo = 1e-9;
  double hi = std::max(1.0, mean_positive);
  for (int i = 0; i < 200; ++i) {
    const double mid = 0.5 * (lo + hi);
    const double value = mid / (-std::expm1(-mid));
    if (value - mean_positive > 0.0) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return 0.5 * (lo + hi);
}

// P(X <= h) for X ~ Poisson(mean), summed in log space.
double poisson_lower_tail(double mean, uint64_t h) {
  if (!(mean > 0.0)) {
    return 1.0;
  }
  double log_term = -mean; // log P(X = 0)
  double total = std::exp(log_term);
  for (uint64_t k = 1; k <= h; ++k) {
    log_term += std::log(mean) - std::log(static_cast<double>(k));
    total += std::exp(log_term);
    if (total >= 1.0) {
      return 1.0;
    }
  }
  return std::min(1.0, total);
}

struct Loaded {
  psk::SpeciesMarkers markers;
  std::vector<uint32_t> counts; // sample count per marker
};

// Marker statistics of reference `ref`, restricted to `eligible` when non-empty.
MarkerStats compute_stats(const Loaded &data, size_t ref,
                          const std::vector<char> &eligible, double bases,
                          const CallOptions &options, double survival) {
  MarkerStats s;
  uint64_t positive_sum = 0;
  const auto &m = data.markers;
  const uint64_t genome_bases = m.refs[ref].bases;
  const uint64_t begin = m.ref_begin[ref];
  const uint64_t end = m.ref_begin[ref + 1];
  for (uint64_t j = begin; j < end; ++j) {
    const uint32_t i = m.index[j];
    if (!eligible.empty() && !eligible[i]) {
      continue;
    }
    ++s.markers;
    const uint32_t c = data.counts[i];
    if (c > 0) {
      ++s.positive;
      positive_sum += c;
      if (c >= 2) {
        ++s.repeated;
      }
    }
  }
  if (s.markers == 0) {
    return s;
  }
  s.containment = static_cast<double>(s.positive) / static_cast<double>(s.markers);
  const double effective_survival = survival > 0.0 ? std::min(1.0, survival) : 1.0;
  if (genome_bases > 0) {
    s.lambda_bases = bases / static_cast<double>(genome_bases);
    s.expected_hits = static_cast<double>(s.markers) *
                      (-std::expm1(-s.lambda_bases * effective_survival));
  }
  s.rho = s.expected_hits > 0.0 ? static_cast<double>(s.positive) / s.expected_hits
                                : 0.0;
  // the lower tail only matters when the hits fall short of the expectation
  s.tail_p = s.rho < options.tau ? poisson_lower_tail(s.expected_hits, s.positive) : 1.0;
  if (s.positive > 0) {
    s.mean_positive =
        static_cast<double>(positive_sum) / static_cast<double>(s.positive);
  }
  if (s.positive >= options.min_positive && s.repeated >= options.min_repeated) {
    const double lambda = ztp_lambda(s.mean_positive);
    if (lambda > 0.0) {
      s.lambda_ztp = lambda;
      if (s.mean_positive >= options.min_mean_positive) {
        // The multiplicity estimate of the coverage is inflated at low
        // coverage by markers shared with abundant relatives, which would
        // make a divergent low-abundance strain look like retained markers
        // are missing. The reads assigned to the species bound the coverage
        // from below; an absence call is judged on the smaller estimate,
        // positive evidence for an otherwise unreportable species on the
        // multiplicity estimate alone.
        s.retention_strict = s.containment / (-std::expm1(-lambda));
        double effective_lambda = lambda;
        if (s.lambda_bases > 0.0) {
          effective_lambda =
              std::min(lambda, s.lambda_bases * effective_survival);
        }
        s.retention = s.containment / (-std::expm1(-effective_lambda));
      }
    }
  }
  return s;
}

struct BestRef {
  bool found{false};
  size_t ref{0};
  MarkerStats stats;
};

// Best reference by containment among those keeping >= min_markers eligible markers.
BestRef best_reference(const Loaded &data, const std::vector<char> &eligible,
                       double bases, const CallOptions &options,
                       double survival) {
  BestRef best;
  for (size_t r = 0; r < data.markers.refs.size(); ++r) {
    if (data.markers.ref_key_count(r) < options.min_markers) {
      continue; // cannot reach min_markers even before exclusion
    }
    MarkerStats stats = compute_stats(data, r, eligible, bases, options, survival);
    if (stats.markers < options.min_markers) {
      continue;
    }
    if (!best.found || stats.containment > best.stats.containment ||
        (stats.containment == best.stats.containment &&
         stats.markers > best.stats.markers)) {
      best.found = true;
      best.ref = r;
      best.stats = stats;
    }
  }
  return best;
}

// key -> index (into the exposure list) of the species that claimed it
using ClaimMap = robin_hood::unordered_flat_map<uint64_t, uint32_t>;

// Claims the markers of `ref` for `owner`; earlier claims win.
void claim_reference(const Loaded &data, size_t ref, uint32_t owner,
                     ClaimMap &claimed, std::vector<char> &eligible) {
  const auto &m = data.markers;
  for (uint64_t j = m.ref_begin[ref]; j < m.ref_begin[ref + 1]; ++j) {
    const uint32_t i = m.index[j];
    eligible[i] = 0;
    claimed.try_emplace(m.keys[i], owner);
  }
}

void lookup_counts(const SampleSketch &sample, Loaded &data) {
  const auto &keys = data.markers.keys;
  data.counts.assign(keys.size(), 0);
  auto pos = sample.keys.begin();
  for (size_t i = 0; i < keys.size(); ++i) {
    pos = std::lower_bound(pos, sample.keys.end(), keys[i]);
    if (pos == sample.keys.end()) {
      break;
    }
    if (*pos == keys[i]) {
      data.counts[i] = sample.counts[static_cast<size_t>(pos - sample.keys.begin())];
    }
  }
}

} // namespace

CallResult call_presence(const psk::SketchIndex &index, const SampleSketch &sample,
                         const std::vector<SpeciesExposure> &exposure,
                         const CallOptions &options) {
  CallResult result;
  result.calls.reserve(exposure.size());
  double total_bases = 0.0;
  for (const SpeciesExposure &e : exposure) {
    total_bases += std::max(0.0, e.bases);
  }

  std::vector<Loaded> loaded(exposure.size());
  std::vector<BestRef> full(exposure.size());
  for (size_t i = 0; i < exposure.size(); ++i) {
    const SpeciesExposure &e = exposure[i];
    SpeciesCall call;
    call.species = e.species;
    call.reads = e.reads;
    call.bases = e.bases;
    if (!(e.bases > 0.0) || !index.load_species(e.species, loaded[i].markers)) {
      call.status = "no_sketch";
      result.calls.push_back(std::move(call));
      continue;
    }
    lookup_counts(sample, loaded[i]);
    call.total_markers = loaded[i].markers.keys.size();
    call.references = static_cast<uint32_t>(loaded[i].markers.refs.size());
    full[i] = best_reference(loaded[i], {}, e.bases, options, 1.0);
    if (!full[i].found) {
      call.status = "not_assessable";
    } else {
      call.status = "assessed";
    }
    result.calls.push_back(std::move(call));
  }

  // pass 1: sample k-mer survival from the marker multiplicities of confident species
  {
    std::vector<std::pair<double, double>> calibrators; // (survival, weight)
    auto collect = [&](bool require_fraction) {
      calibrators.clear();
      for (size_t i = 0; i < exposure.size(); ++i) {
        if (!full[i].found) {
          continue;
        }
        const MarkerStats &s = full[i].stats;
        if (s.positive < options.calibration_min_positive ||
            s.repeated < options.calibration_min_repeated ||
            !(s.lambda_ztp > 0.0) || !(s.lambda_bases > 0.0)) {
          continue;
        }
        const double fraction = total_bases > 0.0 ? exposure[i].bases / total_bases : 0.0;
        if (require_fraction && fraction < options.calibration_min_fraction) {
          continue;
        }
        const double survival = std::min(1.0, s.lambda_ztp / s.lambda_bases);
        if (!std::isfinite(survival) || !(survival > 0.0)) {
          continue;
        }
        calibrators.emplace_back(survival, exposure[i].bases);
      }
    };
    collect(true);
    if (calibrators.size() < 3) {
      collect(false);
    }
    result.calibrators = calibrators.size();
    if (!calibrators.empty()) {
      std::sort(calibrators.begin(), calibrators.end());
      double total = 0.0;
      for (const auto &c : calibrators) {
        total += c.second;
      }
      double cumulative = 0.0;
      for (const auto &c : calibrators) {
        cumulative += c.second;
        if (cumulative >= total / 2.0) {
          result.s_hat = c.first;
          break;
        }
      }
    }
  }
  const double survival = result.s_hat > 0.0 ? result.s_hat : 1.0;
  // survival-aware full statistics for the report and the ranking
  for (size_t i = 0; i < exposure.size(); ++i) {
    if (!full[i].found) {
      continue;
    }
    full[i] = best_reference(loaded[i], {}, exposure[i].bases, options, survival);
    result.calls[i].full = full[i].stats;
    result.calls[i].best_ref = loaded[i].markers.refs[full[i].ref].name;
    result.calls[i].best_ref_bases = loaded[i].markers.refs[full[i].ref].bases;
  }

  // pass 2: ranked private-marker calls
  std::vector<size_t> order;
  for (size_t i = 0; i < exposure.size(); ++i) {
    if (full[i].found) {
      order.push_back(i);
    }
  }
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (full[a].stats.containment != full[b].stats.containment) {
      return full[a].stats.containment > full[b].stats.containment;
    }
    return exposure[a].bases > exposure[b].bases;
  });
  ClaimMap claimed;
  result.assessed = order.size();
  // absent when the hits fall far below expectation or the retention is too low
  auto judge = [&](const MarkerStats &s, std::string &reason) -> bool {
    const bool inconsistent = result.s_hat > 0.0 &&
                              s.expected_hits >= options.min_expected_hits &&
                              s.rho < options.tau && s.tail_p < options.tail_alpha;
    const bool low_retention = s.retention >= 0.0 && s.retention < options.min_retention;
    if (inconsistent) {
      reason = low_retention ? "inconsistent+low_retention" : "inconsistent";
    } else if (low_retention) {
      reason = "low_retention";
    }
    return inconsistent || low_retention;
  };
  for (size_t rank = 0; rank < order.size(); ++rank) {
    const size_t i = order[rank];
    SpeciesCall &call = result.calls[i];
    const Loaded &data = loaded[i];
    const size_t n = data.markers.keys.size();
    std::vector<char> priv(n, 1);
    size_t claimed_count = 0;
    robin_hood::unordered_flat_map<uint32_t, size_t> claimants;
    if (!claimed.empty()) {
      for (size_t j = 0; j < n; ++j) {
        const auto it = claimed.find(data.markers.keys[j]);
        if (it != claimed.end()) {
          priv[j] = 0;
          ++claimed_count;
          ++claimants[it->second];
        }
      }
    }
    call.rank = static_cast<int>(rank);
    call.claimed_fraction = n > 0 ? static_cast<double>(claimed_count) / n : 0.0;
    // species that took most of this species' markers
    size_t dominant = exposure.size();
    size_t dominant_count = 0;
    for (const auto &[owner, count] : claimants) {
      if (count > dominant_count) {
        dominant_count = count;
        dominant = owner;
      }
    }
    if (dominant < exposure.size()) {
      call.claimant = exposure[dominant].species;
      call.claimant_containment = full[dominant].stats.containment;
      call.claimant_retention = full[dominant].stats.retention;
    }
    const BestRef best =
        best_reference(data, priv, exposure[i].bases, options, survival);
    bool absent = false;
    if (!best.found) {
      absent = true;
      call.reason = "redundant";
    } else {
      call.priv = best.stats;
      call.best_ref = data.markers.refs[best.ref].name;
      call.best_ref_bases = data.markers.refs[best.ref].bases;
      absent = judge(best.stats, call.reason);
    }
    if (absent && dominant < exposure.size()) {
      // tie: the claimant does not fit the sample better; the classifier's call stands
      std::string ignored;
      const bool own_consistent = !judge(full[i].stats, ignored);
      const double own = full[i].stats.retention;
      const double theirs = full[dominant].stats.retention;
      if (own_consistent && own >= 0.0 && theirs >= 0.0 &&
          theirs - own < options.tie_retention_delta) {
        absent = false;
        call.reason = "tie";
      }
    }
    if (absent) {
      call.status = "absent";
      result.absent.insert(call.species);
      continue;
    }
    call.status = "present";
    result.present.insert(call.species);
    if (!best.found) {
      continue; // tie on fully claimed markers: nothing left to claim
    }
    // the best genome claims its markers, then further well-covered strains
    claim_reference(data, best.ref, static_cast<uint32_t>(i), claimed, priv);
    call.strains = 1;
    while (call.strains < options.max_strains) {
      const BestRef next =
          best_reference(data, priv, exposure[i].bases, options, survival);
      if (!next.found) {
        break;
      }
      const MarkerStats &s = next.stats;
      const bool covered = s.containment >= options.strain_min_containment &&
                           s.positive >= options.min_positive;
      const bool inside = s.retention < 0.0 || s.retention >= options.min_retention;
      if (!covered || !inside) {
        break;
      }
      claim_reference(data, next.ref, static_cast<uint32_t>(i), claimed, priv);
      ++call.strains;
    }
  }
  // species without a judgeable reference stay reportable (not absent)
  for (SpeciesCall &call : result.calls) {
    if (call.status == "assessed") {
      call.status = "not_assessable";
    }
  }
  // positive genome evidence, judged on private markers
  std::unordered_map<uint32_t, size_t> call_index;
  for (size_t i = 0; i < result.calls.size(); ++i) {
    call_index[result.calls[i].species] = i;
  }
  for (SpeciesCall &call : result.calls) {
    if (call.status != "present" || !call.reason.empty()) {
      continue;
    }
    const MarkerStats &s = call.priv;
    if (s.positive < options.min_positive || s.rho < options.evidence_min_rho) {
      continue;
    }
    bool fits = false;
    if (s.retention_strict >= 0.0) {
      fits = s.retention_strict >= options.evidence_min_retention;
    } else if (call.full.containment > 0.0) {
      fits = s.containment / call.full.containment >= options.evidence_min_ratio;
    }
    call.evidence = fits;
  }
  // a tie claimant with fewer reads is the same genome under another name
  for (const SpeciesCall &call : result.calls) {
    if (call.reason != "tie" || call.claimant == 0) {
      continue;
    }
    const auto it = call_index.find(call.claimant);
    if (it != call_index.end() && result.calls[it->second].reads < call.reads) {
      result.calls[it->second].evidence = false;
    }
  }
  for (const SpeciesCall &call : result.calls) {
    if (call.evidence) {
      result.evidence.insert(call.species);
    }
  }
  return result;
}

void write_call_table(const std::string &path, const CallResult &result,
                      const NcbiTaxdump *taxdump) {
  const auto parent = std::filesystem::path(path).parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(path, std::ios::out | std::ios::binary);
  if (!out) {
    throw std::runtime_error("failed to open presence call table: " + path);
  }
  out << "species_taxid\tname\tstatus\treason\tevidence\trank\treads\tassigned_bases"
      << "\ttotal_markers\treferences\tstrains\tclaimed_fraction\tclaimant\tclaimant_containment\tclaimant_retention\tbest_ref\tbest_ref_bases"
      << "\tfull_markers\tfull_hits\tfull_containment\tfull_expected_hits\tfull_rho\tfull_retention"
      << "\tprivate_markers\tprivate_hits\tprivate_containment\tprivate_expected_hits\tprivate_rho\tprivate_retention"
      << "\tlambda_bases\tprivate_lambda_ztp\tprivate_mean_positive\tprivate_tail_p\ts_hat\n";
  auto name_of = [&](uint32_t species) -> std::string {
    if (taxdump != nullptr && taxdump->enabled() &&
        species < taxdump->scientific_name.size()) {
      return taxdump->scientific_name[species];
    }
    return {};
  };
  auto fmt = [](double v) -> std::string {
    if (!std::isfinite(v)) {
      return "nan";
    }
    std::ostringstream oss;
    oss << std::setprecision(6) << v;
    return oss.str();
  };
  std::vector<const SpeciesCall *> rows;
  rows.reserve(result.calls.size());
  for (const SpeciesCall &call : result.calls) {
    rows.push_back(&call);
  }
  std::sort(rows.begin(), rows.end(), [](const SpeciesCall *a, const SpeciesCall *b) {
    return a->bases > b->bases;
  });
  for (const SpeciesCall *call : rows) {
    out << call->species << '\t' << name_of(call->species) << '\t' << call->status << '\t'
        << call->reason << '\t' << (call->evidence ? 1 : 0) << '\t' << call->rank << '\t'
        << fmt(call->reads) << '\t'
        << fmt(call->bases) << '\t' << call->total_markers << '\t'
        << call->references << '\t' << call->strains << '\t'
        << fmt(call->claimed_fraction) << '\t' << call->claimant << '\t'
        << fmt(call->claimant_containment) << '\t'
        << (call->claimant_retention < 0.0 ? "nan" : fmt(call->claimant_retention)) << '\t'
        << call->best_ref << '\t'
        << call->best_ref_bases << '\t' << call->full.markers << '\t'
        << call->full.positive << '\t' << fmt(call->full.containment) << '\t'
        << fmt(call->full.expected_hits) << '\t' << fmt(call->full.rho) << '\t'
        << (call->full.retention < 0.0 ? "nan" : fmt(call->full.retention)) << '\t'
        << call->priv.markers << '\t' << call->priv.positive << '\t'
        << fmt(call->priv.containment) << '\t' << fmt(call->priv.expected_hits) << '\t'
        << fmt(call->priv.rho) << '\t'
        << (call->priv.retention < 0.0 ? "nan" : fmt(call->priv.retention)) << '\t'
        << fmt(call->full.lambda_bases) << '\t'
        << (call->priv.lambda_ztp < 0.0 ? "nan" : fmt(call->priv.lambda_ztp)) << '\t'
        << fmt(call->priv.mean_positive) << '\t' << fmt(call->priv.tail_p) << '\t'
        << fmt(result.s_hat) << '\n';
  }
}

} // namespace ChimeraClassify::presence_call
