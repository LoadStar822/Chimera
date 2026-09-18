#pragma once

#include <utils/PresenceSketch.hpp>

#include <cstdint>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace ChimeraClassify {
struct NcbiTaxdump;
}

namespace ChimeraClassify::presence_call {

// Sampled k-mers of the whole read set (with multiplicities).
struct SampleSketch {
  std::vector<uint64_t> keys;   // sorted unique
  std::vector<uint32_t> counts; // occurrences per key
  uint64_t sampled_kmers{0};
  uint64_t sequences{0};
  uint64_t bases{0};
};

// Collects per-thread hash buffers during streaming classification.
class SampleSketchCollector {
public:
  explicit SampleSketchCollector(const chimera::presence_sketch::Params &params)
      : params_(params) {}
  const chimera::presence_sketch::Params &params() const { return params_; }
  // Called by each worker thread with its local buffer once it is done.
  void absorb(std::vector<uint64_t> &&hashes, uint64_t sequences, uint64_t bases);
  SampleSketch finalize();

private:
  chimera::presence_sketch::Params params_;
  std::mutex mutex_;
  std::vector<std::vector<uint64_t>> parts_;
  uint64_t sequences_{0};
  uint64_t bases_{0};
};

struct SpeciesExposure {
  uint32_t species{0};
  double reads{0.0};
  double bases{0.0};
};

struct CallOptions {
  // consistency: observed private hits H vs. H_exp = M (1 - exp(-lambda s))
  double tau{0.10};                 // absent when H / H_exp < tau ...
  double min_expected_hits{5.0};    // ... and H_exp >= this ...
  double tail_alpha{1e-3};          // ... and P(Poisson(H_exp) <= H) < alpha
  // retention: markers retained after correcting for the species' own coverage
  double min_retention{0.20};
  uint32_t min_positive{30};        // retention needs H >= this ...
  uint32_t min_repeated{5};         // ... this many markers seen twice ...
  double min_mean_positive{1.2};    // ... and a mean multiplicity that pins lambda
  uint32_t min_markers{100};        // reference needs this many (private) markers to be judged
  // further strains of a present species claim their private markers too
  uint32_t max_strains{8};
  double strain_min_containment{0.10};
  // tie: stay present when the claimant's retention is not better by delta
  double tie_retention_delta{0.05};
  // evidence: private markers hit consistently with the reads on a fitting genome
  double evidence_min_rho{0.5};
  double evidence_min_retention{0.5};
  double evidence_min_ratio{0.6};
  // survival calibration over confident species
  uint32_t calibration_min_positive{200};
  uint32_t calibration_min_repeated{20};
  double calibration_min_fraction{0.01};
};

struct MarkerStats {
  uint64_t markers{0};       // M
  uint64_t positive{0};      // H
  uint64_t repeated{0};      // markers with count >= 2
  double containment{0.0};   // H / M
  double lambda_bases{0.0};  // assigned bases / G
  double expected_hits{0.0}; // M (1 - exp(-lambda_bases * survival))
  double rho{0.0};           // H / H_exp
  double tail_p{1.0};        // P(Poisson(H_exp) <= H)
  double mean_positive{0.0}; // mean multiplicity of positive markers
  double lambda_ztp{-1.0};   // effective coverage from multiplicities (<0: n/a)
  double retention{-1.0};    // containment / (1 - exp(-lambda_ztp)) (<0: n/a)
};

struct SpeciesCall {
  uint32_t species{0};
  double reads{0.0};
  double bases{0.0};
  std::string status;       // no_sketch | not_assessable | present | absent
  std::string reason;       // redundant | inconsistent | low_retention | tie | (empty)
  bool evidence{false};     // present with positive genome evidence (see CallOptions)
  int rank{-1};
  double claimed_fraction{0.0};
  uint32_t claimant{0};              // species that claimed most of the markers (0: none)
  double claimant_containment{0.0};  // its best-genome containment
  double claimant_retention{-1.0};   // its best-genome retention (<0: n/a)
  uint64_t total_markers{0};
  uint32_t references{0};   // genomes of the species in the sketch
  uint32_t strains{0};      // genomes that claimed markers (present species)
  std::string best_ref;
  uint64_t best_ref_bases{0};
  MarkerStats full;    // best reference, all own markers
  MarkerStats priv;    // best reference, private markers
};

struct CallResult {
  std::vector<SpeciesCall> calls;
  std::unordered_set<uint32_t> present;
  std::unordered_set<uint32_t> absent;
  std::unordered_set<uint32_t> evidence; // present species with positive genome evidence
  double s_hat{-1.0};       // per-k-mer survival estimate (<0: not calibrated)
  size_t calibrators{0};
  size_t assessed{0};
};

CallResult call_presence(const chimera::presence_sketch::SketchIndex &index,
                         const SampleSketch &sample,
                         const std::vector<SpeciesExposure> &exposure,
                         const CallOptions &options);

void write_call_table(const std::string &path, const CallResult &result,
                      const NcbiTaxdump *taxdump);

} // namespace ChimeraClassify::presence_call
