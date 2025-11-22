/*
 * Lightweight metadata container for the coverage-based presence model.
 *
 * Stored inside the .imcf archive so classify can reuse per-genome signature
 * statistics without re-scanning the database.
 */
#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

namespace chimera::presence {

struct HashFreqStats {
  uint32_t df_high_threshold{std::numeric_limits<uint32_t>::max()};
  uint32_t df_max_observed{0};
  uint64_t nonzero_counters{0};

  template <class Archive> void serialize(Archive &archive) {
    archive(df_high_threshold, df_max_observed, nonzero_counters);
  }
};

struct HashFrequencyModel {
  uint32_t depth{0};
  uint32_t width{0};
  double quantile{0.0};
  HashFreqStats stats{};
  uint64_t total_hashes{0};
  uint64_t filtered_hashes{0};
  std::vector<uint32_t> counters;

  bool enabled() const {
    return depth > 0 && width > 0 && !counters.empty();
  }

  template <class Archive> void serialize(Archive &archive) {
    archive(depth, width, quantile, stats, total_hashes, filtered_hashes,
            counters);
  }
};

struct CoverageEntry {
  std::string taxid;
  uint64_t unique_signatures{0};
  uint64_t total_signatures{0};
  uint64_t genome_length{0};
  double unique_density{0.0};
  double expected_unique_per_ref_read{0.0};

  template <class Archive> void serialize(Archive &archive) {
    archive(taxid, unique_signatures, total_signatures, genome_length,
            unique_density, expected_unique_per_ref_read);
  }
};

struct CoverageEntryV1 {
  std::string taxid;
  uint64_t unique_signatures{0};
  uint64_t total_signatures{0};

  template <class Archive> void serialize(Archive &archive) {
    archive(taxid, unique_signatures, total_signatures);
  }
};

struct CoverageMeta {
  uint32_t unique_deg_threshold{1};
  uint16_t ref_read_length{150};
  uint16_t effective_span{0};
  HashFrequencyModel freq_model;
  std::vector<CoverageEntry> entries;

  template <class Archive> void serialize(Archive &archive) {
    archive(unique_deg_threshold, ref_read_length, effective_span, freq_model,
            entries);
  }

  bool empty() const { return entries.empty(); }
};

} // namespace chimera::presence
