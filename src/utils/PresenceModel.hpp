/*
 * Lightweight metadata container for the coverage-based presence model.
 *
 * Stored inside the .imcf archive so classify can reuse per-genome signature
 * statistics without re-scanning the database.
 */
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

namespace chimera::presence {

struct CoverageEntry {
  std::string taxid;
  uint64_t unique_signatures{0};
  uint64_t total_signatures{0};

  template <class Archive> void serialize(Archive &archive) {
    archive(taxid, unique_signatures, total_signatures);
  }
};

struct CoverageMeta {
  uint32_t version{1};
  uint32_t unique_deg_threshold{1};
  std::vector<CoverageEntry> entries;

  template <class Archive> void serialize(Archive &archive) {
    archive(version, unique_deg_threshold, entries);
  }

  bool empty() const { return entries.empty(); }
};

} // namespace chimera::presence

