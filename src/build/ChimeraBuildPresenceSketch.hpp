#pragma once

#include "buildConfig.hpp"

#include <cstdint>
#include <filesystem>
#include <string>

namespace ChimeraBuild {

struct PresenceSketchBuildOptions {
  std::filesystem::path input_file;   // build input: <fasta path> <taxid> per line
  std::filesystem::path output_path;  // sketch.psk destination
  std::string taxonomy_dir;           // directory with nodes.dmp (may be empty -> env)
  uint64_t scaled{1000};
  uint32_t max_refs{0};               // 0 = keep every genome (recommended)
  uint16_t threads{1};
  bool verbose{true};
};

struct PresenceSketchBuildStats {
  uint64_t genomes{0};
  uint64_t unreadable_genomes{0};
  uint64_t species{0};
  uint64_t species_capped{0};   // species with more genomes than max_refs
  uint64_t references{0};       // genomes kept in the sketch
  uint64_t keys{0};
  uint64_t bases{0};
  double seconds{0.0};
};

PresenceSketchBuildStats
build_presence_sketch(const PresenceSketchBuildOptions &options);

// Default taxonomy directory for a database when none is given explicitly.
std::string default_presence_taxonomy_dir(const std::filesystem::path &db_path);

} // namespace ChimeraBuild
