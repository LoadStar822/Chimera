#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <vector>

namespace chimera::local_resolution {

struct ArtifactStamp {
  std::filesystem::path relative_path;
  uint64_t size{};
  int64_t mtime{};
};

struct BuildManifest {
  bool local_available{};
  ArtifactStamp core;
  ArtifactStamp local_index;
  ArtifactStamp rep_metadata;
  ArtifactStamp shard_manifest;
  std::vector<ArtifactStamp> shards;
  uint32_t k{};
  uint32_t w{};
  uint32_t targets_per_species{};
  bool presence_available{};
  ArtifactStamp presence_sketch;
  // Species all of whose input sequences are local targets.
  bool complete_species_available{};
  ArtifactStamp complete_species;
};

std::filesystem::path core_archive_path_for(const std::filesystem::path &db_path);
std::filesystem::path manifest_path_for_core(const std::filesystem::path &core_path);

void write_manifest(const std::filesystem::path &core_path,
                    bool local_available,
                    const std::filesystem::path &local_index_path,
                    const std::filesystem::path &rep_metadata_path,
                    const std::filesystem::path &shard_manifest_path,
                    const std::filesystem::path &complete_species_path,
                    uint32_t k,
                    uint32_t w,
                    uint32_t targets_per_species);

// Parses the manifest without checking the artifacts.
std::optional<BuildManifest>
load_manifest_for_db(const std::filesystem::path &db_path);

// Parses the manifest and checks the core and local resolution artifacts.
std::optional<BuildManifest>
load_and_verify_manifest_for_db(const std::filesystem::path &db_path);

// Registers (or replaces) the presence sketch row; false when the database has
// no manifest.
bool stamp_presence_sketch(const std::filesystem::path &core_path,
                           const std::filesystem::path &sketch_path);

std::filesystem::path materialize_manifest_path(
    const std::filesystem::path &core_path,
    const ArtifactStamp &stamp);

} // namespace chimera::local_resolution
