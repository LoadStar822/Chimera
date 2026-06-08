#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace chimera::local_resolution {

struct TargetRep {
  uint32_t genus{};
  uint32_t species{};
  std::string target_name;
  uint32_t target_len{};
  uint32_t anchor_count{};
  uint64_t anchor_byte_offset{};
  uint64_t anchor_byte_size{};
  uint32_t rep_rank{};
};

struct SpeciesEntry {
  uint32_t species{};
  uint32_t genus{};
  uint64_t first_target{};
  uint32_t target_count{};
};

class RepMetadata {
public:
  static RepMetadata open(const std::filesystem::path &path);

  std::optional<SpeciesEntry> find_species(uint32_t species) const;
  std::vector<TargetRep> load_targets(uint32_t species,
                                      size_t max_targets) const;
  std::unordered_map<uint32_t, std::vector<TargetRep>>
  load_targets_many(const std::vector<uint32_t> &species_ids,
                    size_t max_targets) const;

private:
  std::filesystem::path path_;
  uint64_t target_rows_offset_{};
  uint64_t string_data_offset_{};
  std::unordered_map<uint32_t, SpeciesEntry> by_species_;
};

void write_rep_metadata(const std::filesystem::path &path,
                        std::vector<TargetRep> rows);

} // namespace chimera::local_resolution
