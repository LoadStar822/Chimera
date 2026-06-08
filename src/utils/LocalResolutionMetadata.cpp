#include "LocalResolutionMetadata.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace chimera::local_resolution {
namespace {

constexpr char kMagic[] = "CHIMERA_LRMETA_V1";
constexpr uint32_t kVersion = 2;
constexpr uint32_t kMinSupportedVersion = 2;

struct HeaderDisk {
  char magic[sizeof(kMagic)]{};
  uint32_t version{};
  uint64_t species_count{};
  uint64_t target_count{};
  uint64_t species_offset{};
  uint64_t target_offset{};
  uint64_t string_offset{};
};

struct SpeciesDisk {
  uint32_t species{};
  uint32_t genus{};
  uint64_t first_target{};
  uint32_t target_count{};
  uint32_t reserved{};
};

struct TargetDisk {
  uint32_t genus{};
  uint32_t species{};
  uint32_t target_len{};
  uint32_t anchor_count{};
  uint64_t anchor_byte_offset{};
  uint64_t anchor_byte_size{};
  uint32_t rep_rank{};
  uint32_t target_name_len{};
  uint64_t target_name_offset{};
};

template <class T> void write_raw(std::ostream &out, const T &value) {
  out.write(reinterpret_cast<const char *>(&value), sizeof(T));
  if (!out.good()) {
    throw std::runtime_error("failed to write local resolution metadata");
  }
}

template <class T> T read_raw(std::istream &in) {
  T value{};
  in.read(reinterpret_cast<char *>(&value), sizeof(T));
  if (!in.good()) {
    throw std::runtime_error("failed to read local resolution metadata");
  }
  return value;
}

bool target_order(const TargetRep &lhs, const TargetRep &rhs) {
  if (lhs.species != rhs.species) {
    return lhs.species < rhs.species;
  }
  if (lhs.rep_rank != rhs.rep_rank) {
    return lhs.rep_rank < rhs.rep_rank;
  }
  if (lhs.anchor_count != rhs.anchor_count) {
    return lhs.anchor_count > rhs.anchor_count;
  }
  return lhs.target_name < rhs.target_name;
}

HeaderDisk read_header(std::istream &in) {
  HeaderDisk header = read_raw<HeaderDisk>(in);
  if (std::memcmp(header.magic, kMagic, sizeof(kMagic)) != 0 ||
      header.version < kMinSupportedVersion || header.version > kVersion) {
    throw std::runtime_error("invalid local resolution metadata");
  }
  return header;
}

} // namespace

void write_rep_metadata(const std::filesystem::path &path,
                        std::vector<TargetRep> rows) {
  std::sort(rows.begin(), rows.end(), target_order);

  std::vector<SpeciesDisk> species_rows;
  species_rows.reserve(rows.size());
  for (size_t i = 0; i < rows.size();) {
    const uint32_t species = rows[i].species;
    const uint32_t genus = rows[i].genus;
    size_t j = i + 1;
    while (j < rows.size() && rows[j].species == species) {
      ++j;
    }
    if (j - i > std::numeric_limits<uint32_t>::max()) {
      throw std::runtime_error("too many local targets for one species");
    }
    species_rows.push_back(SpeciesDisk{
        species,
        genus,
        static_cast<uint64_t>(i),
        static_cast<uint32_t>(j - i),
        0,
    });
    i = j;
  }

  HeaderDisk header{};
  std::memcpy(header.magic, kMagic, sizeof(kMagic));
  header.version = kVersion;
  header.species_count = species_rows.size();
  header.target_count = rows.size();
  header.species_offset = sizeof(HeaderDisk);
  header.target_offset =
      header.species_offset + species_rows.size() * sizeof(SpeciesDisk);
  header.string_offset =
      header.target_offset + rows.size() * sizeof(TargetDisk);

  std::filesystem::create_directories(path.parent_path().empty()
                                          ? std::filesystem::path(".")
                                          : path.parent_path());
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("failed to open local resolution metadata: " +
                             path.string());
  }

  write_raw(out, header);
  for (const auto &row : species_rows) {
    write_raw(out, row);
  }

  uint64_t string_offset = 0;
  for (const auto &row : rows) {
    if (row.target_name.size() > std::numeric_limits<uint32_t>::max()) {
      throw std::runtime_error("local target name is too large");
    }
    TargetDisk disk{};
    disk.genus = row.genus;
    disk.species = row.species;
    disk.target_len = row.target_len;
    disk.anchor_count = row.anchor_count;
    disk.anchor_byte_offset = row.anchor_byte_offset;
    disk.anchor_byte_size = row.anchor_byte_size;
    disk.rep_rank = row.rep_rank;
    disk.target_name_len = static_cast<uint32_t>(row.target_name.size());
    disk.target_name_offset = string_offset;
    write_raw(out, disk);
    string_offset += disk.target_name_len;
  }

  for (const auto &row : rows) {
    out.write(row.target_name.data(),
              static_cast<std::streamsize>(row.target_name.size()));
    if (!out.good()) {
      throw std::runtime_error("failed to write local resolution strings");
    }
  }
}

RepMetadata RepMetadata::open(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("failed to open local resolution metadata: " +
                             path.string());
  }
  const HeaderDisk header = read_header(in);
  RepMetadata metadata;
  metadata.path_ = path;
  metadata.target_rows_offset_ = header.target_offset;
  metadata.string_data_offset_ = header.string_offset;
  metadata.by_species_.reserve(static_cast<size_t>(header.species_count));

  in.seekg(static_cast<std::streamoff>(header.species_offset));
  if (!in.good()) {
    throw std::runtime_error("failed to seek local resolution species index");
  }
  for (uint64_t i = 0; i < header.species_count; ++i) {
    const SpeciesDisk disk = read_raw<SpeciesDisk>(in);
    metadata.by_species_.emplace(
        disk.species,
        SpeciesEntry{disk.species, disk.genus, disk.first_target,
                     disk.target_count});
  }
  return metadata;
}

std::optional<SpeciesEntry> RepMetadata::find_species(uint32_t species) const {
  const auto found = by_species_.find(species);
  if (found == by_species_.end()) {
    return std::nullopt;
  }
  return found->second;
}

std::vector<TargetRep> RepMetadata::load_targets(uint32_t species,
                                                 size_t max_targets) const {
  auto rows_by_species = load_targets_many({species}, max_targets);
  auto found = rows_by_species.find(species);
  if (found == rows_by_species.end()) {
    return {};
  }
  return std::move(found->second);
}

std::unordered_map<uint32_t, std::vector<TargetRep>>
RepMetadata::load_targets_many(const std::vector<uint32_t> &species_ids,
                               size_t max_targets) const {
  std::unordered_map<uint32_t, std::vector<TargetRep>> result;
  if (max_targets == 0 || species_ids.empty()) {
    return result;
  }

  struct Request {
    uint32_t species{};
    uint64_t first_target{};
    size_t count{};
  };
  std::vector<Request> requests;
  requests.reserve(species_ids.size());
  std::unordered_set<uint32_t> seen;
  seen.reserve(species_ids.size());
  for (uint32_t species : species_ids) {
    if (!seen.insert(species).second) {
      continue;
    }
    const auto found = by_species_.find(species);
    if (found == by_species_.end() || found->second.target_count == 0) {
      continue;
    }
    const SpeciesEntry entry = found->second;
    const size_t count = std::min<size_t>(entry.target_count, max_targets);
    if (count == 0) {
      continue;
    }
    requests.push_back(Request{species, entry.first_target, count});
    result.emplace(species, std::vector<TargetRep>{});
    result[species].reserve(count);
  }
  if (requests.empty()) {
    return result;
  }
  std::sort(requests.begin(), requests.end(), [](const auto &lhs,
                                                 const auto &rhs) {
    return lhs.first_target < rhs.first_target;
  });

  std::ifstream in(path_, std::ios::binary);
  if (!in) {
    throw std::runtime_error("failed to open local resolution metadata: " +
                             path_.string());
  }
  for (const auto &request : requests) {
    std::vector<TargetDisk> disks(request.count);
    const uint64_t byte_offset =
        target_rows_offset_ + request.first_target * sizeof(TargetDisk);
    in.seekg(static_cast<std::streamoff>(byte_offset));
    if (!in.good()) {
      throw std::runtime_error("failed to seek local resolution target rows");
    }
    for (auto &disk : disks) {
      disk = read_raw<TargetDisk>(in);
    }
    auto &rows = result[request.species];
    for (const auto &disk : disks) {
      std::string name(disk.target_name_len, '\0');
      if (disk.target_name_len > 0) {
        const uint64_t name_offset =
            string_data_offset_ + disk.target_name_offset;
        in.seekg(static_cast<std::streamoff>(name_offset));
        if (!in.good()) {
          throw std::runtime_error("failed to seek local resolution strings");
        }
        in.read(name.data(), static_cast<std::streamsize>(name.size()));
        if (!in.good()) {
          throw std::runtime_error(
              "failed to read local resolution target name");
        }
      }
      rows.push_back(TargetRep{disk.genus,
                               disk.species,
                               std::move(name),
                               disk.target_len,
                               disk.anchor_count,
                               disk.anchor_byte_offset,
                               disk.anchor_byte_size,
                               disk.rep_rank});
    }
  }
  return result;
}

} // namespace chimera::local_resolution
