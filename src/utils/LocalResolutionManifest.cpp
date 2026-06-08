#include "LocalResolutionManifest.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace chimera::local_resolution {
namespace {

constexpr const char *kManifestMagic = "chimera_local_resolution_manifest_v1";

std::vector<std::string> split_tab(const std::string &line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) {
    fields.push_back(std::move(field));
  }
  return fields;
}

int64_t artifact_mtime(const std::filesystem::path &path) {
  return std::filesystem::last_write_time(path).time_since_epoch().count();
}

ArtifactStamp stamp_artifact(const std::filesystem::path &base_dir,
                             const std::filesystem::path &path) {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("local resolution artifact is missing: " +
                             path.string());
  }
  ArtifactStamp stamp;
  stamp.relative_path = std::filesystem::relative(path, base_dir);
  stamp.size = static_cast<uint64_t>(std::filesystem::file_size(path));
  stamp.mtime = artifact_mtime(path);
  return stamp;
}

std::filesystem::path resolve_shard_path_from_manifest(
    const std::filesystem::path &shard_manifest_path,
    const std::string &raw_path) {
  std::filesystem::path path(raw_path);
  if (path.is_absolute()) {
    return path;
  }
  const auto parent = shard_manifest_path.parent_path().empty()
                          ? std::filesystem::path(".")
                          : shard_manifest_path.parent_path();
  const auto shard_dir = parent / shard_manifest_path.stem().string();
  const std::array<std::filesystem::path, 3> candidates{
      parent / path, shard_dir / path, shard_dir / path.filename()};
  for (const auto &candidate : candidates) {
    if (std::filesystem::exists(candidate)) {
      return candidate;
    }
  }
  return shard_dir / path.filename();
}

std::vector<ArtifactStamp>
stamp_shards(const std::filesystem::path &base_dir,
             const std::filesystem::path &shard_manifest_path) {
  std::ifstream in(shard_manifest_path);
  if (!in) {
    throw std::runtime_error("failed to open local resolution shard manifest: " +
                             shard_manifest_path.string());
  }
  std::vector<ArtifactStamp> shards;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    const auto fields = split_tab(line);
    if (fields.empty() || fields[0] == "genus") {
      continue;
    }
    if (fields.size() < 2) {
      throw std::runtime_error("invalid local resolution shard manifest row");
    }
    shards.push_back(stamp_artifact(
        base_dir, resolve_shard_path_from_manifest(shard_manifest_path,
                                                   fields[1])));
  }
  return shards;
}

void write_stamp(std::ostream &out,
                 const char *name,
                 const ArtifactStamp &stamp) {
  out << name << '\t' << stamp.relative_path.generic_string() << '\t'
      << stamp.size << '\t' << stamp.mtime << '\n';
}

ArtifactStamp parse_stamp(const std::vector<std::string> &fields,
                          const std::string &label) {
  if (fields.size() != 4) {
    throw std::runtime_error("invalid local resolution manifest row: " + label);
  }
  ArtifactStamp stamp;
  stamp.relative_path = fields[1];
  stamp.size = static_cast<uint64_t>(std::stoull(fields[2]));
  stamp.mtime = static_cast<int64_t>(std::stoll(fields[3]));
  return stamp;
}

void verify_stamp(const std::filesystem::path &core_path,
                  const ArtifactStamp &stamp,
                  const std::string &label) {
  const std::filesystem::path path =
      materialize_manifest_path(core_path, stamp);
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("local resolution " + label +
                             " is missing: " + path.string());
  }
  const uint64_t size = static_cast<uint64_t>(std::filesystem::file_size(path));
  if (size != stamp.size) {
    throw std::runtime_error("local resolution " + label +
                             " does not match the database manifest: " +
                             path.string());
  }
}

uint32_t parse_u32_field(const std::unordered_map<std::string, std::string> &kv,
                         const std::string &key) {
  const auto found = kv.find(key);
  if (found == kv.end()) {
    return 0;
  }
  return static_cast<uint32_t>(std::stoul(found->second));
}

} // namespace

std::filesystem::path
core_archive_path_for(const std::filesystem::path &db_path) {
  if (std::filesystem::exists(db_path) &&
      std::filesystem::is_directory(db_path)) {
    return db_path / "core.imcf";
  }
  if (std::filesystem::exists(db_path) &&
      db_path.extension() == ".imcf") {
    return db_path;
  }
  std::filesystem::path with_extension = db_path;
  with_extension.replace_extension(".imcf");
  if (std::filesystem::exists(with_extension)) {
    return with_extension;
  }
  if (std::filesystem::exists(db_path)) {
    return db_path;
  }
  return with_extension;
}

std::filesystem::path
manifest_path_for_core(const std::filesystem::path &core_path) {
  if (core_path.filename() == "core.imcf") {
    return core_path.parent_path() / "manifest.tsv";
  }
  return core_path.parent_path() /
         (core_path.stem().string() + ".local_manifest.tsv");
}

void write_manifest(const std::filesystem::path &core_path,
                    bool local_available,
                    const std::filesystem::path &local_index_path,
                    const std::filesystem::path &rep_metadata_path,
                    const std::filesystem::path &shard_manifest_path,
                    uint32_t k,
                    uint32_t w,
                    uint32_t targets_per_species) {
  if (!std::filesystem::exists(core_path)) {
    throw std::runtime_error("cannot write local resolution manifest without core archive: " +
                             core_path.string());
  }
  const std::filesystem::path manifest_path = manifest_path_for_core(core_path);
  const std::filesystem::path base_dir =
      core_path.parent_path().empty() ? std::filesystem::path(".")
                                      : core_path.parent_path();
  std::ofstream out(manifest_path, std::ios::trunc);
  if (!out) {
    throw std::runtime_error("failed to open local resolution manifest: " +
                             manifest_path.string());
  }
  out << "format\t" << kManifestMagic << '\n';
  out << "local_available\t" << (local_available ? 1 : 0) << '\n';
  out << "k\t" << k << '\n';
  out << "w\t" << w << '\n';
  out << "targets_per_species\t" << targets_per_species << '\n';
  write_stamp(out, "core", stamp_artifact(base_dir, core_path));
  if (local_available) {
    write_stamp(out, "local_index", stamp_artifact(base_dir, local_index_path));
    write_stamp(out, "rep_metadata",
                stamp_artifact(base_dir, rep_metadata_path));
    write_stamp(out, "shard_manifest",
                stamp_artifact(base_dir, shard_manifest_path));
    const auto shards = stamp_shards(base_dir, shard_manifest_path);
    for (const auto &shard : shards) {
      write_stamp(out, "shard", shard);
    }
  }
}

std::optional<BuildManifest>
load_and_verify_manifest_for_db(const std::filesystem::path &db_path) {
  const std::filesystem::path core_path = core_archive_path_for(db_path);
  const std::filesystem::path manifest_path = manifest_path_for_core(core_path);
  if (!std::filesystem::exists(manifest_path)) {
    return std::nullopt;
  }

  std::ifstream in(manifest_path);
  if (!in) {
    throw std::runtime_error("failed to open local resolution manifest: " +
                             manifest_path.string());
  }
  BuildManifest manifest;
  std::unordered_map<std::string, std::string> scalars;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    const auto fields = split_tab(line);
    if (fields.size() < 2) {
      throw std::runtime_error("invalid local resolution manifest row");
    }
    const std::string &key = fields[0];
    if (key == "format") {
      if (fields.size() != 2 || fields[1] != kManifestMagic) {
        throw std::runtime_error("invalid local resolution manifest format");
      }
    } else if (key == "core") {
      manifest.core = parse_stamp(fields, key);
    } else if (key == "local_index") {
      manifest.local_index = parse_stamp(fields, key);
    } else if (key == "rep_metadata") {
      manifest.rep_metadata = parse_stamp(fields, key);
    } else if (key == "shard_manifest") {
      manifest.shard_manifest = parse_stamp(fields, key);
    } else if (key == "shard") {
      manifest.shards.push_back(parse_stamp(fields, key));
    } else if (fields.size() == 2) {
      scalars.emplace(key, fields[1]);
    } else {
      throw std::runtime_error("invalid local resolution manifest row: " + key);
    }
  }

  manifest.local_available = scalars["local_available"] == "1";
  manifest.k = parse_u32_field(scalars, "k");
  manifest.w = parse_u32_field(scalars, "w");
  manifest.targets_per_species =
      parse_u32_field(scalars, "targets_per_species");

  verify_stamp(core_path, manifest.core, "core archive");
  if (manifest.local_available) {
    verify_stamp(core_path, manifest.local_index, "index");
    verify_stamp(core_path, manifest.rep_metadata, "metadata");
    verify_stamp(core_path, manifest.shard_manifest, "shard manifest");
    for (size_t i = 0; i < manifest.shards.size(); ++i) {
      verify_stamp(core_path, manifest.shards[i],
                   "shard " + std::to_string(i));
    }
  }
  return manifest;
}

std::filesystem::path
materialize_manifest_path(const std::filesystem::path &core_path,
                          const ArtifactStamp &stamp) {
  if (stamp.relative_path.is_absolute()) {
    return stamp.relative_path;
  }
  const std::filesystem::path base_dir =
      core_path.parent_path().empty() ? std::filesystem::path(".")
                                      : core_path.parent_path();
  return base_dir / stamp.relative_path;
}

} // namespace chimera::local_resolution
