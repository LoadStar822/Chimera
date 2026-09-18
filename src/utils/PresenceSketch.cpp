#include "PresenceSketch.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <stdexcept>

namespace chimera::presence_sketch {

namespace {

constexpr char kMagic[8] = {'C', 'H', 'M', 'P', 'S', 'K', '0', '2'};
constexpr uint32_t kFormatVersion = 2;
constexpr size_t kHeaderBytes = 64;

inline uint64_t mix64(uint64_t x, uint64_t seed) {
  x ^= seed;
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ull;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebull;
  x ^= x >> 31;
  return x;
}

template <typename T> void write_pod(std::ostream &os, const T &value) {
  os.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

template <typename T> bool read_pod(std::istream &is, T &value) {
  is.read(reinterpret_cast<char *>(&value), sizeof(T));
  return static_cast<bool>(is);
}

} // namespace

void sample_hashes(const std::vector<seqan3::dna4> &seq, const Params &params,
                   std::vector<uint64_t> &out) {
  const uint32_t k = params.k;
  if (k < 2 || k > 32 || seq.size() < k) {
    return;
  }
  const uint64_t threshold = params.threshold();
  const uint64_t mask = (k == 32) ? ~0ull : ((1ull << (2 * k)) - 1ull);
  const uint32_t rc_shift = 2 * (k - 1);
  uint64_t fwd = 0;
  uint64_t rc = 0;
  uint32_t filled = 0;
  for (const seqan3::dna4 base : seq) {
    const uint64_t r = static_cast<uint64_t>(seqan3::to_rank(base)) & 3ull;
    fwd = ((fwd << 2) | r) & mask;
    rc = (rc >> 2) | ((3ull - r) << rc_shift);
    if (++filled < k) {
      continue;
    }
    const uint64_t canonical = std::min(fwd, rc);
    const uint64_t h = mix64(canonical, params.seed);
    if (h <= threshold) {
      out.push_back(h);
    }
  }
}

void sort_unique(std::vector<uint64_t> &hashes) {
  std::sort(hashes.begin(), hashes.end());
  hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
}

// ----------------------------------------------------------------------------
// Writer

SketchWriter::SketchWriter(const std::filesystem::path &path,
                           const Params &params, uint32_t max_refs,
                           uint64_t input_genomes)
    : path_(path), params_(params), max_refs_(max_refs),
      input_genomes_(input_genomes) {
  const auto parent = path_.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  out_.open(path_, std::ios::binary | std::ios::trunc);
  if (!out_) {
    throw std::runtime_error("failed to open presence sketch for writing: " +
                             path_.string());
  }
  write_header();
}

SketchWriter::~SketchWriter() {
  if (!finished_) {
    try {
      finish();
    } catch (...) {
    }
  }
}

void SketchWriter::write_header() {
  out_.seekp(0);
  out_.write(kMagic, sizeof(kMagic));
  write_pod(out_, kFormatVersion);
  write_pod(out_, params_.k);
  write_pod(out_, params_.scaled);
  write_pod(out_, params_.seed);
  write_pod(out_, max_refs_);
  const uint32_t species_count = static_cast<uint32_t>(directory_.size());
  write_pod(out_, species_count);
  const uint64_t directory_offset = 0; // patched in finish()
  write_pod(out_, directory_offset);
  write_pod(out_, input_genomes_);
  const uint64_t reserved = 0;
  write_pod(out_, reserved);
  if (static_cast<size_t>(out_.tellp()) != kHeaderBytes) {
    throw std::runtime_error("presence sketch header size mismatch");
  }
}

void SketchWriter::add_species(const SpeciesMarkers &markers) {
  if (markers.refs.empty() ||
      (max_refs_ != 0 && markers.refs.size() > max_refs_)) {
    throw std::runtime_error("presence sketch: invalid reference count for species " +
                             std::to_string(markers.species));
  }
  if (markers.ref_begin.size() != markers.refs.size() + 1 ||
      markers.ref_begin.back() != markers.index.size()) {
    throw std::runtime_error("presence sketch: inconsistent reference index for species " +
                             std::to_string(markers.species));
  }
  for (size_t r = 0; r < markers.refs.size(); ++r) {
    if (markers.ref_key_count(r) != markers.refs[r].sketch_size) {
      throw std::runtime_error("presence sketch: sketch size mismatch for species " +
                               std::to_string(markers.species));
    }
  }
  std::lock_guard<std::mutex> lock(mutex_);
  if (finished_) {
    throw std::runtime_error("presence sketch writer already finished");
  }
  DirEntry entry;
  entry.species = markers.species;
  entry.ref_count = static_cast<uint32_t>(markers.refs.size());
  entry.key_count = markers.keys.size();
  entry.index_count = markers.index.size();
  entry.refs_offset = static_cast<uint64_t>(out_.tellp());
  for (const RefInfo &ref : markers.refs) {
    write_pod(out_, ref.bases);
    write_pod(out_, ref.sketch_size);
    const uint32_t name_len = static_cast<uint32_t>(ref.name.size());
    write_pod(out_, name_len);
    out_.write(ref.name.data(), static_cast<std::streamsize>(name_len));
  }
  entry.keys_offset = static_cast<uint64_t>(out_.tellp());
  out_.write(reinterpret_cast<const char *>(markers.keys.data()),
             static_cast<std::streamsize>(markers.keys.size() * sizeof(uint64_t)));
  out_.write(reinterpret_cast<const char *>(markers.index.data()),
             static_cast<std::streamsize>(markers.index.size() * sizeof(uint32_t)));
  if (!out_) {
    throw std::runtime_error("failed while writing presence sketch: " + path_.string());
  }
  keys_written_ += markers.keys.size();
  refs_written_ += markers.refs.size();
  directory_.push_back(entry);
}

void SketchWriter::finish() {
  std::lock_guard<std::mutex> lock(mutex_);
  if (finished_) {
    return;
  }
  finished_ = true;
  std::sort(directory_.begin(), directory_.end(),
            [](const DirEntry &a, const DirEntry &b) { return a.species < b.species; });
  const uint64_t directory_offset = static_cast<uint64_t>(out_.tellp());
  for (const DirEntry &entry : directory_) {
    write_pod(out_, entry.species);
    write_pod(out_, entry.ref_count);
    write_pod(out_, entry.key_count);
    write_pod(out_, entry.index_count);
    write_pod(out_, entry.refs_offset);
    write_pod(out_, entry.keys_offset);
  }
  // patch header: species_count and directory_offset
  out_.seekp(static_cast<std::streamoff>(8 + 4 + 4 + 8 + 8 + 4));
  const uint32_t species_count = static_cast<uint32_t>(directory_.size());
  write_pod(out_, species_count);
  write_pod(out_, directory_offset);
  out_.flush();
  out_.close();
  if (!out_) {
    throw std::runtime_error("failed to finalize presence sketch: " + path_.string());
  }
}

// ----------------------------------------------------------------------------
// Reader

std::optional<SketchIndex> SketchIndex::open(const std::filesystem::path &path,
                                             std::string *error) {
  auto fail = [&](const std::string &message) -> std::optional<SketchIndex> {
    if (error != nullptr) {
      *error = message;
    }
    return std::nullopt;
  };
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return fail("cannot open presence sketch: " + path.string());
  }
  char magic[8];
  in.read(magic, sizeof(magic));
  if (!in || std::memcmp(magic, kMagic, sizeof(kMagic)) != 0) {
    return fail("not a presence sketch: " + path.string());
  }
  SketchIndex index;
  index.path_ = path;
  uint32_t version = 0;
  uint32_t species_count = 0;
  uint64_t directory_offset = 0;
  uint64_t reserved = 0;
  if (!read_pod(in, version) || !read_pod(in, index.params_.k) ||
      !read_pod(in, index.params_.scaled) || !read_pod(in, index.params_.seed) ||
      !read_pod(in, index.max_refs_) || !read_pod(in, species_count) ||
      !read_pod(in, directory_offset) || !read_pod(in, index.input_genomes_) ||
      !read_pod(in, reserved)) {
    return fail("truncated presence sketch header: " + path.string());
  }
  if (version != kFormatVersion) {
    return fail("unsupported presence sketch version " + std::to_string(version));
  }
  if (directory_offset == 0) {
    return fail("presence sketch was not finalized: " + path.string());
  }
  in.seekg(static_cast<std::streamoff>(directory_offset));
  index.directory_.reserve(species_count);
  for (uint32_t i = 0; i < species_count; ++i) {
    uint32_t species = 0;
    DirEntry entry;
    if (!read_pod(in, species) || !read_pod(in, entry.ref_count) ||
        !read_pod(in, entry.key_count) || !read_pod(in, entry.index_count) ||
        !read_pod(in, entry.refs_offset) || !read_pod(in, entry.keys_offset)) {
      return fail("truncated presence sketch directory: " + path.string());
    }
    index.directory_.emplace(species, entry);
  }
  return index;
}

bool SketchIndex::load_species(uint32_t species, SpeciesMarkers &markers) const {
  const auto it = directory_.find(species);
  if (it == directory_.end()) {
    return false;
  }
  const DirEntry &entry = it->second;
  std::ifstream in(path_, std::ios::binary);
  if (!in) {
    throw std::runtime_error("cannot reopen presence sketch: " + path_.string());
  }
  markers.species = species;
  markers.refs.clear();
  markers.refs.reserve(entry.ref_count);
  markers.ref_begin.assign(1, 0);
  markers.ref_begin.reserve(entry.ref_count + 1);
  in.seekg(static_cast<std::streamoff>(entry.refs_offset));
  for (uint32_t i = 0; i < entry.ref_count; ++i) {
    RefInfo ref;
    uint32_t name_len = 0;
    if (!read_pod(in, ref.bases) || !read_pod(in, ref.sketch_size) ||
        !read_pod(in, name_len) || name_len > (1u << 16)) {
      throw std::runtime_error("corrupt presence sketch reference record");
    }
    ref.name.assign(name_len, '\0');
    in.read(ref.name.data(), static_cast<std::streamsize>(name_len));
    markers.ref_begin.push_back(markers.ref_begin.back() + ref.sketch_size);
    markers.refs.push_back(std::move(ref));
  }
  if (markers.ref_begin.back() != entry.index_count) {
    throw std::runtime_error("corrupt presence sketch reference index");
  }
  markers.keys.resize(entry.key_count);
  markers.index.resize(entry.index_count);
  in.seekg(static_cast<std::streamoff>(entry.keys_offset));
  in.read(reinterpret_cast<char *>(markers.keys.data()),
          static_cast<std::streamsize>(entry.key_count * sizeof(uint64_t)));
  in.read(reinterpret_cast<char *>(markers.index.data()),
          static_cast<std::streamsize>(entry.index_count * sizeof(uint32_t)));
  if (!in) {
    throw std::runtime_error("corrupt presence sketch species block");
  }
  return true;
}

std::filesystem::path default_sketch_path_for_db(const std::filesystem::path &db) {
  if (std::filesystem::is_directory(db)) {
    return db / "presence" / "sketch.psk";
  }
  std::filesystem::path core = db;
  if (core.extension() != ".imcf") {
    std::filesystem::path with_extension = core;
    with_extension.replace_extension(".imcf");
    if (std::filesystem::exists(with_extension)) {
      core = with_extension;
    }
  }
  const std::filesystem::path parent =
      core.parent_path().empty() ? std::filesystem::path(".") : core.parent_path();
  return parent / (core.stem().string() + ".presence") / "sketch.psk";
}

} // namespace chimera::presence_sketch
