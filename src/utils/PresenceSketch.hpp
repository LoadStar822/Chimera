/*
 * Genome presence sketch: FracMinHash-sampled canonical 31-mers per species.
 *
 * Build side writes one sidecar file next to the database:
 *   <db>/presence/sketch.psk            (directory databases)
 *   <db stem>.presence/sketch.psk       (legacy single-file databases)
 * Every genome of the build input is kept (Sylph-style): a species block holds
 * the sorted union of its genome sketches plus, for each genome, the indices
 * of its own k-mers in that union.  Classify side samples the reads with the
 * same hash function and compares the observed marker hits of every reported
 * species against what its assigned bases predict (see ChimeraPresenceCall.hpp).
 */
#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace chimera::presence_sketch {

constexpr uint32_t kDefaultK = 31;
constexpr uint64_t kDefaultScaled = 1000;
constexpr uint32_t kDefaultMaxRefs = 0; // 0 = keep every genome
constexpr uint64_t kDefaultSeed = 0x9E3779B97F4A7C15ull;

struct Params {
  uint32_t k{kDefaultK};
  uint64_t scaled{kDefaultScaled};
  uint64_t seed{kDefaultSeed};

  // FracMinHash keeps hashes <= threshold, i.e. a 1/scaled fraction.
  uint64_t threshold() const {
    return scaled <= 1 ? ~0ull : (~0ull) / scaled;
  }
};

// Appends the sampled canonical k-mer hashes of `seq` to `out`.
void sample_hashes(const std::vector<seqan3::dna4> &seq, const Params &params,
                   std::vector<uint64_t> &out);

// Sorts and deduplicates a hash vector in place.
void sort_unique(std::vector<uint64_t> &hashes);

struct RefInfo {
  uint64_t bases{0};        // genome length used as G in coverage estimates
  uint32_t sketch_size{0};  // sampled distinct k-mers of this reference
  std::string name;         // source label (file stem / accession)
};

struct SpeciesMarkers {
  uint32_t species{0};
  std::vector<RefInfo> refs;
  std::vector<uint64_t> keys;       // sorted union of the reference sketches
  // CSR layout: the k-mers of refs[r] are keys[index[j]] for
  // j in [ref_begin[r], ref_begin[r + 1]).
  std::vector<uint64_t> ref_begin;  // refs.size() + 1 entries
  std::vector<uint32_t> index;      // sum of sketch sizes entries

  size_t ref_key_count(size_t r) const {
    return static_cast<size_t>(ref_begin[r + 1] - ref_begin[r]);
  }
};

// Thread-safe appender used by the build.
class SketchWriter {
public:
  SketchWriter(const std::filesystem::path &path, const Params &params,
               uint32_t max_refs, uint64_t input_genomes);
  ~SketchWriter();
  void add_species(const SpeciesMarkers &markers);
  void finish();
  size_t species_written() const { return directory_.size(); }
  uint64_t keys_written() const { return keys_written_; }
  uint64_t refs_written() const { return refs_written_; }

private:
  struct DirEntry {
    uint32_t species{0};
    uint32_t ref_count{0};
    uint64_t key_count{0};
    uint64_t index_count{0};
    uint64_t refs_offset{0};
    uint64_t keys_offset{0};
  };
  void write_header();
  std::filesystem::path path_;
  std::filesystem::path partial_path_;
  Params params_;
  uint32_t max_refs_{0};
  uint64_t input_genomes_{0};
  std::ofstream out_;
  std::mutex mutex_;
  std::vector<DirEntry> directory_;
  uint64_t keys_written_{0};
  uint64_t refs_written_{0};
  bool finished_{false};
};

// Random-access reader used by classify; loads one species at a time.
class SketchIndex {
public:
  static std::optional<SketchIndex> open(const std::filesystem::path &path,
                                         std::string *error = nullptr);
  const Params &params() const { return params_; }
  uint32_t max_refs() const { return max_refs_; }
  size_t species_count() const { return directory_.size(); }
  uint64_t input_genomes() const { return input_genomes_; }
  bool has_species(uint32_t species) const {
    return directory_.count(species) != 0;
  }
  // Returns false when the species is not in the sketch.
  bool load_species(uint32_t species, SpeciesMarkers &markers) const;
  const std::filesystem::path &path() const { return path_; }

private:
  struct DirEntry {
    uint32_t ref_count{0};
    uint64_t key_count{0};
    uint64_t index_count{0};
    uint64_t refs_offset{0};
    uint64_t keys_offset{0};
  };
  std::filesystem::path path_;
  Params params_;
  uint32_t max_refs_{0};
  uint64_t input_genomes_{0};
  std::unordered_map<uint32_t, DirEntry> directory_;
};

// Conventional sidecar location for a database path (directory or .imcf file).
std::filesystem::path default_sketch_path_for_db(const std::filesystem::path &db);

} // namespace chimera::presence_sketch
