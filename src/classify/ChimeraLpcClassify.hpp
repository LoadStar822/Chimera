#pragma once

#include "classifyConfig.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace ChimeraClassify {

// Minimizer keys (2k-bit canonical k-mers, k <= 16) seen in the sample.
// Filled during the main classify pass so local resolution needs no read
// pass of its own and can filter shard anchors with one exact bit test.
class SampleKeyBitset {
public:
  SampleKeyBitset() = default;
  SampleKeyBitset(uint32_t k, uint32_t w);

  bool enabled() const { return !words_.empty(); }
  uint32_t k() const { return k_; }
  uint32_t w() const { return w_; }
  void add(const std::vector<seqan3::dna4> &sequence); // thread-safe
  bool test(uint64_t key) const;
  const std::vector<uint64_t> &words() const { return words_; }
  std::vector<uint64_t> take() { return std::move(words_); }

private:
  uint32_t k_{0};
  uint32_t w_{0};
  std::vector<uint64_t> words_;
};

// Bit per read ordinal; set bits are reads that skip chaining.
struct ReadBitset {
  std::vector<uint64_t> words;

  explicit ReadBitset(uint64_t reads = 0) : words((reads + 63) / 64, 0) {}
  void set(uint64_t ordinal); // thread-safe
  bool test(uint64_t ordinal) const {
    return ordinal / 64 < words.size() &&
           ((words[ordinal / 64] >> (ordinal % 64)) & 1ULL) != 0;
  }
};

struct LocalResolutionCandidate {
  uint32_t taxid{};
  uint32_t score{};
  uint32_t support{}; // representatives chaining nearly as well as the best
};

struct LocalResolutionReadCall {
  uint64_t read_ordinal{};
  std::vector<LocalResolutionCandidate> candidates;
};

struct LocalResolutionReadCallView {
  uint64_t read_ordinal{};
  std::span<const LocalResolutionCandidate> candidates;
};

struct LocalResolutionCallStore {
  std::vector<uint64_t> offsets{0};
  std::vector<LocalResolutionCandidate> candidates;

  uint64_t read_count() const {
    return offsets.empty() ? 0 : static_cast<uint64_t>(offsets.size() - 1);
  }

  bool contains(uint64_t ordinal) const {
    return ordinal + 1 < offsets.size();
  }

  LocalResolutionReadCall materialize(uint64_t ordinal) const {
    LocalResolutionReadCall call;
    call.read_ordinal = ordinal;
    if (!contains(ordinal)) {
      return call;
    }
    const uint64_t begin = offsets[ordinal];
    const uint64_t end = offsets[ordinal + 1];
    if (end > begin) {
      call.candidates.assign(candidates.begin() + static_cast<std::ptrdiff_t>(begin),
                             candidates.begin() + static_cast<std::ptrdiff_t>(end));
    }
    return call;
  }

  LocalResolutionReadCallView view(uint64_t ordinal) const {
    LocalResolutionReadCallView call;
    call.read_ordinal = ordinal;
    if (!contains(ordinal)) {
      return call;
    }
    const uint64_t begin = offsets[ordinal];
    const uint64_t end = offsets[ordinal + 1];
    if (end > begin) {
      const auto begin_offset = static_cast<std::ptrdiff_t>(begin);
      call.candidates = std::span<const LocalResolutionCandidate>(
          candidates.data() + begin_offset,
          static_cast<size_t>(end - begin));
    }
    return call;
  }
};

struct LocalResolutionStats {
  uint64_t reads{};
  uint64_t query_hashes{};
  uint64_t target_filter{};
  uint64_t target_routes{};
  uint64_t core_candidate_reads{};
  uint64_t scanned_shards{};
  uint64_t skipped_shards{};
  uint64_t selected_targets{};
  uint64_t skipped_targets{};
  uint64_t direct_targets{};
  uint64_t target_anchor_records_scanned{};
  uint64_t target_anchor_bytes_read{};
  uint64_t target_anchor_records_matched{};
  uint64_t target_hash_prefilter_rejects{};
  uint64_t direct_load_batches{};
  uint64_t pread_calls{};
  uint64_t pread_bytes{};
  uint64_t raw_chain_records{};
  uint64_t kept_chain_records{};
  uint64_t index_hash_keys{};
  uint64_t overflow_hash_keys{};
  uint64_t dropped_broad_keys{};
  uint64_t dropped_broad_records{};
  uint64_t local_hits{};
  uint64_t local_absent{};
  uint64_t skipped_reads{};
  uint64_t probe_reads{};
  uint64_t probe_chained{};
  uint64_t probe_agree{};
  bool trust_revoked{false};
  uint32_t threads{};
  uint8_t k{};
  uint16_t w{};
  double read_seconds{};
  double ref_seconds{};
  double target_io_seconds{};
  double target_read_seconds{};
  double target_filter_seconds{};
  double target_collect_seconds{};
  double posting_merge_seconds{};
  double index_finalize_seconds{};
  double chain_seconds{};
};

struct LocalResolutionResult {
  LocalResolutionCallStore calls;
  LocalResolutionStats stats;
};

struct LocalResolutionTarget {
  uint32_t genus{};
  uint32_t species{};
  uint32_t target_len{};
  uint32_t anchor_count{};
  uint64_t anchor_byte_offset{};
  uint64_t anchor_byte_size{};
  std::string target_name;
};

// Trusted reads are chained until enough of them have been compared with the
// core call; trust is then revoked for the sample when the two disagree too
// often, otherwise the remaining trusted reads are skipped.
struct TrustProbe {
  std::unordered_map<uint64_t, uint32_t> core_species; // probe ordinal -> species
  uint64_t last_ordinal{};
  uint32_t target_chained{2000};
  uint32_t min_decidable{50};
  double min_agreement{0.7};
};

struct LocalResolutionRequest {
  std::vector<std::string> read_files;
  bool paired{false};
  std::string index_file;
  std::string shard_manifest_file;
  std::vector<LocalResolutionTarget> targets;
  uint32_t diag_bin{};
  uint32_t max_occ{};
  uint32_t min_chain{};
  // A target's chains must span min_coverage of the read or min_coverage_span
  // bases, whichever is less.
  double min_coverage{0.5};
  uint32_t min_coverage_span{300};
  uint32_t threads{};
  const SampleKeyBitset *sample_keys{nullptr}; // replaces the read hash pass
  const ReadBitset *skip_reads{nullptr};       // trusted reads
  const TrustProbe *trust_probe{nullptr};      // when null, skip_reads is final
};

LocalResolutionResult run_local_resolution_engine(
    const LocalResolutionRequest &request);

} // namespace ChimeraClassify
