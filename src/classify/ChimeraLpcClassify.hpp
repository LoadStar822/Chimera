#pragma once

#include "classifyConfig.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace ChimeraClassify {

struct LocalResolutionCandidate {
  std::string taxid;
  uint32_t score{};
};

struct LocalResolutionReadCall {
  std::string read_id;
  std::vector<LocalResolutionCandidate> candidates;
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
  std::vector<LocalResolutionReadCall> reads;
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

struct LocalResolutionRequest {
  std::vector<std::string> read_files;
  std::string index_file;
  std::string shard_manifest_file;
  std::vector<LocalResolutionTarget> targets;
  uint32_t diag_bin{};
  uint32_t max_occ{};
  uint32_t min_chain{};
  uint32_t threads{};
};

LocalResolutionResult run_local_resolution_engine(
    const LocalResolutionRequest &request);

} // namespace ChimeraClassify
