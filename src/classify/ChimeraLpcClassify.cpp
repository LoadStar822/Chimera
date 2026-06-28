#include "ChimeraLpcClassify.hpp"

#include <utils/NativeBoundedIndex.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <dna4_traits.hpp>

#include <fcntl.h>
#include <unistd.h>

namespace {

constexpr size_t kMaxChainTokensPerRead = 256;
constexpr size_t kLocalResolutionReadBatchSize = 4096;

struct PendingRead {
  uint64_t ordinal{};
  uint32_t length{};
  std::vector<seqan3::dna4> sequence;
};

struct ReadRecord {
  uint64_t ordinal{};
  uint32_t anchor_count{};
  std::vector<chimera::native_bounded::Anchor> anchors;
  std::vector<uint32_t> anchor_qids;
  uint32_t length{};
};

struct TargetRecord {
  std::string name;
  uint32_t len{};
  uint32_t species{};
};

struct Posting {
  uint32_t tid_strand{};
  uint32_t pos{};

  Posting() = default;
  Posting(uint32_t tid, uint32_t position, uint8_t strand)
      : tid_strand((tid & kTidMask) |
                   (static_cast<uint32_t>(strand != 0) << 31)),
        pos(position) {}

  uint32_t tid() const { return tid_strand & kTidMask; }
  uint8_t strand() const { return static_cast<uint8_t>(tid_strand >> 31); }

private:
  static constexpr uint32_t kTidMask = 0x7fffffffU;
};
static_assert(sizeof(Posting) == 8, "local posting must stay packed");

struct MatchedPosting {
  uint32_t qid{};
  Posting posting;
};

struct PostingSpan {
  uint32_t offset{};
  uint32_t count{};
};

struct CompactPostingIndex {
  std::vector<Posting> postings;
  std::vector<PostingSpan> spans;
  std::vector<uint8_t> overflow;
};

struct CompactPostingBuild {
  CompactPostingIndex index;
  std::vector<uint32_t> cursor;
};

struct ChainKey {
  uint32_t tid{};
  int32_t diag{};
  uint8_t orient{};

  bool operator==(const ChainKey &other) const {
    return tid == other.tid && diag == other.diag && orient == other.orient;
  }
};

struct ChainKeyHash {
  size_t operator()(const ChainKey &key) const {
    uint64_t x = key.tid;
    x = x * 1315423911u + static_cast<uint32_t>(key.diag);
    x = x * 2654435761u + key.orient;
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    return static_cast<size_t>(x);
  }
};

struct ChainStats {
  uint32_t count{};
  uint32_t q_min{std::numeric_limits<uint32_t>::max()};
  uint32_t q_max{};
  uint32_t t_min{std::numeric_limits<uint32_t>::max()};
  uint32_t t_max{};
};

struct HitOut {
  uint32_t tid{};
  uint32_t score{};
  uint32_t count{};
  uint32_t q_start{};
  uint32_t q_end{};
  uint32_t t_start{};
  uint32_t t_end{};
  uint8_t orient{};
  int32_t diag{};
};

struct TaxonScore {
  uint32_t taxid{};
  uint32_t score{};
};

struct TargetFilter {
  struct Route {
    uint32_t genus{};
    uint32_t species{};
    uint32_t len{};
    uint32_t anchor_count{};
    uint64_t anchor_byte_offset{};
    uint64_t anchor_byte_size{};
    bool direct{};
  };

  std::unordered_set<std::string> names;
  std::unordered_map<std::string, Route> route_by_name;
  std::vector<std::string> ordered_names;
};

struct ShardEntry {
  std::filesystem::path path;
};

struct LoadStats {
  uint64_t selected_targets{};
  uint64_t skipped_targets{};
  uint64_t scanned_shards{};
  uint64_t skipped_shards{};
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
  double target_io_seconds{};
  double target_read_seconds{};
  double target_filter_seconds{};
  double target_collect_seconds{};
  double posting_merge_seconds{};
  double index_finalize_seconds{};
};

struct DirectTargetWork {
  std::filesystem::path path;
  std::string name;
  TargetFilter::Route route;
  uint32_t tid{};
  uint32_t k{};
  uint32_t index_version{};
};

struct DirectTargetLoad {
  TargetRecord target;
  std::vector<MatchedPosting> matches;
  uint64_t raw_chain_records{};
  uint64_t target_anchor_records_scanned{};
  uint64_t target_anchor_bytes_read{};
  uint64_t target_hash_prefilter_rejects{};
  uint64_t pread_calls{};
  uint64_t pread_bytes{};
  double target_read_seconds{};
  double target_filter_seconds{};
  uint32_t tid{};
};

struct DirectTargetPlan {
  std::vector<DirectTargetWork> work_items;
  std::vector<TargetRecord> targets;
  uint64_t scanned_shards{};
  uint64_t skipped_shards{};
};

class DirectShardReader {
public:
  DirectShardReader() = default;
  DirectShardReader(const DirectShardReader &) = delete;
  DirectShardReader &operator=(const DirectShardReader &) = delete;

  ~DirectShardReader() {
    for (const auto &[_, fd] : fds_) {
      if (fd >= 0) {
        ::close(fd);
      }
    }
  }

  size_t read_exact(const std::filesystem::path &path, uint64_t offset,
                    char *data, size_t bytes) {
    int fd = fd_for(path);
    size_t done = 0;
    size_t calls = 0;
    while (done < bytes) {
      ++calls;
      const ssize_t n =
          ::pread(fd, data + done, bytes - done,
                  static_cast<off_t>(offset + done));
      if (n < 0) {
        if (errno == EINTR) {
          continue;
        }
        throw std::runtime_error("failed to read local resolution shard: " +
                                 path.string());
      }
      if (n == 0) {
        throw std::runtime_error("truncated local resolution anchor range: " +
                                 path.string());
      }
      done += static_cast<size_t>(n);
    }
    return calls;
  }

private:
  int fd_for(const std::filesystem::path &path) {
    const std::string key = path.string();
    const auto found = fds_.find(key);
    if (found != fds_.end()) {
      return found->second;
    }
    const int fd = ::open(key.c_str(), O_RDONLY);
    if (fd < 0) {
      throw std::runtime_error("failed to open local resolution shard: " +
                               key);
    }
    fds_.emplace(key, fd);
    return fd;
  }

  std::unordered_map<std::string, int> fds_;
};

constexpr uint64_t kDirectTargetLoadBatchAnchorCap = 64'000'000ULL;
constexpr size_t kDirectReadBufferBytes = 4ULL << 20;

class QueryHashIndex {
public:
  static constexpr uint32_t kNotFound =
      std::numeric_limits<uint32_t>::max();

  uint32_t insert_or_get(uint64_t key) {
    if (keys_.empty() || (dense_keys_.size() + 1) * 2 >= keys_.size()) {
      rehash(keys_.empty() ? 4096 : keys_.size() * 2);
    }
    return insert_no_grow(key);
  }

  uint32_t find_id(uint64_t key) const {
    if (!maybe_contains(key)) {
      return kNotFound;
    }
    return find_id_in_table(key);
  }

  bool maybe_contains(uint64_t key) const {
    if (prefilter_.empty()) {
      return false;
    }
    const size_t bit = static_cast<size_t>(mix(key)) & kPrefilterMask;
    return (prefilter_[bit >> 6] & (uint64_t{1} << (bit & 63))) != 0;
  }

  uint32_t find_id_in_table(uint64_t key) const {
    if (keys_.empty()) {
      return kNotFound;
    }
    const size_t mask = keys_.size() - 1;
    size_t slot = static_cast<size_t>(mix(key)) & mask;
    while (ids_[slot] != kNotFound) {
      if (keys_[slot] == key) {
        return ids_[slot];
      }
      slot = (slot + 1) & mask;
    }
    return kNotFound;
  }

  bool contains(uint64_t key) const { return find_id(key) != kNotFound; }

  size_t size() const { return dense_keys_.size(); }

  uint64_t hash_at(uint32_t id) const { return dense_keys_.at(id); }

  void reserve_slots(size_t requested_slots) {
    if (requested_slots <= keys_.size()) {
      return;
    }
    rehash(requested_slots);
  }

private:
  static uint64_t mix(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
  }

  static size_t round_capacity(size_t requested) {
    size_t capacity = 1;
    while (capacity < requested) {
      capacity <<= 1;
    }
    return std::max<size_t>(capacity, 16);
  }

  void rehash(size_t requested) {
    const size_t new_capacity = round_capacity(requested);
    std::vector<uint64_t> new_keys(new_capacity, 0);
    std::vector<uint32_t> new_ids(new_capacity, kNotFound);
    const size_t mask = new_capacity - 1;
    for (uint32_t id = 0; id < dense_keys_.size(); ++id) {
      const uint64_t key = dense_keys_[id];
      size_t slot = static_cast<size_t>(mix(key)) & mask;
      while (new_ids[slot] != kNotFound) {
        slot = (slot + 1) & mask;
      }
      new_keys[slot] = key;
      new_ids[slot] = id;
    }
    keys_.swap(new_keys);
    ids_.swap(new_ids);
  }

  uint32_t insert_no_grow(uint64_t key) {
    const size_t mask = keys_.size() - 1;
    size_t slot = static_cast<size_t>(mix(key)) & mask;
    while (ids_[slot] != kNotFound) {
      if (keys_[slot] == key) {
        return ids_[slot];
      }
      slot = (slot + 1) & mask;
    }
    const uint32_t id = static_cast<uint32_t>(dense_keys_.size());
    keys_[slot] = key;
    ids_[slot] = id;
    dense_keys_.push_back(key);
    set_prefilter_bit(key);
    return id;
  }

  void insert_existing_no_grow(uint64_t key, uint32_t id) {
    const size_t mask = keys_.size() - 1;
    size_t slot = static_cast<size_t>(mix(key)) & mask;
    while (ids_[slot] != kNotFound) {
      slot = (slot + 1) & mask;
    }
    keys_[slot] = key;
    ids_[slot] = id;
  }

  void set_prefilter_bit(uint64_t key) {
    if (prefilter_.empty()) {
      prefilter_.assign(kPrefilterBitCount / 64, 0);
    }
    const size_t bit = static_cast<size_t>(mix(key)) & kPrefilterMask;
    prefilter_[bit >> 6] |= uint64_t{1} << (bit & 63);
  }

  static constexpr size_t kPrefilterBitCount = size_t{1} << 24;
  static constexpr size_t kPrefilterMask = kPrefilterBitCount - 1;

  std::vector<uint64_t> keys_;
  std::vector<uint32_t> ids_;
  std::vector<uint64_t> dense_keys_;
  std::vector<uint64_t> prefilter_;
};

std::string trim_copy(std::string s) {
  auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
  while (!s.empty() && is_space(static_cast<unsigned char>(s.front()))) {
    s.erase(s.begin());
  }
  while (!s.empty() && is_space(static_cast<unsigned char>(s.back()))) {
    s.pop_back();
  }
  return s;
}

std::string first_token(const std::string &s) {
  size_t n = 0;
  while (n < s.size() && !std::isspace(static_cast<unsigned char>(s[n]))) {
    ++n;
  }
  return s.substr(0, n);
}

uint32_t parse_u32_or_zero(const std::string &s) {
  try {
    size_t consumed = 0;
    const auto value = std::stoul(s, &consumed);
    if (consumed != s.size() ||
        value > std::numeric_limits<uint32_t>::max()) {
      return 0;
    }
    return static_cast<uint32_t>(value);
  } catch (...) {
    return 0;
  }
}

uint64_t parse_u64_or_zero(const std::string &s) {
  try {
    size_t consumed = 0;
    const auto value = std::stoull(s, &consumed);
    if (consumed != s.size()) {
      return 0;
    }
    return value;
  } catch (...) {
    return 0;
  }
}

std::vector<std::string> split_tab_fields(const std::string &line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) {
    fields.push_back(trim_copy(field));
  }
  return fields;
}

ChimeraClassify::LocalResolutionReadCall
make_local_resolution_call(const ReadRecord &read,
                           const std::vector<TaxonScore> &scores) {
  ChimeraClassify::LocalResolutionReadCall call;
  call.read_ordinal = read.ordinal;
  call.candidates.reserve(scores.size());
  for (const auto &score : scores) {
    call.candidates.push_back(
        ChimeraClassify::LocalResolutionCandidate{score.taxid, score.score});
  }
  return call;
}

std::filesystem::path resolve_shard_entry_path(
    const std::filesystem::path &manifest_path, const std::string &raw_path) {
  std::filesystem::path path(raw_path);
  if (path.is_absolute()) {
    return path;
  }

  const auto shard_dir =
      manifest_path.parent_path() / manifest_path.stem().string();
  const std::array<std::filesystem::path, 3> candidates{
      manifest_path.parent_path() / path,
      shard_dir / path,
      shard_dir / path.filename()};
  for (const auto &candidate : candidates) {
    if (std::filesystem::exists(candidate)) {
      return candidate;
    }
  }

  return shard_dir / path.filename();
}

TargetFilter load_target_filter(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot open LPC target list: " + path);
  }
  TargetFilter filter;
  std::string line;
  while (std::getline(in, line)) {
    line = trim_copy(line);
    if (!line.empty()) {
      std::vector<std::string> fields;
      std::stringstream ss(line);
      std::string field;
      while (std::getline(ss, field, '\t')) {
        fields.push_back(trim_copy(field));
      }
      if (fields.size() >= 6) {
        const uint32_t genus = parse_u32_or_zero(fields[0]);
        const uint32_t species = parse_u32_or_zero(fields[1]);
        const uint32_t len = parse_u32_or_zero(fields[2]);
        const uint32_t anchor_count = parse_u32_or_zero(fields[3]);
        const uint64_t anchor_byte_offset = parse_u64_or_zero(fields[4]);
        uint64_t anchor_byte_size = 0;
        std::string name;
        if (fields.size() >= 7) {
          anchor_byte_size = parse_u64_or_zero(fields[5]);
          name = fields[6];
        } else {
          name = fields[5];
        }
        if (genus != 0 && species != 0 && len != 0 && anchor_count != 0 &&
            !name.empty()) {
          if (filter.names.insert(name).second) {
            filter.ordered_names.push_back(name);
          }
          filter.route_by_name[name] = TargetFilter::Route{
              genus, species, len, anchor_count, anchor_byte_offset,
              anchor_byte_size, true};
          continue;
        }
      } else if (fields.size() == 2) {
        const uint32_t genus = parse_u32_or_zero(fields[0]);
        const std::string &name = fields[1];
        if (genus != 0 && !name.empty()) {
          if (filter.names.insert(name).second) {
            filter.ordered_names.push_back(name);
          }
          filter.route_by_name[name] =
              TargetFilter::Route{genus, 0, 0, 0, 0, 0, false};
          continue;
        }
      }
      if (filter.names.insert(line).second) {
        filter.ordered_names.push_back(line);
      }
    }
  }
  if (filter.names.empty()) {
    throw std::runtime_error("LPC target list is empty: " + path);
  }
  return filter;
}

TargetFilter target_filter_from_targets(
    const std::vector<ChimeraClassify::LocalResolutionTarget> &targets) {
  TargetFilter filter;
  for (const auto &target : targets) {
    if (target.genus == 0 || target.species == 0 || target.target_len == 0 ||
        target.anchor_count == 0 || target.target_name.empty()) {
      continue;
    }
    if (filter.names.insert(target.target_name).second) {
      filter.ordered_names.push_back(target.target_name);
    }
    filter.route_by_name[target.target_name] = TargetFilter::Route{
        target.genus, target.species, target.target_len, target.anchor_count,
        target.anchor_byte_offset, target.anchor_byte_size, true};
  }
  if (filter.names.empty()) {
    throw std::runtime_error("Local resolution target panel is empty");
  }
  return filter;
}

std::unordered_map<uint32_t, ShardEntry>
load_shard_manifest(const std::filesystem::path &manifest_path) {
  std::unordered_map<uint32_t, ShardEntry> out;
  std::ifstream in(manifest_path);
  if (!in) {
    return out;
  }
  std::string line;
  if (!std::getline(in, line)) {
    return out;
  }
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string genus;
    std::string path;
    if (!std::getline(ss, genus, '\t') || !std::getline(ss, path, '\t')) {
      continue;
    }
    genus = trim_copy(genus);
    path = trim_copy(path);
    const uint32_t genus_id = parse_u32_or_zero(genus);
    if (genus_id != 0 && !path.empty()) {
      out.emplace(genus_id,
                  ShardEntry{resolve_shard_entry_path(manifest_path, path)});
    }
  }
  return out;
}

std::vector<chimera::native_bounded::Anchor> select_chain_anchors(
    const std::vector<chimera::native_bounded::Anchor> &anchors) {
  const size_t token_budget = std::min(anchors.size(), kMaxChainTokensPerRead);
  std::vector<chimera::native_bounded::Anchor> selected;
  selected.reserve(token_budget);
  for (size_t i = 0; i < token_budget; ++i) {
    const size_t token_index =
        token_budget == anchors.size()
            ? i
            : (i * anchors.size()) / token_budget;
    selected.push_back(anchors[token_index]);
  }
  return selected;
}

void insert_query_hashes(ReadRecord &read,
                         QueryHashIndex &query_hashes) {
  read.anchor_qids.clear();
  read.anchor_qids.reserve(read.anchors.size());
  for (const auto &anchor : read.anchors) {
    read.anchor_qids.push_back(query_hashes.insert_or_get(anchor.hash));
  }
}

void assign_query_hashes(ReadRecord &read, const QueryHashIndex &query_hashes) {
  read.anchor_qids.clear();
  read.anchor_qids.reserve(read.anchors.size());
  for (const auto &anchor : read.anchors) {
    const uint32_t qid = query_hashes.find_id(anchor.hash);
    if (qid == QueryHashIndex::kNotFound) {
      throw std::runtime_error(
          "local resolution query hash missing from first pass");
    }
    read.anchor_qids.push_back(qid);
  }
}

std::vector<ReadRecord>
extract_read_records(const std::vector<PendingRead> &pending_reads, uint32_t k,
                     uint32_t w, uint32_t threads) {
  std::vector<ReadRecord> reads(pending_reads.size());
  if (pending_reads.empty()) {
    return reads;
  }
  const uint32_t worker_count = std::max<uint32_t>(
      1, std::min<uint32_t>(
             threads == 0 ? 1 : threads,
             static_cast<uint32_t>(std::max<size_t>(1, pending_reads.size()))));
  std::atomic<size_t> next_read{0};
  auto worker = [&]() {
    while (true) {
      const size_t idx = next_read.fetch_add(1, std::memory_order_relaxed);
      if (idx >= pending_reads.size()) {
        break;
      }
      const auto &pending = pending_reads[idx];
      ReadRecord read;
      read.ordinal = pending.ordinal;
      read.length = pending.length;
      auto anchors =
          chimera::native_bounded::extract_minimizers(pending.sequence, k, w);
      read.anchor_count = static_cast<uint32_t>(
          std::min<size_t>(anchors.size(), std::numeric_limits<uint32_t>::max()));
      read.anchors = select_chain_anchors(anchors);
      reads[idx] = std::move(read);
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (uint32_t i = 0; i < worker_count; ++i) {
    workers.emplace_back(worker);
  }
  for (auto &worker_thread : workers) {
    worker_thread.join();
  }
  return reads;
}

ReadRecord make_read_record_from_pending(const PendingRead &pending, uint32_t k,
                                         uint32_t w) {
  ReadRecord read;
  read.ordinal = pending.ordinal;
  read.length = pending.length;
  auto anchors =
      chimera::native_bounded::extract_minimizers(pending.sequence, k, w);
  read.anchor_count = static_cast<uint32_t>(
      std::min<size_t>(anchors.size(), std::numeric_limits<uint32_t>::max()));
  read.anchors = select_chain_anchors(anchors);
  return read;
}

template <typename BatchConsumer>
uint64_t for_each_pending_read_batch(const std::vector<std::string> &paths,
                                     size_t batch_size,
                                     BatchConsumer &&consume_batch) {
  std::vector<PendingRead> pending;
  pending.reserve(batch_size);
  uint64_t ordinal = 0;
  auto flush = [&]() {
    if (pending.empty()) {
      return;
    }
    consume_batch(pending);
    std::vector<PendingRead> fresh;
    fresh.reserve(batch_size);
    pending = std::move(fresh);
  };

  for (const auto &path : paths) {
    seqan3::sequence_file_input<raptor::dna4_traits,
                                seqan3::fields<seqan3::field::seq>>
        input{path};
    for (auto &record : input) {
      PendingRead pending_read;
      pending_read.ordinal = ordinal++;
      pending_read.length = static_cast<uint32_t>(
          std::min<size_t>(record.sequence().size(),
                           std::numeric_limits<uint32_t>::max()));
      pending_read.sequence = std::move(record.sequence());
      pending.push_back(std::move(pending_read));
      if (pending.size() >= batch_size) {
        flush();
      }
    }
  }
  flush();
  return ordinal;
}

uint64_t collect_query_hashes_from_reads(
    const std::vector<std::string> &paths, uint32_t k, uint32_t w,
    QueryHashIndex &query_hashes, uint32_t threads) {
  return for_each_pending_read_batch(
      paths, kLocalResolutionReadBatchSize,
      [&](const std::vector<PendingRead> &pending) {
        auto reads = extract_read_records(pending, k, w, threads);
        for (auto &read : reads) {
          insert_query_hashes(read, query_hashes);
        }
      });
}

void append_local_resolution_scores(
    ChimeraClassify::LocalResolutionCallStore &store,
    const std::vector<TaxonScore> &scores) {
  for (const auto &score : scores) {
    store.candidates.push_back(
        ChimeraClassify::LocalResolutionCandidate{score.taxid, score.score});
  }
  store.offsets.push_back(static_cast<uint64_t>(store.candidates.size()));
}

std::vector<ReadRecord> load_reads(
    const std::vector<std::string> &paths, uint32_t k, uint32_t w,
    QueryHashIndex &query_hashes, uint32_t threads) {
  struct PendingRead {
    uint32_t length{};
    std::vector<seqan3::dna4> sequence;
  };

  std::vector<PendingRead> pending_reads;
  for (const auto &path : paths) {
    seqan3::sequence_file_input<
        raptor::dna4_traits,
        seqan3::fields<seqan3::field::seq>>
        input{path};
    for (auto &record : input) {
      PendingRead pending;
      pending.length = static_cast<uint32_t>(
          std::min<size_t>(record.sequence().size(),
                           std::numeric_limits<uint32_t>::max()));
      pending.sequence = std::move(record.sequence());
      pending_reads.push_back(std::move(pending));
    }
  }

  std::vector<ReadRecord> reads(pending_reads.size());
  const uint32_t worker_count = std::max<uint32_t>(
      1, std::min<uint32_t>(
             threads == 0 ? 1 : threads,
             static_cast<uint32_t>(std::max<size_t>(1, pending_reads.size()))));
  std::atomic<size_t> next_read{0};
  auto worker = [&]() {
    while (true) {
      const size_t idx = next_read.fetch_add(1, std::memory_order_relaxed);
      if (idx >= pending_reads.size()) {
        break;
      }
      auto &pending = pending_reads[idx];
      ReadRecord read;
      read.length = pending.length;
      auto anchors = chimera::native_bounded::extract_minimizers(
          pending.sequence, k, w);
      read.anchor_count = static_cast<uint32_t>(
          std::min<size_t>(anchors.size(), std::numeric_limits<uint32_t>::max()));
      read.anchors = select_chain_anchors(anchors);
      reads[idx] = std::move(read);
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (uint32_t i = 0; i < worker_count; ++i) {
    workers.emplace_back(worker);
  }
  for (auto &worker_thread : workers) {
    worker_thread.join();
  }

  std::vector<PendingRead>().swap(pending_reads);

  size_t selected_query_tokens = 0;
  for (const auto &read : reads) {
    selected_query_tokens += read.anchors.size();
  }
  query_hashes.reserve_slots(selected_query_tokens);
  for (auto &read : reads) {
    insert_query_hashes(read, query_hashes);
  }
  return reads;
}

void add_target_records(
    const std::filesystem::path &path,
    const chimera::native_bounded::IndexMeta &meta,
    const TargetFilter &target_filter,
    const QueryHashIndex &query_hashes, std::vector<TargetRecord> &targets,
    std::vector<MatchedPosting> &matches, LoadStats &stats) {
  for (const auto &source : meta.targets) {
    if (!target_filter.names.contains(source.name)) {
      ++stats.skipped_targets;
      continue;
    }
    const uint32_t tid = static_cast<uint32_t>(targets.size());
    targets.push_back(TargetRecord{source.name, source.len, source.species});
    ++stats.selected_targets;
    const auto read_started = std::chrono::steady_clock::now();
    const auto anchors =
        chimera::native_bounded::read_anchor_range(path, meta, source);
    const auto read_finished = std::chrono::steady_clock::now();
    stats.target_read_seconds +=
        std::chrono::duration<double>(read_finished - read_started).count();
    stats.target_anchor_records_scanned += anchors.size();
    stats.target_anchor_bytes_read +=
        static_cast<uint64_t>(anchors.size()) *
        chimera::native_bounded::anchor_record_bytes(meta.k);
    const auto filter_started = std::chrono::steady_clock::now();
    for (const auto &anchor : anchors) {
      if (!query_hashes.maybe_contains(anchor.hash)) {
        ++stats.target_hash_prefilter_rejects;
        continue;
      }
      const uint32_t qid = query_hashes.find_id_in_table(anchor.hash);
      if (qid == QueryHashIndex::kNotFound) {
        continue;
      }
      ++stats.raw_chain_records;
      ++stats.target_anchor_records_matched;
      matches.push_back(MatchedPosting{qid, Posting{tid, anchor.pos,
                                                    anchor.strand}});
    }
    const auto filter_finished = std::chrono::steady_clock::now();
    stats.target_filter_seconds +=
        std::chrono::duration<double>(filter_finished - filter_started).count();
  }
}

bool read_local_varuint(const char *&cursor, const char *end,
                        uint64_t &value) {
  value = 0;
  uint32_t shift = 0;
  while (cursor < end && shift <= 63) {
    const auto byte = static_cast<uint8_t>(*cursor++);
    value |= static_cast<uint64_t>(byte & 0x7fU) << shift;
    if ((byte & 0x80U) == 0) {
      return true;
    }
    shift += 7;
  }
  return false;
}

uint64_t read_local_bits(const char *data, uint64_t bit_offset,
                         uint32_t bits) {
  uint64_t value = 0;
  uint32_t written = 0;
  uint32_t remaining = bits;
  while (remaining > 0) {
    const uint64_t byte_index = bit_offset / 8ULL;
    const uint32_t shift = static_cast<uint32_t>(bit_offset % 8ULL);
    const uint32_t take = std::min<uint32_t>(remaining, 8U - shift);
    const uint64_t mask = (1ULL << take) - 1ULL;
    const auto byte = static_cast<uint8_t>(data[byte_index]);
    value |= ((static_cast<uint64_t>(byte >> shift) & mask) << written);
    bit_offset += take;
    written += take;
    remaining -= take;
  }
  return value;
}

uint64_t local_bitpacked_key_bytes(uint32_t anchor_count, uint32_t k) {
  const uint32_t key_bits = static_cast<uint32_t>(2U * k);
  return (static_cast<uint64_t>(anchor_count) * key_bits + 7ULL) / 8ULL;
}

template <typename AnchorConsumer>
void for_each_encoded_anchor(const char *data, size_t size,
                             uint32_t anchor_count, uint32_t k,
                             uint32_t version, bool fixed_record_format,
                             AnchorConsumer &&consume_anchor) {
  const char *cursor = data;
  const char *end = data + size;
  if (fixed_record_format) {
    for (uint32_t i = 0; i < anchor_count; ++i) {
      chimera::native_bounded::Anchor anchor;
      if (k <= 16) {
        if (static_cast<size_t>(end - cursor) < sizeof(uint32_t)) {
          throw std::runtime_error("truncated native bounded anchor key");
        }
        uint32_t key = 0;
        std::memcpy(&key, cursor, sizeof(key));
        cursor += sizeof(key);
        anchor.hash = key;
      } else {
        if (static_cast<size_t>(end - cursor) < sizeof(anchor.hash)) {
          throw std::runtime_error("truncated native bounded anchor key");
        }
        std::memcpy(&anchor.hash, cursor, sizeof(anchor.hash));
        cursor += sizeof(anchor.hash);
      }
      if (static_cast<size_t>(end - cursor) <
          sizeof(anchor.pos) + sizeof(anchor.strand)) {
        throw std::runtime_error("truncated native bounded anchor record");
      }
      std::memcpy(&anchor.pos, cursor, sizeof(anchor.pos));
      cursor += sizeof(anchor.pos);
      std::memcpy(&anchor.strand, cursor, sizeof(anchor.strand));
      cursor += sizeof(anchor.strand);
      consume_anchor(anchor);
    }
    if (cursor != end) {
      throw std::runtime_error("native bounded anchor block has trailing bytes");
    }
    return;
  }

  const uint32_t key_bits = static_cast<uint32_t>(2U * k);
  const uint64_t key_bytes =
      version >= 4 && k <= 16 ? local_bitpacked_key_bytes(anchor_count, k) : 0;
  if (key_bytes > size) {
    throw std::runtime_error("truncated native bounded anchor keys");
  }
  cursor = version >= 4 && k <= 16 ? data + key_bytes : data;
  uint32_t previous_pos = 0;
  for (uint32_t i = 0; i < anchor_count; ++i) {
    chimera::native_bounded::Anchor anchor;
    if (version >= 4 && k <= 16) {
      anchor.hash = read_local_bits(data, static_cast<uint64_t>(i) * key_bits,
                                    key_bits);
    } else if (k <= 16) {
      if (static_cast<size_t>(end - cursor) < sizeof(uint32_t)) {
        throw std::runtime_error("truncated native bounded anchor key");
      }
      uint32_t key = 0;
      std::memcpy(&key, cursor, sizeof(key));
      cursor += sizeof(key);
      anchor.hash = key;
    } else {
      if (static_cast<size_t>(end - cursor) < sizeof(anchor.hash)) {
        throw std::runtime_error("truncated native bounded anchor key");
      }
      std::memcpy(&anchor.hash, cursor, sizeof(anchor.hash));
      cursor += sizeof(anchor.hash);
    }

    uint64_t packed = 0;
    if (!read_local_varuint(cursor, end, packed)) {
      throw std::runtime_error("truncated native bounded anchor position");
    }
    const uint64_t delta = packed >> 1;
    if (delta > std::numeric_limits<uint32_t>::max()) {
      throw std::runtime_error("native bounded anchor position delta is too large");
    }
    if (i == 0) {
      anchor.pos = static_cast<uint32_t>(delta);
    } else {
      if (delta >
          static_cast<uint64_t>(std::numeric_limits<uint32_t>::max() -
                                previous_pos)) {
        throw std::runtime_error("native bounded anchor position overflow");
      }
      anchor.pos = previous_pos + static_cast<uint32_t>(delta);
    }
    anchor.strand = static_cast<uint8_t>(packed & 1U);
    previous_pos = anchor.pos;
    consume_anchor(anchor);
  }
  if (cursor != end) {
    throw std::runtime_error("native bounded anchor block has trailing bytes");
  }
}

template <typename MatchConsumer>
DirectTargetLoad scan_direct_target(
    DirectShardReader &reader, const DirectTargetWork &work,
    const QueryHashIndex &query_hashes, MatchConsumer &&consume_match) {
  DirectTargetLoad out;
  out.target = TargetRecord{work.name, work.route.len, work.route.species};
  const size_t anchor_record_bytes =
      chimera::native_bounded::anchor_record_bytes(work.k);
  const uint64_t target_byte_size =
      work.route.anchor_byte_size == 0
          ? static_cast<uint64_t>(work.route.anchor_count) *
                anchor_record_bytes
          : work.route.anchor_byte_size;
  thread_local std::vector<char> buffer;
  buffer.resize(static_cast<size_t>(target_byte_size));
  uint64_t byte_offset = work.route.anchor_byte_offset;
  const auto read_started = std::chrono::steady_clock::now();
  if (target_byte_size > 0) {
    out.pread_calls += reader.read_exact(work.path, byte_offset, buffer.data(),
                                         target_byte_size);
  }
  const auto read_finished = std::chrono::steady_clock::now();
  out.target_read_seconds +=
      std::chrono::duration<double>(read_finished - read_started).count();
  out.pread_bytes += target_byte_size;
  out.target_anchor_records_scanned += work.route.anchor_count;
  out.target_anchor_bytes_read += target_byte_size;
  const auto filter_started = std::chrono::steady_clock::now();
  for_each_encoded_anchor(
      buffer.data(), static_cast<size_t>(target_byte_size),
      work.route.anchor_count, work.k, work.index_version,
      work.route.anchor_byte_size == 0, [&](const auto &anchor) {
    if (!query_hashes.maybe_contains(anchor.hash)) {
      ++out.target_hash_prefilter_rejects;
      return;
    }
    const uint32_t qid = query_hashes.find_id_in_table(anchor.hash);
    if (qid == QueryHashIndex::kNotFound) {
      return;
    }
    ++out.raw_chain_records;
    consume_match(out, qid, Posting{work.tid, anchor.pos, anchor.strand});
  });
  const auto filter_finished = std::chrono::steady_clock::now();
  out.target_filter_seconds +=
      std::chrono::duration<double>(filter_finished - filter_started).count();
  out.tid = work.tid;
  return out;
}

DirectTargetLoad load_direct_target(
    DirectShardReader &reader,
    const DirectTargetWork &work,
    const QueryHashIndex &query_hashes) {
  return scan_direct_target(
      reader, work, query_hashes,
      [](DirectTargetLoad &out, uint32_t qid, const Posting &posting) {
        out.matches.push_back(MatchedPosting{qid, posting});
      });
}

std::vector<DirectTargetLoad> load_direct_targets_parallel(
    const std::vector<DirectTargetWork> &work_items, size_t begin, size_t end,
    const QueryHashIndex &query_hashes, uint32_t threads) {
  if (begin > end || end > work_items.size()) {
    throw std::runtime_error("invalid direct LPC target batch range");
  }
  const size_t batch_size = end - begin;
  std::vector<DirectTargetLoad> loaded(batch_size);
  if (batch_size == 0) {
    return loaded;
  }
  const uint32_t worker_count = std::max<uint32_t>(
      1, std::min<uint32_t>(threads == 0 ? 1 : threads,
                            static_cast<uint32_t>(batch_size)));
  std::atomic<size_t> next{0};
  std::atomic<bool> stop{false};
  std::exception_ptr error;
  std::mutex error_mutex;
  auto worker = [&]() {
    DirectShardReader reader;
    while (!stop.load(std::memory_order_relaxed)) {
      const size_t local_idx = next.fetch_add(1, std::memory_order_relaxed);
      if (local_idx >= batch_size) {
        break;
      }
      const size_t idx = begin + local_idx;
      try {
        loaded[local_idx] =
            load_direct_target(reader, work_items[idx], query_hashes);
      } catch (...) {
        std::lock_guard<std::mutex> lock(error_mutex);
        if (error == nullptr) {
          error = std::current_exception();
          stop.store(true, std::memory_order_relaxed);
        }
        break;
      }
    }
  };
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (uint32_t i = 0; i < worker_count; ++i) {
    workers.emplace_back(worker);
  }
  for (auto &thread : workers) {
    thread.join();
  }
  if (error != nullptr) {
    std::rethrow_exception(error);
  }
  return loaded;
}

std::vector<DirectTargetLoad> count_direct_targets_parallel(
    const std::vector<DirectTargetWork> &work_items,
    const QueryHashIndex &query_hashes,
    std::vector<std::atomic<uint32_t>> &counts, uint32_t threads) {
  std::vector<DirectTargetLoad> loaded(work_items.size());
  if (work_items.empty()) {
    return loaded;
  }
  const uint32_t worker_count = std::max<uint32_t>(
      1, std::min<uint32_t>(threads == 0 ? 1 : threads,
                            static_cast<uint32_t>(work_items.size())));
  std::atomic<size_t> next{0};
  std::atomic<bool> stop{false};
  std::exception_ptr error;
  std::mutex error_mutex;
  auto worker = [&]() {
    DirectShardReader reader;
    while (!stop.load(std::memory_order_relaxed)) {
      const size_t idx = next.fetch_add(1, std::memory_order_relaxed);
      if (idx >= work_items.size()) {
        break;
      }
      try {
        loaded[idx] = scan_direct_target(
            reader, work_items[idx], query_hashes,
            [&](DirectTargetLoad &, uint32_t qid, const Posting &) {
              counts[qid].fetch_add(1, std::memory_order_relaxed);
            });
      } catch (...) {
        std::lock_guard<std::mutex> lock(error_mutex);
        if (error == nullptr) {
          error = std::current_exception();
          stop.store(true, std::memory_order_relaxed);
        }
        break;
      }
    }
  };
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (uint32_t i = 0; i < worker_count; ++i) {
    workers.emplace_back(worker);
  }
  for (auto &thread : workers) {
    thread.join();
  }
  if (error != nullptr) {
    std::rethrow_exception(error);
  }
  return loaded;
}

void append_loaded_matches_preserving_order(
    std::vector<MatchedPosting> &matches,
    std::vector<DirectTargetLoad> &loaded,
    const std::vector<size_t> &match_offsets,
    size_t batch_matches,
    uint32_t threads) {
  if (batch_matches == 0) {
    return;
  }
  const size_t base = matches.size();
  matches.resize(base + batch_matches);
  const uint32_t worker_count = std::max<uint32_t>(
      1, std::min<uint32_t>(threads == 0 ? 1 : threads,
                            static_cast<uint32_t>(loaded.size())));
  std::atomic<size_t> next{0};
  auto worker = [&]() {
    while (true) {
      const size_t idx = next.fetch_add(1, std::memory_order_relaxed);
      if (idx >= loaded.size()) {
        break;
      }
      auto &source = loaded[idx].matches;
      if (source.empty()) {
        continue;
      }
      std::memcpy(matches.data() + base + match_offsets[idx], source.data(),
                  source.size() * sizeof(MatchedPosting));
    }
  };
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (uint32_t i = 0; i < worker_count; ++i) {
    workers.emplace_back(worker);
  }
  for (auto &thread : workers) {
    thread.join();
  }
}

size_t direct_target_batch_end(
    const std::vector<DirectTargetWork> &work_items, size_t begin) {
  uint64_t anchors = 0;
  size_t end = begin;
  while (end < work_items.size()) {
    const uint64_t target_anchors = work_items[end].route.anchor_count;
    if (end > begin &&
        anchors + target_anchors > kDirectTargetLoadBatchAnchorCap) {
      break;
    }
    anchors += target_anchors;
    ++end;
  }
  return end;
}

bool build_direct_target_plan(
    const TargetFilter &target_filter,
    const chimera::native_bounded::IndexMeta &root_meta,
    const std::unordered_map<uint32_t, ShardEntry> &shards,
    DirectTargetPlan &plan) {
  if (shards.empty() ||
      target_filter.route_by_name.size() != target_filter.names.size()) {
    return false;
  }
  for (const auto &[_, route] : target_filter.route_by_name) {
    if (!route.direct) {
      return false;
    }
  }

  struct RoutedDirectTarget {
    std::string name;
    TargetFilter::Route route;
    uint32_t tid{};
  };
  std::unordered_map<uint32_t, std::vector<RoutedDirectTarget>> by_genus;
  by_genus.reserve(target_filter.route_by_name.size());
  uint32_t tid = 0;
  for (const auto &name : target_filter.ordered_names) {
    const auto route_found = target_filter.route_by_name.find(name);
    if (route_found == target_filter.route_by_name.end()) {
      throw std::runtime_error("LPC target route missing ordered target");
    }
    by_genus[route_found->second.genus].push_back(
        {name, route_found->second, tid});
    ++tid;
  }

  std::vector<uint32_t> genera;
  genera.reserve(by_genus.size());
  for (const auto &[genus, _] : by_genus) {
    genera.push_back(genus);
  }
  std::sort(genera.begin(), genera.end());
  plan.skipped_shards =
      shards.size() > genera.size() ? shards.size() - genera.size() : 0;
  plan.targets.resize(target_filter.route_by_name.size());
  plan.work_items.reserve(target_filter.route_by_name.size());
  for (const uint32_t genus : genera) {
    const auto shard_found = shards.find(genus);
    if (shard_found == shards.end()) {
      throw std::runtime_error("LPC target route references missing shard");
    }
    ++plan.scanned_shards;
    auto &group = by_genus[genus];
    std::sort(group.begin(), group.end(),
              [](const RoutedDirectTarget &lhs,
                 const RoutedDirectTarget &rhs) {
      if (lhs.route.anchor_byte_offset != rhs.route.anchor_byte_offset) {
        return lhs.route.anchor_byte_offset < rhs.route.anchor_byte_offset;
      }
      return lhs.name < rhs.name;
    });
    for (const auto &entry : group) {
      plan.targets[entry.tid] =
          TargetRecord{entry.name, entry.route.len, entry.route.species};
      plan.work_items.push_back(DirectTargetWork{
          shard_found->second.path, entry.name, entry.route, entry.tid,
          root_meta.k, root_meta.version});
    }
  }
  if (plan.work_items.empty()) {
    throw std::runtime_error("LPC selected no target references");
  }
  for (const auto &target : plan.targets) {
    if (target.name.empty()) {
      throw std::runtime_error("direct LPC target order invariant violated");
    }
  }
  return true;
}

void add_direct_physical_stats(const DirectTargetLoad &loaded,
                               LoadStats &stats) {
  stats.target_anchor_records_scanned +=
      loaded.target_anchor_records_scanned;
  stats.target_anchor_bytes_read += loaded.target_anchor_bytes_read;
  stats.target_hash_prefilter_rejects +=
      loaded.target_hash_prefilter_rejects;
  stats.pread_calls += loaded.pread_calls;
  stats.pread_bytes += loaded.pread_bytes;
  stats.target_read_seconds += loaded.target_read_seconds;
  stats.target_filter_seconds += loaded.target_filter_seconds;
}

void add_direct_logical_match_stats(const DirectTargetLoad &loaded,
                                    LoadStats &stats) {
  stats.target_anchor_records_matched += loaded.raw_chain_records;
  stats.raw_chain_records += loaded.raw_chain_records;
}

std::vector<TargetRecord> load_lpc_targets(
    const std::filesystem::path &index_path,
    const std::filesystem::path &shard_manifest_path,
    const TargetFilter &target_filter,
    const QueryHashIndex &query_hashes, std::vector<MatchedPosting> &matches,
    LoadStats &stats,
    uint32_t threads) {
  std::vector<TargetRecord> targets;
  const auto root_meta = chimera::native_bounded::read_index_header(index_path);
  const auto shards = load_shard_manifest(shard_manifest_path);
  if (!shards.empty()) {
    const bool complete_routes =
        target_filter.route_by_name.size() == target_filter.names.size();
    bool direct_routes = complete_routes;
    if (direct_routes) {
      for (const auto &[_, route] : target_filter.route_by_name) {
        if (!route.direct) {
          direct_routes = false;
          break;
        }
      }
    }
    if (direct_routes) {
      struct RoutedDirectTarget {
        std::string name;
        TargetFilter::Route route;
        uint32_t tid{};
      };
      std::unordered_map<uint32_t,
                         std::vector<RoutedDirectTarget>>
          by_genus;
      by_genus.reserve(target_filter.route_by_name.size());
      uint32_t tid = 0;
      for (const auto &name : target_filter.ordered_names) {
        const auto route_found = target_filter.route_by_name.find(name);
        if (route_found == target_filter.route_by_name.end()) {
          throw std::runtime_error("LPC target route missing ordered target");
        }
        by_genus[route_found->second.genus].push_back(
            {name, route_found->second, tid});
        ++tid;
      }
      std::vector<uint32_t> genera;
      genera.reserve(by_genus.size());
      for (const auto &[genus, _] : by_genus) {
        genera.push_back(genus);
      }
      std::sort(genera.begin(), genera.end());
      stats.skipped_shards =
          shards.size() > genera.size() ? shards.size() - genera.size() : 0;
      std::vector<DirectTargetWork> work_items;
      work_items.reserve(target_filter.route_by_name.size());
      for (const uint32_t genus : genera) {
        const auto shard_found = shards.find(genus);
        if (shard_found == shards.end()) {
          throw std::runtime_error("LPC target route references missing shard");
        }
        ++stats.scanned_shards;
        auto &group = by_genus[genus];
        std::sort(group.begin(), group.end(),
                  [](const RoutedDirectTarget &lhs,
                     const RoutedDirectTarget &rhs) {
          if (lhs.route.anchor_byte_offset != rhs.route.anchor_byte_offset) {
            return lhs.route.anchor_byte_offset < rhs.route.anchor_byte_offset;
          }
          return lhs.name < rhs.name;
        });
        for (const auto &entry : group) {
          work_items.push_back(DirectTargetWork{
              shard_found->second.path, entry.name, entry.route, entry.tid,
              root_meta.k});
        }
      }
      if (work_items.empty()) {
        throw std::runtime_error("LPC selected no target references");
      }
      targets.resize(work_items.size());
      for (size_t begin = 0; begin < work_items.size();) {
        const size_t end = direct_target_batch_end(work_items, begin);
        const auto load_started = std::chrono::steady_clock::now();
        auto loaded = load_direct_targets_parallel(work_items, begin, end,
                                                   query_hashes, threads);
        const auto load_finished = std::chrono::steady_clock::now();
        stats.target_io_seconds +=
            std::chrono::duration<double>(load_finished - load_started)
                .count();
        const auto collect_started = std::chrono::steady_clock::now();
        size_t batch_matches = 0;
        std::vector<size_t> match_offsets(loaded.size() + 1, 0);
        for (size_t loaded_idx = 0; loaded_idx < loaded.size(); ++loaded_idx) {
          auto &loaded_target = loaded[loaded_idx];
          if (loaded_target.tid >= targets.size()) {
            throw std::runtime_error(
                "direct LPC target order invariant violated");
          }
          targets[loaded_target.tid] = std::move(loaded_target.target);
          ++stats.selected_targets;
          ++stats.direct_targets;
          stats.target_anchor_records_scanned +=
              loaded_target.target_anchor_records_scanned;
          stats.target_anchor_bytes_read +=
              loaded_target.target_anchor_bytes_read;
          stats.target_anchor_records_matched +=
              loaded_target.raw_chain_records;
          stats.target_hash_prefilter_rejects +=
              loaded_target.target_hash_prefilter_rejects;
          stats.pread_calls += loaded_target.pread_calls;
          stats.pread_bytes += loaded_target.pread_bytes;
          stats.target_read_seconds += loaded_target.target_read_seconds;
          stats.target_filter_seconds += loaded_target.target_filter_seconds;
          stats.raw_chain_records += loaded_target.raw_chain_records;
          batch_matches += loaded_target.matches.size();
          match_offsets[loaded_idx + 1] = batch_matches;
        }
        append_loaded_matches_preserving_order(matches, loaded, match_offsets,
                                               batch_matches, threads);
        const auto collect_finished = std::chrono::steady_clock::now();
        stats.target_collect_seconds +=
            std::chrono::duration<double>(collect_finished - collect_started)
                .count();
        begin = end;
        ++stats.direct_load_batches;
      }
      for (const auto &target : targets) {
        if (target.name.empty()) {
          throw std::runtime_error("direct LPC target order invariant violated");
        }
      }
      return targets;
    }

    std::vector<uint32_t> genera;
    genera.reserve(shards.size());
    if (complete_routes) {
      for (const auto &[_, route] : target_filter.route_by_name) {
        genera.push_back(route.genus);
      }
      std::sort(genera.begin(), genera.end());
      genera.erase(std::unique(genera.begin(), genera.end()), genera.end());
      stats.skipped_shards =
          shards.size() > genera.size() ? shards.size() - genera.size() : 0;
    } else {
      for (const auto &[genus, _] : shards) {
        genera.push_back(genus);
      }
      std::sort(genera.begin(), genera.end());
    }
    for (const uint32_t genus : genera) {
      const auto shard_found = shards.find(genus);
      if (shard_found == shards.end()) {
        continue;
      }
      ++stats.scanned_shards;
      const auto meta =
          chimera::native_bounded::read_index_meta(shard_found->second.path);
      add_target_records(shard_found->second.path, meta, target_filter,
                         query_hashes, targets, matches, stats);
    }
  } else {
    ++stats.scanned_shards;
    const auto meta = chimera::native_bounded::read_index_meta(index_path);
    add_target_records(index_path, meta, target_filter, query_hashes, targets,
                       matches, stats);
  }
  if (targets.empty()) {
    throw std::runtime_error("LPC selected no target references");
  }
  return targets;
}

CompactPostingIndex build_compact_posting_index(
    const std::vector<MatchedPosting> &matches, size_t query_hash_count,
    uint32_t max_occ, LoadStats &stats) {
  if (query_hash_count >
      static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
    throw std::runtime_error("too many local query hashes for compact index");
  }
  CompactPostingIndex out;
  out.spans.assign(query_hash_count, PostingSpan{});
  out.overflow.assign(query_hash_count, 0);
  std::vector<uint32_t> counts(query_hash_count, 0);
  for (const auto &match : matches) {
    ++counts[match.qid];
  }

  uint64_t kept = 0;
  for (uint32_t qid = 0; qid < counts.size(); ++qid) {
    const uint32_t count = counts[qid];
    if (count == 0) {
      continue;
    }
    if (count > max_occ) {
      out.overflow[qid] = 1;
      ++stats.dropped_broad_keys;
      ++stats.overflow_hash_keys;
      stats.dropped_broad_records += count;
      continue;
    }
    if (kept + count >
        static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
      throw std::runtime_error(
          "too many local postings for compact posting index");
    }
    out.spans[qid] = PostingSpan{static_cast<uint32_t>(kept), count};
    kept += count;
    ++stats.index_hash_keys;
  }

  out.postings.resize(static_cast<size_t>(kept));
  for (uint32_t qid = 0; qid < counts.size(); ++qid) {
    counts[qid] = out.spans[qid].offset;
  }
  for (const auto &match : matches) {
    if (out.overflow[match.qid] != 0) {
      continue;
    }
    out.postings[counts[match.qid]++] = match.posting;
  }
  stats.raw_chain_records = matches.size();
  stats.kept_chain_records = out.postings.size();
  return out;
}

CompactPostingBuild build_compact_posting_index_from_counts(
    const std::vector<std::atomic<uint32_t>> &counts, uint32_t max_occ,
    LoadStats &stats) {
  CompactPostingBuild build;
  const size_t query_hash_count = counts.size();
  if (query_hash_count >
      static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
    throw std::runtime_error("too many local query hashes for compact index");
  }
  build.index.spans.assign(query_hash_count, PostingSpan{});
  build.index.overflow.assign(query_hash_count, 0);
  build.cursor.assign(query_hash_count, 0);

  uint64_t kept = 0;
  for (uint32_t qid = 0; qid < query_hash_count; ++qid) {
    const uint32_t count = counts[qid].load(std::memory_order_relaxed);
    if (count == 0) {
      continue;
    }
    if (count > max_occ) {
      build.index.overflow[qid] = 1;
      ++stats.dropped_broad_keys;
      ++stats.overflow_hash_keys;
      stats.dropped_broad_records += count;
      continue;
    }
    if (kept + count >
        static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
      throw std::runtime_error(
          "too many local postings for compact posting index");
    }
    const uint32_t offset = static_cast<uint32_t>(kept);
    build.index.spans[qid] = PostingSpan{offset, count};
    build.cursor[qid] = offset;
    kept += count;
    ++stats.index_hash_keys;
  }

  build.index.postings.resize(static_cast<size_t>(kept));
  stats.kept_chain_records = build.index.postings.size();
  return build;
}

void fill_compact_postings_from_loaded(
    CompactPostingIndex &index, std::vector<uint32_t> &cursor,
    const std::vector<DirectTargetLoad> &loaded) {
  for (const auto &loaded_target : loaded) {
    for (const auto &match : loaded_target.matches) {
      if (index.overflow[match.qid] != 0) {
        continue;
      }
      index.postings[cursor[match.qid]++] = match.posting;
    }
  }
}

bool try_load_direct_lpc_targets_compact(
    const std::filesystem::path &index_path,
    const std::filesystem::path &shard_manifest_path,
    const TargetFilter &target_filter, const QueryHashIndex &query_hashes,
    uint32_t max_occ, uint32_t threads, std::vector<TargetRecord> &targets,
    CompactPostingIndex &index, LoadStats &stats) {
  const auto root_meta = chimera::native_bounded::read_index_header(index_path);
  const auto shards = load_shard_manifest(shard_manifest_path);
  DirectTargetPlan plan;
  if (!build_direct_target_plan(target_filter, root_meta, shards, plan)) {
    return false;
  }

  targets = std::move(plan.targets);
  stats.scanned_shards += plan.scanned_shards;
  stats.skipped_shards += plan.skipped_shards;
  stats.selected_targets += plan.work_items.size();
  stats.direct_targets += plan.work_items.size();

  std::vector<std::atomic<uint32_t>> counts(query_hashes.size());
  for (auto &count : counts) {
    count.store(0, std::memory_order_relaxed);
  }

  const auto count_started = std::chrono::steady_clock::now();
  auto counted =
      count_direct_targets_parallel(plan.work_items, query_hashes, counts,
                                    threads);
  const auto count_finished = std::chrono::steady_clock::now();
  stats.target_io_seconds +=
      std::chrono::duration<double>(count_finished - count_started).count();
  ++stats.direct_load_batches;
  for (const auto &loaded_target : counted) {
    add_direct_physical_stats(loaded_target, stats);
    add_direct_logical_match_stats(loaded_target, stats);
  }

  const auto index_started = std::chrono::steady_clock::now();
  auto build = build_compact_posting_index_from_counts(counts, max_occ, stats);
  std::vector<std::atomic<uint32_t>>().swap(counts);
  const auto index_built = std::chrono::steady_clock::now();
  stats.posting_merge_seconds +=
      std::chrono::duration<double>(index_built - index_started).count();

  for (size_t begin = 0; begin < plan.work_items.size();) {
    const size_t end = direct_target_batch_end(plan.work_items, begin);
    const auto load_started = std::chrono::steady_clock::now();
    auto loaded = load_direct_targets_parallel(plan.work_items, begin, end,
                                               query_hashes, threads);
    const auto load_finished = std::chrono::steady_clock::now();
    stats.target_io_seconds +=
        std::chrono::duration<double>(load_finished - load_started).count();
    ++stats.direct_load_batches;
    for (const auto &loaded_target : loaded) {
      add_direct_physical_stats(loaded_target, stats);
    }

    const auto fill_started = std::chrono::steady_clock::now();
    fill_compact_postings_from_loaded(build.index, build.cursor, loaded);
    const auto fill_finished = std::chrono::steady_clock::now();
    stats.target_collect_seconds +=
        std::chrono::duration<double>(fill_finished - fill_started).count();
    begin = end;
  }
  std::vector<uint32_t>().swap(build.cursor);
  index = std::move(build.index);
  return true;
}

void sort_chain_hits(std::vector<HitOut> &hits) {
  std::sort(hits.begin(), hits.end(), [](const HitOut &a, const HitOut &b) {
    if (a.score != b.score) {
      return a.score > b.score;
    }
    if (a.count != b.count) {
      return a.count > b.count;
    }
    if (a.tid != b.tid) {
      return a.tid < b.tid;
    }
    if (a.q_start != b.q_start) {
      return a.q_start < b.q_start;
    }
    if (a.q_end != b.q_end) {
      return a.q_end < b.q_end;
    }
    if (a.t_start != b.t_start) {
      return a.t_start < b.t_start;
    }
    if (a.t_end != b.t_end) {
      return a.t_end < b.t_end;
    }
    return a.orient < b.orient;
  });
}

void chain_read_fast_into(const ReadRecord &read, const CompactPostingIndex &index,
                          const std::vector<TargetRecord> &targets,
                          int diag_bin, uint32_t min_chain,
                          std::vector<HitOut> &hits) {
  if (read.anchor_qids.size() != read.anchors.size()) {
    throw std::runtime_error("local query token id invariant violated");
  }
  hits.clear();
  thread_local std::unordered_map<ChainKey, ChainStats, ChainKeyHash> chains;
  chains.clear();
  const size_t desired_chain_capacity = read.anchors.size() * 16 + 1;
  if (chains.bucket_count() < desired_chain_capacity) {
    chains.reserve(desired_chain_capacity);
  }
  for (size_t token_idx = 0; token_idx < read.anchors.size(); ++token_idx) {
    const auto &query = read.anchors[token_idx];
    const uint32_t qid = read.anchor_qids[token_idx];
    if (index.overflow[qid] != 0) {
      continue;
    }
    const auto span = index.spans[qid];
    for (uint32_t i = span.offset; i < span.offset + span.count; ++i) {
      const auto &posting = index.postings[i];
      const uint8_t orient = query.strand ^ posting.strand();
      const int32_t raw_diag =
          orient == 0 ? static_cast<int32_t>(posting.pos) -
                            static_cast<int32_t>(query.pos)
                      : static_cast<int32_t>(posting.pos) +
                            static_cast<int32_t>(query.pos);
      const int32_t diag = raw_diag / std::max(1, diag_bin);
      auto &stats = chains[ChainKey{posting.tid(), diag, orient}];
      ++stats.count;
      stats.q_min = std::min(stats.q_min, query.pos);
      stats.q_max = std::max(stats.q_max, query.pos);
      stats.t_min = std::min(stats.t_min, posting.pos);
      stats.t_max = std::max(stats.t_max, posting.pos);
    }
  }

  hits.reserve(chains.size());
  for (const auto &[key, stats] : chains) {
    if (stats.count < min_chain) {
      continue;
    }
    const uint32_t q_span =
        stats.q_max >= stats.q_min ? (stats.q_max - stats.q_min + 1) : 1;
    const uint32_t t_span =
        stats.t_max >= stats.t_min ? (stats.t_max - stats.t_min + 1) : 1;
    const uint64_t score64 =
        static_cast<uint64_t>(stats.count) * 50ULL +
        static_cast<uint64_t>(std::min(q_span, t_span)) / 20ULL;
    const uint32_t score = static_cast<uint32_t>(
        std::min<uint64_t>(score64, std::numeric_limits<uint32_t>::max()));
    hits.push_back(HitOut{
        key.tid,
        score,
        stats.count,
        stats.q_min,
        std::min<uint32_t>(read.length, stats.q_max + 15),
        stats.t_min,
        std::min<uint32_t>(targets[key.tid].len, stats.t_max + 15),
        key.orient,
        key.diag});
  }
  sort_chain_hits(hits);
}

std::vector<TaxonScore> aggregate_species_scores(
    const std::vector<HitOut> &hits, const std::vector<TargetRecord> &targets) {
  thread_local std::unordered_map<uint32_t, uint64_t> mass_by_species;
  mass_by_species.clear();
  if (mass_by_species.bucket_count() < hits.size()) {
    mass_by_species.reserve(hits.size());
  }
  for (const auto &hit : hits) {
    const auto &target = targets[hit.tid];
    if (target.species == 0) {
      continue;
    }
    mass_by_species[target.species] += hit.score;
  }
  std::vector<TaxonScore> scores;
  scores.reserve(mass_by_species.size());
  for (const auto &[taxid, score] : mass_by_species) {
    scores.push_back({taxid,
                      static_cast<uint32_t>(std::min<uint64_t>(
                          score, std::numeric_limits<uint32_t>::max()))});
  }
	  std::sort(scores.begin(), scores.end(),
	            [](const TaxonScore &a, const TaxonScore &b) {
	              if (a.score != b.score) {
	                return a.score > b.score;
	              }
	              return std::to_string(a.taxid) < std::to_string(b.taxid);
	            });
  return scores;
}

void chain_reads_to_call_store(
    const std::vector<std::string> &read_files, uint32_t k, uint32_t w,
    const QueryHashIndex &query_hashes, const CompactPostingIndex &index,
    const std::vector<TargetRecord> &targets, int diag_bin, uint32_t min_chain,
    uint32_t threads, ChimeraClassify::LocalResolutionCallStore &store,
    std::atomic<uint64_t> &local_hits, std::atomic<uint64_t> &local_absent) {
  for_each_pending_read_batch(
      read_files, kLocalResolutionReadBatchSize,
      [&](const std::vector<PendingRead> &pending) {
        std::vector<std::vector<TaxonScore>> batch_scores(pending.size());
        const uint32_t worker_count = std::max<uint32_t>(
            1, std::min<uint32_t>(
                   threads == 0 ? 1 : threads,
                   static_cast<uint32_t>(std::max<size_t>(1, pending.size()))));
        std::atomic<size_t> next{0};
        auto worker = [&]() {
          std::vector<HitOut> hit_scratch;
          while (true) {
            const size_t idx = next.fetch_add(1, std::memory_order_relaxed);
            if (idx >= pending.size()) {
              break;
            }
            auto read = make_read_record_from_pending(pending[idx], k, w);
            assign_query_hashes(read, query_hashes);
            chain_read_fast_into(read, index, targets, diag_bin, min_chain,
                                 hit_scratch);
            batch_scores[idx] = aggregate_species_scores(hit_scratch, targets);
          }
        };
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        for (uint32_t i = 0; i < worker_count; ++i) {
          workers.emplace_back(worker);
        }
        for (auto &worker_thread : workers) {
          worker_thread.join();
        }
        for (size_t i = 0; i < pending.size(); ++i) {
          if (pending[i].ordinal != store.read_count()) {
            throw std::runtime_error(
                "local resolution read ordinal stream is not contiguous");
          }
          if (batch_scores[i].empty()) {
            local_absent.fetch_add(1, std::memory_order_relaxed);
          } else {
            local_hits.fetch_add(1, std::memory_order_relaxed);
          }
          append_local_resolution_scores(store, batch_scores[i]);
        }
      });
}

} // namespace

namespace ChimeraClassify {

namespace {

LocalResolutionResult
run_local_resolution_engine_impl(const std::vector<std::string> &read_files,
                                 const std::filesystem::path &index_path,
                                 const std::filesystem::path &shard_manifest_path,
                                 const TargetFilter &target_filter,
                                 uint32_t diag_bin, uint32_t max_occ,
                                 uint32_t min_chain, uint32_t threads) {
  if (read_files.empty()) {
    throw std::runtime_error("Local resolution route requires --single input");
  }
  if (index_path.empty()) {
    throw std::runtime_error("Local resolution route requires an index file");
  }
  if (shard_manifest_path.empty()) {
    throw std::runtime_error(
        "Local resolution route requires a shard manifest file");
  }

	  const auto started = std::chrono::steady_clock::now();
	  const auto root_meta = chimera::native_bounded::read_index_header(index_path);
	  QueryHashIndex query_hashes;
	  const uint64_t read_count = collect_query_hashes_from_reads(
	      read_files, root_meta.k, root_meta.w, query_hashes, threads);
	  const auto reads_loaded = std::chrono::steady_clock::now();

  LoadStats stats;
  std::vector<TargetRecord> targets;
  CompactPostingIndex index;
  const bool used_direct_compact = try_load_direct_lpc_targets_compact(
      index_path, shard_manifest_path, target_filter, query_hashes, max_occ,
      threads, targets, index, stats);
  if (!used_direct_compact) {
    std::vector<MatchedPosting> matched_postings;
    targets = load_lpc_targets(index_path, shard_manifest_path, target_filter,
                               query_hashes, matched_postings, stats, threads);
    const auto index_finalize_started = std::chrono::steady_clock::now();
    index = build_compact_posting_index(matched_postings, query_hashes.size(),
                                        max_occ, stats);
    std::vector<MatchedPosting>().swap(matched_postings);
    const auto index_finalized = std::chrono::steady_clock::now();
    stats.posting_merge_seconds =
        std::chrono::duration<double>(index_finalized - index_finalize_started)
            .count();
  }
  const auto refs_loaded = std::chrono::steady_clock::now();

	  LocalResolutionResult result;
	  result.calls.offsets.reserve(static_cast<size_t>(
	      std::min<uint64_t>(read_count + 1,
	                         static_cast<uint64_t>(
	                             std::numeric_limits<size_t>::max()))));
	  result.calls.offsets.clear();
	  result.calls.offsets.push_back(0);
	  std::atomic<uint64_t> local_hits{0};
	  std::atomic<uint64_t> local_absent{0};
	  chain_reads_to_call_store(read_files, root_meta.k, root_meta.w, query_hashes,
	                            index, targets, diag_bin, min_chain, threads,
	                            result.calls, local_hits, local_absent);
	  if (result.calls.read_count() != read_count) {
	    throw std::runtime_error(
	        "local resolution call store read count mismatch");
	  }
	  const auto chained = std::chrono::steady_clock::now();

	  result.stats.reads = read_count;
  result.stats.query_hashes = query_hashes.size();
  result.stats.target_filter = target_filter.names.size();
  result.stats.target_routes = target_filter.route_by_name.size();
  result.stats.scanned_shards = stats.scanned_shards;
  result.stats.skipped_shards = stats.skipped_shards;
  result.stats.selected_targets = stats.selected_targets;
  result.stats.skipped_targets = stats.skipped_targets;
  result.stats.direct_targets = stats.direct_targets;
  result.stats.target_anchor_records_scanned =
      stats.target_anchor_records_scanned;
  result.stats.target_anchor_bytes_read = stats.target_anchor_bytes_read;
  result.stats.target_anchor_records_matched =
      stats.target_anchor_records_matched;
  result.stats.target_hash_prefilter_rejects =
      stats.target_hash_prefilter_rejects;
  result.stats.direct_load_batches = stats.direct_load_batches;
  result.stats.pread_calls = stats.pread_calls;
  result.stats.pread_bytes = stats.pread_bytes;
  result.stats.raw_chain_records = stats.raw_chain_records;
  result.stats.kept_chain_records = stats.kept_chain_records;
  result.stats.index_hash_keys = stats.index_hash_keys;
  result.stats.overflow_hash_keys = stats.overflow_hash_keys;
  result.stats.dropped_broad_keys = stats.dropped_broad_keys;
  result.stats.dropped_broad_records = stats.dropped_broad_records;
  result.stats.local_hits = local_hits.load(std::memory_order_relaxed);
  result.stats.local_absent = local_absent.load(std::memory_order_relaxed);
	  result.stats.threads = std::max<uint32_t>(1, threads == 0 ? 1 : threads);
  result.stats.k = root_meta.k;
  result.stats.w = root_meta.w;
  result.stats.read_seconds =
      std::chrono::duration<double>(reads_loaded - started).count();
  result.stats.ref_seconds =
      std::chrono::duration<double>(refs_loaded - reads_loaded).count();
  result.stats.target_io_seconds = stats.target_io_seconds;
  result.stats.target_read_seconds = stats.target_read_seconds;
  result.stats.target_filter_seconds = stats.target_filter_seconds;
  result.stats.target_collect_seconds = stats.target_collect_seconds;
  result.stats.posting_merge_seconds = stats.posting_merge_seconds;
  result.stats.index_finalize_seconds = stats.index_finalize_seconds;
  result.stats.chain_seconds =
      std::chrono::duration<double>(chained - refs_loaded).count();
  return result;
}

} // namespace

LocalResolutionResult
run_local_resolution_engine(const LocalResolutionRequest &request) {
  const auto target_filter = target_filter_from_targets(request.targets);
  return run_local_resolution_engine_impl(
      request.read_files, std::filesystem::path(request.index_file),
      std::filesystem::path(request.shard_manifest_file), target_filter,
      request.diag_bin, request.max_occ, request.min_chain, request.threads);
}

} // namespace ChimeraClassify
