#include "NativeBoundedIndex.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <deque>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace chimera::native_bounded {
namespace {

constexpr char kMagic[] = "CHIMERA_NBCIDX_V1";
constexpr uint32_t kIndexVersion = 4;
constexpr uint32_t kMinSupportedIndexVersion = 2;

uint64_t mix64(uint64_t x) {
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}

int base_bits(char c) {
  switch (std::toupper(static_cast<unsigned char>(c))) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return -1;
  }
}

template <class T> void write_raw(std::ostream &out, const T &value) {
  out.write(reinterpret_cast<const char *>(&value), sizeof(T));
  if (!out.good()) {
    throw std::runtime_error("failed to write native bounded index");
  }
}

template <class T> T read_raw(std::istream &in) {
  T value{};
  in.read(reinterpret_cast<char *>(&value), sizeof(T));
  if (!in.good()) {
    throw std::runtime_error("failed to read native bounded index");
  }
  return value;
}

void append_varuint(std::vector<char> &out, uint64_t value) {
  while (value >= 0x80) {
    out.push_back(static_cast<char>((value & 0x7fU) | 0x80U));
    value >>= 7;
  }
  out.push_back(static_cast<char>(value));
}

uint32_t varuint_size(uint64_t value) {
  uint32_t size = 1;
  while (value >= 0x80) {
    value >>= 7;
    ++size;
  }
  return size;
}

bool read_varuint(const char *&cursor, const char *end, uint64_t &value) {
  value = 0;
  uint32_t shift = 0;
  while (cursor < end && shift <= 63) {
    const uint8_t byte = static_cast<uint8_t>(*cursor++);
    value |= static_cast<uint64_t>(byte & 0x7fU) << shift;
    if ((byte & 0x80U) == 0) {
      return true;
    }
    shift += 7;
  }
  return false;
}

uint32_t key_bits_for_k(uint32_t k) {
  return k <= 16 ? static_cast<uint32_t>(2U * k) : 64U;
}

uint64_t bitpacked_key_bytes(uint32_t anchor_count, uint32_t k) {
  const uint32_t key_bits = key_bits_for_k(k);
  return (static_cast<uint64_t>(anchor_count) * key_bits + 7ULL) / 8ULL;
}

void write_bits(std::vector<char> &out, uint64_t bit_offset, uint64_t value,
                uint32_t bits) {
  uint32_t remaining = bits;
  while (remaining > 0) {
    const uint64_t byte_index = bit_offset / 8ULL;
    const uint32_t shift = static_cast<uint32_t>(bit_offset % 8ULL);
    const uint32_t take = std::min<uint32_t>(remaining, 8U - shift);
    const uint64_t mask = (1ULL << take) - 1ULL;
    out[static_cast<size_t>(byte_index)] |=
        static_cast<char>((value & mask) << shift);
    value >>= take;
    bit_offset += take;
    remaining -= take;
  }
}

uint64_t read_bits(const char *data, uint64_t bit_offset, uint32_t bits) {
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

void write_string(std::ostream &out, const std::string &value) {
  if (value.size() > std::numeric_limits<uint32_t>::max()) {
    throw std::runtime_error("native bounded string is too large");
  }
  const uint32_t len = static_cast<uint32_t>(value.size());
  write_raw(out, len);
  out.write(value.data(), static_cast<std::streamsize>(value.size()));
  if (!out.good()) {
    throw std::runtime_error("failed to write native bounded string");
  }
}

std::string read_string(std::istream &in) {
  const uint32_t len = read_raw<uint32_t>(in);
  std::string value(len, '\0');
  if (len > 0) {
    in.read(value.data(), static_cast<std::streamsize>(len));
    if (!in.good()) {
      throw std::runtime_error("failed to read native bounded string");
    }
  }
  return value;
}

template <typename BaseBits>
std::vector<Anchor> extract_minimizers_buffered(size_t seq_size, int k, int w,
                                                BaseBits base_bits_at) {
  struct Candidate {
    uint64_t rank{};
    Anchor anchor;
  };
  std::vector<Candidate> kmers;
  const uint64_t mask = (1ULL << (2 * k)) - 1ULL;
  uint64_t fwd = 0;
  uint64_t rev = 0;
  int valid = 0;
  for (size_t i = 0; i < seq_size; ++i) {
    const int b = base_bits_at(i);
    if (b < 0) {
      fwd = 0;
      rev = 0;
      valid = 0;
      continue;
    }
    fwd = ((fwd << 2) | static_cast<uint64_t>(b)) & mask;
    rev = (rev >> 2) | (static_cast<uint64_t>(3 - b) << (2 * (k - 1)));
    ++valid;
    if (valid < k) {
      continue;
    }
    const uint32_t pos = static_cast<uint32_t>(i + 1 - k);
    const uint8_t strand = fwd <= rev ? 0 : 1;
    const uint64_t canonical = fwd <= rev ? fwd : rev;
    kmers.push_back(Candidate{mix64(canonical),
                              Anchor{canonical, pos, strand}});
  }

  std::vector<Anchor> out;
  std::deque<size_t> dq;
  uint32_t last_pos = std::numeric_limits<uint32_t>::max();
  for (size_t i = 0; i < kmers.size(); ++i) {
    while (!dq.empty() && dq.front() + static_cast<size_t>(w) <= i) {
      dq.pop_front();
    }
    while (!dq.empty()) {
      const Candidate &old = kmers[dq.back()];
      const Candidate &now = kmers[i];
      if (old.rank < now.rank ||
          (old.rank == now.rank && old.anchor.pos <= now.anchor.pos)) {
        break;
      }
      dq.pop_back();
    }
    dq.push_back(i);
    if (i + 1 >= static_cast<size_t>(w)) {
      const Anchor &m = kmers[dq.front()].anchor;
      if (m.pos != last_pos) {
        out.push_back(m);
        last_pos = m.pos;
      }
    }
  }
  return out;
}

template <typename BaseBits>
std::vector<Anchor> extract_minimizers_streamed(size_t seq_size, int k, int w,
                                                BaseBits base_bits_at) {
  std::vector<Anchor> out;
  out.reserve(seq_size / static_cast<size_t>(w) + 1);
  const uint64_t mask = (1ULL << (2 * k)) - 1ULL;
  uint64_t fwd = 0;
  uint64_t rev = 0;
  int valid = 0;
  uint64_t kmer_index = 0;
  uint32_t last_pos = std::numeric_limits<uint32_t>::max();
  struct WindowEntry {
    uint64_t index{};
    uint64_t rank{};
    Anchor anchor;
  };
  std::deque<WindowEntry> window;
  for (size_t i = 0; i < seq_size; ++i) {
    const int b = base_bits_at(i);
    if (b < 0) {
      fwd = 0;
      rev = 0;
      valid = 0;
      continue;
    }
    fwd = ((fwd << 2) | static_cast<uint64_t>(b)) & mask;
    rev = (rev >> 2) | (static_cast<uint64_t>(3 - b) << (2 * (k - 1)));
    ++valid;
    if (valid < k) {
      continue;
    }
    const uint32_t pos = static_cast<uint32_t>(i + 1 - k);
    const uint8_t strand = fwd <= rev ? 0 : 1;
    const uint64_t canonical = fwd <= rev ? fwd : rev;
    const uint64_t rank = mix64(canonical);
    const Anchor current{canonical, pos, strand};
    while (!window.empty() &&
           window.front().index + static_cast<uint64_t>(w) <= kmer_index) {
      window.pop_front();
    }
    while (!window.empty()) {
      const WindowEntry &old = window.back();
      if (old.rank < rank ||
          (old.rank == rank && old.anchor.pos <= current.pos)) {
        break;
      }
      window.pop_back();
    }
    window.push_back(WindowEntry{kmer_index, rank, current});
    if (kmer_index + 1 >= static_cast<uint64_t>(w)) {
      const Anchor &m = window.front().anchor;
      if (m.pos != last_pos) {
        out.push_back(m);
        last_pos = m.pos;
      }
    }
    ++kmer_index;
  }
  return out;
}

template <typename BaseBits>
std::vector<Anchor> extract_minimizers_impl(size_t seq_size, int k, int w,
                                            BaseBits base_bits_at) {
  if (k <= 0 || w <= 0 || seq_size < static_cast<size_t>(k)) {
    return {};
  }
  if (k > 31) {
    throw std::runtime_error("native bounded k must be <= 31");
  }
  constexpr size_t kStreamedSequenceThreshold = 64 * 1024;
  if (seq_size < kStreamedSequenceThreshold) {
    return extract_minimizers_buffered(seq_size, k, w, base_bits_at);
  }
  return extract_minimizers_streamed(seq_size, k, w, base_bits_at);
}

} // namespace

std::vector<Anchor> extract_minimizers(std::string_view seq, int k, int w) {
  return extract_minimizers_impl(seq.size(), k, w, [&](size_t i) {
    return base_bits(seq[i]);
  });
}

std::vector<Anchor> extract_minimizers(const std::vector<seqan3::dna4> &seq,
                                       int k, int w) {
  return extract_minimizers_impl(seq.size(), k, w, [&](size_t i) {
    return static_cast<int>(seq[i].to_rank());
  });
}

uint64_t anchor_record_bytes(uint32_t k) {
  return k <= 16 ? 9ULL : 13ULL;
}

std::vector<char> encode_anchor_block(const std::vector<Anchor> &anchors,
                                      uint32_t k) {
  std::vector<char> out;
  const uint32_t key_bits = key_bits_for_k(k);
  if (k <= 16) {
    out.assign(static_cast<size_t>(
                   bitpacked_key_bytes(static_cast<uint32_t>(anchors.size()),
                                       k)),
               0);
  }
  out.reserve(static_cast<size_t>(encoded_anchor_bytes(anchors, k)));
  uint32_t previous_pos = 0;
  bool first = true;
  for (size_t i = 0; i < anchors.size(); ++i) {
    const auto &anchor = anchors[i];
    if (k <= 16) {
      const uint64_t max_key =
          key_bits == 64 ? std::numeric_limits<uint64_t>::max()
                         : ((1ULL << key_bits) - 1ULL);
      if (anchor.hash > max_key) {
        throw std::runtime_error("native bounded anchor key does not fit uint32");
      }
      write_bits(out, static_cast<uint64_t>(i) * key_bits, anchor.hash,
                 key_bits);
    } else {
      const auto *ptr = reinterpret_cast<const char *>(&anchor.hash);
      out.insert(out.end(), ptr, ptr + sizeof(anchor.hash));
    }
    if (!first && anchor.pos < previous_pos) {
      throw std::runtime_error("native bounded anchors are not position sorted");
    }
    const uint32_t delta = first ? anchor.pos : anchor.pos - previous_pos;
    const uint64_t packed =
        (static_cast<uint64_t>(delta) << 1) | static_cast<uint64_t>(anchor.strand & 1U);
    append_varuint(out, packed);
    previous_pos = anchor.pos;
    first = false;
  }
  return out;
}

uint64_t encoded_anchor_bytes(const std::vector<Anchor> &anchors, uint32_t k) {
  uint64_t bytes = k <= 16
                       ? bitpacked_key_bytes(
                             static_cast<uint32_t>(anchors.size()), k)
                       : static_cast<uint64_t>(anchors.size()) *
                             sizeof(uint64_t);
  uint32_t previous_pos = 0;
  bool first = true;
  for (const auto &anchor : anchors) {
    if (!first && anchor.pos < previous_pos) {
      throw std::runtime_error("native bounded anchors are not position sorted");
    }
    const uint32_t delta = first ? anchor.pos : anchor.pos - previous_pos;
    const uint64_t packed =
        (static_cast<uint64_t>(delta) << 1) |
        static_cast<uint64_t>(anchor.strand & 1U);
    bytes += varuint_size(packed);
    previous_pos = anchor.pos;
    first = false;
  }
  return bytes;
}

std::vector<Anchor> decode_anchor_block(const char *data, size_t size,
                                        uint32_t anchor_count, uint32_t k,
                                        uint32_t version) {
  std::vector<Anchor> anchors;
  anchors.reserve(anchor_count);
  const uint32_t key_bits = key_bits_for_k(k);
  const uint64_t key_bytes =
      version >= 4 && k <= 16 ? bitpacked_key_bytes(anchor_count, k) : 0;
  if (key_bytes > size) {
    throw std::runtime_error("truncated native bounded anchor keys");
  }
  const char *cursor =
      version >= 4 && k <= 16 ? data + key_bytes : data;
  const char *end = data + size;
  uint32_t previous_pos = 0;
  for (uint32_t i = 0; i < anchor_count; ++i) {
    Anchor anchor;
    if (version >= 4 && k <= 16) {
      anchor.hash =
          read_bits(data, static_cast<uint64_t>(i) * key_bits, key_bits);
    } else if (k <= 16) {
      if (static_cast<size_t>(end - cursor) < sizeof(uint32_t)) {
        throw std::runtime_error("truncated native bounded anchor key");
      }
      uint32_t key = 0;
      std::memcpy(&key, cursor, sizeof(key));
      cursor += sizeof(key);
      anchor.hash = key;
    } else {
      if (static_cast<size_t>(end - cursor) < sizeof(uint64_t)) {
        throw std::runtime_error("truncated native bounded anchor key");
      }
      std::memcpy(&anchor.hash, cursor, sizeof(anchor.hash));
      cursor += sizeof(anchor.hash);
    }
    uint64_t packed = 0;
    if (!read_varuint(cursor, end, packed)) {
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
    anchors.push_back(anchor);
  }
  if (cursor != end) {
    throw std::runtime_error("native bounded anchor block has trailing bytes");
  }
  return anchors;
}

void write_anchor(std::ostream &out, const Anchor &anchor, uint32_t k) {
  if (k <= 16) {
    if (anchor.hash > std::numeric_limits<uint32_t>::max()) {
      throw std::runtime_error("native bounded anchor key does not fit uint32");
    }
    const uint32_t key = static_cast<uint32_t>(anchor.hash);
    write_raw(out, key);
  } else {
    write_raw(out, anchor.hash);
  }
  write_raw(out, anchor.pos);
  write_raw(out, anchor.strand);
}

bool read_anchor(std::istream &in, Anchor &anchor, uint32_t k) {
  if (k <= 16) {
    uint32_t key = 0;
    in.read(reinterpret_cast<char *>(&key), sizeof(key));
    if (!in.good()) {
      return false;
    }
    anchor.hash = key;
  } else {
    in.read(reinterpret_cast<char *>(&anchor.hash), sizeof(anchor.hash));
    if (!in.good()) {
      return false;
    }
  }
  in.read(reinterpret_cast<char *>(&anchor.pos), sizeof(anchor.pos));
  in.read(reinterpret_cast<char *>(&anchor.strand), sizeof(anchor.strand));
  if (!in.good()) {
    throw std::runtime_error("truncated native bounded anchor record");
  }
  return true;
}

uint64_t write_index(const std::filesystem::path &path, uint32_t k, uint32_t w,
                     const std::vector<TargetMeta> &targets,
                     const std::filesystem::path &anchor_spool_path,
                     uint64_t anchor_count) {
  return write_index_from_anchor_parts(path, k, w, targets, {anchor_spool_path},
                                       anchor_count);
}

uint64_t write_index_header(const std::filesystem::path &path, uint32_t k,
                            uint32_t w,
                            const std::vector<TargetMeta> &targets,
                            uint64_t anchor_count) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("failed to open native bounded index: " +
                             path.string());
  }
  out.write(kMagic, sizeof(kMagic));
  write_raw(out, kIndexVersion);
  write_raw(out, k);
  write_raw(out, w);
  write_raw(out, static_cast<uint64_t>(targets.size()));
  write_raw(out, anchor_count);
  for (const auto &target : targets) {
    write_string(out, target.name);
    write_raw(out, target.species);
    write_raw(out, target.len);
    write_raw(out, target.anchor_offset);
    write_raw(out, target.anchor_count);
    write_raw(out, target.anchor_bytes);
  }
  const uint64_t anchor_data_offset = static_cast<uint64_t>(out.tellp());
  if (!out.good()) {
    throw std::runtime_error("failed to finalize native bounded index header");
  }
  return anchor_data_offset;
}

uint64_t write_index_from_anchor_parts(
    const std::filesystem::path &path, uint32_t k, uint32_t w,
    const std::vector<TargetMeta> &targets,
    const std::vector<std::filesystem::path> &anchor_part_paths,
    uint64_t anchor_count) {
  const uint64_t anchor_data_offset =
      write_index_header(path, k, w, targets, anchor_count);
  std::ofstream out(path, std::ios::binary | std::ios::app);
  if (!out) {
    throw std::runtime_error("failed to open native bounded index: " +
                             path.string());
  }

  std::vector<char> buffer(1 << 20);
  for (const auto &anchor_part_path : anchor_part_paths) {
    std::ifstream anchors(anchor_part_path, std::ios::binary);
    if (!anchors) {
      throw std::runtime_error("failed to open native bounded anchor spool: " +
                               anchor_part_path.string());
    }
    while (anchors) {
      anchors.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
      const std::streamsize got = anchors.gcount();
      if (got > 0) {
        out.write(buffer.data(), got);
      }
    }
  }
  if (!out.good()) {
    throw std::runtime_error("failed to finalize native bounded index");
  }
  return anchor_data_offset;
}

void write_index_metadata_only(const std::filesystem::path &path, uint32_t k,
                               uint32_t w,
                               const std::vector<TargetMeta> &targets) {
  (void)targets;
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("failed to open native bounded index: " +
                             path.string());
  }
  out.write(kMagic, sizeof(kMagic));
  write_raw(out, kIndexVersion);
  write_raw(out, k);
  write_raw(out, w);
  write_raw(out, uint64_t{0});
  write_raw(out, uint64_t{0});
  if (!out.good()) {
    throw std::runtime_error("failed to finalize native bounded index");
  }
}

IndexMeta read_index_meta(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("failed to open native bounded index: " +
                             path.string());
  }
  char magic[sizeof(kMagic)]{};
  in.read(magic, sizeof(magic));
  if (!in.good() || std::memcmp(magic, kMagic, sizeof(kMagic)) != 0) {
    throw std::runtime_error("invalid native bounded index magic: " +
                             path.string());
  }
  const uint32_t version = read_raw<uint32_t>(in);
  if (version < kMinSupportedIndexVersion || version > kIndexVersion) {
    throw std::runtime_error("unsupported native bounded index version");
  }
  IndexMeta meta;
  meta.version = version;
  meta.k = read_raw<uint32_t>(in);
  meta.w = read_raw<uint32_t>(in);
  const uint64_t target_count = read_raw<uint64_t>(in);
  (void)read_raw<uint64_t>(in); // total anchor count
  if (target_count > static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
    throw std::runtime_error("native bounded index has too many targets");
  }
  meta.targets.reserve(static_cast<size_t>(target_count));
  for (uint64_t i = 0; i < target_count; ++i) {
    TargetMeta target;
    target.name = read_string(in);
    target.species = read_raw<uint32_t>(in);
    target.len = read_raw<uint32_t>(in);
    target.anchor_offset = read_raw<uint64_t>(in);
    target.anchor_count = read_raw<uint32_t>(in);
    if (version >= 3) {
      target.anchor_bytes = read_raw<uint64_t>(in);
    } else {
      target.anchor_bytes = static_cast<uint64_t>(target.anchor_count) *
                            anchor_record_bytes(meta.k);
    }
    meta.targets.push_back(std::move(target));
  }
  meta.anchor_data_offset = static_cast<uint64_t>(in.tellg());
  return meta;
}

IndexMeta read_index_header(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("failed to open native bounded index: " +
                             path.string());
  }
  char magic[sizeof(kMagic)]{};
  in.read(magic, sizeof(magic));
  if (!in.good() || std::memcmp(magic, kMagic, sizeof(kMagic)) != 0) {
    throw std::runtime_error("invalid native bounded index magic: " +
                             path.string());
  }
  const uint32_t version = read_raw<uint32_t>(in);
  if (version < kMinSupportedIndexVersion || version > kIndexVersion) {
    throw std::runtime_error("unsupported native bounded index version");
  }
  IndexMeta meta;
  meta.version = version;
  meta.k = read_raw<uint32_t>(in);
  meta.w = read_raw<uint32_t>(in);
  return meta;
}

std::vector<Anchor> read_anchor_range(const std::filesystem::path &path,
                                      const IndexMeta &meta,
                                      const TargetMeta &target) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("failed to open native bounded index: " +
                             path.string());
  }
  const uint64_t byte_offset =
      meta.anchor_data_offset +
      (meta.version >= 3
           ? target.anchor_offset
           : target.anchor_offset * anchor_record_bytes(meta.k));
  in.seekg(static_cast<std::streamoff>(byte_offset), std::ios::beg);
  if (!in.good()) {
    throw std::runtime_error("failed to seek native bounded anchors");
  }
  if (meta.version >= 3) {
    std::vector<char> block(static_cast<size_t>(target.anchor_bytes));
    if (!block.empty()) {
      in.read(block.data(), static_cast<std::streamsize>(block.size()));
      if (!in.good()) {
        throw std::runtime_error("truncated native bounded anchor block");
      }
    }
    return decode_anchor_block(block.data(), block.size(), target.anchor_count,
                               meta.k, meta.version);
  }
  std::vector<Anchor> anchors;
  anchors.reserve(target.anchor_count);
  for (uint32_t i = 0; i < target.anchor_count; ++i) {
    Anchor anchor;
    if (!read_anchor(in, anchor, meta.k)) {
      throw std::runtime_error("truncated native bounded anchor range");
    }
    anchors.push_back(anchor);
  }
  return anchors;
}

std::vector<Anchor> read_anchor_summary(const std::filesystem::path &path,
                                        const IndexMeta &meta,
                                        const TargetMeta &target,
                                        uint32_t max_anchors) {
  if (max_anchors == 0 || target.anchor_count <= max_anchors) {
    return read_anchor_range(path, meta, target);
  }

  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("failed to open native bounded index: " +
                             path.string());
  }
  const uint64_t base_offset =
      meta.anchor_data_offset +
      (meta.version >= 3
           ? target.anchor_offset
           : target.anchor_offset * anchor_record_bytes(meta.k));
  if (meta.version >= 3) {
    auto anchors = read_anchor_range(path, meta, target);
    std::vector<Anchor> summary;
    summary.reserve(max_anchors);
    uint32_t last_idx = std::numeric_limits<uint32_t>::max();
    for (uint32_t i = 0; i < max_anchors; ++i) {
      uint32_t idx = static_cast<uint32_t>(
          (static_cast<uint64_t>(i) * anchors.size()) / max_anchors);
      if (idx >= anchors.size()) {
        idx = static_cast<uint32_t>(anchors.size() - 1);
      }
      if (idx == last_idx) {
        continue;
      }
      last_idx = idx;
      summary.push_back(anchors[idx]);
    }
    return summary;
  }
  std::vector<Anchor> anchors;
  anchors.reserve(max_anchors);
  uint32_t last_idx = std::numeric_limits<uint32_t>::max();
  for (uint32_t i = 0; i < max_anchors; ++i) {
    uint32_t idx = static_cast<uint32_t>(
        (static_cast<uint64_t>(i) * target.anchor_count) / max_anchors);
    if (idx >= target.anchor_count) {
      idx = target.anchor_count - 1;
    }
    if (idx == last_idx) {
      continue;
    }
    last_idx = idx;
    const uint64_t byte_offset =
        base_offset + static_cast<uint64_t>(idx) * anchor_record_bytes(meta.k);
    in.seekg(static_cast<std::streamoff>(byte_offset), std::ios::beg);
    if (!in.good()) {
      throw std::runtime_error("failed to seek native bounded anchor summary");
    }
    Anchor anchor;
    if (!read_anchor(in, anchor, meta.k)) {
      throw std::runtime_error("truncated native bounded anchor summary");
    }
    anchors.push_back(anchor);
  }
  return anchors;
}

} // namespace chimera::native_bounded
