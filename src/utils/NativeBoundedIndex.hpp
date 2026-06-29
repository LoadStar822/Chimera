#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <string>
#include <string_view>
#include <vector>

namespace chimera::native_bounded {

struct Anchor {
  uint64_t hash{};
  uint32_t pos{};
  uint8_t strand{};
};

struct TargetMeta {
  std::string name;
  uint32_t species{};
  uint32_t len{};
  uint64_t anchor_offset{};
  uint32_t anchor_count{};
  uint64_t anchor_bytes{};
};

struct IndexMeta {
  uint32_t version{};
  uint32_t k{};
  uint32_t w{};
  uint64_t anchor_data_offset{};
  std::vector<TargetMeta> targets;
};

std::vector<Anchor> extract_minimizers(std::string_view seq, int k, int w);
std::vector<Anchor> extract_minimizers(const std::vector<seqan3::dna4> &seq,
                                       int k, int w);

uint64_t anchor_record_bytes(uint32_t k);
uint64_t encoded_anchor_bytes(const std::vector<Anchor> &anchors, uint32_t k);
void encode_anchor_block_into(const std::vector<Anchor> &anchors, uint32_t k,
                              std::vector<char> &out);
void encode_anchor_block_into(const std::vector<Anchor> &anchors, uint32_t k,
                              uint64_t reserve_bytes,
                              std::vector<char> &out);
std::vector<char> encode_anchor_block(const std::vector<Anchor> &anchors,
                                      uint32_t k);
std::vector<Anchor> decode_anchor_block(const char *data, size_t size,
                                        uint32_t anchor_count, uint32_t k,
                                        uint32_t version);

uint64_t write_index(const std::filesystem::path &path, uint32_t k, uint32_t w,
                     const std::vector<TargetMeta> &targets,
                     const std::filesystem::path &anchor_spool_path,
                     uint64_t anchor_count);

uint64_t write_index_from_anchor_parts(
    const std::filesystem::path &path, uint32_t k, uint32_t w,
    const std::vector<TargetMeta> &targets,
    const std::vector<std::filesystem::path> &anchor_part_paths,
    uint64_t anchor_count);

uint64_t write_index_header(const std::filesystem::path &path, uint32_t k,
                            uint32_t w,
                            const std::vector<TargetMeta> &targets,
                            uint64_t anchor_count);

void write_index_metadata_only(const std::filesystem::path &path, uint32_t k,
                               uint32_t w,
                               const std::vector<TargetMeta> &targets);

IndexMeta read_index_meta(const std::filesystem::path &path);

IndexMeta read_index_header(const std::filesystem::path &path);

std::vector<Anchor> read_anchor_range(const std::filesystem::path &path,
                                      const IndexMeta &meta,
                                      const TargetMeta &target);

std::vector<Anchor> read_anchor_summary(const std::filesystem::path &path,
                                        const IndexMeta &meta,
                                        const TargetMeta &target,
                                        uint32_t max_anchors);

void write_anchor(std::ostream &out, const Anchor &anchor, uint32_t k);

bool read_anchor(std::istream &in, Anchor &anchor, uint32_t k);

} // namespace chimera::native_bounded
