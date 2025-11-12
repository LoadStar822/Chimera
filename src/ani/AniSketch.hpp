#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace chimera::ani {

struct AniSketchHeader {
  char magic[6]{'C','H','A','N','I','\0'};
  uint32_t version{0};
  uint16_t k{0};
  uint8_t method{0};
  uint32_t scale{1};
  uint64_t seed64{0};
  uint32_t num_taxa{0};
  uint64_t index_offset{0};
};

struct TaxonSketchRef {
  std::string taxid;
  const uint8_t *data{nullptr};
  uint32_t count{0};
};

class AniSketchReader {
public:
  explicit AniSketchReader(const std::filesystem::path &path);
  ~AniSketchReader();
  AniSketchReader(AniSketchReader &&other) noexcept;
  AniSketchReader &operator=(AniSketchReader &&other) noexcept;
  AniSketchReader(const AniSketchReader &) = delete;
  AniSketchReader &operator=(const AniSketchReader &) = delete;

  bool ok() const { return ok_; }
  const AniSketchHeader &header() const { return header_; }
  const std::filesystem::path &path() const { return path_; }

  void for_each(const std::function<void(const TaxonSketchRef &)> &fn) const;
  const TaxonSketchRef *find(std::string_view taxid) const;

private:
  bool parse();
  void close();

  std::filesystem::path path_;
  int fd_{-1};
  size_t size_{0};
  uint8_t *map_{nullptr};
  AniSketchHeader header_{};
  std::vector<TaxonSketchRef> records_;
  std::unordered_map<std::string, size_t> index_;
  bool ok_{false};
};

struct AniStats {
  uint32_t q_size{0};
  uint32_t t_size{0};
  uint32_t intersect{0};
  double containment{0.0};
  double jaccard{0.0};
  double mashD{1.0};
  double ani{0.0};
};

std::vector<uint64_t> build_query_sketch(const std::vector<uint64_t> &raw,
                                         uint32_t scale);

AniStats score(const std::vector<uint64_t> &query,
               const TaxonSketchRef &target,
               uint16_t k);

} // namespace chimera::ani
