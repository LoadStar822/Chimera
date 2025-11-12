#include <ani/AniSketch.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

#if defined(__unix__) || defined(__APPLE__) || defined(__linux__)
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#else
#error "ANI sketch reader currently requires POSIX mmap"
#endif

namespace {

template <typename T>
bool read_primitive(const uint8_t *&ptr, const uint8_t *end, T &out) {
  if (ptr + sizeof(T) > end) {
    return false;
  }
  std::memcpy(&out, ptr, sizeof(T));
  ptr += sizeof(T);
  return true;
}

inline uint64_t load64(const uint8_t *base, uint32_t idx) {
  uint64_t value;
  std::memcpy(&value, base + static_cast<size_t>(idx) * sizeof(uint64_t),
              sizeof(uint64_t));
  return value;
}

} // namespace

namespace chimera::ani {

AniSketchReader::AniSketchReader(const std::filesystem::path &path)
    : path_(path) {
  fd_ = ::open(path.c_str(), O_RDONLY);
  if (fd_ < 0) {
    return;
  }
  struct stat st {
  };
  if (fstat(fd_, &st) != 0) {
    close();
    return;
  }
  if (st.st_size <= 0) {
    close();
    return;
  }
  size_ = static_cast<size_t>(st.st_size);
  map_ = static_cast<uint8_t *>(::mmap(nullptr, size_, PROT_READ, MAP_PRIVATE,
                                       fd_, 0));
  if (map_ == MAP_FAILED) {
    map_ = nullptr;
    close();
    return;
  }
  ::close(fd_);
  fd_ = -1;
  ok_ = parse();
  if (!ok_) {
    close();
  }
}

AniSketchReader::~AniSketchReader() { close(); }

AniSketchReader::AniSketchReader(AniSketchReader &&other) noexcept {
  *this = std::move(other);
}

AniSketchReader &AniSketchReader::operator=(AniSketchReader &&other) noexcept {
  if (this == &other) {
    return *this;
  }
  close();
  path_ = std::move(other.path_);
  fd_ = other.fd_;
  other.fd_ = -1;
  size_ = other.size_;
  other.size_ = 0;
  map_ = other.map_;
  other.map_ = nullptr;
  header_ = other.header_;
  records_ = std::move(other.records_);
  index_ = std::move(other.index_);
  ok_ = other.ok_;
  other.ok_ = false;
  return *this;
}

void AniSketchReader::close() {
  if (map_) {
    ::munmap(map_, size_);
    map_ = nullptr;
  }
  if (fd_ >= 0) {
    ::close(fd_);
    fd_ = -1;
  }
  size_ = 0;
  records_.clear();
  index_.clear();
  ok_ = false;
}

bool AniSketchReader::parse() {
  if (!map_ || size_ < 32) {
    return false;
  }
  const uint8_t *ptr = map_;
  const uint8_t *end = map_ + size_;
  std::memcpy(header_.magic, ptr, sizeof(header_.magic));
  ptr += sizeof(header_.magic);
  const char expected[6] = {'C', 'H', 'A', 'N', 'I', '\0'};
  if (std::memcmp(header_.magic, expected, sizeof(expected)) != 0) {
    return false;
  }
  if (!read_primitive(ptr, end, header_.version))
    return false;
  if (!read_primitive(ptr, end, header_.k))
    return false;
  if (!read_primitive(ptr, end, header_.method))
    return false;
  if (!read_primitive(ptr, end, header_.scale))
    return false;
  if (!read_primitive(ptr, end, header_.seed64))
    return false;
  if (!read_primitive(ptr, end, header_.num_taxa))
    return false;
  if (!read_primitive(ptr, end, header_.index_offset)) {
    header_.index_offset = 0;
  }

  records_.clear();
  records_.reserve(header_.num_taxa);
  for (uint32_t i = 0; i < header_.num_taxa; ++i) {
    uint16_t len = 0;
    if (!read_primitive(ptr, end, len)) {
      return false;
    }
    if (ptr + len > end) {
      return false;
    }
    std::string taxid(reinterpret_cast<const char *>(ptr), len);
    ptr += len;
    uint32_t count = 0;
    if (!read_primitive(ptr, end, count)) {
      return false;
    }
    size_t bytes = static_cast<size_t>(count) * sizeof(uint64_t);
    if (ptr + bytes > end) {
      return false;
    }
    TaxonSketchRef ref;
    ref.taxid = std::move(taxid);
    ref.data = ptr;
    ref.count = count;
    records_.push_back(std::move(ref));
    ptr += bytes;
  }

  index_.clear();
  index_.reserve(records_.size());
  for (size_t i = 0; i < records_.size(); ++i) {
    index_[records_[i].taxid] = i;
  }
  return true;
}

void AniSketchReader::for_each(
    const std::function<void(const TaxonSketchRef &)> &fn) const {
  if (!ok_) {
    return;
  }
  for (const auto &rec : records_) {
    fn(rec);
  }
}

const TaxonSketchRef *AniSketchReader::find(std::string_view taxid) const {
  if (!ok_) {
    return nullptr;
  }
  auto it = index_.find(std::string(taxid));
  if (it == index_.end()) {
    return nullptr;
  }
  return &records_[it->second];
}

std::vector<uint64_t> build_query_sketch(const std::vector<uint64_t> &raw,
                                         uint32_t scale) {
  std::vector<uint64_t> kept;
  uint32_t effectiveScale = (scale == 0) ? 1u : scale;
  if (raw.empty()) {
    return kept;
  }
  const uint64_t threshold =
      std::numeric_limits<uint64_t>::max() /
      static_cast<uint64_t>(effectiveScale);
  kept.reserve(raw.size() / std::max<uint32_t>(1u, effectiveScale));
  for (uint64_t value : raw) {
    if (value <= threshold) {
      kept.push_back(value);
    }
  }
  if (kept.empty()) {
    return kept;
  }
  std::sort(kept.begin(), kept.end());
  kept.erase(std::unique(kept.begin(), kept.end()), kept.end());
  return kept;
}

AniStats score(const std::vector<uint64_t> &query,
               const TaxonSketchRef &target,
               uint16_t k) {
  AniStats stats;
  stats.q_size = static_cast<uint32_t>(query.size());
  stats.t_size = target.count;
  if (stats.q_size == 0 || stats.t_size == 0) {
    stats.mashD = 1.0;
    stats.ani = 0.0;
    return stats;
  }
  size_t i = 0;
  size_t j = 0;
  while (i < query.size() && j < target.count) {
    uint64_t qv = query[i];
    uint64_t tv = load64(target.data, static_cast<uint32_t>(j));
    if (qv == tv) {
      ++stats.intersect;
      ++i;
      ++j;
    } else if (qv < tv) {
      ++i;
    } else {
      ++j;
    }
  }
  if (stats.intersect == 0) {
    stats.mashD = 1.0;
    stats.ani = 0.0;
    return stats;
  }
  stats.containment = static_cast<double>(stats.intersect) /
                      static_cast<double>(stats.q_size);
  double denom = static_cast<double>(stats.q_size) +
                 static_cast<double>(stats.t_size) -
                 static_cast<double>(stats.intersect);
  stats.jaccard = (denom > 0.0)
                      ? static_cast<double>(stats.intersect) / denom
                      : 0.0;
  uint16_t effK = (k == 0) ? 1 : k;
  if (stats.jaccard > 0.0) {
    double ratio = (2.0 * stats.jaccard) / (1.0 + stats.jaccard);
    ratio = std::clamp(ratio, 1e-12, 1.0);
    stats.mashD = -std::log(ratio) / static_cast<double>(effK);
    stats.ani = std::max(0.0, 1.0 - stats.mashD);
  } else {
    stats.mashD = 1.0;
    stats.ani = 0.0;
  }
  return stats;
}

} // namespace chimera::ani
