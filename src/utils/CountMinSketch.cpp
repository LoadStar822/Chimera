#include <utils/CountMinSketch.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

CountMinSketch::CountMinSketch(uint32_t depth, uint32_t width)
    : d_(depth), w_(width),
      size_(static_cast<size_t>(depth) * static_cast<size_t>(width)),
      index_mask_(width > 0 ? static_cast<uint64_t>(width - 1u) : 0u),
      use_mask_(width > 0 && (width & (width - 1u)) == 0u) {
  if (d_ == 0 || w_ == 0) {
    throw std::invalid_argument("CountMinSketch dimensions must be non-zero");
  }
  atomic_table_ = std::make_unique<std::atomic<uint32_t>[]>(size_);
  for (size_t i = 0; i < size_; ++i) {
    atomic_table_[i].store(0u, std::memory_order_relaxed);
  }
}

CountMinSketch::CountMinSketch(uint32_t depth, uint32_t width,
                               const std::vector<uint32_t> &counters)
    : d_(depth), w_(width),
      size_(static_cast<size_t>(depth) * static_cast<size_t>(width)),
      index_mask_(width > 0 ? static_cast<uint64_t>(width - 1u) : 0u),
      use_mask_(width > 0 && (width & (width - 1u)) == 0u) {
  if (d_ == 0 || w_ == 0) {
    throw std::invalid_argument("CountMinSketch dimensions must be non-zero");
  }
  if (counters.size() != size_) {
    throw std::invalid_argument("CountMinSketch counters size mismatch");
  }
  plain_table_ = std::make_unique<uint32_t[]>(size_);
  std::memcpy(plain_table_.get(), counters.data(), size_ * sizeof(uint32_t));
}

uint32_t CountMinSketch::column_index(uint32_t row, uint64_t key,
                                      uint32_t width) {
  const uint64_t salt =
      0x9e3779b97f4a7c15ULL * (static_cast<uint64_t>(row) + 1ULL);
  const uint64_t mixed = key ^ salt;
  if (width > 0 && (width & (width - 1u)) == 0u) {
    return static_cast<uint32_t>(mixed & static_cast<uint64_t>(width - 1u));
  }
  return static_cast<uint32_t>(mixed % static_cast<uint64_t>(width));
}

void CountMinSketch::add(uint64_t key) {
  if (!atomic_table_) {
    throw std::logic_error("CountMinSketch::add called after freeze");
  }
  if (d_ == 4) {
    increment_counter(index(0, key));
    increment_counter(static_cast<size_t>(w_) + index(1, key));
    increment_counter(static_cast<size_t>(w_) * 2u + index(2, key));
    increment_counter(static_cast<size_t>(w_) * 3u + index(3, key));
    return;
  }
  for (uint32_t row = 0; row < d_; ++row) {
    const uint32_t column = index(row, key);
    const size_t offset = static_cast<size_t>(row) * static_cast<size_t>(w_) + column;
    increment_counter(offset);
  }
}

uint32_t CountMinSketch::estimate(uint64_t key) const {
  if (d_ == 4) {
    const uint32_t v0 = load_counter(index(0, key));
    const uint32_t v1 = load_counter(static_cast<size_t>(w_) + index(1, key));
    const uint32_t v2 = load_counter(static_cast<size_t>(w_) * 2u + index(2, key));
    const uint32_t v3 = load_counter(static_cast<size_t>(w_) * 3u + index(3, key));
    return std::min(std::min(v0, v1), std::min(v2, v3));
  }
  uint32_t estimate_value = std::numeric_limits<uint32_t>::max();
  for (uint32_t row = 0; row < d_; ++row) {
    const uint32_t column = index(row, key);
    const size_t offset = static_cast<size_t>(row) * static_cast<size_t>(w_) + column;
    const uint32_t value = load_counter(offset);
    estimate_value = std::min(estimate_value, value);
  }
  if (estimate_value == std::numeric_limits<uint32_t>::max()) {
    return 0;
  }
  return estimate_value;
}

void CountMinSketch::freeze() {
  if (plain_table_) {
    return;
  }
  plain_table_ = std::make_unique<uint32_t[]>(size_);
  for (size_t i = 0; i < size_; ++i) {
    plain_table_[i] = atomic_table_[i].load(std::memory_order_relaxed);
  }
  atomic_table_.reset();
}

void CountMinSketch::exportCounts(std::vector<uint32_t> &out) const {
  out.resize(size_);
  for (size_t i = 0; i < size_; ++i) {
    out[i] = load_counter(i);
  }
}

uint32_t CountMinSketch::index(uint32_t row, uint64_t key) const {
  if (use_mask_) {
    return static_cast<uint32_t>(
        (key ^ (0x9e3779b97f4a7c15ULL * (static_cast<uint64_t>(row) + 1ULL))) &
        index_mask_);
  }
  return column_index(row, key, w_);
}

uint32_t CountMinSketch::load_counter(size_t offset) const {
  if (plain_table_) {
    return plain_table_[offset];
  }
  return atomic_table_[offset].load(std::memory_order_relaxed);
}

void CountMinSketch::increment_counter(size_t offset) {
  atomic_table_[offset].fetch_add(1u, std::memory_order_relaxed);
}
