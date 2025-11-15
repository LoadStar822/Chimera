#include <utils/CountMinSketch.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

CountMinSketch::CountMinSketch(uint32_t depth, uint32_t width)
    : d_(depth), w_(width), size_(static_cast<size_t>(depth) * static_cast<size_t>(width)) {
  if (d_ == 0 || w_ == 0) {
    throw std::invalid_argument("CountMinSketch dimensions must be non-zero");
  }
  table_ = std::make_unique<std::atomic<uint32_t>[]>(size_);
  for (size_t i = 0; i < size_; ++i) {
    table_[i].store(0u, std::memory_order_relaxed);
  }
}

void CountMinSketch::add(uint64_t key) {
  for (uint32_t row = 0; row < d_; ++row) {
    const uint32_t column = index(row, key);
    const size_t offset = static_cast<size_t>(row) * static_cast<size_t>(w_) + column;
    table_[offset].fetch_add(1u, std::memory_order_relaxed);
  }
}

uint32_t CountMinSketch::estimate(uint64_t key) const {
  uint32_t estimate_value = std::numeric_limits<uint32_t>::max();
  for (uint32_t row = 0; row < d_; ++row) {
    const uint32_t column = index(row, key);
    const size_t offset = static_cast<size_t>(row) * static_cast<size_t>(w_) + column;
    const uint32_t value = table_[offset].load(std::memory_order_relaxed);
    estimate_value = std::min(estimate_value, value);
  }
  if (estimate_value == std::numeric_limits<uint32_t>::max()) {
    return 0;
  }
  return estimate_value;
}

uint32_t CountMinSketch::index(uint32_t row, uint64_t key) const {
  const uint64_t salt = 0x9e3779b97f4a7c15ULL * (static_cast<uint64_t>(row) + 1ULL);
  const uint64_t mixed = key ^ salt;
  return static_cast<uint32_t>(mixed % static_cast<uint64_t>(w_));
}
