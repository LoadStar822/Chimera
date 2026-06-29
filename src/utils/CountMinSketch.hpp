#pragma once

#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>

class CountMinSketch {
public:
  CountMinSketch(uint32_t depth, uint32_t width);
  CountMinSketch(uint32_t depth, uint32_t width,
                 const std::vector<uint32_t> &counters);

  static uint32_t column_index(uint32_t row, uint64_t key, uint32_t width);

  void add(uint64_t key);
  void add_many(const std::vector<uint64_t> &keys);
  uint32_t estimate(uint64_t key) const;
  void freeze();

  uint32_t depth() const { return d_; }
  uint32_t width() const { return w_; }

  void exportCounts(std::vector<uint32_t> &out) const;

private:
  uint32_t d_;
  uint32_t w_;
  size_t size_;
  uint64_t index_mask_{0};
  bool use_mask_{false};
  std::unique_ptr<std::atomic<uint32_t>[]> atomic_table_;
  std::unique_ptr<uint32_t[]> plain_table_;

  uint32_t index(uint32_t row, uint64_t key) const;
  uint32_t load_counter(size_t offset) const;
  void increment_counter(size_t offset);
};
