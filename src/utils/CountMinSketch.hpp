#pragma once

#include <atomic>
#include <cstdint>
#include <memory>

class CountMinSketch {
public:
  CountMinSketch(uint32_t depth, uint32_t width);

  void add(uint64_t key);
  uint32_t estimate(uint64_t key) const;

  uint32_t depth() const { return d_; }
  uint32_t width() const { return w_; }

  template <typename F>
  void forEachCounter(F&& f) const {
    for (size_t i = 0; i < size_; ++i) {
      f(table_[i].load(std::memory_order_relaxed));
    }
  }

private:
  uint32_t d_;
  uint32_t w_;
  size_t size_;
  std::unique_ptr<std::atomic<uint32_t>[]> table_;

  inline uint32_t index(uint32_t row, uint64_t key) const;
};
