#pragma once

#include <charconv>
#include <cstdint>
#include <string>

namespace chimera::utils {

inline bool try_parse_u32(const std::string &s, uint32_t &out) {
  if (s.empty()) {
    return false;
  }
  uint32_t value = 0;
  const char *begin = s.data();
  const char *end = s.data() + s.size();
  const auto [ptr, ec] = std::from_chars(begin, end, value, 10);
  if (ec != std::errc() || ptr != end) {
    return false;
  }
  out = value;
  return true;
}

} // namespace chimera::utils
