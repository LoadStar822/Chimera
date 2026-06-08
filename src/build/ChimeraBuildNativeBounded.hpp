#pragma once

#include "ChimeraBuildCommon.hpp"

#include <filesystem>

namespace ChimeraBuild {

struct NativeBoundedBuildStats {
  uint64_t targets{0};
  uint64_t sequences{0};
  uint64_t bp{0};
  uint64_t anchors{0};
  uint64_t activeRecords{0};
  uint64_t representativeRecords{0};
  double count_seconds{0.0};
  double selection_seconds{0.0};
  double layout_seconds{0.0};
  double write_seconds{0.0};
};

NativeBoundedBuildStats build_native_bounded_index(
    const BuildConfig &config,
    const robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    const std::filesystem::path &outputPath);

} // namespace ChimeraBuild
