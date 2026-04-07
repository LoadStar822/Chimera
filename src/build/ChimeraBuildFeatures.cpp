#include "ChimeraBuildCommon.hpp"

#include <stdexcept>

namespace ChimeraBuild {

uint64_t mix64(uint64_t value) { return XXH3_64bits(&value, sizeof(value)); }

void print_build_time(long long milliseconds) {
  long long total_seconds = milliseconds / 1000;
  long long seconds = total_seconds % 60;
  long long total_minutes = total_seconds / 60;
  long long minutes = total_minutes % 60;
  long long hours = total_minutes / 60;

  if (hours > 0) {
    std::cout << hours << "h " << minutes << "min " << seconds << "s "
              << milliseconds % 1000 << "ms" << std::endl;
  } else if (minutes > 0) {
    std::cout << minutes << "min " << seconds << "s "
              << milliseconds % 1000 << "ms" << std::endl;
  } else {
    std::cout << seconds << "s " << milliseconds % 1000 << "ms" << std::endl;
  }
}

chimera::feature::Params make_feature_params(const BuildConfig &config,
                                             uint64_t &seed_out) {
  chimera::feature::Params params{};
  if (!chimera::feature::strobemer_available()) {
    throw std::runtime_error(
        "This Chimera build does not include strobemer support.");
  }
  params.strobe.k = config.strobemer_k;
  params.strobe.order = config.strobemer_order;
  params.strobe.w_min = config.strobemer_w_min;
  params.strobe.w_max = config.strobemer_w_max;
  seed_out = ChimeraBuild::adjust_seed(params.strobe.k);
  params.strobe.seed = seed_out;
  params.strobe.canonical = true;
  return params;
}

} // namespace ChimeraBuild
