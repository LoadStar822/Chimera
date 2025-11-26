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

chimera::feature::Method resolve_feature_method(const BuildConfig &config) {
  if (config.feature == "syncmer")
    return chimera::feature::Method::Syncmer;
  if (config.feature == "strobemer")
    return chimera::feature::Method::Strobemer;
  return chimera::feature::Method::Strobemer;
}

chimera::feature::Params make_feature_params(const BuildConfig &config,
                                             chimera::feature::Method &selected,
                                             uint64_t &seed_out) {
  selected = resolve_feature_method(config);
#ifndef CHIMERA_HAS_STROBEMERS
  if (selected == chimera::feature::Method::Strobemer) {
    selected = chimera::feature::Method::Syncmer;
  }
#endif
  chimera::feature::Params params{};
  if (selected == chimera::feature::Method::Strobemer) {
    params.method = chimera::feature::Method::Strobemer;
    params.strobe.k = config.strobemer_k;
    params.strobe.order = config.strobemer_order;
    params.strobe.w_min = config.strobemer_w_min;
    params.strobe.w_max = config.strobemer_w_max;
    seed_out = ChimeraBuild::adjust_seed(params.strobe.k);
    params.strobe.seed = seed_out;
    params.strobe.canonical = true;
  } else {
    params.method = chimera::feature::Method::Syncmer;
    params.sync.k = config.kmer_size;
    params.sync.s = config.smer_size;
    params.sync.pos = config.syncmer_position;
    seed_out = ChimeraBuild::adjust_seed(config.kmer_size);
    params.sync.seed = seed_out;
    params.sync.canonical = true;
  }
  return params;
}

} // namespace ChimeraBuild
