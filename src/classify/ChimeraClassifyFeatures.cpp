#include "ChimeraClassifyCommon.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace ChimeraClassify {

MarginDecision decide_high_conf(size_t best, size_t second, double eff_eval) {
  MarginDecision dc;
  dc.margin = static_cast<double>(best) - static_cast<double>(second);
  dc.ratio =
      second ? static_cast<double>(best) / static_cast<double>(second)
             : static_cast<double>(best);
  double base = std::max(1.0, std::sqrt(std::max(0.0, eff_eval)));
  dc.need = static_cast<size_t>(std::floor(0.35 * base + 1.0));
  if (best >= 3 && best >= second + dc.need && dc.ratio >= 1.25) {
    dc.accept = true;
  }
  return dc;
}

chimera::feature::Params prepare_feature_params_for_classify(
    const ChimeraBuild::IMCFConfig &imcfConfig,
    FeatureMethod method, size_t &feature_min_len) {
  chimera::feature::Params params{};
  if (method == FeatureMethod::Strobemer) {
    if (imcfConfig.strobeK == 0) {
      throw std::runtime_error(
          "The IMCF database is missing strobemer parameters and cannot run in strobemer mode.");
    }
    params.method = FeatureMethod::Strobemer;
    params.strobe.k = imcfConfig.strobeK;
    params.strobe.order = imcfConfig.strobeOrder;
    params.strobe.w_min = imcfConfig.strobeWmin;
    params.strobe.w_max = imcfConfig.strobeWmax;
    params.strobe.seed = imcfConfig.seed64;
    params.strobe.canonical = true;
  } else {
    params.method = FeatureMethod::Syncmer;
    params.sync.k = imcfConfig.kmerSize;
    params.sync.s = imcfConfig.smerSize;
    params.sync.pos = imcfConfig.syncmerPosition;
    params.sync.seed = imcfConfig.seed64;
    params.sync.canonical = true;
  }
  feature_min_len = chimera::feature::min_required_length(params);
  return params;
}

} // namespace ChimeraClassify
