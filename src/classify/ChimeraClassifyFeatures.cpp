#include "ChimeraClassifyCommon.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace ChimeraClassify {

FeatureMethod parse_feature_method_string(std::string feature) {
  std::transform(feature.begin(), feature.end(), feature.begin(),
                 [](unsigned char ch) {
                   return static_cast<char>(std::tolower(ch));
                 });
  if (feature == "syncmer")
    return FeatureMethod::Syncmer;
  if (feature == "strobemer")
    return FeatureMethod::Strobemer;
  return FeatureMethod::Auto;
}

std::string feature_method_to_string(FeatureMethod method) {
  switch (method) {
  case FeatureMethod::Syncmer:
    return "syncmer";
  case FeatureMethod::Strobemer:
    return "strobemer";
  case FeatureMethod::Auto:
  default:
    return "auto";
  }
}

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
    const ChimeraBuild::IMCFConfig &imcfConfig, ClassifyConfig &config,
    FeatureMethod method, size_t &feature_min_len) {
  chimera::feature::Params params{};
  if (method == FeatureMethod::Strobemer) {
    if (imcfConfig.strobeK == 0) {
      throw std::runtime_error(
          "IMCF 数据库缺少 strobemer 参数，无法以 strobemer 模式分类。");
    }
    auto ensure_match = [](auto &cfg_value, auto reference,
                           const char *name) {
      if (cfg_value == 0) {
        cfg_value = reference;
      } else if (cfg_value != reference) {
        std::ostringstream oss;
        oss << "分类参数 " << name << "=" << static_cast<uint32_t>(cfg_value)
            << " 与数据库设置 " << static_cast<uint32_t>(reference)
            << " 不一致，请调整后重试。";
        throw std::runtime_error(oss.str());
      }
    };
    ensure_match(config.strobemer_k, imcfConfig.strobeK, "--strobe-k");
    ensure_match(config.strobemer_order, imcfConfig.strobeOrder,
                 "--strobe-order");
    ensure_match(config.strobemer_w_min, imcfConfig.strobeWmin,
                 "--strobe-w-min");
    ensure_match(config.strobemer_w_max, imcfConfig.strobeWmax,
                 "--strobe-w-max");

    params.method = FeatureMethod::Strobemer;
    params.strobe.k = config.strobemer_k;
    params.strobe.order = config.strobemer_order;
    params.strobe.w_min = config.strobemer_w_min;
    params.strobe.w_max = config.strobemer_w_max;
    params.strobe.seed = imcfConfig.seed64;
    params.strobe.canonical = true;
  } else {
    if (config.strobemer_k != 0 || config.strobemer_order != 0 ||
        config.strobemer_w_min != 0 || config.strobemer_w_max != 0) {
      throw std::runtime_error(
          "分类参数 --strobe-* 仅在 feature=strobemer 时可用。");
    }
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
