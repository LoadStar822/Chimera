/*
 * -----------------------------------------------------------------------------
 * Filename:      ChimeraClassify.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Description:
 *  Entry for Chimera classify pipeline (split from legacy monolithic file)
 * -----------------------------------------------------------------------------
 */
#include "ChimeraClassifyCommon.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <thread>

namespace ChimeraClassify {

void run(ClassifyConfig config) {
  if (config.threads == 0) {
    unsigned int hardwareThreads = std::thread::hardware_concurrency();
    if (hardwareThreads == 0) {
      hardwareThreads = 1;
    }
    const auto maxThreads =
        static_cast<unsigned int>(std::numeric_limits<uint16_t>::max());
    if (hardwareThreads > maxThreads) {
      hardwareThreads = maxThreads;
    }
    config.threads = static_cast<uint16_t>(hardwareThreads);
  }

  if (!config.em) {
    config.em = true;
  }
  if (config.verbose) {
    std::cout << config << std::endl;
  }
  omp_set_num_threads(config.threads);
  auto TotalclassifyStart = std::chrono::high_resolution_clock::now();
  auto readStart = std::chrono::high_resolution_clock::now();
  std::cout << "Reading input files..." << std::endl;

  FileInfo fileInfo;
  seqan3::contrib::bgzf_thread_count = config.threads;
  moodycamel::ConcurrentQueue<batchReads> readQueue;
  std::vector<classifyResult> classifyResults;
  std::unordered_map<std::string, double> classWeights;
  bool posteriorModelUsed = false;
  long long rebuildActiveMs = 0;

  std::atomic<bool> producer_done{false};
  auto readEnd = readStart;
  std::thread producer([&]() {
    parseReads(readQueue, config, fileInfo);
    readEnd = std::chrono::high_resolution_clock::now();
    if (config.verbose) {
      auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
          readEnd - readStart);
      std::cout << "\nRead time: ";
      print_classify_time(readDuration.count());
    }
    producer_done.store(true, std::memory_order_release);
  });

  std::vector<std::vector<std::string>> indexToTaxid;
  chimera::imcf::InterleavedMergedCuckooFilter imcf;
  ChimeraBuild::IMCFConfig imcfConfig;
  chimera::presence::CoverageMeta coverageMeta;
  auto indexStatus =
      loadFilter(config.dbFile, imcf, imcfConfig, indexToTaxid, &coverageMeta);
  WeightingContext weightCtx;
  std::unique_ptr<CountMinSketch> freqSketch;
  if (coverageMeta.freq_model.enabled()) {
    try {
      freqSketch = std::make_unique<CountMinSketch>(
          coverageMeta.freq_model.depth, coverageMeta.freq_model.width,
          coverageMeta.freq_model.counters);
      weightCtx.freqSketch = freqSketch.get();
      weightCtx.freqStats = coverageMeta.freq_model.stats;
      weightCtx.freqQuantile = coverageMeta.freq_model.quantile;
    } catch (const std::exception &ex) {
      std::cerr << "Warning: 无法加载频率 Sketch，回退到旧权重模型: " << ex.what()
                << std::endl;
    }
  }
  auto normalize_kind = [](std::string &value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) {
                     return static_cast<char>(std::tolower(ch));
                   });
  };
  std::string resolvedKind = imcfConfig.taxonomyKind;
  if (resolvedKind.empty()) {
    resolvedKind = "ncbi";
  } else {
    normalize_kind(resolvedKind);
  }
  if (config.taxonomyKind == "auto" || config.taxonomyKind.empty()) {
    config.taxonomyKind = resolvedKind;
  } else {
    std::string requestedKind = config.taxonomyKind;
    normalize_kind(requestedKind);
    if (requestedKind != resolvedKind) {
      throw std::runtime_error("数据库 taxonomy_kind (" + resolvedKind +
                               ") 与分类请求 (" + requestedKind +
                               ") 不匹配，请检查参数 –-taxonomy-kind。");
    }
    config.taxonomyKind = requestedKind;
  }
  std::string resolvedVersion = imcfConfig.taxonomyVersion;
  if (resolvedVersion.empty()) {
    resolvedVersion = resolvedKind == "gtdb" ? "gtdb-auto" : "ncbi-taxdump";
  }
  if (config.taxonomyVersion == "auto" || config.taxonomyVersion.empty()) {
    config.taxonomyVersion = resolvedVersion;
  } else {
    if (config.taxonomyVersion != resolvedVersion) {
      std::ostringstream oss;
      oss << "数据库 taxonomy_version (" << resolvedVersion
          << ") 与分类请求 (" << config.taxonomyVersion
          << ") 不一致，请检查参数 –-taxonomy-version。";
      throw std::runtime_error(oss.str());
    }
  }
  if (indexStatus.builtActive) {
    rebuildActiveMs = indexStatus.activeMs;
    std::cout << "IMCF index: active-group list missing, rebuilding in memory ("
              << rebuildActiveMs << " ms)" << std::endl;
  }
  FeatureMethod db_method =
      (imcfConfig.featureMethod == 1) ? FeatureMethod::Strobemer
                                      : FeatureMethod::Syncmer;
  FeatureMethod user_method = parse_feature_method_string(config.feature);
  FeatureMethod desired_method = user_method;
  std::string db_method_str = feature_method_to_string(db_method);
  std::string user_method_str = feature_method_to_string(user_method);
  uint8_t desired_k = config.strobemer_k;
  uint8_t desired_order = config.strobemer_order;
  uint16_t desired_w_min = config.strobemer_w_min;
  uint16_t desired_w_max = config.strobemer_w_max;

  if (user_method == FeatureMethod::Auto) {
    ReadStats stats = sample_read_stats(config);
    double avg_len =
        (stats.count == 0)
            ? 0.0
            : static_cast<double>(stats.total_len) /
                  static_cast<double>(stats.count);
    if (stats.count == 0) {
      desired_method = db_method;
    } else {
      size_t representative_len = static_cast<size_t>(std::llround(avg_len));
      if (representative_len == 0)
        representative_len = stats.max_len;
      chimera::feature::Params suggested =
          chimera::feature::auto_params_from_readlen(representative_len);
      if (suggested.method == chimera::feature::Method::Strobemer) {
        size_t required = chimera::feature::min_required_length(suggested);
        if (stats.max_len < required) {
          suggested.method = chimera::feature::Method::Syncmer;
        }
      }
      desired_method = static_cast<FeatureMethod>(suggested.method);
      if (suggested.method == chimera::feature::Method::Strobemer) {
        desired_k = suggested.strobe.k;
        desired_order = suggested.strobe.order;
        desired_w_min = suggested.strobe.w_min;
        desired_w_max = suggested.strobe.w_max;
      } else {
        desired_k = 0;
        desired_order = 0;
        desired_w_min = 0;
        desired_w_max = 0;
      }
    }
  }

  FeatureMethod final_method = desired_method;
  if (final_method == FeatureMethod::Strobemer) {
    if (db_method != FeatureMethod::Strobemer) {
      std::cout << "[info] 输入 reads 建议使用 strobemer，但数据库为 syncmer，自动改用 syncmer。"
                << std::endl;
      final_method = db_method;
    } else {
      if (imcfConfig.strobeK == 0) {
        throw std::runtime_error("IMCF 数据库缺少 strobemer 参数，无法分类。");
      }
      if (user_method == FeatureMethod::Auto) {
        if (desired_k != 0 &&
            (desired_k != imcfConfig.strobeK ||
             desired_order != imcfConfig.strobeOrder ||
             desired_w_min != imcfConfig.strobeWmin ||
             desired_w_max != imcfConfig.strobeWmax)) {
          std::cout
              << "[warn] 自动参数建议 strobemer(k=" << static_cast<int>(desired_k)
              << ", order=" << static_cast<int>(desired_order) << ", w=["
              << desired_w_min << ',' << desired_w_max
              << "])，但数据库为 strobemer(k="
              << static_cast<int>(imcfConfig.strobeK)
              << ", order=" << static_cast<int>(imcfConfig.strobeOrder)
              << ", w=[" << imcfConfig.strobeWmin << ','
              << imcfConfig.strobeWmax
              << "])，将沿用数据库参数。如需匹配自动建议，请重新构建数据库。"
              << std::endl;
        }
      }
      config.strobemer_k = imcfConfig.strobeK;
      config.strobemer_order = imcfConfig.strobeOrder;
      config.strobemer_w_min = imcfConfig.strobeWmin;
      config.strobemer_w_max = imcfConfig.strobeWmax;
    }
  } else {
    if (db_method == FeatureMethod::Strobemer) {
      std::cout << "[warn] 自动模式建议使用 syncmer，但数据库是 strobemer，将沿用数据库参数。"
                << std::endl;
      final_method = db_method;
      config.strobemer_k = imcfConfig.strobeK;
      config.strobemer_order = imcfConfig.strobeOrder;
      config.strobemer_w_min = imcfConfig.strobeWmin;
      config.strobemer_w_max = imcfConfig.strobeWmax;
    } else {
      config.strobemer_k = 0;
      config.strobemer_order = 0;
      config.strobemer_w_min = 0;
      config.strobemer_w_max = 0;
    }
  }

  if (final_method != db_method) {
    std::ostringstream oss;
    oss << "数据库构建特征方法为 " << feature_method_to_string(db_method)
        << "，但分类配置要求 " << feature_method_to_string(final_method)
        << "，请调整 --feature 参数或重新构建数据库。";
    throw std::runtime_error(oss.str());
  }

  std::string final_method_str = feature_method_to_string(final_method);
  config.feature = final_method_str;

  if (final_method == FeatureMethod::Strobemer &&
      !chimera::feature::strobemer_available()) {
    throw std::runtime_error("当前 Chimera 构建未启用 strobemer 支持，无法加载使用 strobemer 的数据库，请重新编译或改用 syncmer 数据库。");
  }

  size_t feature_min_len = 0;
  chimera::feature::Params feature_params =
      prepare_feature_params_for_classify(imcfConfig, config, final_method,
                                          feature_min_len);

  if (config.verbose) {
    if (final_method == FeatureMethod::Strobemer) {
      std::cout << "Feature method: strobemer (k="
                << static_cast<int>(config.strobemer_k)
                << ", order=" << static_cast<int>(config.strobemer_order)
                << ", w=[" << config.strobemer_w_min << ','
                << config.strobemer_w_max
                << "], seed=" << static_cast<unsigned long long>(imcfConfig.seed64)
                << ")" << std::endl;
    } else {
      std::cout << "Feature method: syncmer (k="
                << static_cast<int>(imcfConfig.kmerSize)
                << ", s=" << imcfConfig.smerSize
                << ", pos=" << imcfConfig.syncmerPosition
                << ", seed=" << static_cast<unsigned long long>(imcfConfig.seed64)
                << ")" << std::endl;
    }
  }

  const TaxDict tax = build_tax_dict(indexToTaxid);
  if (config.decoy_reps == 0 && !(config.presence_noise > 0.0) &&
      config.verbose) {
    std::cerr << "Warning: decoy_reps=0，覆盖模型的噪声将退回默认 μ=1e-4；建议 decoy_reps>=1"
              << std::endl;
  }
  PresenceSummary presenceSummary(static_cast<size_t>(config.decoy_reps));
  PresenceSummary *presencePtr = &presenceSummary;
  uint64_t presenceSeed = 0;
  std::hash<std::string> hasher;
  presenceSeed = hasher(config.outputFile);
  presenceSeed ^= (hasher(config.dbFile) << 1);
  presenceSeed ^= static_cast<uint64_t>(imcfConfig.fpSalt);

  auto classifyStart = std::chrono::high_resolution_clock::now();
  std::cout << "Classifying sequences by imcf (feature=" << config.feature
            << ")..." << std::endl;
  classify_streaming(imcfConfig, readQueue, config, imcf, indexToTaxid, tax,
                     classifyResults, fileInfo, producer_done, feature_params,
                     feature_min_len, weightCtx, presencePtr, presenceSeed);
  auto classifyEnd = std::chrono::high_resolution_clock::now();
  auto classifyDuration =
      std::chrono::duration_cast<std::chrono::milliseconds>(classifyEnd -
                                                            classifyStart);
  producer.join();
  if (config.verbose) {
    std::cout << "Classify time: ";
    print_classify_time(classifyDuration.count());
    if (classifyDuration.count() > 0 && fileInfo.sequenceNum > 0) {
      double readsPerSec =
          static_cast<double>(fileInfo.sequenceNum) /
          (static_cast<double>(classifyDuration.count()) / 1000.0);
      std::cout << "平均分类速度: " << std::fixed << std::setprecision(1)
                << readsPerSec << " reads/s" << std::defaultfloat << std::endl;
    }
  }
  if (fileInfo.sequenceNum > 0) {
    fileInfo.avgLen = fileInfo.bpLength / fileInfo.sequenceNum;
  }
  if (config.verbose && fileInfo.sequenceNum > 0) {
    size_t min_print =
        (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength)
            ? 0
            : fileInfo.minLen;
    std::cout << "Read length stats: min=" << min_print
              << ", max=" << fileInfo.maxLen << ", avg=" << fileInfo.avgLen
              << std::endl;
  }

  PresenceDecision presenceDecision =
      evaluate_presence_coverage(presenceSummary, tax, config, coverageMeta,
                                 fileInfo.sequenceNum, fileInfo.avgLen);
  if (config.verbose) {
    auto oldFlags = std::cout.flags();
    auto oldPrecision = std::cout.precision();
    std::cout << "Presence caller (coverage): tests=" << presenceDecision.tested
              << ", accepted=" << presenceDecision.acceptedCount
              << ", mu=" << std::scientific << std::setprecision(3)
              << presenceDecision.noiseMu << ", pi=" << config.presence_pi
              << ", tau=" << config.presence_tau << std::defaultfloat
              << std::endl;
    std::cout.flags(oldFlags);
    std::cout.precision(oldPrecision);
  }
  PresenceFilterStats filterStats{};
  if (!config.em) {
    filterStats =
        apply_presence_filter(presenceDecision, tax, classifyResults, fileInfo);
    if (config.verbose &&
        (filterStats.trimmedAssignments > 0 ||
         filterStats.forcedUnclassified > 0)) {
      std::cout << "Presence filter: trimmed " << filterStats.trimmedAssignments
                << " assignments, forced " << filterStats.forcedUnclassified
                << " reads to unclassified" << std::endl;
    }
  } else if (config.verbose) {
    std::cout
        << "[Info] EM enabled: Skipping pre-EM presence filter to preserve "
           "shared candidates."
        << std::endl;
  }

  std::unordered_map<std::string, double> emPriorScale;
  if (!coverageMeta.entries.empty()) {
    const uint16_t span =
        (coverageMeta.effective_span > 0) ? coverageMeta.effective_span
                                          : static_cast<uint16_t>(1);
    const size_t read_len_for_prior =
        (fileInfo.avgLen > 0) ? fileInfo.avgLen
                              : static_cast<size_t>(coverageMeta.ref_read_length);
    const double window_current =
        std::max<int64_t>(1, static_cast<int64_t>(read_len_for_prior) -
                                 static_cast<int64_t>(span) + 1);
    const double window_ref =
        std::max<int64_t>(1, static_cast<int64_t>(coverageMeta.ref_read_length) -
                                 static_cast<int64_t>(span) + 1);
    for (const auto &entry : coverageMeta.entries) {
      double expected = 0.0;
      if (entry.unique_density > 0.0) {
        expected = entry.unique_density * window_current;
      } else if (entry.expected_unique_per_ref_read > 0.0 &&
                 window_ref > 0.0) {
        expected =
            entry.expected_unique_per_ref_read * (window_current / window_ref);
      }
      double prior_mix = expected;
      auto it_tid = tax.str2id.find(entry.taxid);
      if (it_tid != tax.str2id.end()) {
        auto it_post = presenceDecision.posteriors.find(it_tid->second);
        if (it_post != presenceDecision.posteriors.end()) {
          double pres = std::max(0.0, it_post->second);
          if (prior_mix > 0.0) {
            prior_mix = 0.5 * prior_mix + 0.5 * pres;
          } else {
            prior_mix = pres;
          }
        }
      }
      if (prior_mix > 0.0) {
        emPriorScale[entry.taxid] = prior_mix;
      }
    }
  }

  if (config.em) {
    auto EMstart = std::chrono::high_resolution_clock::now();
    std::cout << "Running EM algorithm..." << std::endl;
    EMOptions options;
    options.temp = 1.05;
    options.prior_strength = 0.0; // 关闭先验强度，避免剪枝被“回血”
    options.coexist_penalty = 0.0;
    auto [posterior, weights] =
        EMAlgorithm(classifyResults, config.emIter, config.emThreshold, options,
                    emPriorScale.empty() ? nullptr : &emPriorScale);
    classifyResults = std::move(posterior);
    classWeights = std::move(weights);
    posteriorModelUsed = true;
    auto EMend = std::chrono::high_resolution_clock::now();
    auto EMduration =
        std::chrono::duration_cast<std::chrono::milliseconds>(EMend - EMstart);
    if (config.verbose) {
      std::cout << "EM time: ";
      print_classify_time(EMduration.count());
    }
  }
  if (posteriorModelUsed) {
    DecisionConfig decisionConfig;
    decisionConfig.posterior_threshold = config.post_thres;
    decisionConfig.min_class_weight = config.post_pi_min;

    postEmDecision(classifyResults, decisionConfig, classWeights);
    fileInfo.classifiedNum = 0;
    fileInfo.unclassifiedNum = 0;
    for (const auto &result : classifyResults) {
      if (!result.taxidCount.empty() &&
          result.taxidCount.front().first == "unclassified") {
        ++fileInfo.unclassifiedNum;
      } else {
        ++fileInfo.classifiedNum;
      }
    }
  }

  auto accumulate_rejects = [&]() {
    fileInfo.rejectReasons.clear();
    fileInfo.rejectByTaxid.clear();
    for (const auto &res : classifyResults) {
      if (res.reject_reason.empty()) {
        continue;
      }
      fileInfo.rejectReasons[res.reject_reason] += 1;
      std::string hint = res.best_taxid_hint;
      if (hint.empty() && !res.taxidCount.empty() &&
          res.taxidCount.front().first != "unclassified") {
        hint = res.taxidCount.front().first;
      }
      if (!hint.empty()) {
        fileInfo.rejectByTaxid[hint][res.reject_reason] += 1;
      }
    }
  };
  accumulate_rejects();

  auto print_rejects = [&]() {
    if (fileInfo.rejectReasons.empty()) {
      return;
    }
    std::cout << "Reject breakdown (reads):" << std::endl;
    std::vector<std::pair<std::string, size_t>> reasons;
    reasons.reserve(fileInfo.rejectReasons.size());
    for (const auto &kv : fileInfo.rejectReasons) {
      reasons.emplace_back(kv.first, kv.second);
    }
    std::sort(reasons.begin(), reasons.end(),
              [](auto &a, auto &b) { return a.second > b.second; });
    size_t topR = std::min<size_t>(8, reasons.size());
    for (size_t i = 0; i < topR; ++i) {
      std::cout << "  - " << reasons[i].first << ": " << reasons[i].second
                << std::endl;
    }
    std::cout << "  total rejected: " << fileInfo.unclassifiedNum << std::endl;
    std::vector<std::pair<std::string, size_t>> taxa;
    for (const auto &kv : fileInfo.rejectByTaxid) {
      size_t sum = 0;
      for (const auto &r : kv.second)
        sum += r.second;
      taxa.emplace_back(kv.first, sum);
    }
    std::sort(taxa.begin(), taxa.end(),
              [](auto &a, auto &b) { return a.second > b.second; });
    size_t topT = std::min<size_t>(5, taxa.size());
    for (size_t i = 0; i < topT; ++i) {
      std::cout << "  taxid " << taxa[i].first << ": " << taxa[i].second
                << " rejected (";
      const auto &mp = fileInfo.rejectByTaxid[taxa[i].first];
      size_t shown = 0;
      for (auto it = mp.begin(); it != mp.end() && shown < 3; ++it, ++shown) {
        std::cout << it->first << "=" << it->second;
        if (shown + 1 < 3 && std::next(it) != mp.end())
          std::cout << ", ";
      }
      std::cout << ")" << std::endl;
    }
  };
  if (config.verbose) {
    print_rejects();
  }

  auto saveStart = std::chrono::high_resolution_clock::now();
  std::cout << "Saving classification results..." << std::endl;
  saveResult(classifyResults, config);
  auto saveEnd = std::chrono::high_resolution_clock::now();
  auto saveDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
      saveEnd - saveStart);
  if (config.verbose) {
    std::cout << "\nSave time: ";
    print_classify_time(saveDuration.count());
    std::cout << "Total sequences: " << fileInfo.sequenceNum << std::endl;

    const auto format_percentage = [](size_t part, size_t total) {
      std::ostringstream oss;
      if (total == 0) {
        oss << "N/A";
        return oss.str();
      }
      oss.setf(std::ios::fixed);
      oss << std::setprecision(2)
          << static_cast<double>(part) * 100.0 /
                 static_cast<double>(total)
          << '%';
      return oss.str();
    };

    std::cout << "Classified sequences: " << fileInfo.classifiedNum << " ("
              << format_percentage(fileInfo.classifiedNum, fileInfo.sequenceNum)
              << ")" << std::endl;
    std::cout << "Unclassified sequences: " << fileInfo.unclassifiedNum << " ("
              << format_percentage(fileInfo.unclassifiedNum, fileInfo.sequenceNum)
              << ")" << std::endl;
    if (rebuildActiveMs > 0) {
      std::cout << "Index rebuild summary:" << std::endl;
      if (rebuildActiveMs > 0) {
        std::cout << "  Active index: ";
        print_classify_time(rebuildActiveMs);
      }
    }
  }

  auto TotalclassifyEnd = std::chrono::high_resolution_clock::now();
  auto TotalclassifyDuration =
      std::chrono::duration_cast<std::chrono::milliseconds>(TotalclassifyEnd -
                                                            TotalclassifyStart);

  if (config.verbose) {
    std::cout << "\nTotal classify time: ";
    print_classify_time(TotalclassifyDuration.count());
  }
}

} // namespace ChimeraClassify
