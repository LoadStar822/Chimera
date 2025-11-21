/*
 * -----------------------------------------------------------------------------
 * Filename:      ChimeraClassify.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-08-09
 *
 * Last Modified: 2024-11-18
 *
 * Description:
 *  Classify functions for Chimera
 *
 * Version:
 *  1.4
 * -----------------------------------------------------------------------------
 */
#include "ChimeraClassify.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <numeric>
#include <queue>
#include <sstream>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <span>

#include <utils/FeatureHasher.hpp>
#include <utils/PresenceModel.hpp>

namespace ChimeraClassify {

struct PresenceSummary;
struct PresenceDecision;

namespace {

using FeatureMethod = chimera::feature::Method;

constexpr size_t kInvalidLength = std::numeric_limits<size_t>::max();

FeatureMethod parse_feature_method_string(std::string feature)
{
  std::transform(feature.begin(), feature.end(), feature.begin(),
                 [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
  if (feature == "syncmer")
    return FeatureMethod::Syncmer;
  if (feature == "strobemer")
    return FeatureMethod::Strobemer;
  return FeatureMethod::Auto;
}

std::string feature_method_to_string(FeatureMethod method)
{
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

struct MarginDecision {
  bool accept{false};
  size_t need{0};
  double margin{0.0};
  double ratio{0.0};
};

static inline MarginDecision decide_high_conf(size_t best, size_t second,
                                              double eff_eval)
{
  MarginDecision dc;
  dc.margin = static_cast<double>(best) - static_cast<double>(second);
  dc.ratio = second ? static_cast<double>(best) / static_cast<double>(second)
                    : static_cast<double>(best);
  double base = std::max(1.0, std::sqrt(std::max(0.0, eff_eval)));
  dc.need = static_cast<size_t>(std::floor(0.35 * base + 1.0));
  if (best >= 3 && best >= second + dc.need && dc.ratio >= 1.25) {
    dc.accept = true;
  }
  return dc;
}

struct ReadStats {
  size_t count{0};
  size_t total_len{0};
  size_t min_len{kInvalidLength};
  size_t max_len{0};

  void update(size_t len) {
    if (len == 0)
      return;
    ++count;
    total_len += len;
    if (len < min_len)
      min_len = len;
    if (len > max_len)
      max_len = len;
  }
};

ReadStats sample_read_stats(const ClassifyConfig &config, size_t max_reads = 20000)
{
  ReadStats stats;
  if (!config.singleFiles.empty()) {
    for (const auto &file : config.singleFiles) {
      try {
        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin{file};
        for (auto &record : fin) {
          stats.update(record.sequence().size());
          if (stats.count >= max_reads)
            break;
        }
      } catch (const std::exception &e) {
        std::cerr << "Warning: 采样读取文件 " << file << " 失败: " << e.what() << std::endl;
      }
      if (stats.count >= max_reads)
        break;
    }
  } else if (!config.pairedFiles.empty()) {
    for (size_t i = 0; i + 1 < config.pairedFiles.size(); i += 2) {
      try {
        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin1{config.pairedFiles[i]};
        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin2{config.pairedFiles[i + 1]};
        auto it1 = fin1.begin();
        auto it2 = fin2.begin();
        for (; it1 != fin1.end() && it2 != fin2.end(); ++it1, ++it2) {
          stats.update((*it1).sequence().size());
          stats.update((*it2).sequence().size());
          if (stats.count >= max_reads)
            break;
        }
      } catch (const std::exception &e) {
        std::cerr << "Warning: 采样读取文件 " << config.pairedFiles[i]
                  << " 或 " << config.pairedFiles[i + 1]
                  << " 失败: " << e.what() << std::endl;
      }
      if (stats.count >= max_reads)
        break;
    }
  }
  if (stats.count == 0) {
    stats.min_len = 0;
    stats.max_len = 0;
  }
  return stats;
}

chimera::feature::Params prepare_feature_params_for_classify(
    const ChimeraBuild::IMCFConfig &imcfConfig,
    ClassifyConfig &config,
    FeatureMethod method,
    size_t &feature_min_len)
{
  chimera::feature::Params params{};
  if (method == FeatureMethod::Strobemer)
  {
    if (imcfConfig.strobeK == 0)
    {
      throw std::runtime_error("IMCF 数据库缺少 strobemer 参数，无法以 strobemer 模式分类。");
    }
    auto ensure_match = [](auto &cfg_value, auto reference, const char *name) {
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
    ensure_match(config.strobemer_order, imcfConfig.strobeOrder, "--strobe-order");
    ensure_match(config.strobemer_w_min, imcfConfig.strobeWmin, "--strobe-w-min");
    ensure_match(config.strobemer_w_max, imcfConfig.strobeWmax, "--strobe-w-max");

    params.method = FeatureMethod::Strobemer;
    params.strobe.k = config.strobemer_k;
    params.strobe.order = config.strobemer_order;
    params.strobe.w_min = config.strobemer_w_min;
    params.strobe.w_max = config.strobemer_w_max;
    params.strobe.seed = imcfConfig.seed64;
    params.strobe.canonical = true;
  }
  else
  {
    if (config.strobemer_k != 0 || config.strobemer_order != 0 ||
        config.strobemer_w_min != 0 || config.strobemer_w_max != 0)
    {
      throw std::runtime_error("分类参数 --strobe-* 仅在 feature=strobemer 时可用。");
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

} // namespace
// --- fast taxid dictionary ---
struct TaxDict {
  std::vector<std::vector<uint32_t>> idx2id;  // [bin][species] -> tid_id
  std::vector<std::string> id2str;            // tid_id -> taxid string
  std::vector<std::vector<uint32_t>> tid2bin; // tid_id -> 所在 bin 列表（去重）
  robin_hood::unordered_flat_map<std::string, uint32_t> str2id; // taxid -> tid_id
};

static TaxDict
build_tax_dict(const std::vector<std::vector<std::string>> &idx2tax) {
  robin_hood::unordered_flat_map<std::string, uint32_t> dict;
  dict.reserve(1ull << 20);

  TaxDict td;
  td.idx2id.resize(idx2tax.size());
  for (size_t b = 0; b < idx2tax.size(); ++b) {
    td.idx2id[b].resize(idx2tax[b].size());
    for (size_t s = 0; s < idx2tax[b].size(); ++s) {
      const std::string &t = idx2tax[b][s];
      auto it = dict.find(t);
      uint32_t id;
      if (it == dict.end()) {
        id = static_cast<uint32_t>(td.id2str.size());
        dict.emplace(t, id);
        td.id2str.push_back(t);
        td.tid2bin.emplace_back();
      } else {
        id = it->second;
      }
      td.idx2id[b][s] = id;
      td.tid2bin[id].push_back(static_cast<uint32_t>(b));
    }
  }
  for (auto &bins : td.tid2bin) {
    std::sort(bins.begin(), bins.end());
    bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
  }
  td.str2id = std::move(dict);
  return td;
}
/**
 * @brief Print the time taken for classification in a human-readable format.
 *
 * This function prints the time taken for classification in a human-readable
 * format. It calculates the seconds, minutes, and hours from the milliseconds.
 * The function outputs different formats based on the length of time.
 *
 * @param milliseconds The time taken for classification in milliseconds.
 */
void print_classify_time(long long milliseconds) {
  // Calculate seconds, minutes, and hours
  long long total_seconds = milliseconds / 1000;
  long long seconds = total_seconds % 60;
  long long total_minutes = total_seconds / 60;
  long long minutes = total_minutes % 60;
  long long hours = total_minutes / 60;

  // Output different formats based on the length of time
  if (hours > 0) {
    std::cout << hours << "h " << minutes << "min " << seconds << "s "
              << milliseconds % 1000 << "ms" << std::endl;
  } else if (minutes > 0) {
    std::cout << minutes << "min " << seconds << "s " << milliseconds % 1000
              << "ms" << std::endl;
  } else {
    std::cout << seconds << "s " << milliseconds % 1000 << "ms" << std::endl;
  }
}

/**
 * @brief Parse the reads from input files and store them in a queue.
 *
 * This function reads the sequences from input files and creates batchReads
 * objects. It processes single files in parallel or paired files in parallel
 * based on the configuration. The sequences are read in chunks and stored in
 * batchReads objects. The batchReads objects are then pushed to the global
 * queue for further processing. The number of sequences for each thread is
 * counted and accumulated to the fileInfo.sequenceNum.
 *
 * @param readQueue The queue to store the batchReads objects.
 * @param config The configuration for parsing the reads.
 * @param fileInfo The information about the files and sequences.
 */
void parseReads(moodycamel::ConcurrentQueue<batchReads> &readQueue,
                ClassifyConfig config, FileInfo &fileInfo) {
  size_t totalSequences = 0;
  size_t totalFiles = 0;

  if (!config.singleFiles.empty()) {
#pragma omp parallel
    {
      std::vector<batchReads> localBatches;
      size_t localSeq = 0;
      size_t localFiles = 0;

#pragma omp for schedule(dynamic)
      for (size_t i = 0; i < config.singleFiles.size(); ++i) {
        const auto &file = config.singleFiles[i];
        ++localFiles;

        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin{file};

        for (auto &&rec : fin | seqan3::views::chunk(config.batchSize)) {
          batchReads batch;
          for (auto &&r : rec) {
            batch.ids.emplace_back(std::move(r.id()));
            batch.seqs.emplace_back(std::move(r.sequence()));
          }
          localSeq += batch.ids.size();
          localBatches.emplace_back(std::move(batch));
        }
      }

      for (auto &batch : localBatches) {
        readQueue.enqueue(std::move(batch));
      }

#pragma omp atomic
      totalSequences += localSeq;
#pragma omp atomic
      totalFiles += localFiles;
    }
  } else if (!config.pairedFiles.empty()) {
    if (config.pairedFiles.size() % 2 != 0) {
      throw std::runtime_error("Paired input requires an even number of files");
    }

#pragma omp parallel
    {
      std::vector<batchReads> localBatches;
      size_t localSeq = 0;
      size_t localFiles = 0;

#pragma omp for schedule(dynamic)
      for (size_t i = 0; i < config.pairedFiles.size(); i += 2) {
        localFiles += 2;

        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin1{config.pairedFiles[i]};
        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin2{config.pairedFiles[i + 1]};

        auto range1 = fin1 | seqan3::views::chunk(config.batchSize);
        auto range2 = fin2 | seqan3::views::chunk(config.batchSize);

        auto it1 = range1.begin();
        auto end1 = range1.end();
        auto it2 = range2.begin();
        auto end2 = range2.end();

        for (; it1 != end1 && it2 != end2; ++it1, ++it2) {
          batchReads batch;
          for (auto &&r : *it1) {
            batch.ids.emplace_back(std::move(r.id()));
            batch.seqs.emplace_back(std::move(r.sequence()));
          }
          for (auto &&r : *it2) {
            batch.seqs2.emplace_back(std::move(r.sequence()));
          }

          localSeq += batch.ids.size();
          localBatches.emplace_back(std::move(batch));
        }
      }

      for (auto &batch : localBatches) {
        readQueue.enqueue(std::move(batch));
      }

#pragma omp atomic
      totalSequences += localSeq;
#pragma omp atomic
      totalFiles += localFiles;
    }
  } else {
    throw std::runtime_error("No input files specified");
  }

  fileInfo.sequenceNum += totalSequences;
  fileInfo.fileNum += totalFiles;

  if (config.verbose) {
    std::cout << "Number of files: " << fileInfo.fileNum << std::endl;
    std::cout << "Number of sequences: " << fileInfo.sequenceNum << std::endl
              << std::endl;
  }
}

struct IMCFIndexStatus {
  bool builtActive{false};
  long long activeMs{0};
};

/**
 * @brief 加载 IMCF 主体与索引文件。
 *
 * 该函数会按顺序读取 `.imcf` 主档案以及可选的 `.imcf.idx`
 * 索引文件；若索引缺失将自动重建并记录耗时，用于后续提示用户。
 */
IMCFIndexStatus
loadFilter(const std::string &input_file,
           chimera::imcf::InterleavedMergedCuckooFilter &imcf,
           ChimeraBuild::IMCFConfig &imcfConfig,
           std::vector<std::vector<std::string>> &indexToTaxid,
           chimera::presence::CoverageMeta *coverageMeta = nullptr) {
  namespace fs = std::filesystem;
  using Clock = std::chrono::steady_clock;

  fs::path archivePath{input_file};
  if (!fs::exists(archivePath)) {
    archivePath = fs::path{input_file + ".imcf"};
  }
  if (!fs::exists(archivePath)) {
    throw std::runtime_error("无法找到 IMCF 主档案: " + input_file);
  }

  fs::path indexBase = archivePath;
  if (indexBase.extension() == ".imcf") {
    indexBase.replace_extension("");
  }

  std::ifstream is(archivePath, std::ios::binary);
  if (!is.is_open()) {
    throw std::runtime_error("无法打开 IMCF 主档案: " + archivePath.string());
  }

  cereal::BinaryInputArchive archive(is);
  try {
    archive(imcf);
  } catch (const cereal::Exception &exc) {
    throw std::runtime_error(std::string("加载 IMCF 主档案失败: ") + exc.what());
  }
  try {
    archive(indexToTaxid);
  } catch (const cereal::Exception &exc) {
    throw std::runtime_error(std::string("加载 IMCF taxid 索引失败: ") + exc.what());
  }
  try {
    archive(imcfConfig);
  } catch (const cereal::Exception &exc) {
    throw std::runtime_error(std::string("加载 IMCF 配置失败: ") + exc.what());
  }
  if (coverageMeta) {
    try {
      archive(*coverageMeta);
    } catch (const cereal::Exception &) {
      coverageMeta->entries.clear();
      coverageMeta->version = 0;
      coverageMeta->unique_deg_threshold = 1;
    }
  } else {
    // consume optional coverage meta if present, ignore failures for backward
    // compatibility
    try {
      chimera::presence::CoverageMeta tmp;
      archive(tmp);
    } catch (const cereal::Exception &) {
    }
  }
  is.close();

  if (imcfConfig.hashVersion != ChimeraBuild::IMCFConfig::CurrentHashVersion) {
    throw std::runtime_error(
        "IMCF 数据库版本不兼容：hash_version=" +
        std::to_string(imcfConfig.hashVersion) +
        "，请使用当前 Chimera 重新构建数据库。");
  }
  if (imcfConfig.seed64 == 0) {
    throw std::runtime_error(
        "IMCF 数据库缺少 syncmer 种子信息，请重新执行 Chimera build。"
    );
  }
  if (imcfConfig.fpSalt != ChimeraBuild::IMCFConfig::DefaultFingerprintSalt) {
    throw std::runtime_error(
        "IMCF 数据库指纹盐值与当前实现不一致，检测到 fp_salt=" +
        std::to_string(imcfConfig.fpSalt) +
        "，请重新构建数据库或升级程序。");
  }
  auto timed = [](auto &&fn) {
    auto start = Clock::now();
    bool ok = fn();
    auto end = Clock::now();
    long long elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    return std::pair<bool, long long>{ok, elapsed};
  };

  IMCFIndexStatus status;

  std::string activePath = indexBase.string() + ".imcf.idx";
  auto [activeLoaded, activeMs] =
      timed([&]() { return imcf.loadActiveIndex(activePath); });
  if (!activeLoaded || !imcf.hasActiveIndex()) {
    auto rebuild = timed([&]() {
      imcf.buildActiveGroups();
      return true;
    });
    status.builtActive = true;
    status.activeMs = rebuild.second;
  } else {
    status.activeMs = activeMs;
  }

  return status;
}

/**
 * @brief Save the classification results to an output file.
 *
 * This function saves the classification results to an output file in TSV
 * format. It takes the classification results and the configuration for saving
 * the results. The function opens the output file and writes the classification
 * results to the file. Each line in the file contains the sequence ID followed
 * by the taxid and count pairs. The function closes the file after writing the
 * results.
 *
 * @param classifyResults The classification results to save.
 * @param config The configuration for saving the results.
 */
void saveResult(std::vector<classifyResult> classifyResults,
                ClassifyConfig config) {
  // Ensure the output file has a .tsv extension
  std::string outputFile = config.outputFile;
  if (std::filesystem::path(outputFile).extension() != ".tsv") {
    outputFile += ".tsv";
  }

  // Open the output file
  std::ofstream os(outputFile, std::ios::out);

  // Check if the file is successfully opened
  if (!os.is_open()) {
    throw std::runtime_error("Failed to open file: " + config.outputFile);
  }

  // Write the classification results to the output file
  for (const auto &result : classifyResults) {
    os << result.id << '\t';
    bool handled = false;
    if (!result.taxidCount.empty() &&
        result.taxidCount.front().first == "unclassified") {
      os << "unclassified";
      handled = true;
    }
    if (!handled) {
      for (const auto &[taxid, count] : result.taxidCount) {
        if (taxid == "unclassified") {
          os << taxid;
          continue;
        }
        os << taxid << ':' << count << '\t';
      }
    }
    os << '\n';
  }

  // Close the file
  os.close();
}

void postEmDecision(std::vector<classifyResult> &results,
                    const DecisionConfig &decisionConfig,
                    const std::unordered_map<std::string, double> &classWeights) {
  constexpr const char *kUnclassified = "unclassified";

  auto format_val = [](double value) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(4) << value;
    return oss.str();
  };

  auto prune_by_global_pi = [&](classifyResult &res, double pi_min) {
    if (pi_min <= 0.0 || res.posteriors.empty() || classWeights.empty()) {
      return;
    }
    std::vector<std::pair<std::string, double>> kept;
    kept.reserve(res.posteriors.size());
    double sum = 0.0;
    for (const auto &kv : res.posteriors) {
      auto weight_it = classWeights.find(kv.first);
      double w = (weight_it != classWeights.end()) ? weight_it->second : 0.0;
      if (w >= pi_min) {
        kept.push_back(kv);
        sum += kv.second;
      }
    }
    if (kept.empty()) {
      res.posteriors.clear();
      return;
    }
    if (sum > 0.0) {
      for (auto &kv : kept) {
        kv.second /= sum;
      }
    }
    res.posteriors.swap(kept);
  };

  for (auto &result : results) {
    // 先用全局权重剪枝，减少长尾“假阳性类”对后验的稀释
    if (decisionConfig.min_class_weight > 0.0) {
      prune_by_global_pi(result, decisionConfig.min_class_weight);
    }
    if (result.posteriors.empty()) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(kUnclassified, 1);
      continue;
    }

    auto posterior = result.posteriors;
    std::sort(posterior.begin(), posterior.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    const auto &top = posterior.front();
    double top_score = top.second;

    // 动态阈值：短读段必须更“自信”才放行
    double dyn_post = decisionConfig.posterior_threshold;
    double evalWeight = result.evaluated;
    if (evalWeight < 24.0) {
      dyn_post = std::max(dyn_post, 0.60);
    } else if (evalWeight < 48.0) {
      dyn_post = std::max(dyn_post, 0.56);
    } else {
      dyn_post = std::max(dyn_post, 0.52);
    }
    double class_weight = 0.0;
    bool weight_ok = true;
    if (!classWeights.empty()) {
      auto weight_it = classWeights.find(top.first);
      if (weight_it != classWeights.end()) {
        class_weight = weight_it->second;
        weight_ok = (class_weight >= decisionConfig.min_class_weight);
      }
    }
    bool pass = weight_ok && (top_score >= dyn_post);

    result.posteriors = std::move(posterior);

    if (pass) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(top.first, 0);
      continue;
    }

    result.taxidCount.clear();
    result.taxidCount.emplace_back(kUnclassified, 1);
  }
}

struct GroupHeat {
  std::vector<uint32_t> score;
  uint32_t decay_shift = 5;   // divide by 32
  uint32_t decay_period = 64; // decay every 64 sequences
  uint32_t counter = 0;

  void ensure(size_t bins) {
    if (score.size() < bins) {
      score.resize(bins, 0);
    }
  }

  void decay_if_needed() {
    if (decay_period == 0) {
      return;
    }
    ++counter;
    if (counter >= decay_period) {
      counter = 0;
      for (auto &v : score) {
        v -= (v >> decay_shift);
      }
    }
  }

  void boost(uint32_t bin, uint32_t delta) {
    if (bin >= score.size()) {
      score.resize(static_cast<size_t>(bin) + 1, 0);
    }
    uint64_t next =
        static_cast<uint64_t>(score[bin]) + static_cast<uint64_t>(delta);
    score[bin] = static_cast<uint32_t>(
        std::min<uint64_t>(next, std::numeric_limits<uint32_t>::max()));
  }
};

static inline uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  x ^= (x >> 31);
  return x;
}

struct PresenceStats {
  double score{0.0};
  double uniqueScore{0.0};
  double hits{0.0};
  uint64_t uniqueHits{0};
  std::vector<double> decoys;
};

struct PresenceAccumulator {
  size_t decoyReps{0};
  robin_hood::unordered_flat_map<uint32_t, PresenceStats> stats;

  explicit PresenceAccumulator(size_t reps = 0) : decoyReps(reps) {}

  PresenceStats &touch(uint32_t tid) {
    auto [it, inserted] = stats.try_emplace(tid);
    if (decoyReps > 0 && it->second.decoys.size() != decoyReps) {
      it->second.decoys.resize(decoyReps, 0.0);
    }
    (void)inserted;
    return it->second;
  }

  void add_target(uint32_t tid, double weight, bool uniqueEdge) {
    PresenceStats &entry = touch(tid);
    entry.score += weight;
    entry.hits += 1.0;
    if (uniqueEdge) {
      entry.uniqueHits += 1;
      entry.uniqueScore += weight;
    }
  }

  void add_decoy(size_t rep, uint32_t tid, double weight) {
    if (decoyReps == 0 || rep >= decoyReps) {
      return;
    }
    PresenceStats &entry = touch(tid);
    entry.decoys[rep] += weight;
  }
};

struct PresenceSummary {
  size_t decoyReps{0};
  robin_hood::unordered_flat_map<uint32_t, PresenceStats> stats;

  explicit PresenceSummary(size_t reps = 0) : decoyReps(reps) {}

  void merge(const PresenceAccumulator &acc) {
    for (const auto &[tid, entry] : acc.stats) {
      auto &dst = stats[tid];
      dst.score += entry.score;
      dst.uniqueScore += entry.uniqueScore;
      dst.hits += entry.hits;
      dst.uniqueHits += entry.uniqueHits;
      if (decoyReps > 0) {
        if (dst.decoys.size() != decoyReps) {
          dst.decoys.resize(decoyReps, 0.0);
        }
        size_t copy = std::min(decoyReps, entry.decoys.size());
        for (size_t i = 0; i < copy; ++i) {
          dst.decoys[i] += entry.decoys[i];
        }
      }
    }
  }
};

struct PresenceDecision {
  std::unordered_set<uint32_t> accepted;
  robin_hood::unordered_flat_map<uint32_t, double> qValues;
  robin_hood::unordered_flat_map<uint32_t, double> posteriors;
  robin_hood::unordered_flat_map<uint32_t, double> logPosteriors;
  robin_hood::unordered_flat_map<uint32_t, double> lambdaHats;
  double threshold{1.0};
  double noiseMu{0.0};
  double priorPi{0.0};
  size_t tested{0};
  size_t acceptedCount{0};
  size_t decoyPositives{0};
};

static PresenceDecision evaluate_presence_tdFDR(const PresenceSummary &summary,
                                                const TaxDict &tax,
                                                const ClassifyConfig &config);

inline void processSequence(
    const std::vector<uint64_t> &hashs1, size_t readLen,
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    ClassifyConfig &config, GroupHeat &heat,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const std::string &id,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    PresenceAccumulator *presenceAcc, uint64_t decoySeed) {
  // Calculate the number of hash values and determine the threshold for
  // classification
  size_t hashNum = hashs1.size();
  auto xor_reduce = [](const std::vector<uint64_t> &vals) {
    uint64_t acc = 0;
    for (uint64_t v : vals) {
      acc ^= v;
    }
    return acc;
  };

  const size_t binNumAll = indexToTaxid.size();
  heat.ensure(binNumAll);
  const bool presenceEnabled = (presenceAcc != nullptr);
  const size_t presenceDecoyReps = presenceEnabled ? presenceAcc->decoyReps : 0;
  const double exclusiveGamma = std::max(0.0, config.exclusive_gamma);

  size_t targetSample = hashNum / 4;
  targetSample = std::clamp<size_t>(targetSample, 16, 96);
  const size_t sampleBudget = std::min<size_t>(targetSample, hashs1.size());
  std::vector<uint64_t> sampleVals;
  if (sampleBudget > 0) {
    sampleVals.reserve(sampleBudget);
    size_t step = std::max<size_t>(1, hashs1.size() / sampleBudget);
    for (size_t i = 0; i < hashs1.size() && sampleVals.size() < sampleBudget;
         i += step) {
      sampleVals.push_back(hashs1[i]);
    }
    if (sampleVals.size() < sampleBudget && !hashs1.empty()) {
      sampleVals.push_back(hashs1.back());
    }
  }

  std::vector<std::vector<uint32_t>> sampleCount;
  std::vector<std::pair<uint32_t, uint16_t>> touchedS;
  touchedS.reserve(64);
  robin_hood::unordered_flat_map<uint32_t, uint32_t> sampleBinScore;
  uint64_t coarseTotal = 0;
  uint64_t sampleCovered = 0;
  std::vector<uint64_t> deferredEval;
  if (hashs1.size() > 64) {
    deferredEval.reserve(hashs1.size() - 64);
  }

  if (sampleBudget > 0) {
    imcf.bulkCount_sparse(sampleVals, sampleCount, &touchedS);
    sampleBinScore.reserve(touchedS.size());
    for (auto [bi, sp] : touchedS) {
      uint32_t contrib = sampleCount[bi][sp];
      if (contrib == 0) {
        continue;
      }
      coarseTotal += contrib;
      sampleBinScore[bi] += contrib;
    }
  }

  size_t sqrtBins = std::max<size_t>(
      1, static_cast<size_t>(std::sqrt(static_cast<double>(binNumAll))));
  const size_t candidateCap =
      std::min<size_t>(binNumAll, static_cast<size_t>(256));
  std::vector<uint32_t> topBins;
  bool fallback_full = (coarseTotal == 0);

  robin_hood::unordered_flat_set<uint32_t> candidateSet;
  candidateSet.reserve(256);
  robin_hood::unordered_flat_set<uint32_t> lowDegPreserve;
  lowDegPreserve.reserve(64);

  if (!sampleVals.empty()) {
    std::vector<uint32_t> routed;
    routed.reserve(16);
    for (auto v : sampleVals) {
      routed.clear();
      imcf.route(v, routed);
      if (routed.empty()) {
        continue;
      }
      std::sort(routed.begin(), routed.end());
      routed.erase(std::unique(routed.begin(), routed.end()), routed.end());
      if (routed.size() <= 2) {
        for (uint32_t b : routed) {
          if (b < binNumAll) {
            candidateSet.insert(b);
            lowDegPreserve.insert(b);
          }
        }
      }
    }
  }

  constexpr double coverageTarget = 0.92;
  if (coarseTotal > 0 && !sampleBinScore.empty()) {
    std::vector<std::pair<uint32_t, uint32_t>> ranked;
    ranked.reserve(sampleBinScore.size());
    for (const auto &kv : sampleBinScore) {
      ranked.emplace_back(kv.first, kv.second);
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    uint64_t goal = static_cast<uint64_t>(
        std::ceil(static_cast<double>(coarseTotal) * coverageTarget));
    uint64_t covered = 0;
    for (const auto &[bin, score] : ranked) {
      if (bin >= binNumAll) {
        continue;
      }
      candidateSet.insert(bin);
      covered += score;
      if (goal > 0 && covered >= goal) {
        break;
      }
    }
    sampleCovered = covered;
    fallback_full = false;
  }

  if (!fallback_full) {
    for (const auto &[bin, _] : sampleBinScore) {
      if (bin < binNumAll) {
        candidateSet.insert(bin);
      }
    }
  }

  if (!candidateSet.empty()) {
    std::vector<std::pair<uint32_t, uint64_t>> weighted;
    weighted.reserve(candidateSet.size());
    for (uint32_t bin : candidateSet) {
      uint64_t weight = 0;
      if (auto it = sampleBinScore.find(bin); it != sampleBinScore.end()) {
        weight += it->second * 4ull;
      }
      if (bin < heat.score.size()) {
        weight += heat.score[bin];
      }
      if (lowDegPreserve.find(bin) != lowDegPreserve.end()) {
        weight += (1ull << 32);
      }
      weighted.emplace_back(bin, weight);
    }
    std::sort(weighted.begin(), weighted.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    topBins.reserve(std::min(candidateCap, weighted.size()));
    for (const auto &[bin, _] : weighted) {
      topBins.push_back(bin);
      if (topBins.size() >= candidateCap) {
        break;
      }
    }
  }

  if (topBins.empty()) {
    fallback_full = true;
    topBins.resize(binNumAll);
    std::iota(topBins.begin(), topBins.end(), 0u);
  }

  if (!fallback_full) {
    std::sort(topBins.begin(), topBins.end());
    topBins.erase(std::unique(topBins.begin(), topBins.end()), topBins.end());
  }

  bool candidateEmpty = fallback_full || topBins.empty();
  if (!fallback_full) {
    for (auto bin : topBins) {
      uint32_t delta = 1;
      if (auto it = sampleBinScore.find(bin); it != sampleBinScore.end()) {
        delta = std::max<uint32_t>(delta, it->second);
      }
      heat.boost(bin, delta);
    }
  }
  heat.decay_if_needed();

  if (!fallback_full && topBins.size() != binNumAll) {
    std::sort(topBins.begin(), topBins.end());
    topBins.erase(std::unique(topBins.begin(), topBins.end()), topBins.end());
  }

  robin_hood::unordered_flat_map<uint32_t, double> tidScore;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> uniqueHits;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> consistencyHits;
  tidScore.reserve(128);
  uniqueHits.reserve(128);
  consistencyHits.reserve(128);

  robin_hood::unordered_flat_map<uint32_t, uint32_t> binHitCount;
  binHitCount.reserve(128);

  std::vector<uint32_t> minimizerTids;
  minimizerTids.reserve(16);
  std::vector<uint32_t> routedBinsBuf;
  routedBinsBuf.reserve(16);

  double eff_eval = 0.0;
  size_t n_eval = 0;

  const std::vector<uint32_t> *activeSubset =
      (fallback_full || topBins.size() == binNumAll) ? nullptr : &topBins;

  auto bucket_degree = [](size_t d) -> size_t {
    if (d <= 1) {
      return 0;
    }
    if (d <= 3) {
      return 1;
    }
    if (d <= 7) {
      return 2;
    }
    if (d <= 15) {
      return 3;
    }
    if (d <= 31) {
      return 4;
    }
    return 5;
  };

  auto evaluate_minimizer = [&](uint64_t value,
                                const std::vector<uint32_t> *subset) -> double {
    minimizerTids.clear();
    auto emit = [&](uint32_t bin, uint16_t sp) {
      if (bin >= tax.idx2id.size()) {
        return;
      }
      const auto &speciesVec = tax.idx2id[bin];
      if (sp >= speciesVec.size()) {
        return;
      }
      uint32_t tid = speciesVec[sp];
      if (tid >= tax.id2str.size()) {
        return;
      }
      minimizerTids.push_back(tid);
      ++binHitCount[bin];
    };

    if (!subset) {
      imcf.bulkContain_events(value, emit);
    } else {
      imcf.bulkContain_events_subset(value, *subset, emit);
    }

    if (minimizerTids.empty()) {
      return 0.0;
    }

    std::sort(minimizerTids.begin(), minimizerTids.end());
    minimizerTids.erase(
        std::unique(minimizerTids.begin(), minimizerTids.end()),
        minimizerTids.end());

    size_t deg = minimizerTids.size();
    if (deg == 0) {
      return 0.0;
    }
    double exclusivityWeight = 1.0;
    if (exclusiveGamma > 0.0 && deg > 0) {
      exclusivityWeight = std::pow(static_cast<double>(deg), -exclusiveGamma);
    }

    double totalBins = subset ? static_cast<double>(subset->size())
                              : static_cast<double>(binNumAll);
    if (totalBins <= 0.0) {
      totalBins = static_cast<double>(binNumAll);
    }

    routedBinsBuf.clear();
    imcf.route(value, routedBinsBuf);
    if (!routedBinsBuf.empty()) {
      std::sort(routedBinsBuf.begin(), routedBinsBuf.end());
      routedBinsBuf.erase(
          std::unique(routedBinsBuf.begin(), routedBinsBuf.end()),
          routedBinsBuf.end());
    }

    size_t df = routedBinsBuf.empty() ? deg : routedBinsBuf.size();
    double idf = std::log2((totalBins + 1.0) /
                           (static_cast<double>(df) + 1.0));
    idf = std::clamp(idf, 0.5, 5.0);

    double denom = std::log2(2.0 + static_cast<double>(deg));
    double weight = denom > 0.0 ? 1.0 / denom : 1.0;
    double contrib = idf * weight * exclusivityWeight;
    for (uint32_t tid : minimizerTids) {
      tidScore[tid] += contrib;
      ++consistencyHits[tid];
    }
    if (presenceEnabled && !minimizerTids.empty()) {
      bool uniqueEdge = (deg == 1);
      for (uint32_t tid : minimizerTids) {
        presenceAcc->add_target(tid, exclusivityWeight, uniqueEdge);
      }
      if (presenceDecoyReps > 0 && deg > 0) {
        for (size_t rep = 0; rep < presenceDecoyReps; ++rep) {
          uint64_t mix = splitmix64(value ^ (static_cast<uint64_t>(rep + 1) << 17) ^ decoySeed);
          size_t pick = static_cast<size_t>(mix % deg);
          uint32_t decoyTid = minimizerTids[pick];
          presenceAcc->add_decoy(rep, decoyTid, exclusivityWeight);
        }
      }
    }
    if (deg == 1 && !minimizerTids.empty()) {
      ++uniqueHits[minimizerTids.front()];
    }

    return contrib;
  };

  auto recompute_subset_state = [&]() {
    activeSubset = (fallback_full || topBins.size() == binNumAll)
                       ? nullptr
                       : &topBins;
  };

  robin_hood::unordered_flat_set<uint32_t> topBinSet;
  if (!fallback_full) {
    topBinSet.reserve(topBins.size());
    for (auto bin : topBins) {
      topBinSet.insert(bin);
    }
  }

  auto compute_top = [&](uint32_t &bestTid, uint32_t &secondTid, double &best,
                         double &second) {
    best = 0.0;
    second = 0.0;
    bestTid = std::numeric_limits<uint32_t>::max();
    secondTid = std::numeric_limits<uint32_t>::max();
    for (const auto &kv : tidScore) {
      double c = kv.second;
      if (c > best) {
        second = best;
        secondTid = bestTid;
        best = c;
        bestTid = kv.first;
      } else if (c > second) {
        second = c;
        secondTid = kv.first;
      }
    }
  };

  auto compute_ratio = [](double best_score, double second_score) {
    if (second_score > 0.0) {
      return best_score / std::max(second_score, std::numeric_limits<double>::min());
    }
    return best_score <= 0.0 ? 0.0
                             : std::numeric_limits<double>::infinity();
  };

  struct EvidenceStats {
    uint32_t bestTid = std::numeric_limits<uint32_t>::max();
    uint32_t secondTid = std::numeric_limits<uint32_t>::max();
    double best = 0.0;
    double second = 0.0;
    double ratio = 0.0;
    double gap = 0.0;
    size_t uniqueCount = 0;
    double uniqueRatio = 0.0;
    size_t consistency = 0;
  };

  auto collect_stats = [&]() -> EvidenceStats {
    EvidenceStats stats;
    compute_top(stats.bestTid, stats.secondTid, stats.best, stats.second);
    stats.ratio = compute_ratio(stats.best, stats.second);
    stats.gap = stats.best - stats.second;
    if (stats.bestTid != std::numeric_limits<uint32_t>::max()) {
      if (auto it = uniqueHits.find(stats.bestTid); it != uniqueHits.end()) {
        stats.uniqueCount = it->second;
      }
      if (auto it = consistencyHits.find(stats.bestTid);
          it != consistencyHits.end()) {
        stats.consistency = it->second;
      }
    }
    double denom = eff_eval > 0.0 ? eff_eval : 1.0;
    stats.uniqueRatio = static_cast<double>(stats.uniqueCount) / denom;
    return stats;
  };

  auto consistency_ratio = [&](uint32_t tid) -> double {
    if (tid == std::numeric_limits<uint32_t>::max() || n_eval == 0) {
      return 0.0;
    }
    auto it = consistencyHits.find(tid);
    if (it == consistencyHits.end()) {
      return 0.0;
    }
    return static_cast<double>(it->second) /
           static_cast<double>(std::max<size_t>(1, n_eval));
  };

  size_t n0 = std::min<size_t>(64, hashs1.size());
  for (size_t i = 0; i < n0; ++i) {
    eff_eval += evaluate_minimizer(hashs1[i], activeSubset);
  }
  n_eval = n0;

  auto meets_quick = [&](const EvidenceStats &s) {
    if (s.bestTid == std::numeric_limits<uint32_t>::max()) {
      return false;
    }
    double thr_conf_local = std::ceil(config.shotThreshold * eff_eval);
    double gap_need_local = std::max(0.5, eff_eval / 24.0);
    bool strong = (s.best >= thr_conf_local);
    bool unique_ok = (s.uniqueCount >= 3) ||
                     (s.uniqueCount >= 2 && s.uniqueRatio >= 0.12);
    bool stable = (s.gap >= gap_need_local) && (s.ratio >= 1.35);
    size_t bestRoundedLocal = static_cast<size_t>(
        std::max<double>(0.0, std::llround(s.best)));
    size_t secondRoundedLocal = static_cast<size_t>(
        std::max<double>(0.0, std::llround(s.second)));
    size_t needLocal = static_cast<size_t>(std::ceil(thr_conf_local));
    auto dc = decide_high_conf(bestRoundedLocal, secondRoundedLocal, eff_eval);
    bool margin_ok = (bestRoundedLocal >= needLocal) && dc.accept;
    return strong && unique_ok && stable && margin_ok;
  };

  EvidenceStats stats = collect_stats();
  bool highConfPre = meets_quick(stats);

  if (!highConfPre && hashs1.size() > n0) {
    size_t mask = (static_cast<size_t>(1) << 3) - 1; // sample 1/8
    for (size_t i = n0; i < hashs1.size(); ++i) {
      if ((hashs1[i] & mask) != 0) {
        deferredEval.push_back(hashs1[i]);
        continue;
      }
      eff_eval += evaluate_minimizer(hashs1[i], activeSubset);
      ++n_eval;
    }
    stats = collect_stats();
    highConfPre = meets_quick(stats);

    if (!highConfPre && !deferredEval.empty()) {
      for (auto value : deferredEval) {
        eff_eval += evaluate_minimizer(value, activeSubset);
        ++n_eval;
      }
      deferredEval.clear();
      stats = collect_stats();
      highConfPre = meets_quick(stats);
    }
  }

  double thr_conf = std::max(1.0, std::ceil(config.shotThreshold * eff_eval));
  double gap_need = std::max(0.5, eff_eval / 24.0);

  uint32_t bestTid = stats.bestTid;
  double best = stats.best;
  double second = stats.second;
  double best_ratio = stats.ratio;
  double gap = stats.gap;
  size_t uniqueCount = stats.uniqueCount;
  double uniqueRatio = stats.uniqueRatio;

  bool expanded = false;
  if (!fallback_full && bestTid < tax.tid2bin.size()) {
    const auto &shards = tax.tid2bin[bestTid];
    for (uint32_t bin : shards) {
      if (topBinSet.find(bin) == topBinSet.end()) {
        topBins.push_back(bin);
        topBinSet.insert(bin);
        expanded = true;
      }
    }
  }

  if (expanded) {
    std::sort(topBins.begin(), topBins.end());
    topBins.erase(std::unique(topBins.begin(), topBins.end()), topBins.end());
    fallback_full = (topBins.size() == binNumAll);
    recompute_subset_state();

    tidScore.clear();
    binHitCount.clear();
    uniqueHits.clear();
    consistencyHits.clear();
    eff_eval = 0.0;
    n_eval = 0;
    for (auto value : hashs1) {
      eff_eval += evaluate_minimizer(value, activeSubset);
      ++n_eval;
    }
    stats = collect_stats();
    highConfPre = meets_quick(stats);
    thr_conf = std::ceil(config.shotThreshold * eff_eval);
    gap_need = std::max(0.5, eff_eval / 24.0);
    bestTid = stats.bestTid;
    best = stats.best;
    second = stats.second;
    best_ratio = stats.ratio;
    gap = stats.gap;
    uniqueCount = stats.uniqueCount;
    uniqueRatio = stats.uniqueRatio;
  }

  size_t bestRounded = static_cast<size_t>(
      std::max<double>(0.0, std::llround(best)));
  size_t secondRounded = static_cast<size_t>(
      std::max<double>(0.0, std::llround(second)));
  size_t thrConfNeed = static_cast<size_t>(std::ceil(thr_conf));
  auto dc = decide_high_conf(bestRounded, secondRounded, eff_eval);
  bool marginAccept = (bestRounded >= thrConfNeed) && dc.accept;
  highConfPre = highConfPre && marginAccept;

  classifyResult result;
  // 告诉 EM/VEM：这条 read 的有效证据权重（考虑到 deg 与 IDF）
  result.evaluated = eff_eval;
  result.id = id;
  std::pair<std::string, std::size_t> maxCount;
  bool maxCountValid = false;

  size_t effRounded = static_cast<size_t>(
      std::max<double>(1.0, std::llround(eff_eval)));
  if (!maxCountValid && bestRounded > 0 && bestTid < tax.id2str.size()) {
    maxCount =
        std::make_pair(tax.id2str[bestTid], std::min(bestRounded, effRounded));
    maxCountValid = true;
  }

  bool use_em = config.em && !highConfPre;

  double maxEvidence = std::min(best, eff_eval);
  double beta = config.firstFilterBeta;
  if (beta <= 0.0) {
    beta = 0.8;
  }
  beta = std::clamp(beta, 0.0, 1.0);
  // 不再在 EM 路径强行压低 beta；至少保证 0.5 防止极弱证据入场
  if (use_em && beta < 0.50) {
    beta = 0.50;
  }
  size_t thr_beta =
      static_cast<size_t>(std::floor(beta * std::max(0.0, maxEvidence)));
  size_t thr_eval = static_cast<size_t>(std::ceil(
      config.shotThreshold * (config.adaptive_shot ? eff_eval
                                                   : static_cast<double>(hashNum))));
  if (thr_eval == 0) {
    thr_eval = 1;
  }
  if (use_em) {
    double base = config.adaptive_shot ? eff_eval
                                       : static_cast<double>(hashNum);
    // EM 也需要足量证据：soften 到 0.50（原 0.40）
    double softened_ratio = std::min(config.shotThreshold, 0.50);
    size_t em_eval =
        static_cast<size_t>(std::ceil(base * softened_ratio));
    if (em_eval == 0 && base > 0.0) {
      em_eval = 1;
    }
    thr_eval = std::min(thr_eval, em_eval);
  }
  size_t thr_min_eval = 0;
  if (n_eval > 0) {
    thr_min_eval = static_cast<size_t>(
        std::ceil(0.3 * static_cast<double>(n_eval)));
  }
  thr_min_eval = std::max<size_t>(thr_min_eval, 12);
  thr_min_eval = std::max<size_t>(thr_min_eval, 12);

  double TOT = 0.0;
  for (const auto &kv : tidScore) {
    TOT += kv.second;
  }

  size_t M = 0;
  if (!topBins.empty()) {
    for (auto b : topBins) {
      if (b < tax.idx2id.size()) {
        M += tax.idx2id[b].size();
      }
    }
  } else {
    for (size_t b = 0; b < tax.idx2id.size(); ++b) {
      M += tax.idx2id[b].size();
    }
  }
  M = std::max<size_t>(M, 1);

  size_t thr_final = 0;
  if (use_em) {
    size_t thr_beta_eval = std::min(thr_beta, thr_eval);
    size_t min_hits = 1;
    thr_final = std::max({thr_beta_eval, min_hits, thr_min_eval});
  } else {
    thr_final = std::max({thr_eval, thr_beta, thr_min_eval});
  }

  if (highConfPre && bestTid < tax.id2str.size()) {
    size_t bestRoundedDirect = static_cast<size_t>(
        std::max<double>(1.0, std::llround(best)));
    const std::string &taxid = tax.id2str[bestTid];
    result.taxidCount.emplace_back(taxid, bestRoundedDirect);
    maxCount = std::make_pair(taxid, bestRoundedDirect);
    maxCountValid = true;
    
  } else {
    for (const auto &[tid_id, rawScore] : tidScore) {
      size_t rounded = static_cast<size_t>(
          std::max<double>(0.0, std::llround(rawScore)));
      size_t countVal = std::min(rounded, effRounded);
      if (countVal >= thr_final) {
        const std::string &taxid = tax.id2str[tid_id];
        result.taxidCount.emplace_back(taxid, countVal);
        if (countVal == bestRounded) {
          maxCount = std::make_pair(taxid, countVal);
          maxCountValid = true;
        }
      }
    }
  }

  if (result.taxidCount.empty() && use_em && maxCountValid && maxCount.second > 0) {
    result.taxidCount.emplace_back(maxCount);
  }


  size_t baseTopK = config.preEmTopK > 0 ? static_cast<size_t>(config.preEmTopK)
                                         : 32;
  size_t dynamicTopK = baseTopK;
  double tieGapNeed = std::max(2.0, eff_eval / 24.0);
  bool nearTie = (gap < tieGapNeed) || (best_ratio < 1.30);
  bool binOverflow = (!fallback_full && topBins.size() > 40) ||
                     (result.taxidCount.size() > baseTopK);
  if (nearTie || binOverflow) {
    dynamicTopK = std::max<size_t>(64, dynamicTopK);
  } else if (best_ratio >= 2.5 &&
             best >= thr_conf + std::max(1.0, eff_eval / 16.0) &&
             dynamicTopK > 8) {
    dynamicTopK = 8;
  }

  if (!result.taxidCount.empty() && use_em) {
    size_t K = dynamicTopK;
    if (K > 0 && result.taxidCount.size() > K) {
      std::nth_element(
          result.taxidCount.begin(), result.taxidCount.begin() + K,
          result.taxidCount.end(),
          [](const auto &a, const auto &b) { return a.second > b.second; });
      result.taxidCount.resize(K);
    }
  }

  // Update classifyResult based on the number of remaining taxids.
  if (!result.taxidCount.empty()) {

    for (const auto &entry : result.taxidCount) {
      fileInfo.taxidTotalMatches[entry.first] += 1;
    }

    bool isUniqueMapping = (result.taxidCount.size() == 1);
    if (isUniqueMapping) {
      fileInfo.uniqueTaxids.insert(result.taxidCount.front().first);
      const std::string &taxid = result.taxidCount.front().first;
      fileInfo.taxidUniqueMatches[taxid] += 1;
    }

    fileInfo.classifiedNum++;
    if (result.taxidCount.size() == 1) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(maxCount);
    }
  } else {
    fileInfo.unclassifiedNum++;
    result.taxidCount.emplace_back("unclassified", 1);
  }


  // Add the classifyResult to the shared classifyResults vector
  classifyResults.emplace_back(std::move(result));
}

void processBatch(batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
                         std::vector<std::vector<std::string>> &indexToTaxid,
                         const TaxDict &tax, ClassifyConfig &config,
                         chimera::imcf::InterleavedMergedCuckooFilter &imcf,
                         std::vector<classifyResult> &classifyResults,
                         const chimera::feature::Params & feature_params,
                         size_t feature_min_len,
                         FileInfo &fileInfo,
                         GroupHeat &heat,
                         PresenceAccumulator *presenceAcc,
                         uint64_t decoySeed) {
  // Process batch of reads
  std::vector<uint64_t> hashs1;
  hashs1.reserve(2048);
  if (!batch.seqs2.empty()) {
    for (size_t i = 0; i < batch.ids.size(); ++i) {
      hashs1.clear();
      size_t len1 = (i < batch.seqs.size()) ? batch.seqs[i].size() : 0;
      size_t len2 = (i < batch.seqs2.size()) ? batch.seqs2[i].size() : 0;
      size_t readLen = len1 + len2;
      if (readLen > 0) {
        if (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength || readLen < fileInfo.minLen)
          fileInfo.minLen = readLen;
        if (readLen > fileInfo.maxLen)
          fileInfo.maxLen = readLen;
        fileInfo.bpLength += readLen;
      }
      if (i < batch.seqs.size() && batch.seqs[i].size() >= feature_min_len) {
        chimera::feature::compute_hashes_append(batch.seqs[i], feature_params, hashs1);
      }
      if (i < batch.seqs2.size() && batch.seqs2[i].size() >= feature_min_len) {
        chimera::feature::compute_hashes_append(batch.seqs2[i], feature_params, hashs1);
      }
      if (hashs1.size() > 2048) {
        std::sort(hashs1.begin(), hashs1.end());
        hashs1.erase(std::unique(hashs1.begin(), hashs1.end()), hashs1.end());
      }
      processSequence(hashs1, readLen, imcfConfig, indexToTaxid, tax, config,
                      heat, imcf, batch.ids[i], classifyResults, fileInfo,
                      presenceAcc, decoySeed);
    }
  } else {
    // Process single-end reads
    for (size_t i = 0; i < batch.seqs.size(); i++) {
      hashs1.clear();
      size_t readLen = batch.seqs[i].size();
      if (readLen > 0) {
        if (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength || readLen < fileInfo.minLen)
          fileInfo.minLen = readLen;
        if (readLen > fileInfo.maxLen)
          fileInfo.maxLen = readLen;
        fileInfo.bpLength += readLen;
      }
      if (batch.seqs[i].size() >= feature_min_len) {
        chimera::feature::compute_hashes_append(batch.seqs[i], feature_params, hashs1);
      }
      if (hashs1.size() > 2048) {
        std::sort(hashs1.begin(), hashs1.end());
        hashs1.erase(std::unique(hashs1.begin(), hashs1.end()), hashs1.end());
      }
      // Process the hash values for classification
      processSequence(hashs1, readLen, imcfConfig, indexToTaxid, tax, config,
                      heat, imcf, batch.ids[i], classifyResults, fileInfo,
                      presenceAcc, decoySeed);
    }
  }
}

/**
 * @brief Classify reads while streaming batches from producer thread.
 */
void classify_streaming(ChimeraBuild::IMCFConfig &imcfConfig,
                        moodycamel::ConcurrentQueue<batchReads> &readQueue,
                        ClassifyConfig &config,
                        chimera::imcf::InterleavedMergedCuckooFilter &imcf,
                        std::vector<std::vector<std::string>> &indexToTaxid,
                        const TaxDict &tax,
                        std::vector<classifyResult> &classifyResults,
                        FileInfo &fileInfo, std::atomic<bool> &producer_done,
                        const chimera::feature::Params &feature_params,
                        size_t feature_min_len,
                        PresenceSummary *presenceSummary,
                        uint64_t decoySeed) {

#pragma omp parallel
  {
    batchReads batch;
    std::vector<classifyResult> localClassifyResults;
    FileInfo localFileInfo;
    localFileInfo.minLen = kInvalidLength;
    localFileInfo.maxLen = 0;
    localFileInfo.bpLength = 0;
    GroupHeat heat;
    heat.ensure(indexToTaxid.size());
    PresenceAccumulator presenceLocal(presenceSummary ? presenceSummary->decoyReps : 0);
    PresenceAccumulator *presencePtr = presenceSummary ? &presenceLocal : nullptr;

    for (;;) {
      if (readQueue.try_dequeue(batch)) {
        processBatch(batch, imcfConfig, indexToTaxid, tax, config, imcf,
                     localClassifyResults, feature_params, feature_min_len,
                     localFileInfo, heat, presencePtr, decoySeed);
        continue;
      }
      if (producer_done.load(std::memory_order_acquire)) {
        break;
      }
      std::this_thread::yield();
    }

#pragma omp critical
    {
      classifyResults.insert(classifyResults.end(),
                             localClassifyResults.begin(),
                             localClassifyResults.end());

      fileInfo.classifiedNum += localFileInfo.classifiedNum;
      fileInfo.unclassifiedNum += localFileInfo.unclassifiedNum;

      fileInfo.uniqueTaxids.insert(localFileInfo.uniqueTaxids.begin(),
                                   localFileInfo.uniqueTaxids.end());

      for (const auto &[taxid, count] : localFileInfo.taxidTotalMatches) {
        fileInfo.taxidTotalMatches[taxid] += count;
      }

      for (const auto &[taxid, count] : localFileInfo.taxidUniqueMatches) {
        fileInfo.taxidUniqueMatches[taxid] += count;
      }
      size_t localMin = (localFileInfo.minLen == kInvalidLength) ? 0 : localFileInfo.minLen;
      if (localMin > 0 &&
          (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength || localMin < fileInfo.minLen)) {
        fileInfo.minLen = localMin;
      }
      if (localFileInfo.maxLen > fileInfo.maxLen) {
        fileInfo.maxLen = localFileInfo.maxLen;
      }
      fileInfo.bpLength += localFileInfo.bpLength;
      if (presenceSummary) {
        presenceSummary->merge(presenceLocal);
      }
    }
  }

}

/**
 * @brief Classifies reads using the Interleaved Merged Cuckoo Filter (IMCF) and
 * performs post-classification filtering.
 *
 * This function processes batches of reads using the IMCF, classifies them, and
 * then applies three filtering steps (the first step is performed within the
 * `processBatch` function). The classified results are collected, and
 * information about the classified and unclassified reads is updated.
 *
 * @param imcfConfig Configuration parameters for the IMCF, including k-mer size
 * and window size.
 * @param readQueue A concurrent queue holding batches of reads to be
 * classified.
 * @param config Configuration settings for the classification.
 * @param imcf A reference to the Interleaved Merged Cuckoo Filter used for
 * classification.
 * @param indexToTaxid A mapping of indices in the IMCF to taxids for resolving
 * classification results.
 * @param classifyResults A vector to store the classification results for each
 * read.
 * @param fileInfo A structure to store information about the processed reads,
 * including counts of classified and unclassified reads.
 *
 * @details
 * The function processes the input reads in parallel using OpenMP. Each thread
 * dequeues a batch of reads from `readQueue`, processes it using
 * `processBatch`, and stores the results in local variables. These local
 * results are merged into the shared `classifyResults` vector and `fileInfo`
 * structure in a critical section to ensure thread safety.
 *
 * The `processBatch` function performs the actual classification; results are
 * merged without additional post filters.
 */
void classify(ChimeraBuild::IMCFConfig &imcfConfig,
              moodycamel::ConcurrentQueue<batchReads> &readQueue,
              ClassifyConfig &config,
              chimera::imcf::InterleavedMergedCuckooFilter &imcf,
              std::vector<std::vector<std::string>> &indexToTaxid,
              const TaxDict &tax, std::vector<classifyResult> &classifyResults,
              FileInfo &fileInfo,
              const chimera::feature::Params &feature_params,
              size_t feature_min_len) {

  // Parallel processing of batches using OpenMP
#pragma omp parallel
  {
    batchReads batch;
    std::vector<classifyResult> localClassifyResults;
    FileInfo localFileInfo;
    localFileInfo.minLen = kInvalidLength;
    localFileInfo.maxLen = 0;
    localFileInfo.bpLength = 0;
    GroupHeat heat;
    heat.ensure(indexToTaxid.size());
    // Dequeue and process batches until the queue is empty
    while (readQueue.try_dequeue(batch)) {
      processBatch(batch, imcfConfig, indexToTaxid, tax, config, imcf,
                   localClassifyResults, feature_params, feature_min_len,
                   localFileInfo, heat, nullptr, 0);
    }
    // Critical section to merge local results into shared resources
#pragma omp critical
    {
      classifyResults.insert(classifyResults.end(),
                             localClassifyResults.begin(),
                             localClassifyResults.end());

      fileInfo.classifiedNum += localFileInfo.classifiedNum;
      fileInfo.unclassifiedNum += localFileInfo.unclassifiedNum;

      // Merge unique taxids from local results
      fileInfo.uniqueTaxids.insert(localFileInfo.uniqueTaxids.begin(),
                                   localFileInfo.uniqueTaxids.end());

      // Update taxid match counts
      for (const auto &[taxid, count] : localFileInfo.taxidTotalMatches) {
        fileInfo.taxidTotalMatches[taxid] += count;
      }

      for (const auto &[taxid, count] : localFileInfo.taxidUniqueMatches) {
        fileInfo.taxidUniqueMatches[taxid] += count;
      }
      size_t localMin = (localFileInfo.minLen == kInvalidLength) ? 0 : localFileInfo.minLen;
      if (localMin > 0 &&
          (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength || localMin < fileInfo.minLen)) {
        fileInfo.minLen = localMin;
      }
      if (localFileInfo.maxLen > fileInfo.maxLen) {
        fileInfo.maxLen = localFileInfo.maxLen;
      }
      fileInfo.bpLength += localFileInfo.bpLength;
    }
  }

}

struct PresenceFilterStats {
  size_t trimmedAssignments{0};
  size_t forcedUnclassified{0};
};

static PresenceFilterStats apply_presence_filter(
    const PresenceDecision &decision, const TaxDict &tax,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo) {
  PresenceFilterStats stats;
  if (decision.accepted.empty()) {
    return stats;
  }
  auto keepTaxid = [&](const std::string &taxid) {
    if (taxid == "unclassified") {
      return true;
    }
    auto it = tax.str2id.find(taxid);
    if (it == tax.str2id.end()) {
      return true;
    }
    return decision.accepted.find(it->second) != decision.accepted.end();
  };
  for (auto &result : classifyResults) {
    if (result.taxidCount.empty()) {
      continue;
    }
    auto before = result.taxidCount.size();
    result.taxidCount.erase(
        std::remove_if(result.taxidCount.begin(), result.taxidCount.end(),
                       [&](const auto &kv) { return !keepTaxid(kv.first); }),
        result.taxidCount.end());
    stats.trimmedAssignments += before - result.taxidCount.size();
    if (result.taxidCount.empty()) {
      result.taxidCount.emplace_back("unclassified", 1);
      ++stats.forcedUnclassified;
    }
    if (!result.posteriors.empty()) {
      result.posteriors.erase(
          std::remove_if(result.posteriors.begin(), result.posteriors.end(),
                         [&](const auto &kv) { return !keepTaxid(kv.first); }),
          result.posteriors.end());
    }
  }
  for (auto it = fileInfo.uniqueTaxids.begin(); it != fileInfo.uniqueTaxids.end();) {
    auto mapIt = tax.str2id.find(*it);
    if (mapIt != tax.str2id.end() &&
        decision.accepted.find(mapIt->second) == decision.accepted.end()) {
      it = fileInfo.uniqueTaxids.erase(it);
    } else {
      ++it;
    }
  }
  return stats;
}

static PresenceDecision evaluate_presence_coverage(
    const PresenceSummary &summary, const TaxDict &tax,
    const ClassifyConfig &config,
    const chimera::presence::CoverageMeta &meta) {
  PresenceDecision decision;
  decision.threshold = config.presence_tau;
  decision.priorPi = config.presence_pi;
  if (summary.stats.empty()) {
    return decision;
  }

  robin_hood::unordered_flat_map<std::string, uint64_t> uniqueMap;
  uniqueMap.reserve(meta.entries.size());
  for (const auto &entry : meta.entries) {
    uniqueMap.emplace(entry.taxid, entry.unique_signatures);
  }

  std::vector<uint64_t> uniqueCounts(tax.id2str.size(), 0);
  for (size_t i = 0; i < tax.id2str.size(); ++i) {
    auto it = uniqueMap.find(tax.id2str[i]);
    if (it != uniqueMap.end()) {
      uniqueCounts[i] = it->second;
    }
  }

  auto resolve_unique = [&](uint32_t tid) -> double {
    if (tid < uniqueCounts.size() && uniqueCounts[tid] > 0) {
      return static_cast<double>(uniqueCounts[tid]);
    }
    return 0.0;
  };

  double mu = config.presence_noise;
  if (!(mu > 0.0)) {
    double muAccum = 0.0;
    double weightSum = 0.0;
    for (const auto &[tid, stats] : summary.stats) {
      double u = resolve_unique(tid);
      if (u <= 0.0) {
        continue;
      }
      if (!stats.decoys.empty()) {
        double decoyMean =
            std::accumulate(stats.decoys.begin(), stats.decoys.end(), 0.0) /
            static_cast<double>(stats.decoys.size());
        double mu_j = decoyMean / u;
        muAccum += mu_j * u;
        weightSum += u;
      }
    }
    if (weightSum > 0.0) {
      mu = muAccum / weightSum;
    }
  }
  if (!(mu > 0.0)) {
    mu = 1e-4;
  }
  mu = std::max(mu, 1e-8);
  decision.noiseMu = mu;

  double pi = std::clamp(config.presence_pi, 1e-9, 1.0 - 1e-6);
  double logPriorOdds = std::log(pi) - std::log1p(-pi);

  decision.tested = summary.stats.size();
  for (const auto &[tid, stats] : summary.stats) {
    double u = resolve_unique(tid);
    if (u <= 0.0) {
      if (stats.uniqueHits > 0) {
        u = static_cast<double>(stats.uniqueHits);
      } else {
        u = std::max(1.0, stats.hits);
      }
    }
    double u_eff =
        std::max<double>(u, static_cast<double>(config.presence_u_min));
    double C = (stats.uniqueScore > 0.0) ? stats.uniqueScore : stats.score;
    double lambda_hat = std::max(0.0, (C / u_eff) - mu);
    double logBF = 0.0;
    if (mu > 0.0) {
      double ratio = (lambda_hat + mu) / mu;
      if (ratio > 0.0 && C > 0.0) {
        logBF = C * std::log(ratio) - lambda_hat * u_eff;
      } else {
        logBF = -lambda_hat * u_eff;
      }
    }
    double logPosterior = logBF + logPriorOdds;
    decision.logPosteriors[tid] = logPosterior;
    decision.lambdaHats[tid] = lambda_hat;
    double posteriorProb = 0.5;
    if (logPosterior >= 0.0) {
      posteriorProb = 1.0 / (1.0 + std::exp(-logPosterior));
    } else {
      double e = std::exp(logPosterior);
      posteriorProb = e / (1.0 + e);
    }
    decision.qValues[tid] = posteriorProb; // reused为兼容
    decision.posteriors[tid] = posteriorProb;
    if (logPosterior >= config.presence_tau) {
      decision.accepted.insert(tid);
    }
  }
  decision.acceptedCount = decision.accepted.size();
  return decision;
}

/**
 * @brief Run the classification process.
 *
 * This function runs the classification process using the provided
 * configuration. It prints the configuration if the verbose flag is set. The
 * function sets the number of threads for OpenMP and starts the classification
 * process. It measures the time taken for reading, classifying, and saving the
 * results. The function prints the time taken for each step if the verbose flag
 * is set.
 *
 * @param config The configuration for the classification process.
 */
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
        loadFilter(config.dbFile, imcf, imcfConfig, indexToTaxid,
                   &coverageMeta);
    auto normalize_kind = [](std::string &value) {
      std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
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
      std::cout
          << "IMCF index: active-group list missing, rebuilding in memory ("
          << rebuildActiveMs << " ms)" << std::endl;
    }
    FeatureMethod db_method = (imcfConfig.featureMethod == 1)
                                  ? FeatureMethod::Strobemer
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
      double avg_len = (stats.count == 0)
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
        std::cout << "[info] 输入 reads 建议使用 strobemer，但数据库为 syncmer，自动改用 syncmer。" << std::endl;
        final_method = db_method;
      } else {
        if (imcfConfig.strobeK == 0) {
          throw std::runtime_error("IMCF 数据库缺少 strobemer 参数，无法分类。");
        }
        if (user_method == FeatureMethod::Auto) {
          if (desired_k != 0 && (desired_k != imcfConfig.strobeK ||
                                 desired_order != imcfConfig.strobeOrder ||
                                 desired_w_min != imcfConfig.strobeWmin ||
                                 desired_w_max != imcfConfig.strobeWmax)) {
            std::cout << "[warn] 自动参数建议 strobemer(k=" << static_cast<int>(desired_k)
                      << ", order=" << static_cast<int>(desired_order)
                      << ", w=[" << desired_w_min << ',' << desired_w_max
                      << "])，但数据库为 strobemer(k=" << static_cast<int>(imcfConfig.strobeK)
                      << ", order=" << static_cast<int>(imcfConfig.strobeOrder)
                      << ", w=[" << imcfConfig.strobeWmin << ',' << imcfConfig.strobeWmax
                      << "])，将沿用数据库参数。如需匹配自动建议，请重新构建数据库。" << std::endl;
          }
        }
        config.strobemer_k = imcfConfig.strobeK;
        config.strobemer_order = imcfConfig.strobeOrder;
        config.strobemer_w_min = imcfConfig.strobeWmin;
        config.strobemer_w_max = imcfConfig.strobeWmax;
      }
    } else { // final_method syncmer
      if (db_method == FeatureMethod::Strobemer) {
        std::cout << "[warn] 自动模式建议使用 syncmer，但数据库是 strobemer，将沿用数据库参数。" << std::endl;
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

    if (final_method == FeatureMethod::Strobemer && !chimera::feature::strobemer_available()) {
      throw std::runtime_error("当前 Chimera 构建未启用 strobemer 支持，无法加载使用 strobemer 的数据库，请重新编译或改用 syncmer 数据库。");
    }

    size_t feature_min_len = 0;
    chimera::feature::Params feature_params =
        prepare_feature_params_for_classify(imcfConfig, config, final_method, feature_min_len);

    

    if (config.verbose) {
      if (final_method == FeatureMethod::Strobemer) {
        std::cout << "Feature method: strobemer (k="
                  << static_cast<int>(config.strobemer_k)
                  << ", order=" << static_cast<int>(config.strobemer_order)
                  << ", w=[" << config.strobemer_w_min << ',' << config.strobemer_w_max
                  << "], seed=" << static_cast<unsigned long long>(imcfConfig.seed64)
                  << ")" << std::endl;
      } else {
        std::cout << "Feature method: syncmer (k=" << static_cast<int>(imcfConfig.kmerSize)
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
    std::cout << "Classifying sequences by imcf (feature=" << config.feature << ")..." << std::endl;
    classify_streaming(imcfConfig, readQueue, config, imcf, indexToTaxid, tax,
                       classifyResults, fileInfo, producer_done,
                       feature_params, feature_min_len, presencePtr,
                       presenceSeed);
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
                  << readsPerSec << " reads/s" << std::defaultfloat
                  << std::endl;
      }
    }
    if (fileInfo.sequenceNum > 0) {
      fileInfo.avgLen = fileInfo.bpLength / fileInfo.sequenceNum;
    }
    if (config.verbose && fileInfo.sequenceNum > 0) {
      size_t min_print = (fileInfo.minLen == 0 || fileInfo.minLen == kInvalidLength) ? 0 : fileInfo.minLen;
      std::cout << "Read length stats: min=" << min_print
                << ", max=" << fileInfo.maxLen
                << ", avg=" << fileInfo.avgLen << std::endl;
    }
    

    PresenceDecision presenceDecision =
        evaluate_presence_coverage(presenceSummary, tax, config, coverageMeta);
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
    auto filterStats = apply_presence_filter(presenceDecision, tax,
                                             classifyResults, fileInfo);
    if (config.verbose && (filterStats.trimmedAssignments > 0 ||
                           filterStats.forcedUnclassified > 0)) {
      std::cout << "Presence filter: trimmed " << filterStats.trimmedAssignments
                << " assignments, forced " << filterStats.forcedUnclassified
                << " reads to unclassified" << std::endl;
    }

  if (config.em) {
    auto EMstart = std::chrono::high_resolution_clock::now();
    std::cout << "Running EM algorithm..." << std::endl;
    EMOptions options;
    options.temp = 1.10;
    options.prior_strength = 0.25;
    options.coexist_penalty = 0.20;
    auto [posterior, weights] = EMAlgorithm(classifyResults, config.emIter,
                                            config.emThreshold, options);
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
          << static_cast<double>(part) * 100.0 / static_cast<double>(total)
          << '%';
      return oss.str();
    };

    std::cout << "Classified sequences: " << fileInfo.classifiedNum << " ("
              << format_percentage(fileInfo.classifiedNum, fileInfo.sequenceNum)
              << ")" << std::endl;
    std::cout << "Unclassified sequences: " << fileInfo.unclassifiedNum << " ("
              << format_percentage(fileInfo.unclassifiedNum,
                                   fileInfo.sequenceNum)
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
