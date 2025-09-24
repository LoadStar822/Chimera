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
#include <filesystem>
#include <iomanip>
#include <limits>
#include <mutex>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace ChimeraClassify {

namespace {
class ProgressTracker {
public:
  ProgressTracker(bool enabled, std::string stage, size_t min_step,
                  double min_interval_seconds)
      : enabled_(enabled && min_step != 0), stage_(std::move(stage)),
        min_step_(std::max<size_t>(min_step, 1)),
        interval_(std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::duration<double>(
                std::max(0.0, min_interval_seconds)))),
        start_time_(std::chrono::steady_clock::now()),
        last_report_time_(start_time_) {}

  bool enabled() const { return enabled_; }

  void add_total(size_t value) {
    if (!enabled_) {
      return;
    }
    total_.fetch_add(value, std::memory_order_relaxed);
  }

  void mark_processed(size_t value) {
    if (!enabled_) {
      return;
    }
    size_t processed =
        processed_.fetch_add(value, std::memory_order_relaxed) + value;
    maybe_report(processed, false);
  }

  void finish() {
    if (!enabled_) {
      return;
    }
    size_t processed = processed_.load(std::memory_order_relaxed);
    size_t total = total_.load(std::memory_order_relaxed);
    if (processed == 0 && total == 0) {
      return;
    }
    maybe_report(processed, true);
  }

private:
  std::string format_rate(double reads_per_second) const {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    if (reads_per_second >= 1'000'000.0) {
      double rate = reads_per_second / 1'000'000.0;
      oss << std::setprecision(rate >= 10.0 ? 0 : 1) << rate << "M reads/s";
    } else if (reads_per_second >= 1'000.0) {
      double rate = reads_per_second / 1'000.0;
      oss << std::setprecision(rate >= 10.0 ? 0 : 1) << rate << "k reads/s";
    } else {
      oss << std::setprecision(reads_per_second >= 10.0 ? 0 : 1)
          << reads_per_second << " reads/s";
    }
    return oss.str();
  }

  std::string format_duration(double seconds) const {
    if (seconds < 0.0) {
      seconds = 0.0;
    }
    long long total_seconds = static_cast<long long>(std::round(seconds));
    long long hours = total_seconds / 3600;
    long long minutes = (total_seconds % 3600) / 60;
    long long secs = total_seconds % 60;

    std::ostringstream oss;
    if (hours > 0) {
      oss << hours << "h";
    }
    if (minutes > 0 || hours > 0) {
      if (hours > 0) {
        oss << ' ';
      }
      oss << minutes << "min";
    }
    if (hours == 0) {
      if (minutes > 0) {
        oss << ' ';
      }
      oss << secs << 's';
    }
    if (oss.str().empty()) {
      return "0s";
    }
    return oss.str();
  }

  void maybe_report(size_t processed, bool force) {
    if (!enabled_) {
      return;
    }
    size_t total = total_.load(std::memory_order_relaxed);
    auto now = std::chrono::steady_clock::now();

    std::unique_lock<std::mutex> lock(mutex_);
    if (!force) {
      bool reached_total = total > 0 && processed >= total;
      bool step_ready = processed >= last_reported_ + min_step_;
      bool interval_ready =
          interval_.count() == 0 || now - last_report_time_ >= interval_;
      if (!reached_total && !step_ready && !interval_ready) {
        return;
      }
    }

    last_report_time_ = now;
    last_reported_ = processed;

    double elapsed_ms = static_cast<double>(
        std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time_)
            .count());
    double rps = elapsed_ms > 0.0 ? processed * 1000.0 / elapsed_ms : 0.0;

    std::ostringstream oss;
    oss << '[' << stage_ << "] 已处理 " << processed;
    if (total > 0) {
      double pct = static_cast<double>(processed) /
                   static_cast<double>(std::max<size_t>(total, 1));
      pct = std::min(pct * 100.0, 100.0);
      oss << '/' << total << " (" << std::fixed << std::setprecision(1) << pct
          << "%)" << std::defaultfloat;
    } else {
      oss << " 条reads";
    }
    if (rps > 0.0) {
      oss << "，速度 " << format_rate(rps);
    }
    if (total > 0 && processed < total && rps > 0.0) {
      double remaining_seconds = static_cast<double>(total - processed) / rps;
      oss << "，剩余约 " << format_duration(remaining_seconds);
    }

    std::string message = oss.str();
    size_t padding = last_line_width_ > message.size()
                         ? last_line_width_ - message.size()
                         : 0;
    last_line_width_ = std::max(last_line_width_, message.size());

    std::cout << '\r' << message << std::string(padding, ' ') << std::flush;
    if (force || (total > 0 && processed >= total)) {
      std::cout << std::endl;
      last_line_width_ = 0;
    }
  }

  bool enabled_;
  std::string stage_;
  size_t min_step_;
  std::chrono::milliseconds interval_;
  std::atomic<size_t> total_{0};
  std::atomic<size_t> processed_{0};
  size_t last_reported_{0};
  std::chrono::steady_clock::time_point start_time_;
  std::chrono::steady_clock::time_point last_report_time_;
  size_t last_line_width_{0};
  std::mutex mutex_;
};
} // namespace
// --- fast taxid dictionary ---
struct TaxDict {
  std::vector<std::vector<uint32_t>> idx2id; // [bin][species] -> tid_id
  std::vector<std::string> id2str;           // tid_id -> taxid string
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
      } else {
        id = it->second;
      }
      td.idx2id[b][s] = id;
    }
  }
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
                ClassifyConfig config, FileInfo &fileInfo,
                ProgressTracker *progress = nullptr) {
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
          if (progress) {
            progress->add_total(batch.ids.size());
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

          if (progress) {
            progress->add_total(batch.ids.size());
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
  bool builtRouter{false};
  long long activeMs{0};
  long long routerMs{0};
};

/**
 * @brief 加载 IMCF 主体与索引文件。
 *
 * 该函数会按顺序读取 `.imcf` 主档案以及可选的 `.imcf.idx` 与 `.imcf.rtr`
 * 索引文件；若索引缺失将自动重建并记录耗时，用于后续提示用户。
 */
IMCFIndexStatus
loadFilter(const std::string &input_file,
           chimera::imcf::InterleavedMergedCuckooFilter &imcf,
           ChimeraBuild::IMCFConfig &imcfConfig,
           std::vector<std::vector<std::string>> &indexToTaxid) {
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
  archive(imcf);
  archive(indexToTaxid);
  archive(imcfConfig);
  is.close();

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

  std::string routerPath = indexBase.string() + ".imcf.rtr";
  auto [routerLoaded, routerMs] =
      timed([&]() { return imcf.loadRouterIndex(routerPath); });
  if (!routerLoaded || !imcf.hasRouterIndex()) {
    auto rebuild = timed([&]() {
      imcf.buildRouterIndex();
      return true;
    });
    status.builtRouter = true;
    status.routerMs = rebuild.second;
  } else {
    status.routerMs = routerMs;
  }

  return status;
}

/**
 * Adjust the seed value based on the kmer size.
 *
 * @param kmer_size The size of the kmer.
 * @param seed The seed value to adjust (default: 0x8F3F73B5CF1C9ADEULL).
 * @return The adjusted seed value.
 */
inline constexpr static uint64_t
adjust_seed(uint8_t const kmer_size,
            uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept {
  unsigned double_k = static_cast<unsigned>(kmer_size) * 2u;
  unsigned shift = double_k >= 64u ? 0u : (64u - double_k);
  return seed >> shift;
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
    } else if (config.output_posterior && !result.posteriors.empty()) {
      auto oldFlags = os.flags();
      auto oldPrecision = os.precision();
      os.setf(std::ios::fixed, std::ios::floatfield);
      os << std::setprecision(4) << result.posteriors.front().first << ':'
         << result.posteriors.front().second;
      if (result.posteriors.size() > 1) {
        os << '\t' << "POST_TOP2=" << result.posteriors[1].first << ':'
           << result.posteriors[1].second;
      }
      os.flags(oldFlags);
      os.precision(oldPrecision);
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
                    const std::unordered_map<std::string, double> &classWeights,
                    std::optional<std::reference_wrapper<LCA>> lca) {
  constexpr const char *kUnclassified = "unclassified";

  for (auto &result : results) {
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
    double runner_up = posterior.size() > 1 ? posterior[1].second : 0.0;

    bool pass = top_score >= decisionConfig.posterior_threshold;
    if (!std::isnan(decisionConfig.margin_ratio)) {
      double ratio = runner_up <= 0.0 ? std::numeric_limits<double>::infinity()
                                      : (top_score / runner_up);
      pass = pass && (ratio >= decisionConfig.margin_ratio);
    } else {
      pass = pass && ((top_score - runner_up) >= decisionConfig.margin_delta);
    }

    double global_weight = 0.0;
    auto weight_it = classWeights.find(top.first);
    if (weight_it != classWeights.end()) {
      global_weight = weight_it->second;
    }
    pass = pass && (global_weight >= decisionConfig.min_class_weight);

    result.posteriors = std::move(posterior);

    if (pass) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(top.first, 0);
      continue;
    }

    if (decisionConfig.use_lca_fallback && lca &&
        result.posteriors.size() > 1) {
      std::vector<std::string> candidates;
      candidates.reserve(result.posteriors.size());
      for (const auto &entry : result.posteriors) {
        candidates.emplace_back(entry.first);
      }
      std::string fallback = lca->get().getLCA(candidates);
      if (!fallback.empty() && fallback != kUnclassified) {
        result.taxidCount.clear();
        result.taxidCount.emplace_back(fallback, 0);
        continue;
      }
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

inline void processSequence(
    const std::vector<size_t> &hashs1, ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    ClassifyConfig &config, GroupHeat &heat,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const std::string &id,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo, LCA &lca,
    ProgressTracker *progress) {
  // Calculate the number of hash values and determine the threshold for
  // classification
  size_t hashNum = hashs1.size();
  auto thr_at_eval = [&](size_t n_eval_now) -> size_t {
    size_t base = config.adaptive_shot ? n_eval_now : hashNum;
    double scaled = static_cast<double>(base) * config.shotThreshold;
    size_t t = static_cast<size_t>(std::ceil(scaled));
    return std::max<size_t>(t, 1);
  };

  const size_t binNumAll = indexToTaxid.size();
  heat.ensure(binNumAll);

  size_t targetSample = hashNum / 4;
  targetSample = std::clamp<size_t>(targetSample, 16, 96);
  const size_t sampleBudget = std::min<size_t>(targetSample, hashs1.size());
  std::vector<size_t> sampleVals;
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
  size_t KmaxLimit = std::min<size_t>(128, binNumAll);
  size_t desiredK = std::clamp<size_t>(sqrtBins, sqrtBins, KmaxLimit);
  std::vector<uint32_t> topBins;
  bool fallback_full = (coarseTotal == 0);

  if (imcf.hasRouterIndex() && !sampleVals.empty()) {
    robin_hood::unordered_flat_map<uint32_t, uint32_t> freq;
    freq.reserve(256);
    std::vector<uint32_t> routed;
    for (auto v : sampleVals) {
      routed.clear();
      imcf.route(v, routed);
      for (uint32_t b : routed) {
        if (b < binNumAll) {
          ++freq[b];
        }
      }
    }
    if (!freq.empty()) {
      std::vector<std::pair<uint32_t, uint32_t>> ranked;
      ranked.reserve(freq.size());
      for (const auto &kv : freq) {
        ranked.emplace_back(kv.first, kv.second);
      }
      std::sort(ranked.begin(), ranked.end(), [](const auto &a, const auto &b) {
        return a.second > b.second;
      });
      size_t take = std::min(desiredK, ranked.size());
      topBins.reserve(take);
      for (size_t i = 0; i < take; ++i) {
        topBins.push_back(ranked[i].first);
      }
      fallback_full = false;
    }
  }

  if (topBins.empty()) {
    size_t adaptiveK = sqrtBins;
    if (coarseTotal > 0 && !sampleBinScore.empty()) {
      std::vector<std::pair<uint32_t, uint32_t>> scoreVec;
      scoreVec.reserve(sampleBinScore.size());
      for (const auto &kv : sampleBinScore) {
        scoreVec.emplace_back(kv.first, kv.second);
      }
      std::sort(
          scoreVec.begin(), scoreVec.end(),
          [](const auto &a, const auto &b) { return a.second > b.second; });
      adaptiveK = scoreVec.size();
      for (size_t i = 1; i < scoreVec.size(); ++i) {
        double ratio = static_cast<double>(scoreVec[i].second) /
                       static_cast<double>(scoreVec[i - 1].second);
        if (ratio <= 0.8) {
          adaptiveK = i;
          break;
        }
      }
      adaptiveK = std::max<size_t>(adaptiveK, sqrtBins);
      adaptiveK = std::min<size_t>(adaptiveK, KmaxLimit);
    }
    adaptiveK = std::clamp<size_t>(adaptiveK, sqrtBins, KmaxLimit);
    size_t K = adaptiveK;

    if (fallback_full) {
      topBins.resize(binNumAll);
      std::iota(topBins.begin(), topBins.end(), 0u);
    } else {
      robin_hood::unordered_flat_set<uint32_t> candidateSet;
      candidateSet.reserve(K * 2 + sampleBinScore.size());
      std::vector<uint32_t> candidateBins;
      candidateBins.reserve(K * 2 + sampleBinScore.size());

      auto push_candidate = [&](uint32_t idx) {
        if (idx >= binNumAll) {
          return;
        }
        if (candidateSet.insert(idx).second) {
          candidateBins.push_back(idx);
        }
      };

      size_t heatTake = std::min<size_t>(
          binNumAll, std::max<size_t>(K / 2, static_cast<size_t>(1)));
      using HeapNode = std::pair<uint32_t, uint32_t>;
      auto cmpNode = [](const HeapNode &a, const HeapNode &b) {
        return a.first > b.first;
      };
      std::priority_queue<HeapNode, std::vector<HeapNode>, decltype(cmpNode)>
          heatHeap(cmpNode);
      for (uint32_t idx = 0; idx < binNumAll; ++idx) {
        uint32_t sc = heat.score[idx];
        if (heatHeap.size() < heatTake) {
          heatHeap.emplace(sc, idx);
        } else if (sc > heatHeap.top().first) {
          heatHeap.pop();
          heatHeap.emplace(sc, idx);
        }
      }
      std::vector<uint32_t> fromHeat;
      fromHeat.reserve(heatHeap.size());
      while (!heatHeap.empty()) {
        fromHeat.push_back(heatHeap.top().second);
        heatHeap.pop();
      }
      std::reverse(fromHeat.begin(), fromHeat.end());
      for (auto idx : fromHeat) {
        push_candidate(idx);
      }

      for (const auto &[bin, score] : sampleBinScore) {
        (void)score;
        push_candidate(bin);
      }

      size_t desired = std::min(K, binNumAll);
      if (candidateBins.size() < desired) {
        size_t need = desired - candidateBins.size();
        if (need > 0) {
          std::priority_queue<HeapNode, std::vector<HeapNode>,
                              decltype(cmpNode)>
              fillHeap(cmpNode);
          for (uint32_t idx = 0; idx < binNumAll; ++idx) {
            if (candidateSet.find(idx) != candidateSet.end()) {
              continue;
            }
            uint32_t sc = heat.score[idx];
            if (fillHeap.size() < need) {
              fillHeap.emplace(sc, idx);
            } else if (sc > fillHeap.top().first) {
              fillHeap.pop();
              fillHeap.emplace(sc, idx);
            }
          }
          std::vector<uint32_t> extra;
          extra.reserve(fillHeap.size());
          while (!fillHeap.empty()) {
            extra.push_back(fillHeap.top().second);
            fillHeap.pop();
          }
          std::reverse(extra.begin(), extra.end());
          for (auto idx : extra) {
            push_candidate(idx);
          }
        }
      }

      desired = std::min(K, binNumAll);
      std::priority_queue<HeapNode, std::vector<HeapNode>, decltype(cmpNode)>
          finalHeap(cmpNode);
      for (auto idx : candidateBins) {
        uint32_t combined = heat.score[idx];
        if (auto it = sampleBinScore.find(idx); it != sampleBinScore.end()) {
          combined += it->second;
        }
        if (finalHeap.size() < desired) {
          finalHeap.emplace(combined, idx);
        } else if (combined > finalHeap.top().first) {
          finalHeap.pop();
          finalHeap.emplace(combined, idx);
        }
      }
      while (!finalHeap.empty()) {
        topBins.push_back(finalHeap.top().second);
        finalHeap.pop();
      }
      std::reverse(topBins.begin(), topBins.end());

      if (topBins.size() < desired) {
        robin_hood::unordered_flat_set<uint32_t> topSet(topBins.begin(),
                                                        topBins.end());
        for (uint32_t idx = 0; idx < binNumAll && topBins.size() < desired;
             ++idx) {
          if (topSet.insert(idx).second) {
            topBins.push_back(idx);
          }
        }
      }

      if (topBins.empty()) {
        topBins.resize(binNumAll);
        std::iota(topBins.begin(), topBins.end(), 0u);
        fallback_full = true;
      }
    }
  }

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

  robin_hood::unordered_flat_map<uint32_t, uint32_t> tidCount;
  tidCount.reserve(128);

  auto emit_to_tid = [&](uint32_t bin, uint16_t sp) {
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
    ++tidCount[tid];
  };

  size_t n_eval = 0;
  size_t n0 = std::min<size_t>(64, hashs1.size());
  for (size_t i = 0; i < n0; ++i) {
    if (fallback_full || topBins.size() == binNumAll) {
      imcf.bulkContain_events(hashs1[i], emit_to_tid);
    } else {
      imcf.bulkContain_events_subset(hashs1[i], topBins, emit_to_tid);
    }
  }
  n_eval = n0;

  auto top2 = [&](uint32_t &bestTid, size_t &best, size_t &second) {
    best = 0;
    second = 0;
    bestTid = std::numeric_limits<uint32_t>::max();
    for (const auto &kv : tidCount) {
      size_t c = kv.second;
      if (c > best) {
        second = best;
        best = c;
        bestTid = kv.first;
      } else if (c > second) {
        second = c;
      }
    }
  };

  [[maybe_unused]] uint32_t bestTid = std::numeric_limits<uint32_t>::max();
  size_t best = 0;
  size_t second = 0;
  top2(bestTid, best, second);

  size_t thr_conf = thr_at_eval(n_eval);
  double best_ratio = std::numeric_limits<double>::infinity();
  if (second > 0) {
    best_ratio = static_cast<double>(best) /
                 static_cast<double>(std::max<size_t>(second, size_t(1)));
  } else if (best == 0) {
    best_ratio = 0.0;
  }
  bool confident = (best >= thr_conf) && (best_ratio >= 1.8);

  if (!confident && hashs1.size() > n0) {
    size_t mask = (static_cast<size_t>(1) << 3) - 1; // sample 1/8
    for (size_t i = n0; i < hashs1.size(); ++i) {
      if ((hashs1[i] & mask) != 0) {
        continue;
      }
      if (fallback_full || topBins.size() == binNumAll) {
        imcf.bulkContain_events(hashs1[i], emit_to_tid);
      } else {
        imcf.bulkContain_events_subset(hashs1[i], topBins, emit_to_tid);
      }
      ++n_eval;
    }
    top2(bestTid, best, second);
  }

  classifyResult result;
  result.id = id;
  std::pair<std::string, std::size_t> maxCount;

  bool use_em = (config.em || config.vem);

  size_t maxBinCount = std::min(best, n_eval);
  double beta = config.firstFilterBeta;
  if (beta <= 0.0) {
    beta = 0.8;
  }
  beta = std::clamp(beta, 0.0, 1.0);
  if (use_em) {
    beta = std::min(beta, 0.45);
  }
  size_t thr_beta =
      static_cast<size_t>(std::floor(beta * static_cast<double>(maxBinCount)));
  size_t thr_eval = thr_at_eval(n_eval);
  if (use_em) {
    size_t base = config.adaptive_shot ? n_eval : hashNum;
    double softened_ratio = std::min(config.shotThreshold, 0.4);
    size_t em_eval = static_cast<size_t>(
        std::ceil(static_cast<double>(base) * softened_ratio));
    if (em_eval == 0 && base > 0) {
      em_eval = 1;
    }
    thr_eval = std::min(thr_eval, em_eval);
  }
  size_t thr_min_eval = (config.min_eval_count > 0) ? config.min_eval_count : 0;

  size_t thr_final = 0;
  if (use_em) {
    size_t adaptive_floor = static_cast<size_t>(
        std::ceil(0.25 * static_cast<double>(std::max<size_t>(n_eval, 1))));
    if (adaptive_floor == 0 && n_eval > 0) {
      adaptive_floor = 1;
    }

    size_t min_eval_gate = 0;
    if (thr_min_eval > 0 && n_eval >= thr_min_eval) {
      min_eval_gate = thr_min_eval;
    }

    size_t thr_cap = std::max(adaptive_floor, min_eval_gate);
    thr_cap = std::min(thr_cap, std::max<size_t>(n_eval, 1));
    thr_final = std::max(thr_beta, thr_cap);
  } else {
    uint64_t TOT = 0;
    for (const auto &kv : tidCount) {
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

    double mu = static_cast<double>(TOT) / static_cast<double>(M);
    double Z = (config.adaptive_fdr ? std::max(0.0, config.fdr_z) : 0.0);
    size_t thr_fdr = (Z > 0.0) ? static_cast<size_t>(std::ceil(
                                     mu + Z * std::sqrt(std::max(mu, 1e-9))))
                               : 1;

    thr_final = std::max({thr_eval, thr_fdr, thr_beta, thr_min_eval});
  }

  for (const auto &[tid_id, rawCount] : tidCount) {
    size_t countVal = std::min<size_t>(rawCount, n_eval);
    if (countVal >= thr_final) {
      const std::string &taxid = tax.id2str[tid_id];
      result.taxidCount.emplace_back(taxid, countVal);
      if (countVal == maxBinCount) {
        maxCount = std::make_pair(taxid, countVal);
      }
    }
  }

  if (!result.taxidCount.empty() && (config.em || config.vem)) {
    size_t K = config.preEmTopK > 0 ? config.preEmTopK : static_cast<size_t>(0);
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
    if (config.lca && result.taxidCount.size() > 1) {
      std::vector<std::string> taxids;
      for (auto &[taxid, count] : result.taxidCount) {
        taxids.push_back(taxid);
      }
      std::string lcaTaxid = lca.getLCA(taxids);
      result.taxidCount.clear();
      result.taxidCount.emplace_back(lcaTaxid, 0);
    } else if (result.taxidCount.size() == 1) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(maxCount);
    }
  } else {
    fileInfo.unclassifiedNum++;
    result.taxidCount.emplace_back("unclassified", 1);
  }

  // Add the classifyResult to the shared classifyResults vector
  if (progress) {
    progress->mark_processed(1);
  }
  classifyResults.emplace_back(std::move(result));
}

/*
 * @brief 从分类用的 taxonomy 文件构建 LCA 结构。
 *
 * taxonomy 文件假定每行提供 `child parent rank name` 四列，以空格分隔；
 * 这里只关心父子关系，因此只解析前两列并写入 LCA。
 */
void buildLCA(LCA &lca, const std::string &taxFile) {
  std::ifstream is(taxFile);
  if (!is.is_open()) {
    throw std::runtime_error("Failed to open file: " + taxFile);
  }

  std::string line;
  while (std::getline(is, line)) {
    std::istringstream iss(line);
    std::string childID;
    std::string parentID;
    std::string rank;
    std::string name;
    if (!(iss >> childID >> parentID >> rank >> name)) {
      continue;
    }
    lca.addEdge(parentID, childID);
  }

  is.close();
  lca.doEulerWalk("1");
}

// 第二阶段筛选：丢弃未出现在 uniqueTaxids 中的税号，防止噪声 taxid 残留。
void secondFilteringStep(std::vector<classifyResult> &classifyResults,
                         const std::unordered_set<std::string> &uniqueTaxids) {
  for (auto &result : classifyResults) {
    result.taxidCount.erase(
        std::remove_if(result.taxidCount.begin(), result.taxidCount.end(),
                       [&uniqueTaxids](const auto &pair) {
                         return uniqueTaxids.find(pair.first) ==
                                uniqueTaxids.end();
                       }),
        result.taxidCount.end());

    if (result.taxidCount.empty()) {
      result.taxidCount.emplace_back("unclassified", 1);
    }
  }
}

/*
 * 第三阶段筛选：识别唯一匹配度过低的 taxid，若大部分 reads 与其他 taxid 重合
 * 则将其替换为重合度更高的 taxid，避免劣质税号污染结果。
 */
void thirdFilteringStep(std::vector<classifyResult> &classifyResults,
                        FileInfo &fileInfo) {
  std::unordered_set<std::string> lowUniqueTaxids;
  for (const auto &[taxid, totalMatches] : fileInfo.taxidTotalMatches) {
    size_t uniqueMatches = fileInfo.taxidUniqueMatches[taxid];
    if (totalMatches == 0) {
      continue;
    }
    if (static_cast<double>(uniqueMatches) < 0.05 * totalMatches) {
      lowUniqueTaxids.insert(taxid);
    }
  }

  robin_hood::unordered_flat_map<std::string, std::unordered_set<std::string>>
      taxidToReads;
  for (const auto &result : classifyResults) {
    for (const auto &[taxid, _] : result.taxidCount) {
      taxidToReads[taxid].insert(result.id);
    }
  }

  for (const auto &taxid : lowUniqueTaxids) {
    const auto &reads = taxidToReads[taxid];
    if (reads.empty()) {
      continue;
    }

    for (const auto &[otherTaxid, otherReads] : taxidToReads) {
      if (taxid == otherTaxid) {
        continue;
      }

      size_t sharedReads = 0;
      for (const auto &readId : reads) {
        if (otherReads.find(readId) != otherReads.end()) {
          ++sharedReads;
        }
      }

      double sharedPercentage =
          static_cast<double>(sharedReads) / static_cast<double>(reads.size());
      if (sharedPercentage < 0.95) {
        continue;
      }

      for (auto &result : classifyResults) {
        if (reads.find(result.id) == reads.end()) {
          continue;
        }
        result.taxidCount.erase(std::remove_if(result.taxidCount.begin(),
                                               result.taxidCount.end(),
                                               [&taxid](const auto &pair) {
                                                 return pair.first == taxid;
                                               }),
                                result.taxidCount.end());
        result.taxidCount.emplace_back(otherTaxid, 0);
      }
      break;
    }
  }
}

/**
 * @brief Processes a batch of reads for classification using the Interleaved
 * Merged Cuckoo Filter (IMCF).
 *
 * This function handles both paired-end and single-end read processing. It
 * generates minimizer hash values from the reads, combines them as needed, and
 * passes the hash values to the `processSequence` function for classification.
 * The results are stored in the `classifyResults` vector and file statistics
 * are updated in `fileInfo`.
 *
 * @param batch A structure containing the batch of reads to be processed,
 * including read sequences and their IDs.
 * @param imcfConfig Configuration settings for the IMCF, including k-mer size
 * and window size.
 * @param indexToTaxid A mapping of indices in the IMCF to the corresponding
 * taxids for classification.
 * @param config Configuration settings for classification, including filtering
 * options.
 * @param imcf The Interleaved Merged Cuckoo Filter used for classification.
 * @param classifyResults A vector to store the classification results for each
 * read.
 * @param minimiser_view A view used to generate minimizer hash values from the
 * read sequences.
 * @param fileInfo A structure to store information about the processed reads,
 * such as the number of classified and unclassified reads.
 * @param lca A reference to an LCA (Lowest Common Ancestor) structure, used if
 * LCA classification is enabled.
 *
 * @details
 * The function first checks if paired-end reads (`batch.seqs2`) are present. If
 * so, it processes the reads in pairs, generating minimizer hash values for
 * both sequences in each pair and combining the results. For single-end reads,
 * it processes each sequence individually.
 *
 * The generated hash values are passed to the `processSequence` function for
 * classification, which updates the classification results and the file
 * information.
 *
 * - If the read length is smaller than the window size specified in
 * `imcfConfig`, the read is skipped.
 * - The function ensures that hash values are generated only for reads that
 * meet the minimum length requirement.
 */
inline void processBatch(batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
                         std::vector<std::vector<std::string>> &indexToTaxid,
                         const TaxDict &tax, ClassifyConfig &config,
                         chimera::imcf::InterleavedMergedCuckooFilter &imcf,
                         std::vector<classifyResult> &classifyResults,
                         const auto &minimiser_view, FileInfo &fileInfo,
                         LCA &lca, GroupHeat &heat,
                         ProgressTracker *progress = nullptr) {
  // Process batch of reads
  std::vector<size_t> hashs1;
  if (!batch.seqs2.empty()) {
    for (size_t i = 0; i < batch.ids.size(); ++i) {
      hashs1.clear();
      if (i < batch.seqs.size() &&
          batch.seqs[i].size() >= imcfConfig.windowSize) {
        hashs1 =
            batch.seqs[i] | minimiser_view | seqan3::ranges::to<std::vector>();
      }
      if (i < batch.seqs2.size() &&
          batch.seqs2[i].size() >= imcfConfig.windowSize) {
        std::vector<size_t> hashs2 =
            batch.seqs2[i] | minimiser_view | seqan3::ranges::to<std::vector>();
        hashs1.insert(hashs1.end(), hashs2.begin(), hashs2.end());
      }
      if (hashs1.size() > 256) {
        std::sort(hashs1.begin(), hashs1.end());
        hashs1.erase(std::unique(hashs1.begin(), hashs1.end()), hashs1.end());
      }
      processSequence(hashs1, imcfConfig, indexToTaxid, tax, config, heat, imcf,
                      batch.ids[i], classifyResults, fileInfo, lca, progress);
    }
  } else {
    // Process single-end reads
    for (size_t i = 0; i < batch.seqs.size(); i++) {
      hashs1.clear();
      if (batch.seqs[i].size() >= imcfConfig.windowSize) {
        // Generate minimizer hash values for the sequence
        hashs1 =
            batch.seqs[i] | minimiser_view | seqan3::ranges::to<std::vector>();
      }
      if (hashs1.size() > 256) {
        std::sort(hashs1.begin(), hashs1.end());
        hashs1.erase(std::unique(hashs1.begin(), hashs1.end()), hashs1.end());
      }
      // Process the hash values for classification
      processSequence(hashs1, imcfConfig, indexToTaxid, tax, config, heat, imcf,
                      batch.ids[i], classifyResults, fileInfo, lca, progress);
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
                        ProgressTracker *progress) {
  auto minimiser_view = seqan3::views::minimiser_hash(
      seqan3::shape{seqan3::ungapped{imcfConfig.kmerSize}},
      seqan3::window_size{imcfConfig.windowSize},
      seqan3::seed{adjust_seed(imcfConfig.kmerSize)});

  LCA lca;
  if (config.lca) {
    buildLCA(lca, config.taxFile);
  }

#pragma omp parallel
  {
    batchReads batch;
    std::vector<classifyResult> localClassifyResults;
    FileInfo localFileInfo;
    GroupHeat heat;
    heat.ensure(indexToTaxid.size());

    for (;;) {
      if (readQueue.try_dequeue(batch)) {
        processBatch(batch, imcfConfig, indexToTaxid, tax, config, imcf,
                     localClassifyResults, minimiser_view, localFileInfo, lca,
                     heat, progress);
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
    }
  }

  if (!(config.em || config.vem) || !config.skip_post_filter) {
    secondFilteringStep(classifyResults, fileInfo.uniqueTaxids);
    thirdFilteringStep(classifyResults, fileInfo);
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
 * @param config Configuration settings for the classification, including LCA
 * settings and taxonomic file paths.
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
 * The `processBatch` function performs the actual classification and the first
 * filtering step. After all batches are processed, the function applies the
 * `secondFilteringStep` and `thirdFilteringStep` functions to refine the
 * classification results.
 *
 * If LCA (Lowest Common Ancestor) classification is enabled in `config`, an LCA
 * structure is built using the `buildLCA` function with the provided taxonomic
 * file.
 *
 * @note
 * - The function assumes that the `processBatch`, `secondFilteringStep`, and
 * `thirdFilteringStep` functions are defined and correctly implemented.
 * - OpenMP is used for parallel processing. Ensure proper thread
 * synchronization and safety when using shared resources.
 */
void classify(ChimeraBuild::IMCFConfig &imcfConfig,
              moodycamel::ConcurrentQueue<batchReads> &readQueue,
              ClassifyConfig &config,
              chimera::imcf::InterleavedMergedCuckooFilter &imcf,
              std::vector<std::vector<std::string>> &indexToTaxid,
              const TaxDict &tax, std::vector<classifyResult> &classifyResults,
              FileInfo &fileInfo, ProgressTracker *progress) {
  // Create a minimiser hash view based on IMCF configuration
  auto minimiser_view = seqan3::views::minimiser_hash(
      seqan3::shape{seqan3::ungapped{imcfConfig.kmerSize}},
      seqan3::window_size{imcfConfig.windowSize},
      seqan3::seed{adjust_seed(imcfConfig.kmerSize)});

  // Initialize LCA structure if LCA mode is enabled in the configuration
  LCA lca;
  if (config.lca) {
    buildLCA(lca, config.taxFile);
  }

  // Parallel processing of batches using OpenMP
#pragma omp parallel
  {
    batchReads batch;
    std::vector<classifyResult> localClassifyResults;
    FileInfo localFileInfo;
    GroupHeat heat;
    heat.ensure(indexToTaxid.size());
    // Dequeue and process batches until the queue is empty
    while (readQueue.try_dequeue(batch)) {
      processBatch(batch, imcfConfig, indexToTaxid, tax, config, imcf,
                   localClassifyResults, minimiser_view, localFileInfo, lca,
                   heat, progress);
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
    }
  }

  if (!(config.em || config.vem) || !config.skip_post_filter) {
    secondFilteringStep(classifyResults, fileInfo.uniqueTaxids);
    thirdFilteringStep(classifyResults, fileInfo);
  }
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

  if (!config.em && !config.vem && !config.lca) {
    config.em = true;
  }
  if (!(config.post_ratio > 0.0)) {
    config.post_ratio = std::numeric_limits<double>::quiet_NaN();
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
  long long rebuildRouterMs = 0;

  std::string progressLabel =
      config.filter.empty() ? "classify" : config.filter;
  ProgressTracker progress(config.verbose && config.progress, progressLabel,
                           config.progressStep, config.progressInterval);
  ProgressTracker *progressPtr = progress.enabled() ? &progress : nullptr;

  if (config.filter == "imcf") {
    std::atomic<bool> producer_done{false};
    auto readEnd = readStart;
    std::thread producer([&]() {
      parseReads(readQueue, config, fileInfo, progressPtr);
      readEnd = std::chrono::high_resolution_clock::now();
      producer_done.store(true, std::memory_order_release);
    });

    std::vector<std::vector<std::string>> indexToTaxid;
    chimera::imcf::InterleavedMergedCuckooFilter imcf;
    ChimeraBuild::IMCFConfig imcfConfig;
    auto indexStatus =
        loadFilter(config.dbFile, imcf, imcfConfig, indexToTaxid);
    if (indexStatus.builtActive) {
      rebuildActiveMs = indexStatus.activeMs;
      std::cout
          << "IMCF index: active-group list missing, rebuilding in memory ("
          << rebuildActiveMs << " ms)" << std::endl;
    }
    if (indexStatus.builtRouter) {
      rebuildRouterMs = indexStatus.routerMs;
      std::cout << "IMCF index: router table missing, rebuilding in memory ("
                << rebuildRouterMs << " ms)" << std::endl;
    }
    const TaxDict tax = build_tax_dict(indexToTaxid);

    auto classifyStart = std::chrono::high_resolution_clock::now();
    std::cout << "Classifying sequences by imcf..." << std::endl;
    classify_streaming(imcfConfig, readQueue, config, imcf, indexToTaxid, tax,
                       classifyResults, fileInfo, producer_done, progressPtr);
    auto classifyEnd = std::chrono::high_resolution_clock::now();
    auto classifyDuration =
        std::chrono::duration_cast<std::chrono::milliseconds>(classifyEnd -
                                                              classifyStart);
    producer.join();
    if (progressPtr) {
      progress.finish();
    }
    auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
        readEnd - readStart);
    if (config.verbose) {
      std::cout << "Read time: ";
      print_classify_time(readDuration.count());
      std::cout << std::endl;
    }
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
  } else {
    throw std::runtime_error("Invalid filter type: " + config.filter +
                             ". 当前仅支持 imcf");
  }

  if (config.em) {
    auto EMstart = std::chrono::high_resolution_clock::now();
    std::cout << "Running EM algorithm..." << std::endl;
    EMOptions options;
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
  if (config.vem) {
    auto VEMstart = std::chrono::high_resolution_clock::now();
    std::cout << "Running VEM algorithm..." << std::endl;
    VEMOptions options;
    auto [posterior, weights] = VEMAlgorithm(classifyResults, config.emIter,
                                             config.emThreshold, options);
    classifyResults = std::move(posterior);
    classWeights = std::move(weights);
    posteriorModelUsed = true;
    auto VEMend = std::chrono::high_resolution_clock::now();
    auto VEMduration = std::chrono::duration_cast<std::chrono::milliseconds>(
        VEMend - VEMstart);
    if (config.verbose) {
      std::cout << "VEM time: ";
      print_classify_time(VEMduration.count());
    }
  }

  if (posteriorModelUsed) {
    DecisionConfig decisionConfig;
    decisionConfig.posterior_threshold = config.post_thres;
    decisionConfig.margin_delta = config.post_margin;
    decisionConfig.margin_ratio = config.post_ratio;
    decisionConfig.min_class_weight = config.post_pi_min;
    decisionConfig.use_lca_fallback = config.lca_fallback;

    std::optional<std::reference_wrapper<LCA>> fallbackLca{};
    LCA lcaInstance;
    if (decisionConfig.use_lca_fallback && !config.taxFile.empty()) {
      buildLCA(lcaInstance, config.taxFile);
      fallbackLca = std::ref(lcaInstance);
    }

    postEmDecision(classifyResults, decisionConfig, classWeights, fallbackLca);
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
    if (config.lca) {
      std::cout << "Total LCA classification: " << fileInfo.lcaNum << " ("
                << format_percentage(fileInfo.lcaNum, fileInfo.classifiedNum)
                << ")" << std::endl;
    }
    if (rebuildActiveMs > 0 || rebuildRouterMs > 0) {
      std::cout << "Index rebuild summary:" << std::endl;
      if (rebuildActiveMs > 0) {
        std::cout << "  Active index: ";
        print_classify_time(rebuildActiveMs);
      }
      if (rebuildRouterMs > 0) {
        std::cout << "  Router index: ";
        print_classify_time(rebuildRouterMs);
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
