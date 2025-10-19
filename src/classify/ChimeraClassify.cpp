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
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <span>

#include <utils/Syncmer.hpp>

namespace ChimeraClassify {

namespace dbg {
struct TraceCfg {
  bool enabled = true;
  bool once = true;
  std::string id_substr;
  std::string file = "chimera_trace.log";
  size_t topk_print = 10;

  static TraceCfg fromEnv() {
    TraceCfg cfg;
    if (const char *v = std::getenv("CHIMERA_TRACE")) {
      cfg.enabled = std::string(v) != "0";
    }
    if (const char *v = std::getenv("CHIMERA_TRACE_ONCE")) {
      cfg.once = std::string(v) != "0";
    }
    if (const char *v = std::getenv("CHIMERA_TRACE_ID")) {
      cfg.id_substr = v;
    }
    if (const char *v = std::getenv("CHIMERA_TRACE_FILE")) {
      cfg.file = v;
    }
    return cfg;
  }
};

struct Counters {
  std::atomic<uint64_t> total_reads{0};
  std::atomic<uint64_t> too_short{0};
  std::atomic<uint64_t> has_minimizers{0};
  std::atomic<uint64_t> imcf_event_zero{0};
  std::atomic<uint64_t> candidate_zero{0};
  std::atomic<uint64_t> entered_em{0};
  std::atomic<uint64_t> passed_post{0};
} g;

static TraceCfg gCfg = TraceCfg::fromEnv();
static std::atomic_flag gDeepUsed = ATOMIC_FLAG_INIT;
static std::mutex gMtx;
static std::unique_ptr<std::ofstream> gOut;

inline bool enabled() { return gCfg.enabled; }

inline void reset() {
  if (!enabled()) {
    return;
  }
  g.total_reads.store(0, std::memory_order_relaxed);
  g.too_short.store(0, std::memory_order_relaxed);
  g.has_minimizers.store(0, std::memory_order_relaxed);
  g.imcf_event_zero.store(0, std::memory_order_relaxed);
  g.candidate_zero.store(0, std::memory_order_relaxed);
  g.entered_em.store(0, std::memory_order_relaxed);
  g.passed_post.store(0, std::memory_order_relaxed);
  gDeepUsed.clear(std::memory_order_release);
}

inline void initSink() {
  if (!enabled()) {
    return;
  }
  std::scoped_lock lk(gMtx);
  if (gOut) {
    return;
  }
  auto sink = std::make_unique<std::ofstream>(gCfg.file,
                                              std::ios::out | std::ios::app);
  if (sink->is_open()) {
    gOut = std::move(sink);
  }
}

inline bool wantDeep(const std::string &id, bool isSuspicious) {
  if (!enabled()) {
    return false;
  }
  if (!gCfg.id_substr.empty() && id.find(gCfg.id_substr) == std::string::npos) {
    return false;
  }
  if (!isSuspicious) {
    return false;
  }
  if (!gCfg.once) {
    return true;
  }
  if (!gDeepUsed.test_and_set()) {
    return true;
  }
  return false;
}

template <class... Args>
inline void log(const char *fmt, Args &&...args) {
  if (!enabled()) {
    return;
  }
  if (!gOut) {
    return;
  }
  std::scoped_lock lk(gMtx);
  int size = std::snprintf(nullptr, 0, fmt, std::forward<Args>(args)...);
  if (size <= 0) {
    return;
  }
  std::vector<char> buf(static_cast<size_t>(size) + 1, '\0');
  std::snprintf(buf.data(), buf.size(), fmt, std::forward<Args>(args)...);
  (*gOut) << buf.data() << '\n';
  gOut->flush();
}

struct TraceRecord {
  size_t readLen = 0;
  size_t minimizerTotal = 0;
  size_t uniqCount = 0;
  size_t sampleCount = 0;
  uint64_t xorAll = 0;
  uint64_t xorSample = 0;
  uint64_t imcfEvents = 0;
  size_t binsWithHits = 0;
  size_t touchedBins = 0;
  size_t candidateSize = 0;
  std::vector<std::pair<uint32_t, uint32_t>> topBinHits;
  std::vector<std::pair<std::string, size_t>> emInputTop;
  std::vector<std::pair<std::string, double>> emFinalTop;
  uint64_t totalTidHits = 0;
  bool candidateEmpty = false;
  bool fallbackFull = false;
  size_t thrConf = 0;
  size_t thrEval = 0;
  size_t thrBeta = 0;
  size_t thrMinEval = 0;
  double shotThreshold = 0.0;
  double firstFilterBeta = 0.0;
  size_t preEmTopK = 0;
  size_t evaluated = 0;
  size_t bestCount = 0;
  size_t secondCount = 0;
  uint32_t bestTid = std::numeric_limits<uint32_t>::max();
  std::string bestTaxid;
  bool useEm = false;
  bool enteredEm = false;
  bool passedPost = false;
  std::string decisionReason;
  bool deepLogged = false;
  bool suspiciousPre = false;
  bool suspiciousPost = false;
  bool emFinalLogged = false;
  double evaluatedWeight = 0.0;
};

inline void dumpTrace(TraceRecord &rec, const std::string &id,
                      const char *stage) {
  if (!enabled()) {
    return;
  }
  const char *stageLabel = stage ? stage : "pre";
  const char *bestTaxid = rec.bestTaxid.empty() ? "N/A" : rec.bestTaxid.c_str();
  log("READ[%s][%s] len=%zu mm_total=%zu mm_uniq=%zu mm_sample=%zu xor_all=%016llx xor_sample=%016llx",
      id.c_str(), stageLabel, rec.readLen, rec.minimizerTotal, rec.uniqCount,
      rec.sampleCount, static_cast<unsigned long long>(rec.xorAll),
      static_cast<unsigned long long>(rec.xorSample));
  log("IMCF: events=%llu bins=%zu touched=%zu candidate=%zu fallback=%d total_tid=%llu",
      static_cast<unsigned long long>(rec.imcfEvents), rec.binsWithHits,
      rec.touchedBins, rec.candidateSize, rec.fallbackFull ? 1 : 0,
      static_cast<unsigned long long>(rec.totalTidHits));
  if (!rec.topBinHits.empty()) {
    size_t idx = 0;
    for (const auto &[bin, count] : rec.topBinHits) {
      log("  bin[%zu]=%u count=%u", idx, bin, count);
      ++idx;
    }
  }
  log("Thresholds: shot=%.3f thr_conf=%zu thr_eval=%zu thr_beta=%zu thr_min_eval=%zu preEmTopK=%zu",
      rec.shotThreshold, rec.thrConf, rec.thrEval, rec.thrBeta, rec.thrMinEval,
      rec.preEmTopK);
  log("Counts: evaluated=%zu weight=%.2f best=%zu second=%zu best_tid=%s",
      rec.evaluated, rec.evaluatedWeight, rec.bestCount, rec.secondCount,
      bestTaxid);
  if (!rec.emInputTop.empty()) {
    size_t idx = 0;
    for (const auto &[taxid, cnt] : rec.emInputTop) {
      log("  EM-in[%zu]=%s:%zu", idx, taxid.c_str(), cnt);
      ++idx;
    }
  }
  if (!rec.emFinalTop.empty()) {
    size_t idx = 0;
    for (const auto &[taxid, val] : rec.emFinalTop) {
      log("  EM-final[%zu]=%s:%.6f", idx, taxid.c_str(), val);
      ++idx;
    }
  }
  log("Decision: passed=%d reason=%s", rec.passedPost ? 1 : 0,
      rec.decisionReason.empty() ? "" : rec.decisionReason.c_str());
  rec.emFinalLogged = !rec.emFinalTop.empty();
}
} // namespace dbg

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
  std::vector<std::vector<uint32_t>> idx2id;  // [bin][species] -> tid_id
  std::vector<std::string> id2str;            // tid_id -> taxid string
  std::vector<std::vector<uint32_t>> tid2bin; // tid_id -> 所在 bin 列表（去重）
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
          readQueue.enqueue(std::move(batch));
        }
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
          readQueue.enqueue(std::move(batch));
        }
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
           ClassifyConfig &expectedConfig) {
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

  uint32_t layoutMagic = 0;
  is.read(reinterpret_cast<char *>(&layoutMagic), sizeof(layoutMagic));
  if (!is.good()) {
    throw std::runtime_error("无法读取 IMCF 布局头部，请重新构建数据库。");
  }
  if (layoutMagic != chimera::imcf::LayoutHeaderV2::Magic) {
    throw std::runtime_error(
        "IMCF 数据库缺少 v2 布局头部，无法加载，请使用当前 Chimera 重新构建数据库。");
  }

  cereal::BinaryInputArchive archive(is);
  archive(imcf);
  archive(indexToTaxid);
  archive(imcfConfig);
  is.close();

  const auto &layout = imcf.layout();
  if (layout.layoutMajor != chimera::imcf::LayoutHeaderV2::CurrentMajor) {
    throw std::runtime_error(
        "IMCF 布局版本不兼容：layout_major=" +
        std::to_string(layout.layoutMajor) +
        "，请重新执行 Chimera build。");
  }

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
  auto to_lower = [](std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return value;
  };
  std::string storedKind = to_lower(imcfConfig.taxonomyKind);
  if (storedKind.empty()) {
    storedKind = "ncbi";
  }
  std::string expectedKind = to_lower(expectedConfig.taxonomyKind);
  const bool kindAuto = expectedKind.empty() || expectedKind == "auto";
  if (kindAuto) {
    expectedConfig.taxonomyKind = storedKind;
  } else if (!storedKind.empty() && storedKind != expectedKind) {
    throw std::runtime_error(
        "IMCF taxonomy_kind 不匹配，数据库记录为 '" + imcfConfig.taxonomyKind +
        "'，当前运行期望 '" + expectedConfig.taxonomyKind +
        "'。请确认使用的数据库与分类参数一致。");
  }

  std::string storedVersion = imcfConfig.taxonomyVersion;
  if (storedVersion.empty()) {
    storedVersion = (storedKind == "gtdb") ? "gtdb-unknown" : "ncbi-taxdump";
    imcfConfig.taxonomyVersion = storedVersion;
  }
  std::string expectedVersion = expectedConfig.taxonomyVersion;
  const bool versionAuto = expectedVersion.empty() || expectedVersion == "auto";
  if (versionAuto) {
    expectedConfig.taxonomyVersion = storedVersion;
  } else if (!storedVersion.empty() && storedVersion != expectedVersion) {
    throw std::runtime_error(
        "IMCF taxonomy_version 不匹配，数据库记录为 '" + imcfConfig.taxonomyVersion +
        "'，当前运行期望 '" + expectedConfig.taxonomyVersion +
        "'。请确认使用相同版本的 taxonomy 数据。");
  }
  if (expectedConfig.verbose && (kindAuto || versionAuto)) {
    std::cout << "Using taxonomy: kind=" << expectedConfig.taxonomyKind
              << ", version=" << expectedConfig.taxonomyVersion << std::endl;
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
    auto trace = result.trace;
    auto appendReason = [&](const std::string &tag) {
      if (!trace || tag.empty()) {
        return;
      }
      if (trace->decisionReason.empty()) {
        trace->decisionReason = tag;
        return;
      }
      if (trace->decisionReason.find(tag) == std::string::npos) {
        trace->decisionReason += "|" + tag;
      }
    };

    if (result.posteriors.empty()) {
      if (trace) {
        trace->passedPost = false;
        trace->suspiciousPost = trace->enteredEm;
        appendReason("posterior_empty");
        if (trace->suspiciousPost && !trace->deepLogged &&
            dbg::wantDeep(result.id, trace->suspiciousPost)) {
          trace->deepLogged = true;
          dbg::dumpTrace(*trace, result.id, "post");
        } else if (trace->deepLogged) {
          dbg::log("Decision[%s][post]: passed=0 reason=%s", result.id.c_str(),
                   trace->decisionReason.c_str());
        }
      }
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
    double ratio = runner_up <= 0.0 ? std::numeric_limits<double>::infinity()
                                    : (top_score / runner_up);
    double delta = top_score - runner_up;

    // 动态阈值：短读段必须更“自信”才放行
    double dyn_post = decisionConfig.posterior_threshold;
    double evalWeight = 0.0;
    if (trace) {
      evalWeight = (trace->evaluatedWeight > 0.0)
                       ? trace->evaluatedWeight
                       : static_cast<double>(trace->evaluated);
    }
    if (evalWeight < 24.0) {
      dyn_post = std::max(dyn_post, 0.60);
    } else if (evalWeight < 48.0) {
      dyn_post = std::max(dyn_post, 0.56);
    } else {
      dyn_post = std::max(dyn_post, 0.52);
    }
    double dyn_ratio = std::isnan(decisionConfig.margin_ratio)
                           ? 1.30
                           : std::max(decisionConfig.margin_ratio, 1.30);
    double dyn_delta = std::max(decisionConfig.margin_delta, 0.03);

    bool threshold_ok = top_score >= dyn_post;
    bool ratio_ok = ratio >= dyn_ratio;
    bool delta_ok = delta >= dyn_delta;
    bool pass_margins = ratio_ok && delta_ok;

    double class_weight = 0.0;
    bool weight_ok = true;
    if (!classWeights.empty()) {
      auto weight_it = classWeights.find(top.first);
      if (weight_it != classWeights.end()) {
        class_weight = weight_it->second;
        weight_ok = (class_weight >= decisionConfig.min_class_weight);
      }
    }
    double weight_guard = std::max(0.7, dyn_post);
    if (top_score >= weight_guard) {
      weight_ok = true;
    }

    if (trace) {
      trace->emFinalTop.clear();
      size_t limit = std::min(posterior.size(), dbg::gCfg.topk_print);
      for (size_t i = 0; i < limit; ++i) {
        trace->emFinalTop.emplace_back(posterior[i].first, posterior[i].second);
      }
      if (trace->deepLogged && !trace->emFinalLogged &&
          !trace->emFinalTop.empty()) {
        size_t idx = 0;
        for (const auto &[taxid, val] : trace->emFinalTop) {
          dbg::log("  EM-final[%zu]=%s:%.6f", idx, taxid.c_str(), val);
          ++idx;
        }
        trace->emFinalLogged = true;
      }
    }

    std::string reason;
    if (!threshold_ok) {
      reason = "below_posterior(" + format_val(top_score) + "<" +
               format_val(dyn_post) + ")";
    } else if (!pass_margins) {
      if (!std::isnan(dyn_ratio) && !ratio_ok) {
        reason = "below_ratio(" + format_val(ratio) + "<" +
                 format_val(dyn_ratio) + ")";
      } else {
        reason = "below_margin(" + format_val(delta) + "<" +
                 format_val(dyn_delta) + ")";
      }
    } else if (!weight_ok) {
      reason = "below_weight(" + format_val(class_weight) + "<" +
               format_val(decisionConfig.min_class_weight) + ")";
    }

    bool pass = threshold_ok && pass_margins && weight_ok;

    result.posteriors = std::move(posterior);

    if (pass) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(top.first, 0);
      if (trace) {
        trace->passedPost = true;
        trace->suspiciousPost = false;
        appendReason("post_pass");
        if (trace->deepLogged) {
          dbg::log("Decision[%s][post]: passed=1 reason=%s", result.id.c_str(),
                   trace->decisionReason.c_str());
        }
      }
      dbg::g.passed_post.fetch_add(1, std::memory_order_relaxed);
      continue;
    }

    bool evidence_ok = false;
    if (decisionConfig.evidence_override && trace) {
      double eval_ref = (trace->evaluatedWeight > 0.0)
                            ? trace->evaluatedWeight
                            : static_cast<double>(trace->evaluated);
      double gap_need = std::max(1.0, eval_ref / 32.0);
      double gap = static_cast<double>(trace->bestCount) -
                   static_cast<double>(trace->secondCount);
      double ratio_pre = (trace->secondCount > 0)
                             ? static_cast<double>(trace->bestCount) /
                                   static_cast<double>(trace->secondCount)
                             : std::numeric_limits<double>::infinity();
      double min_ratio = (eval_ref < 40.0) ? 2.0 : (eval_ref < 80.0 ? 1.6 : 1.4);
      if (trace->bestCount >= trace->thrConf &&
          (ratio_pre >= min_ratio || gap >= gap_need) &&
          !trace->bestTaxid.empty()) {
        evidence_ok = true;
      }
    }

    if (decisionConfig.evidence_override && evidence_ok && trace) {
      const std::string &fallback = trace->bestTaxid;
      result.taxidCount.clear();
      result.taxidCount.emplace_back(fallback, 0);
      if (trace) {
        trace->passedPost = true;
        trace->suspiciousPost = false;
        appendReason("pre_evidence_override");
        if (trace->deepLogged) {
          dbg::log("Decision[%s][post]: passed=1 reason=%s", result.id.c_str(),
                   trace->decisionReason.c_str());
        }
      }
      dbg::g.passed_post.fetch_add(1, std::memory_order_relaxed);
      continue;
    }

    if (reason.empty()) {
      reason = "post_reject";
    }
    appendReason(reason);

    if (trace) {
      trace->passedPost = false;
      trace->suspiciousPost = trace->enteredEm;
      if (trace->suspiciousPost && !trace->deepLogged &&
          dbg::wantDeep(result.id, trace->suspiciousPost)) {
        trace->deepLogged = true;
        dbg::dumpTrace(*trace, result.id, "post");
      } else if (trace->deepLogged) {
        dbg::log("Decision[%s][post]: passed=0 reason=%s", result.id.c_str(),
                 trace->decisionReason.c_str());
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
    const std::vector<uint64_t> &hashs1, size_t readLen,
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid, const TaxDict &tax,
    ClassifyConfig &config, GroupHeat &heat,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf, const std::string &id,
    std::vector<classifyResult> &classifyResults, FileInfo &fileInfo,
    ProgressTracker *progress) {
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

  std::shared_ptr<dbg::TraceRecord> trace;
  if (dbg::enabled()) {
    trace = std::make_shared<dbg::TraceRecord>();
    trace->readLen = readLen;
    trace->minimizerTotal = hashNum;
    trace->shotThreshold = config.shotThreshold;
    trace->firstFilterBeta = config.firstFilterBeta;
    trace->preEmTopK = config.preEmTopK;
    trace->thrMinEval = config.min_eval_count;
    if (!hashs1.empty()) {
      trace->xorAll = xor_reduce(hashs1);
      std::vector<uint64_t> uniqVals = hashs1;
      std::sort(uniqVals.begin(), uniqVals.end());
      uniqVals.erase(std::unique(uniqVals.begin(), uniqVals.end()),
                     uniqVals.end());
      trace->uniqCount = uniqVals.size();
    }
  }

  const size_t binNumAll = indexToTaxid.size();
  heat.ensure(binNumAll);

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

  if (trace) {
    trace->sampleCount = sampleVals.size();
    trace->xorSample = xor_reduce(sampleVals);
    if (!touchedS.empty()) {
      std::unordered_set<uint32_t> touched;
      touched.reserve(touchedS.size());
      for (const auto &[bi, _] : touchedS) {
        touched.insert(bi);
      }
      trace->touchedBins = touched.size();
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
  if (trace) {
    trace->candidateSize = topBins.size();
    trace->fallbackFull = fallback_full;
    trace->candidateEmpty = candidateEmpty;
  }
  if (trace && candidateEmpty) {
    dbg::g.candidate_zero.fetch_add(1, std::memory_order_relaxed);
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

  robin_hood::unordered_flat_map<uint32_t, double> tidScore;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> uniqueHits;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> consistencyHits;
  tidScore.reserve(128);
  uniqueHits.reserve(128);
  consistencyHits.reserve(128);

  uint64_t eventHits = 0;
  robin_hood::unordered_flat_map<uint32_t, uint32_t> binHitCount;
  if (trace) {
    binHitCount.reserve(128);
  }

  std::vector<uint32_t> minimizerTids;
  minimizerTids.reserve(16);
  std::vector<uint32_t> routedBinsBuf;
  routedBinsBuf.reserve(16);

  double eff_eval = 0.0;
  size_t n_eval = 0;

  const std::vector<uint32_t> *activeSubset =
      (fallback_full || topBins.size() == binNumAll) ? nullptr : &topBins;

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
      if (trace) {
        ++eventHits;
        ++binHitCount[bin];
      }
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

    double contrib = idf / std::sqrt(static_cast<double>(deg));
    for (uint32_t tid : minimizerTids) {
      tidScore[tid] += contrib;
      ++consistencyHits[tid];
    }
    if (deg == 1) {
      ++uniqueHits[minimizerTids.front()];
    }

    return idf;
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

  auto compute_top = [&](uint32_t &bestTid, double &best, double &second) {
    best = 0.0;
    second = 0.0;
    bestTid = std::numeric_limits<uint32_t>::max();
    for (const auto &kv : tidScore) {
      double c = kv.second;
      if (c > best) {
        second = best;
        best = c;
        bestTid = kv.first;
      } else if (c > second) {
        second = c;
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
    compute_top(stats.bestTid, stats.best, stats.second);
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
    return strong && unique_ok && stable;
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
    if (trace) {
      binHitCount.clear();
      eventHits = 0;
    }
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

  if (trace) {
    trace->thrConf = static_cast<size_t>(std::ceil(thr_conf));
    trace->imcfEvents = eventHits;
    trace->binsWithHits = binHitCount.size();
    trace->evaluated = n_eval;
    trace->evaluatedWeight = eff_eval;
    trace->bestCount = static_cast<size_t>(std::llround(best));
    trace->secondCount = static_cast<size_t>(std::llround(second));
    trace->bestTid = bestTid;
    if (bestTid < tax.id2str.size()) {
      trace->bestTaxid = tax.id2str[bestTid];
    }
    if (highConfPre) {
      if (trace->decisionReason.empty()) {
        trace->decisionReason = "pre_high_conf";
      } else if (trace->decisionReason.find("pre_high_conf") ==
                 std::string::npos) {
        trace->decisionReason += "|pre_high_conf";
      }
      trace->suspiciousPre = false;
    }
    if (!binHitCount.empty()) {
      std::vector<std::pair<uint32_t, uint32_t>> ranked;
      ranked.reserve(binHitCount.size());
      for (const auto &kv : binHitCount) {
        ranked.emplace_back(kv.first, kv.second);
      }
      std::sort(ranked.begin(), ranked.end(),
                [](const auto &a, const auto &b) { return a.second > b.second; });
      if (ranked.size() > dbg::gCfg.topk_print) {
        ranked.resize(dbg::gCfg.topk_print);
      }
      trace->topBinHits = std::move(ranked);
    }
    bool eventZero = (eventHits == 0);
    if (!hashs1.empty() && eventZero) {
      dbg::g.imcf_event_zero.fetch_add(1, std::memory_order_relaxed);
      trace->suspiciousPre = true;
    }
    if (eventHits > 0 && trace->candidateEmpty) {
      trace->suspiciousPre = true;
    }
  }

  classifyResult result;
  // 告诉 EM/VEM：这条 read 的有效证据权重（考虑到 deg 与 IDF）
  result.evaluated = eff_eval;
  result.id = id;
  result.trace = trace;
  std::pair<std::string, std::size_t> maxCount;
  bool maxCountValid = false;

  size_t bestRounded = static_cast<size_t>(std::max<double>(0.0, std::llround(best)));
  size_t effRounded = static_cast<size_t>(
      std::max<double>(1.0, std::llround(eff_eval)));
  if (!maxCountValid && bestRounded > 0 && bestTid < tax.id2str.size()) {
    maxCount =
        std::make_pair(tax.id2str[bestTid], std::min(bestRounded, effRounded));
    maxCountValid = true;
  }

  bool use_em = (config.em || config.vem) && !highConfPre;

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
  if (config.min_eval_count > 0) {
    thr_min_eval = std::max(thr_min_eval, config.min_eval_count);
  }

  if (trace) {
    trace->thrBeta = thr_beta;
    trace->thrEval = thr_eval;
    trace->thrMinEval = thr_min_eval;
    trace->useEm = use_em;
  }

  double TOT = 0.0;
  for (const auto &kv : tidScore) {
    TOT += kv.second;
  }
  if (trace) {
    trace->totalTidHits = static_cast<uint64_t>(std::llround(TOT));
  }

  if (trace && !tidScore.empty()) {
    std::vector<std::pair<std::string, size_t>> ranked;
    ranked.reserve(tidScore.size());
    for (const auto &[tid_id, raw] : tidScore) {
      if (tid_id < tax.id2str.size()) {
        ranked.emplace_back(
            tax.id2str[tid_id],
            static_cast<size_t>(std::max<double>(0.0, std::llround(raw))));
      }
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    if (ranked.size() > dbg::gCfg.topk_print) {
      ranked.resize(dbg::gCfg.topk_print);
    }
    trace->emInputTop = ranked;
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

  auto compute_thr_fdr = [&](double Z_value) -> size_t {
    if (!config.adaptive_fdr || Z_value <= 0.0) {
      return 1;
    }
    double mu = TOT / static_cast<double>(M);
    return static_cast<size_t>(
        std::ceil(mu + Z_value * std::sqrt(std::max(mu, 1e-9))));
  };

  size_t thr_final = 0;
  if (use_em) {
    size_t thr_fdr = compute_thr_fdr(std::max(0.0, config.fdr_z));
    size_t thr_beta_eval = std::min(thr_beta, thr_eval);
    size_t min_hits = 1;
    // EM 路径也强制最少评估量
    thr_final = std::max({thr_fdr, thr_beta_eval, min_hits, thr_min_eval});
  } else {
    size_t thr_fdr = compute_thr_fdr(std::max(0.0, config.fdr_z));
    thr_final = std::max({thr_eval, thr_fdr, thr_beta, thr_min_eval});
  }

  if (highConfPre && bestTid < tax.id2str.size()) {
    size_t bestRoundedDirect = static_cast<size_t>(
        std::max<double>(1.0, std::llround(best)));
    const std::string &taxid = tax.id2str[bestTid];
    result.taxidCount.emplace_back(taxid, bestRoundedDirect);
    maxCount = std::make_pair(taxid, bestRoundedDirect);
    maxCountValid = true;
    if (trace) {
      trace->passedPost = true;
      trace->suspiciousPost = false;
    }
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

  if (trace) {
    bool enters = use_em && !result.taxidCount.empty() &&
                  !(result.taxidCount.size() == 1 &&
                    result.taxidCount.front().first == "unclassified");
    trace->enteredEm = enters;
    if (enters) {
      dbg::g.entered_em.fetch_add(1, std::memory_order_relaxed);
    }
  }

  size_t dynamicTopK = config.preEmTopK > 0 ? static_cast<size_t>(config.preEmTopK)
                                            : 16;
  bool strong = (best_ratio >= 2.5) &&
                (best >= thr_conf + std::max(1.0, eff_eval / 16.0));
  bool weak = (best < thr_conf) || (best_ratio < 1.20);
  if (strong && dynamicTopK > 8) {
    dynamicTopK = 8;
  } else if (weak && dynamicTopK < 64) {
    dynamicTopK = 64;
  }
  if (trace) {
    trace->preEmTopK = dynamicTopK;
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

  if (trace) {
    if (trace->decisionReason.empty()) {
      if (!hashs1.empty() && trace->imcfEvents == 0) {
        trace->decisionReason = "imcf_no_hit";
      } else if (trace->candidateEmpty) {
        trace->decisionReason = "candidate_empty";
      }
    }
    if (trace->suspiciousPre && !trace->deepLogged &&
        dbg::wantDeep(id, trace->suspiciousPre)) {
      trace->deepLogged = true;
      dbg::dumpTrace(*trace, id, "pre");
    }
  }

  // Add the classifyResult to the shared classifyResults vector
  if (progress) {
    progress->mark_processed(1);
  }
  classifyResults.emplace_back(std::move(result));
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
 * @param syncmer_positions Positions used to determine syncmer offsets.
 * @param syncmer_seed Seed used when hashing syncmers.
 * @param fileInfo A structure to store information about the processed reads,
 * such as the number of classified and unclassified reads.
 *
 * @details
 * The function first checks if paired-end reads (`batch.seqs2`) are present. If
 * so, it processes the reads in pairs, generating syncmer hash values for
 * both sequences in each pair and combining the results. For single-end reads,
 * it processes each sequence individually.
 *
 * The generated hash values are passed to the `processSequence` function for
 * classification, which updates the classification results and the file
 * information.
 *
 * - If the read length is smaller than the k-mer size specified in
 * `imcfConfig`, the read is skipped.
 * - The function ensures that hash values are generated only for reads that
 * meet the minimum length requirement.
 */
inline void processBatch(batchReads batch, ChimeraBuild::IMCFConfig &imcfConfig,
                         std::vector<std::vector<std::string>> &indexToTaxid,
                         const TaxDict &tax, ClassifyConfig &config,
                         chimera::imcf::InterleavedMergedCuckooFilter &imcf,
                         std::vector<classifyResult> &classifyResults,
                         std::span<const size_t> syncmer_positions,
                         uint64_t syncmer_seed,
                         FileInfo &fileInfo,
                         GroupHeat &heat,
                         ProgressTracker *progress = nullptr) {
  // Process batch of reads
  std::vector<uint64_t> hashs1;
  if (!batch.seqs2.empty()) {
    for (size_t i = 0; i < batch.ids.size(); ++i) {
      hashs1.clear();
      size_t len1 = (i < batch.seqs.size()) ? batch.seqs[i].size() : 0;
      size_t len2 = (i < batch.seqs2.size()) ? batch.seqs2[i].size() : 0;
      size_t readLen = len1 + len2;
      if (i < batch.seqs.size() && batch.seqs[i].size() >= imcfConfig.kmerSize) {
        hashs1 = chimera::syncmer::compute_hashes(batch.seqs[i],
                                                  imcfConfig.smerSize,
                                                  imcfConfig.kmerSize,
                                                  syncmer_positions,
                                                  syncmer_seed,
                                                  true);
      }
      if (i < batch.seqs2.size() && batch.seqs2[i].size() >= imcfConfig.kmerSize) {
        std::vector<uint64_t> hashs2 = chimera::syncmer::compute_hashes(
            batch.seqs2[i],
            imcfConfig.smerSize,
            imcfConfig.kmerSize,
            syncmer_positions,
            syncmer_seed,
            true);
        hashs1.insert(hashs1.end(), hashs2.begin(), hashs2.end());
      }
      if (hashs1.size() > 2048) {
        std::sort(hashs1.begin(), hashs1.end());
        hashs1.erase(std::unique(hashs1.begin(), hashs1.end()), hashs1.end());
      }
      if (dbg::enabled()) {
        dbg::g.total_reads.fetch_add(1, std::memory_order_relaxed);
        bool allShort = (len1 < imcfConfig.kmerSize) &&
                        (len2 < imcfConfig.kmerSize);
        if (allShort) {
          dbg::g.too_short.fetch_add(1, std::memory_order_relaxed);
        }
        if (!hashs1.empty()) {
          dbg::g.has_minimizers.fetch_add(1, std::memory_order_relaxed);
        }
      }
      processSequence(hashs1, readLen, imcfConfig, indexToTaxid, tax, config,
                      heat, imcf, batch.ids[i], classifyResults, fileInfo,
                      progress);
    }
  } else {
    // Process single-end reads
    for (size_t i = 0; i < batch.seqs.size(); i++) {
      hashs1.clear();
      size_t readLen = batch.seqs[i].size();
      if (batch.seqs[i].size() >= imcfConfig.kmerSize) {
        // Generate syncmer hash values for the sequence
        hashs1 = chimera::syncmer::compute_hashes(batch.seqs[i],
                                                  imcfConfig.smerSize,
                                                  imcfConfig.kmerSize,
                                                  syncmer_positions,
                                                  syncmer_seed,
                                                  true);
      }
      if (hashs1.size() > 2048) {
        std::sort(hashs1.begin(), hashs1.end());
        hashs1.erase(std::unique(hashs1.begin(), hashs1.end()), hashs1.end());
      }
      if (dbg::enabled()) {
        dbg::g.total_reads.fetch_add(1, std::memory_order_relaxed);
        if (readLen < imcfConfig.kmerSize) {
          dbg::g.too_short.fetch_add(1, std::memory_order_relaxed);
        }
        if (!hashs1.empty()) {
          dbg::g.has_minimizers.fetch_add(1, std::memory_order_relaxed);
        }
      }
      // Process the hash values for classification
      processSequence(hashs1, readLen, imcfConfig, indexToTaxid, tax, config,
                      heat, imcf, batch.ids[i], classifyResults, fileInfo,
                      progress);
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
  if (dbg::enabled()) {
    dbg::reset();
  }
  dbg::initSink();
  if (dbg::enabled()) {
    dbg::log("Trace start: seed64=%llu, k=%u, s=%u, pos=%u, hashVersion=%u, fpSalt=%llu",
             static_cast<unsigned long long>(imcfConfig.seed64),
             static_cast<unsigned>(imcfConfig.kmerSize),
             static_cast<unsigned>(imcfConfig.smerSize),
             static_cast<unsigned>(imcfConfig.syncmerPosition),
             static_cast<unsigned>(imcfConfig.hashVersion),
             static_cast<unsigned long long>(imcfConfig.fpSalt));
  }

  const std::array<size_t, 1> syncmer_positions{ static_cast<size_t>(imcfConfig.syncmerPosition) };
  const std::span<const size_t> syncmer_pos_span(syncmer_positions);

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
                     localClassifyResults, syncmer_pos_span, imcfConfig.seed64,
                     localFileInfo, heat, progress);
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

  if (dbg::enabled()) {
    dbg::log("Summary: reads=%llu, too_short=%llu, has_min=%llu, imcf_event_zero=%llu, candidate_zero=%llu, entered_em=%llu, passed_post=%llu",
             static_cast<unsigned long long>(dbg::g.total_reads.load()),
             static_cast<unsigned long long>(dbg::g.too_short.load()),
             static_cast<unsigned long long>(dbg::g.has_minimizers.load()),
             static_cast<unsigned long long>(dbg::g.imcf_event_zero.load()),
             static_cast<unsigned long long>(dbg::g.candidate_zero.load()),
             static_cast<unsigned long long>(dbg::g.entered_em.load()),
             static_cast<unsigned long long>(dbg::g.passed_post.load()));
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
 * The `processBatch` function performs the actual classification and the first
 * filtering step. After all batches are processed, the function applies the
 * `secondFilteringStep` and `thirdFilteringStep` functions to refine the
 * classification results.
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
  if (dbg::enabled()) {
    dbg::reset();
  }
  dbg::initSink();
  if (dbg::enabled()) {
    dbg::log("Trace start: seed64=%llu, k=%u, s=%u, pos=%u, hashVersion=%u, fpSalt=%llu",
             static_cast<unsigned long long>(imcfConfig.seed64),
             static_cast<unsigned>(imcfConfig.kmerSize),
             static_cast<unsigned>(imcfConfig.smerSize),
             static_cast<unsigned>(imcfConfig.syncmerPosition),
             static_cast<unsigned>(imcfConfig.hashVersion),
             static_cast<unsigned long long>(imcfConfig.fpSalt));
  }

  const std::array<size_t, 1> syncmer_positions{ static_cast<size_t>(imcfConfig.syncmerPosition) };
  const std::span<const size_t> syncmer_pos_span(syncmer_positions);

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
                   localClassifyResults, syncmer_pos_span, imcfConfig.seed64,
                   localFileInfo, heat, progress);
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

  if (dbg::enabled()) {
    dbg::log("Summary: reads=%llu, too_short=%llu, has_min=%llu, imcf_event_zero=%llu, candidate_zero=%llu, entered_em=%llu, passed_post=%llu",
             static_cast<unsigned long long>(dbg::g.total_reads.load()),
             static_cast<unsigned long long>(dbg::g.too_short.load()),
             static_cast<unsigned long long>(dbg::g.has_minimizers.load()),
             static_cast<unsigned long long>(dbg::g.imcf_event_zero.load()),
             static_cast<unsigned long long>(dbg::g.candidate_zero.load()),
             static_cast<unsigned long long>(dbg::g.entered_em.load()),
             static_cast<unsigned long long>(dbg::g.passed_post.load()));
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

  if (!config.em && !config.vem) {
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
        loadFilter(config.dbFile, imcf, imcfConfig, indexToTaxid, config);
    if (indexStatus.builtActive) {
      rebuildActiveMs = indexStatus.activeMs;
      std::cout
          << "IMCF index: active-group list missing, rebuilding in memory ("
          << rebuildActiveMs << " ms)" << std::endl;
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
    decisionConfig.evidence_override = config.evidence_override;

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
