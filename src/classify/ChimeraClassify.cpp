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

#include <utils/Parse.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <thread>

#include <unistd.h>

namespace {

uint64_t current_rss_kb() {
  std::ifstream in("/proc/self/statm");
  uint64_t size_pages = 0;
  uint64_t resident_pages = 0;
  if (!(in >> size_pages >> resident_pages)) {
    return 0;
  }
  long page_size = sysconf(_SC_PAGESIZE);
  if (page_size <= 0) {
    return 0;
  }
  return resident_pages * static_cast<uint64_t>(page_size / 1024);
}

double kb_to_gib(uint64_t kb) { return static_cast<double>(kb) / 1024.0 / 1024.0; }

static std::vector<std::string> split_tab(const std::string &line) {
  std::vector<std::string> parts;
  std::string field;
  std::stringstream ss(line);
  while (std::getline(ss, field, '\t')) {
    parts.push_back(field);
  }
  return parts;
}

static bool try_parse_double(const std::string &s, double &out) {
  try {
    size_t idx = 0;
    out = std::stod(s, &idx);
    return idx > 0;
  } catch (...) {
    return false;
  }
}

static std::unordered_map<std::string, double>
load_weight_map_file(const std::string &path) {
  std::unordered_map<std::string, double> weights;
  std::ifstream is(path);
  if (!is.is_open()) {
    throw std::runtime_error("Failed to open weight map file: " + path);
  }

  std::string line;
  if (!std::getline(is, line)) {
    return weights;
  }

  auto header_cols = split_tab(line);
  auto find_col = [&](const std::string &name) -> int {
    for (size_t i = 0; i < header_cols.size(); ++i) {
      if (header_cols[i] == name) {
        return static_cast<int>(i);
      }
    }
    return -1;
  };

  int w_idx = find_col("number_reads");
  int id_idx = -1;
  const std::vector<std::string> id_candidates = {
      "#anonymous_contig_id", "anonymous_contig_id", "contig_id", "read_id", "id"};
  for (const auto &cand : id_candidates) {
    id_idx = find_col(cand);
    if (id_idx >= 0) {
      break;
    }
  }

  auto add_weight = [&](const std::string &id, double w) {
    if (id.empty()) {
      return;
    }
    if (!(w > 0.0)) {
      w = 1.0;
    }
    weights[id] = w;
  };

  if (w_idx >= 0 && id_idx >= 0) {
    // CAMI-style mapping.tsv with header.
    while (std::getline(is, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }
      auto cols = split_tab(line);
      if (static_cast<int>(cols.size()) <= std::max(w_idx, id_idx)) {
        continue;
      }
      double w = 1.0;
      try_parse_double(cols[w_idx], w);
      add_weight(cols[id_idx], w);
    }
    return weights;
  }

  // Fallback: plain id<TAB>weight (no header). Try to parse the first line too.
  auto parse_two_col = [&](const std::string &ln) {
    if (ln.empty() || ln[0] == '#') {
      return;
    }
    auto cols = split_tab(ln);
    if (cols.size() < 2) {
      return;
    }
    double w = 0.0;
    if (!try_parse_double(cols[1], w)) {
      return;
    }
    add_weight(cols[0], w);
  };

  parse_two_col(line);
  while (std::getline(is, line)) {
    parse_two_col(line);
  }
  return weights;
}

static inline std::string trim_copy(std::string s) {
  auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
  while (!s.empty() && is_space(static_cast<unsigned char>(s.front()))) {
    s.erase(s.begin());
  }
  while (!s.empty() && is_space(static_cast<unsigned char>(s.back()))) {
    s.pop_back();
  }
  return s;
}

static std::unique_ptr<ChimeraClassify::NcbiTaxdump>
maybe_load_ncbi_taxdump(const std::string &taxonomyKind) {
  if (taxonomyKind != "ncbi") {
    return nullptr;
  }
  const char *env_dir = std::getenv("CHIMERA_NCBI_TAXDUMP_DIR");
  std::filesystem::path base =
      env_dir && *env_dir ? env_dir
                          : "/mnt/sda/tianqinzhong/project/SimDataset/taxdump";
  std::filesystem::path nodes = base / "nodes.dmp";
  if (!std::filesystem::exists(nodes)) {
    return nullptr;
  }

  std::ifstream is(nodes);
  if (!is.is_open()) {
    return nullptr;
  }

  auto tax = std::make_unique<ChimeraClassify::NcbiTaxdump>();
  std::string line;
  while (std::getline(is, line)) {
    if (line.empty()) {
      continue;
    }
    // nodes.dmp: tax_id <tab>|<tab> parent_tax_id <tab>|<tab> rank <tab>| ...
    std::vector<std::string> fields;
    fields.reserve(8);
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, '\t')) {
      tok = trim_copy(tok);
      if (tok.empty() || tok == "|") {
        continue;
      }
      fields.push_back(tok);
      if (fields.size() >= 3) {
        // only need the first 3 meaningful fields
        break;
      }
    }
    if (fields.size() < 3) {
      continue;
    }
    uint32_t tid = 0;
    uint32_t parent = 0;
    if (!chimera::utils::try_parse_u32(fields[0], tid) ||
        !chimera::utils::try_parse_u32(fields[1], parent)) {
      continue;
    }
    const std::string &rank = fields[2];
    if (tid >= tax->parent.size()) {
      tax->parent.resize(static_cast<size_t>(tid) + 1, 0);
      tax->is_species.resize(static_cast<size_t>(tid) + 1, 0);
      tax->is_genus.resize(static_cast<size_t>(tid) + 1, 0);
    }
    tax->parent[tid] = parent;
    const bool is_sp = (rank == "species");
    tax->is_species[tid] = is_sp ? 1 : 0;
    const bool is_g = (rank == "genus");
    tax->is_genus[tid] = is_g ? 1 : 0;
  }

  if (!tax->enabled()) {
    return nullptr;
  }
  return tax;
}

} // namespace

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

  omp_set_num_threads(config.threads);

  std::vector<std::vector<std::string>> indexToTaxid;
  chimera::imcf::InterleavedMergedCuckooFilter imcf;
  ChimeraBuild::IMCFConfig imcfConfig;
  chimera::presence::CoverageMeta coverageMeta;
  loadFilter(config.dbFile, imcf, imcfConfig, indexToTaxid, &coverageMeta);
  WeightingContext weightCtx;
  std::unique_ptr<CountMinSketch> freqSketch;
  std::unique_ptr<NcbiTaxdump> ncbiTaxdump;
  if (coverageMeta.freq_model.enabled()) {
    try {
      freqSketch = std::make_unique<CountMinSketch>(
          coverageMeta.freq_model.depth, coverageMeta.freq_model.width,
          coverageMeta.freq_model.counters);
      weightCtx.freqSketch = freqSketch.get();
      weightCtx.freqStats = coverageMeta.freq_model.stats;
      weightCtx.freqQuantile = coverageMeta.freq_model.quantile;
    } catch (const std::exception &) {
    }
  }
  size_t uniq_nonzero = 0;
  if (!coverageMeta.entries.empty()) {
    for (const auto &entry : coverageMeta.entries) {
      if (entry.unique_signatures > 0) {
        ++uniq_nonzero;
      }
    }
  }

  bool freq_trusted = true;
  if (coverageMeta.freq_model.enabled() && !coverageMeta.entries.empty() &&
      uniq_nonzero == 0) {
    freq_trusted = false;
  }
  weightCtx.freq_trusted = freq_trusted;

  std::unordered_map<std::string, double> sampleWeights;
  if (!config.weight_map_file.empty()) {
    sampleWeights = load_weight_map_file(config.weight_map_file);
    if (!sampleWeights.empty()) {
      weightCtx.sampleWeights = &sampleWeights;
    }
  }
  auto normalize_kind = [](std::string &value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) {
                     return static_cast<char>(std::tolower(ch));
                   });
  };
  std::string taxonomyKind = imcfConfig.taxonomyKind;
  if (taxonomyKind.empty()) {
    taxonomyKind = "ncbi";
  }
  normalize_kind(taxonomyKind);

  // Optional NCBI strain/subspecies -> species collapse.
  // This helps avoid candidate saturation by many strain taxids, which can
  // prevent sister species from entering pre-EM/posterior lists.
  ncbiTaxdump = maybe_load_ncbi_taxdump(taxonomyKind);
  if (ncbiTaxdump && ncbiTaxdump->enabled()) {
    weightCtx.ncbiTaxdump = ncbiTaxdump.get();
  }
  FeatureMethod db_method =
      (imcfConfig.featureMethod == 1) ? FeatureMethod::Strobemer
                                      : FeatureMethod::Syncmer;
  FeatureMethod final_method = db_method;
  if (final_method == FeatureMethod::Strobemer && imcfConfig.strobeK == 0) {
    throw std::runtime_error("IMCF 数据库缺少 strobemer 参数，无法分类。");
  }

  if (final_method == FeatureMethod::Strobemer &&
      !chimera::feature::strobemer_available()) {
    throw std::runtime_error("当前 Chimera 构建未启用 strobemer 支持，无法加载使用 strobemer 的数据库，请重新编译或改用 syncmer 数据库。");
  }

  size_t feature_min_len = 0;
  chimera::feature::Params feature_params =
      prepare_feature_params_for_classify(imcfConfig, final_method,
                                          feature_min_len);

  const TaxDict tax = build_tax_dict(indexToTaxid);
  std::vector<uint32_t> tid2speciesRep;
  std::vector<uint32_t> tid2genus;
  if (weightCtx.ncbiTaxdump && weightCtx.ncbiTaxdump->enabled()) {
    tid2speciesRep.resize(tax.id2str.size());
    tid2genus.resize(tax.id2str.size(), 0u);
    for (uint32_t i = 0; i < tid2speciesRep.size(); ++i) {
      tid2speciesRep[i] = i;
    }

    robin_hood::unordered_flat_map<uint32_t, uint32_t> species2rep;
    species2rep.reserve(tax.id2str.size() / 2 + 1);
    for (uint32_t tid_id = 0; tid_id < tax.id2str.size(); ++tid_id) {
      const std::string &taxid = tax.id2str[tid_id];
      uint32_t tid = 0;
      if (!chimera::utils::try_parse_u32(taxid, tid)) {
        continue;
      }
      uint32_t sid = weightCtx.ncbiTaxdump->to_species(tid);
      uint32_t gid = weightCtx.ncbiTaxdump->to_genus(tid);
      tid2genus[tid_id] = gid;
      auto it = species2rep.find(sid);
      if (it == species2rep.end()) {
        species2rep.emplace(sid, tid_id);
        tid2speciesRep[tid_id] = tid_id;
      } else {
        tid2speciesRep[tid_id] = it->second;
      }
    }
    weightCtx.tid2speciesRep = &tid2speciesRep;
    weightCtx.tid2genus = &tid2genus;
  }
  PresenceSummary presenceSummary(config.presence_breadth_bits);
  PresenceSummary *presencePtr = &presenceSummary;

  config.tail_risk_u = 1.0;
  config.tail_risk_s = 1.0;

  constexpr uint32_t kTailRiskProbeReads = 200000;
  constexpr std::size_t kTailRiskTopK = 10;
  constexpr double kTailRiskTopMassAnchor = 0.75;
  constexpr double kTailRiskEffSpeciesAnchor = 32.0;
  constexpr double kTailRiskRAnchor = 1.0 - kTailRiskTopMassAnchor;
  constexpr int kTailRiskBeta = 5; // align with profile ceil(1/PROB_HI_BASE)

  TailRiskProbeStats probeStats;
  if (kTailRiskProbeReads > 0) {
    FileInfo probeInfo;
    AutoClassifyPolicy probePolicy{};
    probePolicy.candidate = derive_candidate_policy(config);
    std::vector<moodycamel::ConcurrentQueue<batchReads>> probeQueues(
        static_cast<size_t>(std::max<uint16_t>(1, config.threads)));
      std::vector<classifyResult> probeResults;
      std::atomic<bool> probe_done{false};
      std::thread probeProducer([&]() {
        parseReads(probeQueues, config, probeInfo,
                 kTailRiskProbeReads);
      probe_done.store(true, std::memory_order_release);
    });

    classify_streaming(imcfConfig, probeQueues, config, imcf, indexToTaxid,
                       tax, probeResults, probeInfo, probe_done, feature_params,
                       feature_min_len, weightCtx, probePolicy, nullptr);
    probeProducer.join();

    if (!probeResults.empty()) {
      std::sort(probeResults.begin(), probeResults.end(),
                [](const classifyResult &a, const classifyResult &b) {
                  return a.id < b.id;
                });
    }

    bool probe_em = false;
    if (!probeResults.empty()) {
      EMOptions options;
      options.temp = 1.05;
      options.prune_ratio = config.em_prune_ratio;
      options.conf_power = config.em_conf_power;
      auto [posterior, weights] =
          EMAlgorithm(probeResults, config.emIter, 0.0, options, nullptr);
      probeResults = std::move(posterior);
      probe_em = true;
    }

    std::unordered_map<std::string, double> species_counts;
    species_counts.reserve(probeResults.size() / 4 + 8);
    double unclassified_weight = 0.0;

    for (const auto &res : probeResults) {
      double weight = (res.sample_weight > 0.0) ? res.sample_weight : 1.0;
      std::string top;
      if (probe_em && !res.posteriors.empty()) {
        double best_prob = -1.0;
        for (const auto &kv : res.posteriors) {
          if (kv.second > best_prob) {
            best_prob = kv.second;
            top = kv.first;
          }
        }
      } else {
        top = res.best_taxid_hint;
        if (top.empty()) {
          double best_score = -1.0;
          for (const auto &kv : res.taxidCount) {
            if (kv.first == "unclassified") {
              continue;
            }
            if (kv.second > best_score) {
              best_score = kv.second;
              top = kv.first;
            }
          }
        }
      }
      if (top.empty() || top == "unclassified") {
        unclassified_weight += weight;
        continue;
      }

      std::string key = top;
      if (weightCtx.ncbiTaxdump && weightCtx.ncbiTaxdump->enabled()) {
        uint32_t tid = 0;
        if (chimera::utils::try_parse_u32(top, tid)) {
          uint32_t sid = weightCtx.ncbiTaxdump->to_species(tid);
          if (sid > 0) {
            key = std::to_string(sid);
          }
        }
      }
      species_counts[key] += weight;
    }

    std::vector<double> counts;
    counts.reserve(species_counts.size());
    for (const auto &kv : species_counts) {
      counts.push_back(kv.second);
    }
    std::size_t unclassified_reads =
        static_cast<std::size_t>(std::llround(unclassified_weight));
    probeStats =
        compute_tail_risk_probe_stats(counts, unclassified_reads, kTailRiskTopK);
    const double tail_risk_r = clamp01(1.0 - probeStats.top_mass);
    config.tail_risk_u = compute_tail_risk_u(tail_risk_r, probeStats.eff_species,
                                             kTailRiskRAnchor,
                                             kTailRiskEffSpeciesAnchor);
    config.tail_risk_s = compute_tail_risk_s(config.tail_risk_u, kTailRiskBeta);
  }

  const double tail_risk_hi = smoothstep(config.tail_risk_s, 0.60, 0.95);
  if (!config.firstFilterBeta_user) {
    config.firstFilterBeta = lerp(0.5, 0.8, tail_risk_hi);
  }
  config.em_conf_power = lerp(1.0, 2.0, tail_risk_hi);
  std::cout << "[classify][auto] tail-risk"
            << " u=" << std::fixed << std::setprecision(4) << config.tail_risk_u
            << " s=" << config.tail_risk_s
            << " top_mass=" << probeStats.top_mass
            << " eff_species=" << probeStats.eff_species
            << " uncls=" << probeStats.unclassified << "\n";
  std::cout << "[classify][auto] firstFilterBeta=" << config.firstFilterBeta
            << " em_conf_power=" << config.em_conf_power
            << " hi_gate=" << tail_risk_hi << "\n";
  std::cout << std::defaultfloat;

  FileInfo fileInfo;
  seqan3::contrib::bgzf_thread_count = config.threads;
  std::vector<moodycamel::ConcurrentQueue<batchReads>> readQueues(
      static_cast<size_t>(std::max<uint16_t>(1, config.threads)));
  std::vector<classifyResult> classifyResults;
  std::unordered_map<std::string, double> classWeights;
  bool posteriorModelUsed = false;
  AutoClassifyPolicy classifyPolicy{};
  classifyPolicy.candidate = derive_candidate_policy(config);

  std::atomic<bool> producer_done{false};
  std::thread producer([&]() {
    parseReads(readQueues, config, fileInfo);
    producer_done.store(true, std::memory_order_release);
  });

  classify_streaming(imcfConfig, readQueues, config, imcf, indexToTaxid, tax,
                     classifyResults, fileInfo, producer_done, feature_params,
                     feature_min_len, weightCtx, classifyPolicy, presencePtr);
  producer.join();
  // Determinism: classify_streaming merges thread-local batches in the order
  // threads finish, which is non-deterministic. Sort by read/contig id before
  // downstream EM/post-processing so repeated runs are stable.
  std::sort(classifyResults.begin(), classifyResults.end(),
            [](const classifyResult &a, const classifyResult &b) {
              return a.id < b.id;
            });
  if (fileInfo.sequenceNum > 0) {
    fileInfo.avgLen = fileInfo.bpLength / fileInfo.sequenceNum;
  }
  const AutoClassifyPolicy autoPolicy = derive_auto_policy(fileInfo, config);

  size_t presenceTotalReads = fileInfo.sequenceNum;
  size_t presenceMeanReadLen = fileInfo.avgLen;

  PresenceDecision presenceDecision = evaluate_presence_coverage(
      presenceSummary, tax, config, coverageMeta, presenceTotalReads,
      presenceMeanReadLen);

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

  // Optional: saturate per-read sample weights to avoid extreme long tails
  // dominating EM in tail-rich datasets (e.g., CAMI-long).
  //
  // Default: auto (enabled only when weights are present and highly skewed).
  {
    const double sat_strength = smoothstep(config.tail_risk_s, 0.60, 0.95);

    std::vector<double> raw;
    raw.reserve(classifyResults.size());
    for (const auto &res : classifyResults) {
      if (res.sample_weight > 0.0) {
        raw.push_back(res.sample_weight);
      }
    }
    if (!raw.empty()) {
      std::sort(raw.begin(), raw.end());
      auto pick = [&](double q) -> double {
        if (raw.empty()) {
          return 0.0;
        }
        q = std::clamp(q, 0.0, 1.0);
        size_t idx = static_cast<size_t>(
            std::floor(q * static_cast<double>(raw.size() - 1)));
        if (idx >= raw.size()) {
          idx = raw.size() - 1;
        }
        return raw[idx];
      };
      const double med = pick(0.50);
      const double p90 = pick(0.90);
      const double p99 = pick(0.99);
      const double mx = raw.back();
      const double denom = std::max(1e-12, med);
      const double ratio_max = mx / denom;
      const double ratio_p99 = p99 / denom;
      // Auto guard: only apply when the tail is clearly extreme.
      const bool skewed =
          (ratio_p99 >= 16.0) || (ratio_max >= 256.0) || (mx >= 1e5);
      if (skewed && sat_strength > 0.0 && (med > 0.0) &&
          (std::log1p(med) > 0.0)) {
        const double scale = med / std::log1p(med);
        for (auto &res : classifyResults) {
          if (res.sample_weight > 0.0) {
            const double sat = scale * std::log1p(res.sample_weight);
            res.sample_weight =
                ((1.0 - sat_strength) * res.sample_weight) + (sat_strength * sat);
          }
        }
      }
    }
  }

  EMOptions options;
  options.temp = 1.05;
  options.prune_ratio = config.em_prune_ratio;
  options.conf_power = config.em_conf_power;
  auto [posterior, weights] =
      EMAlgorithm(classifyResults, config.emIter, 0.0, options,
                  emPriorScale.empty() ? nullptr : &emPriorScale);
  classifyResults = std::move(posterior);
  classWeights = std::move(weights);
  posteriorModelUsed = true;
  if (posteriorModelUsed) {
    DecisionConfig decisionConfig;
    double tuned_post_pi_min = config.post_pi_min;
    double tuned_post_pi_min_weights = tuned_post_pi_min;
    double top10_mass = 0.0;
    double top50_mass = 0.0;
    double tail_base = 0.0;
    double head_heavy_score = 0.0;
    double relax_strength = 0.0;

    if (tuned_post_pi_min > 0.0 && !classWeights.empty()) {
      std::vector<double> weights;
      weights.reserve(classWeights.size());
      double total_mass = 0.0;
      for (const auto &kv : classWeights) {
        const double w = kv.second;
        if (w > 0.0) {
          weights.push_back(w);
          total_mass += w;
        }
      }
      if (total_mass > 0.0 && weights.size() > 1) {
        std::sort(weights.begin(), weights.end(),
                  [](double a, double b) { return a > b; });
        const size_t topK = std::min<size_t>(10, weights.size());
        for (size_t i = 0; i < topK; ++i) {
          top10_mass += weights[i];
        }
        top10_mass /= total_mass;
        const size_t topK50 = std::min<size_t>(50, weights.size());
        for (size_t i = 0; i < topK50; ++i) {
          top50_mass += weights[i];
        }
        top50_mass /= total_mass;

        const double h10 = clamp01((top10_mass - 0.65) / (1.0 - 0.65));
        const double h50 = clamp01((top50_mass - 0.90) / (1.0 - 0.90));
        head_heavy_score = std::max(h10, h50);
        relax_strength =
            clamp01(config.tail_risk_s * (1.0 - head_heavy_score));

        auto tail_mass = [&](double pi_min) -> double {
          if (!(pi_min > 0.0)) {
            return 0.0;
          }
          double tail = 0.0;
          for (size_t i = weights.size(); i-- > 0;) {
            const double w = weights[i];
            if (w >= pi_min) {
              break;
            }
            tail += w;
          }
          return tail / total_mass;
        };

        const double base = tuned_post_pi_min;
        tail_base = tail_mass(base);
        const double tail_trigger = 0.02;
        if (tail_base > tail_trigger && relax_strength > 0.0) {
          double target_tail = 0.01;
          if (fileInfo.avgLen >= 2000) {
            target_tail = 0.005;
          }
          const std::array<double, 4> candidates = {5e-4, 2e-4, 1e-4, 5e-5};
          double tuned_candidate = base;
          for (double c : candidates) {
            if (c > base) {
              continue;
            }
            tuned_candidate = c;
            if (tail_mass(c) <= target_tail) {
              break;
            }
          }
          if (tuned_candidate < base) {
            const double log_base = std::log10(base);
            const double log_tuned = std::log10(tuned_candidate);
            const double log_mix =
                ((1.0 - relax_strength) * log_base) +
                (relax_strength * log_tuned);
            double tuned = std::pow(10.0, log_mix);
            tuned = std::clamp(tuned, tuned_candidate, base);
            tuned_post_pi_min = tuned;
          }
        }
      }
    }
    tuned_post_pi_min_weights = tuned_post_pi_min;

    constexpr double kAutoPostPiMinLoDefault = 1e-4;
    constexpr size_t kAutoPostPiMinL0 = 2000;
    constexpr size_t kAutoPostPiMinL1 = 800;
    const double auto_pi_lo = kAutoPostPiMinLoDefault;

    auto len_tune = ChimeraClassify::tune_post_pi_min_by_avg_len(
        tuned_post_pi_min_weights, auto_pi_lo, fileInfo.avgLen, relax_strength,
        kAutoPostPiMinL0, kAutoPostPiMinL1);
    tuned_post_pi_min = len_tune.tuned;

    const bool short_like = (autoPolicy.regime == ReadRegime::ShortLike);
    const bool fallback_enabled = !autoPolicy.post.disable_fallback;
    const bool enable_selective_reject = autoPolicy.post.enable_selective_reject;

    decisionConfig.min_class_weight = tuned_post_pi_min;
    decisionConfig.fallback_strength =
        fallback_enabled ? relax_strength : 0.0;
    decisionConfig.selective_reject_strength =
        enable_selective_reject ? relax_strength : 0.0;
    decisionConfig.fallback_gap_min = 0.10;
    std::cout << "[classify][auto] post-pi"
              << " s=" << std::fixed << std::setprecision(4)
              << config.tail_risk_s << " head_heavy=" << head_heavy_score
              << " relax=" << relax_strength
              << " base=" << config.post_pi_min
              << " tail_base=" << tail_base
              << " tuned=" << tuned_post_pi_min
              << " len_t=" << len_tune.t
              << " avg_len=" << fileInfo.avgLen
              << " regime=" << (short_like ? "short_like" : "long_like")
              << " fallback_enabled=" << (fallback_enabled ? 1 : 0)
              << " enable_selective_reject="
              << (enable_selective_reject ? 1 : 0)
              << " fallback_strength=" << decisionConfig.fallback_strength
              << " selective_strength="
              << decisionConfig.selective_reject_strength
              << " fallback_gap_min=" << decisionConfig.fallback_gap_min
              << "\n";
    std::cout << std::defaultfloat;

    postEmDecision(classifyResults, decisionConfig, classWeights, tax,
                   &presenceDecision, weightCtx.ncbiTaxdump, autoPolicy.post);
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

  saveResult(classifyResults, config);
}

} // namespace ChimeraClassify
