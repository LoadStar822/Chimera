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
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <thread>

namespace {

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

struct SpoolPreparedCandidate {
  size_t idx{0};
  double log_likelihood{0.0};
};

struct SpoolEMFit {
  std::vector<uint32_t> taxid_list;
  std::unordered_map<uint32_t, size_t> taxid_to_idx;
  std::vector<double> class_pi_vec;
  std::vector<double> posterior_pi_vec;
  std::unordered_map<std::string, double> class_weights;
  size_t read_count{0};
  size_t candidate_count{0};
};

static std::string spool_taxid_to_string(
    uint32_t tid, const ChimeraClassify::TaxDict &tax) {
  if (tid == ChimeraClassify::kSpoolUnclassifiedTid ||
      tid >= tax.id2str.size()) {
    return "unclassified";
  }
  return tax.id2str[tid];
}

static void normalize_distribution(std::vector<double> &dist) {
  double sum = 0.0;
  for (double v : dist) {
    sum += v;
  }
  if (!(sum > 0.0)) {
    const double uniform = dist.empty() ? 0.0 : 1.0 / static_cast<double>(dist.size());
    for (double &v : dist) {
      v = uniform;
    }
    return;
  }
  const double inv_sum = 1.0 / sum;
  for (double &v : dist) {
    v *= inv_sum;
  }
}

static bool prepare_spool_candidates(
    const ChimeraClassify::SpoolReadRecord &record,
    const std::unordered_map<uint32_t, size_t> &taxid_to_idx,
    const ChimeraClassify::EMOptions &options,
    std::vector<SpoolPreparedCandidate> &prepared) {
  prepared.clear();
  if (record.evaluated <= 0.0 || record.candidates.empty()) {
    return false;
  }

  struct TempCandidate {
    size_t idx{0};
    double c{0.0};
  };
  std::vector<TempCandidate> tmp;
  tmp.reserve(record.candidates.size());
  double max_count = 0.0;
  for (const auto &cand : record.candidates) {
    if (cand.tid == ChimeraClassify::kSpoolUnclassifiedTid ||
        !(cand.score > 0.0)) {
      continue;
    }
    auto it = taxid_to_idx.find(cand.tid);
    if (it == taxid_to_idx.end()) {
      continue;
    }
    const double c = cand.score / record.evaluated;
    tmp.push_back(TempCandidate{it->second, c});
    if (c > max_count) {
      max_count = c;
    }
  }
  if (tmp.empty()) {
    return false;
  }

  double denom = max_count + options.eps * static_cast<double>(tmp.size());
  if (denom <= 0.0) {
    denom = options.eps * static_cast<double>(tmp.size());
  }
  prepared.reserve(tmp.size());
  for (const auto &cand : tmp) {
    const double likelihood = (cand.c + options.eps) / denom;
    const double log_likelihood =
        options.temp * std::log(std::max(likelihood, options.eps));
    prepared.push_back(SpoolPreparedCandidate{cand.idx, log_likelihood});
  }
  return !prepared.empty();
}

static SpoolEMFit fit_spool_em(
    const std::vector<std::string> &spoolPaths, const ChimeraClassify::TaxDict &tax,
    size_t maxIter, const ChimeraClassify::EMOptions &options,
    const std::unordered_map<std::string, double> *prior_scale) {
  SpoolEMFit fit;
  std::unordered_set<uint32_t> taxid_set;
  std::unordered_map<uint32_t, double> weighted_evidence_by_tid;

  for (const auto &path : spoolPaths) {
    std::ifstream is(path, std::ios::binary);
    if (!is.is_open()) {
      throw std::runtime_error("Failed to open classify spool: " + path);
    }
    ChimeraClassify::read_spool_header(is, path);
    ChimeraClassify::SpoolReadRecord record;
    while (ChimeraClassify::read_spool_record(is, record)) {
      ++fit.read_count;
      fit.candidate_count += record.candidates.size();
      if (record.evaluated <= 0.0 || record.candidates.empty()) {
        continue;
      }
      double total_count = 0.0;
      for (const auto &cand : record.candidates) {
        if (cand.tid == ChimeraClassify::kSpoolUnclassifiedTid ||
            !(cand.score > 0.0)) {
          continue;
        }
        taxid_set.insert(cand.tid);
        total_count += cand.score;
      }
      if (total_count > 0.0) {
        for (const auto &cand : record.candidates) {
          if (cand.tid != ChimeraClassify::kSpoolUnclassifiedTid &&
              cand.score > 0.0) {
            weighted_evidence_by_tid[cand.tid] += cand.score / total_count;
          }
        }
      }
    }
  }

  fit.taxid_list.reserve(taxid_set.size());
  for (uint32_t tid : taxid_set) {
    fit.taxid_list.push_back(tid);
  }
  if (fit.taxid_list.empty()) {
    return fit;
  }

  const size_t taxid_count = fit.taxid_list.size();
  fit.taxid_to_idx.reserve(taxid_count);
  for (size_t idx = 0; idx < taxid_count; ++idx) {
    fit.taxid_to_idx.emplace(fit.taxid_list[idx], idx);
  }

  std::vector<double> weighted_evidence(taxid_count, 0.0);
  for (size_t idx = 0; idx < taxid_count; ++idx) {
    auto it = weighted_evidence_by_tid.find(fit.taxid_list[idx]);
    if (it != weighted_evidence_by_tid.end()) {
      weighted_evidence[idx] = it->second;
    }
  }

  auto get_scale = [&](size_t idx) -> double {
    if (!prior_scale) {
      return 1.0;
    }
    const std::string taxid = spool_taxid_to_string(fit.taxid_list[idx], tax);
    auto it = prior_scale->find(taxid);
    if (it == prior_scale->end() || !(it->second > 0.0)) {
      return 1.0;
    }
    return std::sqrt(it->second);
  };

  fit.class_pi_vec.assign(taxid_count, 0.0);
  constexpr double kBackgroundPrior = 1e-9;
  for (size_t idx = 0; idx < taxid_count; ++idx) {
    const double scale = get_scale(idx);
    double init_val = weighted_evidence[idx] * scale;
    if (init_val <= 0.0 && prior_scale && scale > 0.0) {
      init_val = 1e-5;
    }
    if (init_val <= 0.0) {
      init_val = kBackgroundPrior;
    }
    fit.class_pi_vec[idx] = init_val;
  }
  normalize_distribution(fit.class_pi_vec);
  const std::vector<double> abundance_prior = fit.class_pi_vec;

#ifdef _OPENMP
  const int num_threads = omp_get_max_threads();
#else
  const int num_threads = 1;
#endif
  std::vector<std::vector<double>> thread_expected(
      num_threads, std::vector<double>(taxid_count, 0.0));
  std::vector<double> expected_counts(taxid_count, 0.0);

  for (size_t iteration = 0; iteration < maxIter; ++iteration) {
    const bool final_iter = (iteration + 1 == maxIter);
    if (final_iter) {
      fit.posterior_pi_vec = fit.class_pi_vec;
    }
    for (auto &local : thread_expected) {
      std::fill(local.begin(), local.end(), 0.0);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      auto &local_expected = thread_expected[tid];
      std::vector<SpoolPreparedCandidate> prepared;
      std::vector<std::pair<size_t, double>> log_components;
      std::vector<double> q_values;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for (size_t path_idx = 0; path_idx < spoolPaths.size(); ++path_idx) {
        const auto &path = spoolPaths[path_idx];
        std::ifstream is(path, std::ios::binary);
        if (!is.is_open()) {
          throw std::runtime_error("Failed to open classify spool: " + path);
        }
        ChimeraClassify::read_spool_header(is, path);
        ChimeraClassify::SpoolReadRecord record;
        while (ChimeraClassify::read_spool_record(is, record)) {
          if (!prepare_spool_candidates(record, fit.taxid_to_idx, options,
                                        prepared)) {
            continue;
          }
          log_components.clear();
          q_values.clear();
          log_components.reserve(prepared.size());
          q_values.reserve(prepared.size());
          double max_log = -std::numeric_limits<double>::infinity();
          for (const auto &cand : prepared) {
            const double prior = std::max(fit.class_pi_vec[cand.idx],
                                          options.eps);
            const double score = std::log(prior) + cand.log_likelihood;
            log_components.emplace_back(cand.idx, score);
            if (score > max_log) {
              max_log = score;
            }
          }
          if (log_components.empty() || !std::isfinite(max_log)) {
            continue;
          }
          double sum_exp = 0.0;
          for (const auto &entry : log_components) {
            sum_exp += std::exp(entry.second - max_log);
          }
          const double normalizer = max_log + std::log(sum_exp);
          double max_post = 0.0;
          for (const auto &entry : log_components) {
            const double q = std::exp(entry.second - normalizer);
            q_values.push_back(q);
            if (q > max_post) {
              max_post = q;
            }
          }
          double conf_w = 1.0;
          if (options.conf_power > 0.0 && max_post > 0.0) {
            conf_w = std::pow(max_post, options.conf_power);
          }
          const double weight =
              (record.evaluated > 0.0 ? record.evaluated : 1.0) * conf_w;
          for (size_t k = 0; k < log_components.size(); ++k) {
            local_expected[log_components[k].first] += weight * q_values[k];
          }
        }
      }
    }

    std::fill(expected_counts.begin(), expected_counts.end(), 0.0);
    double sum_expected = 0.0;
    for (const auto &local : thread_expected) {
      for (size_t idx = 0; idx < taxid_count; ++idx) {
        expected_counts[idx] += local[idx];
      }
    }
    for (double v : expected_counts) {
      sum_expected += v;
    }

    const double prior_mass = std::max(0.0, options.prior_strength);
    const double uniform_mass = options.alpha * static_cast<double>(taxid_count);
    const double pseudo_mass = (prior_mass > 0.0) ? prior_mass : uniform_mass;
    double denominator = pseudo_mass + sum_expected;
    if (denominator <= 0.0) {
      denominator = 1.0;
    }

    double max_expected = 0.0;
    for (double v : expected_counts) {
      if (v > max_expected) {
        max_expected = v;
      }
    }
    const double prune_thres = max_expected * options.prune_ratio;
    for (size_t idx = 0; idx < taxid_count; ++idx) {
      double expected = expected_counts[idx];
      if (expected < prune_thres) {
        expected = 0.0;
      }
      const double prior_component =
          (prior_mass > 0.0) ? prior_mass * abundance_prior[idx]
                             : options.alpha;
      fit.class_pi_vec[idx] = (expected + prior_component) / denominator;
    }
  }

  if (fit.posterior_pi_vec.empty()) {
    fit.posterior_pi_vec = fit.class_pi_vec;
  }
  fit.class_weights.reserve(taxid_count);
  for (size_t idx = 0; idx < taxid_count; ++idx) {
    fit.class_weights.emplace(spool_taxid_to_string(fit.taxid_list[idx], tax),
                              fit.class_pi_vec[idx]);
  }
  return fit;
}

static std::vector<std::pair<std::string, double>>
materialize_spool_posterior(const ChimeraClassify::SpoolReadRecord &record,
                            const SpoolEMFit &fit,
                            const ChimeraClassify::TaxDict &tax,
                            const ChimeraClassify::EMOptions &options) {
  std::vector<std::pair<std::string, double>> posterior;
  std::vector<SpoolPreparedCandidate> prepared;
  if (!prepare_spool_candidates(record, fit.taxid_to_idx, options, prepared)) {
    return posterior;
  }

  std::vector<std::pair<size_t, double>> log_components;
  log_components.reserve(prepared.size());
  double max_log = -std::numeric_limits<double>::infinity();
  for (const auto &cand : prepared) {
    const double prior = std::max(fit.posterior_pi_vec[cand.idx], options.eps);
    const double score = std::log(prior) + cand.log_likelihood;
    log_components.emplace_back(cand.idx, score);
    if (score > max_log) {
      max_log = score;
    }
  }
  if (log_components.empty() || !std::isfinite(max_log)) {
    return posterior;
  }
  double sum_exp = 0.0;
  for (const auto &entry : log_components) {
    sum_exp += std::exp(entry.second - max_log);
  }
  const double normalizer = max_log + std::log(sum_exp);
  posterior.reserve(log_components.size());
  for (const auto &entry : log_components) {
    const double q = std::exp(entry.second - normalizer);
    posterior.emplace_back(
        spool_taxid_to_string(fit.taxid_list[entry.first], tax), q);
  }
  return posterior;
}

static ChimeraClassify::classifyResult
materialize_spool_result(const ChimeraClassify::SpoolReadRecord &record,
                         const SpoolEMFit &fit,
                         const ChimeraClassify::TaxDict &tax,
                         const ChimeraClassify::EMOptions &options) {
  ChimeraClassify::classifyResult result;
  result.id = record.id;
  result.evaluated = record.evaluated;
  result.reject_reason = record.reject_reason;
  if (record.best_taxid_hint != ChimeraClassify::kSpoolUnclassifiedTid) {
    result.best_taxid_hint = spool_taxid_to_string(record.best_taxid_hint, tax);
  }
  result.posteriors = materialize_spool_posterior(record, fit, tax, options);
  return result;
}

static std::string resolve_tsv_output_path(const std::string &path) {
  if (std::filesystem::path(path).extension() == ".tsv") {
    return path;
  }
  return path + ".tsv";
}

static void write_spool_em_results(
    const std::vector<std::string> &spoolPaths, const SpoolEMFit &fit,
    const ChimeraClassify::EMOptions &options,
    const ChimeraClassify::DecisionConfig &decisionConfig,
    const ChimeraClassify::TaxDict &tax,
    const ChimeraClassify::PresenceDecision *presenceDecision,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump,
    const ChimeraClassify::PostDecisionPolicy &postPolicy,
    const ChimeraClassify::ClassifyConfig &config,
    ChimeraClassify::FileInfo &fileInfo) {
  const std::string outputFile = resolve_tsv_output_path(config.outputFile);
  struct PartStats {
    uint64_t classified{0};
    uint64_t unclassified{0};
  };

  const size_t part_count = spoolPaths.size();
  std::vector<std::string> partPaths(part_count);
  for (size_t i = 0; i < part_count; ++i) {
    partPaths[i] = outputFile + ".part." + std::to_string(i) + ".tmp";
  }
  std::vector<PartStats> partStats(part_count);
  std::vector<std::exception_ptr> partErrors(part_count);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::int64_t part_i = 0;
       part_i < static_cast<std::int64_t>(part_count); ++part_i) {
    const size_t part_idx = static_cast<size_t>(part_i);
    try {
      std::ofstream partOs(partPaths[part_idx],
                           std::ios::out | std::ios::binary);
      if (!partOs.is_open()) {
        throw std::runtime_error("Failed to open classify output part: " +
                                 partPaths[part_idx]);
      }
      std::vector<char> partBuffer(1 << 20, '\0');
      partOs.rdbuf()->pubsetbuf(
          partBuffer.data(), static_cast<std::streamsize>(partBuffer.size()));
      std::ostringstream postTopkOss;
      std::vector<ChimeraClassify::classifyResult> chunk;
      chunk.reserve(4096);

      auto flush_chunk = [&]() {
        if (chunk.empty()) {
          return;
        }
        ChimeraClassify::postEmDecision(chunk, decisionConfig,
                                        fit.class_weights, tax, presenceDecision,
                                        ncbiTaxdump, postPolicy);
        for (const auto &result : chunk) {
          if (!result.taxidCount.empty() &&
              result.taxidCount.front().first == "unclassified") {
            ++partStats[part_idx].unclassified;
          } else {
            ++partStats[part_idx].classified;
          }
          ChimeraClassify::writeResultRecord(partOs, result, postTopkOss);
        }
        chunk.clear();
      };

      std::ifstream is(spoolPaths[part_idx], std::ios::binary);
      if (!is.is_open()) {
        throw std::runtime_error("Failed to open classify spool: " +
                                 spoolPaths[part_idx]);
      }
      ChimeraClassify::read_spool_header(is, spoolPaths[part_idx]);
      ChimeraClassify::SpoolReadRecord record;
      while (ChimeraClassify::read_spool_record(is, record)) {
        chunk.push_back(materialize_spool_result(record, fit, tax, options));
        if (chunk.size() >= 4096) {
          flush_chunk();
        }
      }
      flush_chunk();
      partOs.close();
      if (!partOs.good()) {
        throw std::runtime_error("Failed to close classify output part: " +
                                 partPaths[part_idx]);
      }
    } catch (...) {
      partErrors[part_idx] = std::current_exception();
    }
  }

  auto cleanup_parts = [&]() {
    for (const auto &path : partPaths) {
      std::error_code ec;
      std::filesystem::remove(path, ec);
    }
  };
  for (const auto &error : partErrors) {
    if (error) {
      cleanup_parts();
      std::rethrow_exception(error);
    }
  }

  std::ofstream os(outputFile, std::ios::out | std::ios::binary);
  if (!os.is_open()) {
    cleanup_parts();
    throw std::runtime_error("Failed to open file: " + outputFile);
  }
  std::vector<char> outputBuffer(1 << 20, '\0');
  os.rdbuf()->pubsetbuf(outputBuffer.data(),
                        static_cast<std::streamsize>(outputBuffer.size()));
  fileInfo.classifiedNum = 0;
  fileInfo.unclassifiedNum = 0;
  for (size_t i = 0; i < part_count; ++i) {
    std::ifstream partIs(partPaths[i], std::ios::in | std::ios::binary);
    if (!partIs.is_open()) {
      cleanup_parts();
      throw std::runtime_error("Failed to open classify output part: " +
                               partPaths[i]);
    }
    os << partIs.rdbuf();
    if (!os.good()) {
      cleanup_parts();
      throw std::runtime_error("Failed to append classify output part: " +
                               partPaths[i]);
    }
    fileInfo.classifiedNum += partStats[i].classified;
    fileInfo.unclassifiedNum += partStats[i].unclassified;
  }
  os.close();
  cleanup_parts();
}

static std::filesystem::path make_spool_dir(const std::string &outputFile) {
  return std::filesystem::path(outputFile + ".spool");
}

static std::vector<std::string>
make_spool_paths(const std::filesystem::path &spoolDir, size_t workerCount) {
  std::vector<std::string> paths;
  paths.reserve(workerCount);
  for (size_t i = 0; i < workerCount; ++i) {
    paths.push_back((spoolDir / ("worker_" + std::to_string(i) + ".cspool"))
                        .string());
  }
  return paths;
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
  const bool freq_model_enabled = coverageMeta.freq_model.enabled();
  WeightingContext weightCtx;
  std::unique_ptr<CountMinSketch> freqSketch;
  std::unique_ptr<NcbiTaxdump> ncbiTaxdump;
  if (freq_model_enabled) {
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
  if (freq_model_enabled) {
    std::vector<uint32_t>().swap(coverageMeta.freq_model.counters);
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
  if (freq_model_enabled && !coverageMeta.entries.empty() && uniq_nonzero == 0) {
    freq_trusted = false;
  }
  weightCtx.freq_trusted = freq_trusted;
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
  if (imcfConfig.featureMethod != 1) {
    throw std::runtime_error(
        "This Chimera version no longer supports syncmer databases. Rebuild the database with strobemer.");
  }
  if (imcfConfig.strobeK == 0) {
    throw std::runtime_error("The IMCF database is missing strobemer parameters and cannot be used for classification.");
  }

  if (!chimera::feature::strobemer_available()) {
    throw std::runtime_error("This Chimera build does not include strobemer support and cannot load a strobemer database.");
  }

  size_t feature_min_len = 0;
  chimera::feature::Params feature_params =
      prepare_feature_params_for_classify(imcfConfig, feature_min_len);

  TaxDict tax = build_tax_dict(indexToTaxid);
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
  std::vector<std::vector<std::string>>().swap(indexToTaxid);
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
    std::vector<QueueThrottle> probeThrottles(probeQueues.size());
    std::vector<classifyResult> probeResults;
    std::atomic<bool> probe_done{false};
    std::thread probeProducer([&]() {
      parseReads(probeQueues, config, probeInfo, kTailRiskProbeReads,
                 &probeThrottles);
      probe_done.store(true, std::memory_order_release);
    });

    classify_streaming(imcfConfig, probeQueues, config, imcf, tax, probeResults,
                       probeInfo, probe_done, feature_params, feature_min_len,
                       weightCtx, probePolicy, nullptr, &probeThrottles);
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
          EMAlgorithm(std::move(probeResults), config.emIter, 0.0, options,
                      nullptr);
      probeResults = std::move(posterior);
      probe_em = true;
    }

    std::unordered_map<std::string, double> species_counts;
    species_counts.reserve(probeResults.size() / 4 + 8);
    std::size_t unclassified_reads = 0;

    for (const auto &res : probeResults) {
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
        ++unclassified_reads;
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
      species_counts[key] += 1.0;
    }

    std::vector<double> counts;
    counts.reserve(species_counts.size());
    for (const auto &kv : species_counts) {
      counts.push_back(kv.second);
    }
    probeStats = compute_tail_risk_probe_stats(counts, unclassified_reads,
                                               kTailRiskTopK);
    const double tail_risk_r = clamp01(1.0 - probeStats.top_mass);
    config.tail_risk_u = compute_tail_risk_u(
        tail_risk_r, probeStats.eff_species, kTailRiskRAnchor,
        kTailRiskEffSpeciesAnchor);
    config.tail_risk_s = compute_tail_risk_s(config.tail_risk_u, kTailRiskBeta);
    std::vector<classifyResult>().swap(probeResults);
    std::vector<moodycamel::ConcurrentQueue<batchReads>>().swap(probeQueues);
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
  std::vector<QueueThrottle> queueThrottles(readQueues.size());
  std::unordered_map<std::string, double> classWeights;
  bool posteriorModelUsed = false;
  AutoClassifyPolicy classifyPolicy{};
  classifyPolicy.candidate = derive_candidate_policy(config);
  const std::filesystem::path spoolDir = make_spool_dir(config.outputFile);
  std::filesystem::remove_all(spoolDir);
  std::filesystem::create_directories(spoolDir);
  const std::vector<std::string> spoolPaths =
      make_spool_paths(spoolDir, readQueues.size());

  std::atomic<bool> producer_done{false};
  std::thread producer([&]() {
    parseReads(readQueues, config, fileInfo, 0, &queueThrottles);
    producer_done.store(true, std::memory_order_release);
  });

  classify_streaming_spool(imcfConfig, readQueues, config, imcf, tax,
                           spoolPaths, fileInfo, producer_done, feature_params,
                           feature_min_len, weightCtx, classifyPolicy,
                           presencePtr, &queueThrottles);
  producer.join();
  std::vector<moodycamel::ConcurrentQueue<batchReads>>().swap(readQueues);
  if (fileInfo.sequenceNum > 0) {
    fileInfo.avgLen = fileInfo.bpLength / fileInfo.sequenceNum;
  }
  const AutoClassifyPolicy autoPolicy = derive_auto_policy(fileInfo, config);

  size_t presenceTotalReads = fileInfo.sequenceNum;
  size_t presenceMeanReadLen = fileInfo.avgLen;

  PresenceDecision presenceDecision = evaluate_presence_coverage(
      presenceSummary, tax, config, coverageMeta, presenceTotalReads,
      presenceMeanReadLen);
  PresenceSummary().stats.swap(presenceSummary.stats);
  presenceSummary.sketchBits = 0;
  presenceSummary.sketchWords = 0;
  presencePtr = nullptr;

  std::unordered_map<std::string, double> emPriorScale;
  if (!coverageMeta.entries.empty()) {
    const uint16_t span =
        (coverageMeta.effective_span > 0) ? coverageMeta.effective_span
                                          : static_cast<uint16_t>(1);
    const size_t read_len_for_prior =
        (fileInfo.avgLen > 0)
            ? fileInfo.avgLen
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
  std::vector<chimera::presence::CoverageEntry>().swap(coverageMeta.entries);

  EMOptions options;
  options.temp = 1.05;
  options.prune_ratio = config.em_prune_ratio;
  options.conf_power = config.em_conf_power;
  SpoolEMFit speciesFit =
      fit_spool_em(spoolPaths, tax, config.emIter, options,
                   emPriorScale.empty() ? nullptr : &emPriorScale);
  classWeights = speciesFit.class_weights;
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

    write_spool_em_results(spoolPaths, speciesFit, options, decisionConfig,
                           tax, &presenceDecision, weightCtx.ncbiTaxdump,
                           autoPolicy.post, config, fileInfo);
  }
  std::filesystem::remove_all(spoolDir);
}

} // namespace ChimeraClassify
