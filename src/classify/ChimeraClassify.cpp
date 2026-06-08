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
#include "ChimeraClassifyAutoPolicy.hpp"
#include "ChimeraLpcClassify.hpp"
#include "ChimeraClassifyReadout.hpp"

#include <utils/LocalResolutionMetadata.hpp>
#include <utils/LocalResolutionManifest.hpp>
#include <utils/NativeBoundedIndex.hpp>
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
#include <optional>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <unordered_set>

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
  if (env_dir == nullptr || *env_dir == '\0') {
    return nullptr;
  }
  std::filesystem::path base = env_dir;
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

struct SpoolAbundancePreparedCandidate {
  size_t idx{0};
  double score{0.0};
};

struct SpoolAbundanceFit {
  std::vector<uint32_t> taxid_list;
  std::unordered_map<uint32_t, size_t> taxid_to_idx;
  std::vector<double> class_pi_vec;
  std::vector<double> abundance_em_pi_vec;
  std::vector<double> abundance_mass_vec;
  std::vector<double> local_mass_vec;
  std::vector<double> read_support_vec;
  bool abundance_uses_local_mass{false};
  size_t read_count{0};
  size_t candidate_count{0};
};

struct SpoolSampleMixturePreparedCandidate {
  size_t idx{0};
  double likelihood{0.0};
};

struct SpoolSampleMixtureFit {
  std::vector<uint32_t> taxid_list;
  std::unordered_map<uint32_t, size_t> taxid_to_idx;
  std::vector<double> class_pi_vec;
  size_t read_count{0};
  size_t candidate_count{0};
};

struct AbundanceEvidenceOutput {
  double count{0.0};
  double abundance_mass{0.0};
  double primary_mass{0.0};
  double decision_mass{0.0};
  double local_mass{0.0};
  double read_support{0.0};
};

static std::string spool_taxid_to_string(
    uint32_t tid, const ChimeraClassify::TaxDict &tax) {
  if (tid == ChimeraClassify::kSpoolUnclassifiedTid ||
      tid >= tax.id2str.size()) {
    return "unclassified";
  }
  return tax.id2str[tid];
}

static bool env_flag_enabled(const char *name) {
  const char *value = std::getenv(name);
  return value != nullptr && value[0] != '\0' &&
         !(value[0] == '0' && value[1] == '\0');
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

static const std::vector<ChimeraClassify::SpoolCandidate> &
select_spool_abundance_candidates(
    const ChimeraClassify::SpoolReadRecord &record) {
  return !record.abundance_candidates.empty() ? record.abundance_candidates
                                              : record.candidates;
}

static bool prepare_spool_abundance_candidates(
    const ChimeraClassify::SpoolReadRecord &record,
    const std::unordered_map<uint32_t, size_t> &taxid_to_idx,
    const ChimeraClassify::EMOptions &options,
    std::vector<SpoolAbundancePreparedCandidate> &prepared,
    double &read_weight) {
  prepared.clear();
  read_weight = 0.0;
  const auto &candidates = select_spool_abundance_candidates(record);
  if (candidates.empty()) {
    return false;
  }

  prepared.reserve(candidates.size());
  for (const auto &cand : candidates) {
    if (cand.tid == ChimeraClassify::kSpoolUnclassifiedTid ||
        !(cand.score > 0.0)) {
      continue;
    }
    auto it = taxid_to_idx.find(cand.tid);
    if (it == taxid_to_idx.end()) {
      continue;
    }
    const double score = std::max(cand.score, options.eps);
    prepared.push_back(SpoolAbundancePreparedCandidate{it->second, score});
    read_weight += score;
  }
  return read_weight > 0.0 && !prepared.empty();
}

static bool prepare_spool_sample_mixture_candidates(
    const ChimeraClassify::SpoolReadRecord &record,
    const std::unordered_map<uint32_t, size_t> &taxid_to_idx,
    std::vector<SpoolSampleMixturePreparedCandidate> &prepared) {
  prepared.clear();
  if (record.sample_mixture_candidates.empty()) {
    return false;
  }
  double total = 0.0;
  for (const auto &cand : record.sample_mixture_candidates) {
    if (cand.tid != ChimeraClassify::kSpoolUnclassifiedTid &&
        cand.score > 0.0 && taxid_to_idx.find(cand.tid) != taxid_to_idx.end()) {
      total += cand.score;
    }
  }
  if (!(total > 0.0)) {
    return false;
  }
  prepared.reserve(record.sample_mixture_candidates.size());
  const double inv_total = 1.0 / total;
  for (const auto &cand : record.sample_mixture_candidates) {
    if (cand.tid == ChimeraClassify::kSpoolUnclassifiedTid ||
        !(cand.score > 0.0)) {
      continue;
    }
    auto it = taxid_to_idx.find(cand.tid);
    if (it == taxid_to_idx.end()) {
      continue;
    }
    prepared.push_back(
        SpoolSampleMixturePreparedCandidate{it->second, cand.score * inv_total});
  }
  return !prepared.empty();
}

template <typename RecordConsumer>
static void for_each_spool_record(const std::string &path,
                                  RecordConsumer consume_record) {
  std::ifstream is(path, std::ios::binary);
  if (!is.is_open()) {
    throw std::runtime_error("Failed to open classify spool: " + path);
  }
  ChimeraClassify::read_spool_header(is, path);
  ChimeraClassify::SpoolReadRecord record;
  while (ChimeraClassify::read_spool_record(is, record)) {
    consume_record(record);
  }
}

static SpoolEMFit fit_spool_em(
    const std::vector<std::string> &spoolPaths, const ChimeraClassify::TaxDict &tax,
    size_t maxIter, const ChimeraClassify::EMOptions &options,
    const std::unordered_map<std::string, double> *prior_scale) {
  SpoolEMFit fit;
  std::unordered_set<uint32_t> taxid_set;
  std::unordered_map<uint32_t, double> weighted_evidence_by_tid;

  for (const auto &path : spoolPaths) {
    for_each_spool_record(
        path, [&](const ChimeraClassify::SpoolReadRecord &record) {
          ++fit.read_count;
          fit.candidate_count += record.candidates.size();
          if (record.evaluated <= 0.0 || record.candidates.empty()) {
            return;
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
                weighted_evidence_by_tid[cand.tid] +=
                    cand.score / total_count;
              }
            }
          }
        });
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
        for_each_spool_record(
            path, [&](const ChimeraClassify::SpoolReadRecord &record) {
              if (!prepare_spool_candidates(record, fit.taxid_to_idx, options,
                                            prepared)) {
                return;
              }
              log_components.clear();
              q_values.clear();
              log_components.reserve(prepared.size());
              q_values.reserve(prepared.size());
              double max_log = -std::numeric_limits<double>::infinity();
              for (const auto &cand : prepared) {
                const double prior =
                    std::max(fit.class_pi_vec[cand.idx], options.eps);
                const double score = std::log(prior) + cand.log_likelihood;
                log_components.emplace_back(cand.idx, score);
                if (score > max_log) {
                  max_log = score;
                }
              }
              if (log_components.empty() || !std::isfinite(max_log)) {
                return;
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
                local_expected[log_components[k].first] +=
                    weight * q_values[k];
              }
            });
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

static SpoolAbundanceFit fit_spool_abundance_em(
    const std::vector<std::string> &spoolPaths,
    const ChimeraClassify::TaxDict &tax, size_t maxIter,
    const ChimeraClassify::EMOptions &options, double pruneRatio) {
  SpoolAbundanceFit fit;
  std::unordered_set<uint32_t> taxid_set;
  std::unordered_map<uint32_t, double> local_mass_by_tid;
  std::unordered_map<uint32_t, double> read_support_by_tid;

  for (const auto &path : spoolPaths) {
    for_each_spool_record(
        path, [&](const ChimeraClassify::SpoolReadRecord &record) {
          const auto &candidates = select_spool_abundance_candidates(record);
          if (candidates.empty()) {
            return;
          }
          ++fit.read_count;
          fit.candidate_count += candidates.size();
          std::unordered_set<uint32_t> seen_in_read;
          seen_in_read.reserve(candidates.size());
          for (const auto &cand : candidates) {
            if (cand.tid == ChimeraClassify::kSpoolUnclassifiedTid ||
                cand.tid >= tax.id2str.size() || !(cand.score > 0.0)) {
              continue;
            }
            taxid_set.insert(cand.tid);
            local_mass_by_tid[cand.tid] += cand.score;
            seen_in_read.insert(cand.tid);
          }
          for (uint32_t tid : seen_in_read) {
            read_support_by_tid[tid] += 1.0;
          }
        });
  }

  fit.taxid_list.reserve(taxid_set.size());
  for (uint32_t tid : taxid_set) {
    fit.taxid_list.push_back(tid);
  }
  if (fit.taxid_list.empty()) {
    return fit;
  }
  std::sort(fit.taxid_list.begin(), fit.taxid_list.end());

  const size_t taxid_count = fit.taxid_list.size();
  fit.taxid_to_idx.reserve(taxid_count);
  fit.local_mass_vec.assign(taxid_count, 0.0);
  fit.read_support_vec.assign(taxid_count, 0.0);
  for (size_t idx = 0; idx < taxid_count; ++idx) {
    const uint32_t tid = fit.taxid_list[idx];
    fit.taxid_to_idx.emplace(tid, idx);
    auto local_it = local_mass_by_tid.find(tid);
    if (local_it != local_mass_by_tid.end()) {
      fit.local_mass_vec[idx] = local_it->second;
    }
    auto support_it = read_support_by_tid.find(tid);
    if (support_it != read_support_by_tid.end()) {
      fit.read_support_vec[idx] = support_it->second;
    }
  }

  fit.class_pi_vec = fit.local_mass_vec;
  normalize_distribution(fit.class_pi_vec);
  fit.abundance_mass_vec.assign(taxid_count, 0.0);

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
      fit.abundance_em_pi_vec = fit.class_pi_vec;
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
      std::vector<SpoolAbundancePreparedCandidate> prepared;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for (size_t path_idx = 0; path_idx < spoolPaths.size(); ++path_idx) {
        const auto &path = spoolPaths[path_idx];
        for_each_spool_record(
            path, [&](const ChimeraClassify::SpoolReadRecord &record) {
              double read_weight = 0.0;
              if (!prepare_spool_abundance_candidates(
                      record, fit.taxid_to_idx, options, prepared,
                      read_weight)) {
                return;
              }
              double denom = 0.0;
              for (const auto &cand : prepared) {
                denom += std::max(fit.class_pi_vec[cand.idx], options.eps) *
                         cand.score;
              }
              if (!(denom > 0.0)) {
                return;
              }
              for (const auto &cand : prepared) {
                const double prior =
                    std::max(fit.class_pi_vec[cand.idx], options.eps);
                const double q = (prior * cand.score) / denom;
                local_expected[cand.idx] += read_weight * q;
              }
            });
      }
    }

    std::fill(expected_counts.begin(), expected_counts.end(), 0.0);
    double sum_expected = 0.0;
    for (const auto &local : thread_expected) {
      for (size_t idx = 0; idx < taxid_count; ++idx) {
        expected_counts[idx] += local[idx];
      }
    }
    double max_expected = 0.0;
    for (double v : expected_counts) {
      if (v > max_expected) {
        max_expected = v;
      }
    }
    const double prune_thres = max_expected * pruneRatio;
    for (double &v : expected_counts) {
      if (v < prune_thres) {
        v = 0.0;
      }
      sum_expected += v;
    }
    if (final_iter) {
      fit.abundance_mass_vec = expected_counts;
    }

    const double pseudo_mass =
        options.alpha * static_cast<double>(taxid_count);
    double denominator = pseudo_mass + sum_expected;
    if (denominator <= 0.0) {
      denominator = 1.0;
    }
    for (size_t idx = 0; idx < taxid_count; ++idx) {
      fit.class_pi_vec[idx] = (expected_counts[idx] + options.alpha) /
                              denominator;
    }
  }

  double abundance_sum = 0.0;
  for (double v : fit.abundance_mass_vec) {
    abundance_sum += v;
  }
  if (!(abundance_sum > 0.0)) {
    fit.abundance_mass_vec = fit.local_mass_vec;
    fit.abundance_uses_local_mass = true;
  }
  return fit;
}

static SpoolSampleMixtureFit fit_spool_sample_mixture_em(
    const std::vector<std::string> &spoolPaths, size_t maxIter) {
  SpoolSampleMixtureFit fit;
  std::unordered_set<uint32_t> taxid_set;
  std::unordered_map<uint32_t, double> local_evidence_by_tid;

  for (const auto &path : spoolPaths) {
    for_each_spool_record(
        path, [&](const ChimeraClassify::SpoolReadRecord &record) {
          if (record.sample_mixture_candidates.empty()) {
            return;
          }
          ++fit.read_count;
          fit.candidate_count += record.sample_mixture_candidates.size();
          double total = 0.0;
          for (const auto &cand : record.sample_mixture_candidates) {
            if (cand.tid != ChimeraClassify::kSpoolUnclassifiedTid &&
                cand.score > 0.0) {
              total += cand.score;
            }
          }
          if (!(total > 0.0)) {
            return;
          }
          const double inv_total = 1.0 / total;
          for (const auto &cand : record.sample_mixture_candidates) {
            if (cand.tid == ChimeraClassify::kSpoolUnclassifiedTid ||
                !(cand.score > 0.0)) {
              continue;
            }
            taxid_set.insert(cand.tid);
            local_evidence_by_tid[cand.tid] += cand.score * inv_total;
          }
        });
  }

  fit.taxid_list.reserve(taxid_set.size());
  for (uint32_t tid : taxid_set) {
    fit.taxid_list.push_back(tid);
  }
  if (fit.taxid_list.empty()) {
    return fit;
  }
  std::sort(fit.taxid_list.begin(), fit.taxid_list.end());

  const size_t taxid_count = fit.taxid_list.size();
  fit.taxid_to_idx.reserve(taxid_count);
  fit.class_pi_vec.assign(taxid_count, 0.0);
  constexpr double kJeffreysAlpha = 0.5;
  double prior_total = 0.0;
  for (size_t idx = 0; idx < taxid_count; ++idx) {
    const uint32_t tid = fit.taxid_list[idx];
    fit.taxid_to_idx.emplace(tid, idx);
    double value = kJeffreysAlpha;
    auto it = local_evidence_by_tid.find(tid);
    if (it != local_evidence_by_tid.end()) {
      value += it->second;
    }
    fit.class_pi_vec[idx] = value;
    prior_total += value;
  }
  if (!(prior_total > 0.0)) {
    normalize_distribution(fit.class_pi_vec);
  } else {
    const double inv_prior_total = 1.0 / prior_total;
    for (double &v : fit.class_pi_vec) {
      v *= inv_prior_total;
    }
  }

#ifdef _OPENMP
  const int num_threads = omp_get_max_threads();
#else
  const int num_threads = 1;
#endif
  std::vector<std::vector<double>> thread_expected(
      num_threads, std::vector<double>(taxid_count, 0.0));
  std::vector<double> expected_counts(taxid_count, 0.0);
  const size_t iterations = std::max<size_t>(1, maxIter);

  for (size_t iteration = 0; iteration < iterations; ++iteration) {
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
      std::vector<SpoolSampleMixturePreparedCandidate> prepared;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for (size_t path_idx = 0; path_idx < spoolPaths.size(); ++path_idx) {
        const auto &path = spoolPaths[path_idx];
        for_each_spool_record(
            path, [&](const ChimeraClassify::SpoolReadRecord &record) {
              if (!prepare_spool_sample_mixture_candidates(
                      record, fit.taxid_to_idx, prepared)) {
                return;
              }
              double denom = 0.0;
              for (const auto &cand : prepared) {
                denom += cand.likelihood * fit.class_pi_vec[cand.idx];
              }
              if (!(denom > 0.0)) {
                return;
              }
              const double inv_denom = 1.0 / denom;
              for (const auto &cand : prepared) {
                const double q =
                    cand.likelihood * fit.class_pi_vec[cand.idx] * inv_denom;
                local_expected[cand.idx] += q;
              }
            });
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
    const double denominator =
        sum_expected + (kJeffreysAlpha * static_cast<double>(taxid_count));
    if (!(denominator > 0.0)) {
      continue;
    }
    const double inv_denominator = 1.0 / denominator;
    for (size_t idx = 0; idx < taxid_count; ++idx) {
      fit.class_pi_vec[idx] =
          (expected_counts[idx] + kJeffreysAlpha) * inv_denominator;
    }
  }
  return fit;
}

static std::unordered_set<std::string> write_abundance_em_evidence(
    const std::string &evidencePath, const SpoolAbundanceFit &fit,
    const std::unordered_map<std::string, double> &postemPrimaryEvidenceMass,
    const std::unordered_map<std::string, double> &postemDecisionEvidenceMass,
    const ChimeraClassify::TaxDict &tax, double abundanceWeight) {
  std::unordered_map<std::string, AbundanceEvidenceOutput> merged;
  merged.reserve(fit.taxid_list.size() + postemPrimaryEvidenceMass.size() +
                 postemDecisionEvidenceMass.size());

  double abundanceSum = 0.0;
  for (size_t idx = 0; idx < fit.taxid_list.size(); ++idx) {
    if (idx < fit.abundance_mass_vec.size()) {
      abundanceSum += std::max(0.0, fit.abundance_mass_vec[idx]);
    }
  }
  double decisionSum = 0.0;
  for (const auto &[taxid, value] : postemDecisionEvidenceMass) {
    if (!taxid.empty() && taxid != "unclassified") {
      decisionSum += std::max(0.0, value);
    }
  }
  if (!(abundanceSum > 0.0)) {
    abundanceWeight = 0.0;
  }
  if (!(decisionSum > 0.0)) {
    abundanceWeight = 1.0;
  }
  abundanceWeight = ChimeraClassify::clamp01(abundanceWeight);
  const double weightEpsilon =
      std::sqrt(static_cast<double>(std::numeric_limits<float>::epsilon()));
  if (abundanceWeight >= 1.0 - weightEpsilon) {
    abundanceWeight = 1.0;
  } else if (abundanceWeight <= weightEpsilon) {
    abundanceWeight = 0.0;
  }

  for (size_t idx = 0; idx < fit.taxid_list.size(); ++idx) {
    if (idx >= fit.abundance_mass_vec.size() ||
        idx >= fit.local_mass_vec.size() ||
        idx >= fit.read_support_vec.size()) {
      continue;
    }
    const double abundanceMass = fit.abundance_mass_vec[idx];
    const double localMass = fit.local_mass_vec[idx];
    const double readSupport = fit.read_support_vec[idx];
    if (!(abundanceMass > 0.0) && !(localMass > 0.0) &&
        !(readSupport > 0.0)) {
      continue;
    }
    const std::string taxid = spool_taxid_to_string(fit.taxid_list[idx], tax);
    if (taxid.empty() || taxid == "unclassified") {
      continue;
    }
    auto &dst = merged[taxid];
    dst.abundance_mass += abundanceMass;
    dst.local_mass += localMass;
    dst.read_support += readSupport;
  }

  for (const auto &[taxid, value] : postemPrimaryEvidenceMass) {
    if (taxid.empty() || taxid == "unclassified" || !(value > 0.0)) {
      continue;
    }
    merged[taxid].primary_mass += value;
  }
  for (const auto &[taxid, value] : postemDecisionEvidenceMass) {
    if (taxid.empty() || taxid == "unclassified" || !(value > 0.0)) {
      continue;
    }
    merged[taxid].decision_mass += value;
  }

  for (auto &[taxid, evidence] : merged) {
    const double abundanceP =
        abundanceSum > 0.0 ? evidence.abundance_mass / abundanceSum : 0.0;
    const double decisionP =
        decisionSum > 0.0 ? evidence.decision_mass / decisionSum : 0.0;
    evidence.count =
        (abundanceWeight * abundanceP) + ((1.0 - abundanceWeight) * decisionP);
  }

  std::vector<std::pair<std::string, AbundanceEvidenceOutput>> activeOrder;
  activeOrder.reserve(merged.size());
  for (const auto &[taxid, evidence] : merged) {
    if (evidence.count > 0.0) {
      activeOrder.emplace_back(taxid, evidence);
    }
  }
  std::sort(activeOrder.begin(), activeOrder.end(),
            [](const auto &lhs, const auto &rhs) {
              if (lhs.second.count != rhs.second.count) {
                return lhs.second.count > rhs.second.count;
              }
              return lhs.first < rhs.first;
            });

  double countSum = 0.0;
  double countSqSum = 0.0;
  for (const auto &[taxid, evidence] : activeOrder) {
    const double value = std::max(0.0, evidence.count);
    countSum += value;
    countSqSum += value * value;
  }
  std::unordered_set<std::string> activeTaxids;
  activeTaxids.reserve(activeOrder.size());
  if (countSum > 0.0 && countSqSum > 0.0) {
    const double effectiveTaxa = (countSum * countSum) / countSqSum;
    const double residualPower = 2.0 - abundanceWeight;
    const double residualScale =
        std::pow(std::max(1.0, effectiveTaxa), residualPower);
    const double residualLimit = countSum / std::max(1.0, residualScale);
    double cumulative = 0.0;
    for (const auto &[taxid, evidence] : activeOrder) {
      activeTaxids.insert(taxid);
      cumulative += std::max(0.0, evidence.count);
      if ((countSum - cumulative) <= residualLimit) {
        break;
      }
    }
  }
  if (!activeTaxids.empty()) {
    std::vector<std::pair<std::string, double>> residualDecisionOrder;
    residualDecisionOrder.reserve(activeOrder.size());
    double residualDecisionSum = 0.0;
    double residualDecisionSqSum = 0.0;
    for (const auto &[taxid, evidence] : activeOrder) {
      if (activeTaxids.find(taxid) != activeTaxids.end()) {
        continue;
      }
      const double support = std::max(0.0, evidence.decision_mass);
      if (!(support > 0.0)) {
        continue;
      }
      residualDecisionOrder.emplace_back(taxid, support);
      residualDecisionSum += support;
      residualDecisionSqSum += support * support;
    }
    if (residualDecisionSum > 0.0 && residualDecisionSqSum > 0.0) {
      const double residualDecisionEff =
          (residualDecisionSum * residualDecisionSum) / residualDecisionSqSum;
      for (const auto &[taxid, support] : residualDecisionOrder) {
        if ((support * residualDecisionEff) > residualDecisionSum) {
          activeTaxids.insert(taxid);
        }
      }
    }
  }

  std::vector<std::pair<std::string, AbundanceEvidenceOutput>> order;
  order.reserve(activeTaxids.empty() ? activeOrder.size() : activeTaxids.size());
  for (const auto &[taxid, evidence] : activeOrder) {
    if (activeTaxids.empty() || activeTaxids.find(taxid) != activeTaxids.end()) {
      order.emplace_back(taxid, evidence);
    }
  }
  std::sort(order.begin(), order.end(),
            [](const auto &lhs, const auto &rhs) {
              return lhs.first < rhs.first;
            });

  std::ofstream aggregateOut(evidencePath, std::ios::out | std::ios::binary);
  if (!aggregateOut.is_open()) {
    throw std::runtime_error("Failed to open abundance evidence output: " +
                             evidencePath);
  }
  aggregateOut << "taxid\tcount\tabundance_mass\tprimary_mass"
               << "\tdecision_mass\tlocal_mass\tread_support\n";
  aggregateOut << std::setprecision(17);
  for (const auto &[taxid, evidence] : order) {
    aggregateOut << taxid << '\t' << evidence.count << '\t'
                 << evidence.abundance_mass << '\t' << evidence.primary_mass
                 << '\t' << evidence.decision_mass << '\t'
                 << evidence.local_mass << '\t'
                 << evidence.read_support << '\n';
  }
  aggregateOut.close();
  if (!aggregateOut.good()) {
    throw std::runtime_error("Failed to close abundance evidence output: " +
                             evidencePath);
  }
  std::unordered_set<std::string> emittedTaxids;
  emittedTaxids.reserve(order.size());
  for (const auto &entry : order) {
    emittedTaxids.insert(entry.first);
  }
  return emittedTaxids;
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

static std::vector<std::pair<std::string, double>>
materialize_spool_sample_mixture_posterior(
    const ChimeraClassify::SpoolReadRecord &record,
    const SpoolSampleMixtureFit &fit, const ChimeraClassify::TaxDict &tax) {
  std::vector<std::pair<std::string, double>> posterior;
  if (fit.taxid_list.empty() || record.sample_mixture_candidates.empty()) {
    return posterior;
  }
  std::vector<SpoolSampleMixturePreparedCandidate> prepared;
  if (!prepare_spool_sample_mixture_candidates(record, fit.taxid_to_idx,
                                               prepared)) {
    return posterior;
  }
  double denom = 0.0;
  for (const auto &cand : prepared) {
    if (cand.idx < fit.class_pi_vec.size()) {
      denom += cand.likelihood * fit.class_pi_vec[cand.idx];
    }
  }
  if (!(denom > 0.0)) {
    return posterior;
  }
  const double inv_denom = 1.0 / denom;
  posterior.reserve(prepared.size());
  for (const auto &cand : prepared) {
    if (cand.idx >= fit.class_pi_vec.size() ||
        cand.idx >= fit.taxid_list.size()) {
      continue;
    }
    const double q = cand.likelihood * fit.class_pi_vec[cand.idx] * inv_denom;
    posterior.emplace_back(spool_taxid_to_string(fit.taxid_list[cand.idx], tax),
                           q);
  }
  std::sort(posterior.begin(), posterior.end(), [](const auto &lhs,
                                                   const auto &rhs) {
    if (lhs.second != rhs.second) {
      return lhs.second > rhs.second;
    }
    return lhs.first < rhs.first;
  });
  return posterior;
}

static ChimeraClassify::classifyResult
materialize_spool_result(const ChimeraClassify::SpoolReadRecord &record,
                         const SpoolEMFit &fit,
                         const SpoolSampleMixtureFit &sampleMixtureFit,
                         const ChimeraClassify::TaxDict &tax,
                         const ChimeraClassify::EMOptions &options) {
  ChimeraClassify::classifyResult result;
  result.id = record.id;
  result.evaluated = record.evaluated;
  result.reject_reason = record.reject_reason;
  if (record.best_taxid_hint != ChimeraClassify::kSpoolUnclassifiedTid) {
    result.best_taxid_hint = spool_taxid_to_string(record.best_taxid_hint, tax);
  }
  result.abundanceCount.reserve(record.abundance_candidates.size());
  for (const auto &cand : record.abundance_candidates) {
    if (cand.tid != ChimeraClassify::kSpoolUnclassifiedTid &&
        cand.tid < tax.id2str.size() && cand.score > 0.0) {
      result.abundanceCount.emplace_back(spool_taxid_to_string(cand.tid, tax),
                                         cand.score);
    }
  }
  result.posteriors = materialize_spool_posterior(record, fit, tax, options);
  result.sampleMixturePosteriors =
      materialize_spool_sample_mixture_posterior(record, sampleMixtureFit, tax);
  result.sampleMixtureLocalScores.reserve(
      record.sample_mixture_candidates.size());
  for (const auto &cand : record.sample_mixture_candidates) {
    if (cand.tid != ChimeraClassify::kSpoolUnclassifiedTid &&
        cand.tid < tax.id2str.size() && cand.score > 0.0) {
      result.sampleMixtureLocalScores.emplace_back(
          spool_taxid_to_string(cand.tid, tax), cand.score);
    }
    if (cand.score > result.sampleMixtureTopScore) {
      result.sampleMixtureTopScore = cand.score;
    }
  }
  return result;
}

static std::string resolve_tsv_output_path(const std::string &path) {
  if (std::filesystem::path(path).extension() == ".tsv") {
    return path;
  }
  return path + ".tsv";
}

static std::string resolve_evidence_output_path(const std::string &outputFile) {
  std::filesystem::path path(outputFile);
  if (path.filename() == "ChimeraClassify.tsv") {
    return (path.parent_path() / "ChimeraEvidence.tsv").string();
  }
  path.replace_extension(".evidence.tsv");
  return path.string();
}

static std::string resolve_presence_output_path(const std::string &outputFile) {
  std::filesystem::path path(outputFile);
  if (path.filename() == "ChimeraClassify.tsv") {
    return (path.parent_path() / "ChimeraPresence.tsv").string();
  }
  path.replace_extension(".presence.tsv");
  return path.string();
}

static std::string
resolve_local_profile_output_path(const std::string &outputFile) {
  std::filesystem::path path(outputFile);
  if (path.filename() == "ChimeraClassify.tsv") {
    return (path.parent_path() / "ChimeraLocalProfile.json").string();
  }
  path.replace_extension(".local_profile.json");
  return path.string();
}

struct SpoolOutputPartStats {
  uint64_t classified{0};
  uint64_t unclassified{0};
};

using EvidenceAggregateMap = std::unordered_map<std::string, double>;
using LocalResolutionCallMap =
    std::unordered_map<std::string, ChimeraClassify::LocalResolutionReadCall>;

struct LocalResolutionPanel {
  uint32_t k{};
  std::vector<ChimeraClassify::LocalResolutionTarget> targets;
  uint64_t selected_groups{};
  uint64_t selected_species{};
  uint64_t selected_targets{};
  uint64_t selected_anchor_records{};
  uint64_t selected_anchor_bytes{};
  uint64_t selected_species_source_pairs{};
  uint64_t selected_source_extra_targets{};
  uint64_t source_cap1_targets{};
  uint64_t source_cap1_anchor_records{};
  uint64_t source_cap1_anchor_bytes{};
  uint64_t source_cap2_targets{};
  uint64_t source_cap2_anchor_records{};
  uint64_t source_cap2_anchor_bytes{};
  uint32_t targets_per_species{};
  uint32_t max_targets_per_group{};
  uint64_t anchor_byte_budget{};
  uint64_t group_cap_skipped_targets{};
  uint64_t group_cap_skipped_anchor_records{};
  uint64_t group_cap_skipped_anchor_bytes{};
  uint64_t budget_skipped_targets{};
  uint64_t budget_skipped_anchor_records{};
  uint64_t budget_skipped_anchor_bytes{};
};

struct LocalResolutionArtifacts {
  std::filesystem::path index_path;
  std::filesystem::path rep_metadata_path;
  uint32_t k{};
};

static void write_local_resolution_profile_json(
    const std::string &outputFile, const std::string &status,
    double sampleDivergence, double divergenceThreshold,
    const LocalResolutionPanel *panel,
    const ChimeraClassify::LocalResolutionStats *stats,
    double metadataSeconds, double panelSeconds, double engineSeconds) {
  const std::filesystem::path path =
      resolve_local_profile_output_path(outputFile);
  const auto parent = path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(path, std::ios::out | std::ios::binary);
  if (!out) {
    throw std::runtime_error("Failed to open local profile output: " +
                             path.string());
  }

  out << std::setprecision(10);
  out << "{\n";
  out << "  \"status\": \"" << status << "\",\n";
  out << "  \"sample_divergence\": " << sampleDivergence << ",\n";
  out << "  \"divergence_threshold\": " << divergenceThreshold << ",\n";
  out << "  \"metadata_seconds\": " << metadataSeconds << ",\n";
  out << "  \"panel_seconds\": " << panelSeconds << ",\n";
  out << "  \"engine_seconds\": " << engineSeconds << ",\n";
  out << "  \"panel\": {\n";
  out << "    \"selected_groups\": "
      << (panel == nullptr ? 0 : panel->selected_groups) << ",\n";
  out << "    \"selected_species\": "
      << (panel == nullptr ? 0 : panel->selected_species) << ",\n";
  out << "    \"selected_targets\": "
      << (panel == nullptr ? 0 : panel->selected_targets) << ",\n";
  out << "    \"selected_anchor_records\": "
      << (panel == nullptr ? 0 : panel->selected_anchor_records) << ",\n";
  out << "    \"selected_anchor_bytes\": "
      << (panel == nullptr ? 0 : panel->selected_anchor_bytes) << ",\n";
  out << "    \"selected_species_source_pairs\": "
      << (panel == nullptr ? 0 : panel->selected_species_source_pairs)
      << ",\n";
  out << "    \"selected_source_extra_targets\": "
      << (panel == nullptr ? 0 : panel->selected_source_extra_targets)
      << ",\n";
  out << "    \"source_cap1_targets\": "
      << (panel == nullptr ? 0 : panel->source_cap1_targets) << ",\n";
  out << "    \"source_cap1_anchor_records\": "
      << (panel == nullptr ? 0 : panel->source_cap1_anchor_records) << ",\n";
  out << "    \"source_cap1_anchor_bytes\": "
      << (panel == nullptr ? 0 : panel->source_cap1_anchor_bytes) << ",\n";
  out << "    \"source_cap2_targets\": "
      << (panel == nullptr ? 0 : panel->source_cap2_targets) << ",\n";
  out << "    \"source_cap2_anchor_records\": "
      << (panel == nullptr ? 0 : panel->source_cap2_anchor_records) << ",\n";
  out << "    \"source_cap2_anchor_bytes\": "
      << (panel == nullptr ? 0 : panel->source_cap2_anchor_bytes) << ",\n";
  out << "    \"targets_per_species\": "
      << (panel == nullptr ? 0 : panel->targets_per_species) << ",\n";
  out << "    \"max_targets_per_group\": "
      << (panel == nullptr ? 0 : panel->max_targets_per_group) << ",\n";
  out << "    \"anchor_byte_budget\": "
      << (panel == nullptr ? 0 : panel->anchor_byte_budget) << ",\n";
  out << "    \"group_cap_skipped_targets\": "
      << (panel == nullptr ? 0 : panel->group_cap_skipped_targets) << ",\n";
  out << "    \"group_cap_skipped_anchor_records\": "
      << (panel == nullptr ? 0 : panel->group_cap_skipped_anchor_records)
      << ",\n";
  out << "    \"group_cap_skipped_anchor_bytes\": "
      << (panel == nullptr ? 0 : panel->group_cap_skipped_anchor_bytes)
      << ",\n";
  out << "    \"budget_skipped_targets\": "
      << (panel == nullptr ? 0 : panel->budget_skipped_targets) << ",\n";
  out << "    \"budget_skipped_anchor_records\": "
      << (panel == nullptr ? 0 : panel->budget_skipped_anchor_records)
      << ",\n";
  out << "    \"budget_skipped_anchor_bytes\": "
      << (panel == nullptr ? 0 : panel->budget_skipped_anchor_bytes) << "\n";
  out << "  },\n";
  out << "  \"engine\": {\n";
  auto write_stat = [&](const char *name, auto value, bool last = false) {
    out << "    \"" << name << "\": " << value;
    if (!last) {
      out << ',';
    }
    out << '\n';
  };
  if (stats != nullptr) {
    write_stat("reads", stats->reads);
    write_stat("query_hashes", stats->query_hashes);
    write_stat("target_filter", stats->target_filter);
    write_stat("target_routes", stats->target_routes);
    write_stat("core_candidate_reads", stats->core_candidate_reads);
    write_stat("scanned_shards", stats->scanned_shards);
    write_stat("skipped_shards", stats->skipped_shards);
    write_stat("selected_targets", stats->selected_targets);
    write_stat("skipped_targets", stats->skipped_targets);
    write_stat("direct_targets", stats->direct_targets);
    write_stat("target_anchor_records_scanned",
               stats->target_anchor_records_scanned);
    write_stat("target_anchor_bytes_read", stats->target_anchor_bytes_read);
    write_stat("target_anchor_records_matched",
               stats->target_anchor_records_matched);
    write_stat("target_hash_prefilter_rejects",
               stats->target_hash_prefilter_rejects);
    write_stat("direct_load_batches", stats->direct_load_batches);
    write_stat("pread_calls", stats->pread_calls);
    write_stat("pread_bytes", stats->pread_bytes);
    write_stat("raw_chain_records", stats->raw_chain_records);
    write_stat("kept_chain_records", stats->kept_chain_records);
    write_stat("index_hash_keys", stats->index_hash_keys);
    write_stat("overflow_hash_keys", stats->overflow_hash_keys);
    write_stat("dropped_broad_keys", stats->dropped_broad_keys);
    write_stat("dropped_broad_records", stats->dropped_broad_records);
    write_stat("local_hits", stats->local_hits);
    write_stat("local_absent", stats->local_absent);
    write_stat("threads", stats->threads);
    write_stat("k", static_cast<uint32_t>(stats->k));
    write_stat("w", stats->w);
    write_stat("read_seconds", stats->read_seconds);
    write_stat("ref_seconds", stats->ref_seconds);
    write_stat("target_io_seconds", stats->target_io_seconds);
    write_stat("target_read_seconds", stats->target_read_seconds);
    write_stat("target_filter_seconds", stats->target_filter_seconds);
    write_stat("target_collect_seconds", stats->target_collect_seconds);
    write_stat("posting_merge_seconds", stats->posting_merge_seconds);
    write_stat("index_finalize_seconds", stats->index_finalize_seconds);
    write_stat("chain_seconds", stats->chain_seconds, true);
  }
  out << "  }\n";
  out << "}\n";
  if (!out.good()) {
    throw std::runtime_error("Failed to write local profile output: " +
                             path.string());
  }
}

static std::optional<LocalResolutionArtifacts>
resolve_local_resolution_artifacts(const std::string &dbFile) {
  const auto manifest =
      chimera::local_resolution::load_and_verify_manifest_for_db(dbFile);
  if (!manifest.has_value() || !manifest->local_available) {
    return std::nullopt;
  }
  const std::filesystem::path corePath =
      chimera::local_resolution::core_archive_path_for(dbFile);
  return LocalResolutionArtifacts{
      chimera::local_resolution::materialize_manifest_path(
          corePath, manifest->local_index),
      chimera::local_resolution::materialize_manifest_path(
          corePath, manifest->rep_metadata),
      manifest->k,
  };
}

struct LocalResolutionPostTopkScores {
  std::unordered_map<uint32_t, double> species_scores;
  std::unordered_map<uint32_t, uint64_t> species_read_support;
  uint64_t rows{0};
  uint64_t rows_with_post_topk{0};
  uint64_t post_topk_items{0};
  uint64_t usable_post_topk_items{0};
  uint64_t singleton_rows{0};
  uint64_t ambiguous_rows{0};
  double top1_sum{0.0};
};

struct LocalResolutionEligibility {
  bool hard_ambiguity{false};
  double mean_post_items{0.0};
  double singleton_rate{0.0};
  double mean_top1{0.0};
  double ambiguous_rate{0.0};
};

struct LocalResolutionSpeciesPanelScore {
  uint32_t species{0};
  double score{0.0};
  uint64_t read_support{0};
};

static std::string local_resolution_source_key(const std::string &targetName) {
  const size_t pos = targetName.find("|contig=");
  if (pos == std::string::npos) {
    return targetName;
  }
  return targetName.substr(0, pos);
}

static std::vector<chimera::local_resolution::TargetRep>
select_local_resolution_source_targets(
    const std::vector<chimera::local_resolution::TargetRep> &rows,
    size_t max_sources) {
  if (max_sources == 0 || rows.empty()) {
    return {};
  }
  std::vector<chimera::local_resolution::TargetRep> selected;
  selected.reserve(rows.size());
  std::unordered_set<std::string> selectedSources;
  selectedSources.reserve(max_sources);
  for (const auto &row : rows) {
    const std::string source = local_resolution_source_key(row.target_name);
    if (selectedSources.contains(source)) {
      selected.push_back(row);
      continue;
    }
    if (selectedSources.size() >= max_sources) {
      continue;
    }
    selectedSources.insert(source);
    selected.push_back(row);
  }
  return selected;
}

static void finalize_local_resolution_panel_shadow(LocalResolutionPanel &panel) {
  const uint64_t localResolutionAnchorRecordBytes =
      chimera::native_bounded::anchor_record_bytes(panel.k);
  std::unordered_map<std::string, uint32_t> targetsBySpeciesSource;
  targetsBySpeciesSource.reserve(panel.targets.size());
  for (const auto &target : panel.targets) {
    const std::string sourceKey =
        std::to_string(target.species) + "\t" +
        local_resolution_source_key(target.target_name);
    auto &sourceTargetCount = targetsBySpeciesSource[sourceKey];
    ++sourceTargetCount;
    const uint64_t targetAnchorRecords = target.anchor_count;
    const uint64_t targetAnchorBytes =
        target.anchor_byte_size == 0
            ? targetAnchorRecords * localResolutionAnchorRecordBytes
            : target.anchor_byte_size;
    if (sourceTargetCount <= 1) {
      ++panel.source_cap1_targets;
      panel.source_cap1_anchor_records += targetAnchorRecords;
      panel.source_cap1_anchor_bytes += targetAnchorBytes;
    }
    if (sourceTargetCount <= 2) {
      ++panel.source_cap2_targets;
      panel.source_cap2_anchor_records += targetAnchorRecords;
      panel.source_cap2_anchor_bytes += targetAnchorBytes;
    }
  }
  panel.selected_species_source_pairs = targetsBySpeciesSource.size();
  panel.selected_source_extra_targets =
      panel.selected_targets > panel.selected_species_source_pairs
          ? panel.selected_targets - panel.selected_species_source_pairs
          : 0;
}

static LocalResolutionPanel build_local_resolution_panel(
    const std::unordered_map<uint32_t, double> &speciesScores,
    const std::unordered_map<uint32_t, uint64_t> &speciesReadSupport,
    const chimera::local_resolution::RepMetadata &repMetadata,
    const ChimeraClassify::ClassifyConfig &config, uint32_t localK) {
  const uint64_t localResolutionAnchorRecordBytes =
      chimera::native_bounded::anchor_record_bytes(localK);
  std::unordered_map<uint32_t, double> groupScores;
  for (const auto &[species, score] : speciesScores) {
    const auto entry = repMetadata.find_species(species);
    if (!entry.has_value() || entry->target_count == 0) {
      continue;
    }
    groupScores[entry->genus] += score;
  }

  std::vector<std::pair<uint32_t, double>> groups(groupScores.begin(),
                                                  groupScores.end());
  std::sort(groups.begin(), groups.end(), [](const auto &lhs,
                                             const auto &rhs) {
    if (lhs.second != rhs.second) {
      return lhs.second > rhs.second;
    }
    return lhs.first < rhs.first;
  });
  if (groups.size() > config.local_resolution_top_groups) {
    groups.resize(config.local_resolution_top_groups);
  }

  std::unordered_set<uint32_t> selectedGroupSet;
  selectedGroupSet.reserve(groups.size());
  for (const auto &[group, _] : groups) {
    selectedGroupSet.insert(group);
  }

  std::unordered_map<uint32_t, std::vector<LocalResolutionSpeciesPanelScore>>
      speciesByGroup;
  for (const auto &[species, score] : speciesScores) {
    const auto entry = repMetadata.find_species(species);
    if (!entry.has_value() || entry->target_count == 0) {
      continue;
    }
    const uint32_t group = entry->genus;
    if (!selectedGroupSet.contains(group)) {
      continue;
    }
    uint64_t readSupport = 0;
    const auto readSupportIt = speciesReadSupport.find(species);
    if (readSupportIt != speciesReadSupport.end()) {
      readSupport = readSupportIt->second;
    }
    speciesByGroup[group].push_back({species, score, readSupport});
  }

  LocalResolutionPanel panel;
  panel.k = localK;
  panel.targets_per_species = config.local_resolution_targets_per_species;
  panel.max_targets_per_group = config.local_resolution_max_targets_per_group;
  panel.anchor_byte_budget = config.local_resolution_max_anchor_bytes;
  panel.selected_groups = groups.size();
  std::vector<uint32_t> selectedSpecies;
  std::unordered_map<uint32_t, double> selectedSpeciesScores;
  for (const auto &[group, _] : groups) {
    auto speciesRows = speciesByGroup[group];
    std::sort(speciesRows.begin(), speciesRows.end(), [](const auto &lhs,
                                                         const auto &rhs) {
      if (lhs.score != rhs.score) {
        return lhs.score > rhs.score;
      }
      if (lhs.read_support != rhs.read_support) {
        return lhs.read_support > rhs.read_support;
      }
      return lhs.species < rhs.species;
    });
    if (speciesRows.size() > config.local_resolution_species_per_group) {
      speciesRows.resize(config.local_resolution_species_per_group);
    }
    panel.selected_species += speciesRows.size();
    for (const auto &row : speciesRows) {
      selectedSpecies.push_back(row.species);
      selectedSpeciesScores[row.species] = row.score;
    }
  }
  const auto rawTargetsBySpecies =
      repMetadata.load_targets_many(selectedSpecies,
                                    std::numeric_limits<size_t>::max());
  std::unordered_map<uint32_t, std::vector<chimera::local_resolution::TargetRep>>
      targetsBySpecies;
  targetsBySpecies.reserve(rawTargetsBySpecies.size());
  for (const auto &[species, rows] : rawTargetsBySpecies) {
    auto selectedRows = select_local_resolution_source_targets(
        rows, config.local_resolution_targets_per_species);
    if (!selectedRows.empty()) {
      targetsBySpecies.emplace(species, std::move(selectedRows));
    }
  }

  std::unordered_map<uint32_t, uint32_t> selectedTargetsByGroup;
  selectedTargetsByGroup.reserve(selectedGroupSet.size());

  auto admitTarget =
      [&](const chimera::local_resolution::TargetRep &row) {
        const uint64_t targetAnchorRecords = row.anchor_count;
        const uint64_t targetAnchorBytes =
            row.anchor_byte_size == 0
                ? targetAnchorRecords * localResolutionAnchorRecordBytes
                : row.anchor_byte_size;
        const uint32_t currentGroupTargets = selectedTargetsByGroup[row.genus];
        if (currentGroupTargets >= config.local_resolution_max_targets_per_group) {
          ++panel.group_cap_skipped_targets;
          panel.group_cap_skipped_anchor_records += targetAnchorRecords;
          panel.group_cap_skipped_anchor_bytes += targetAnchorBytes;
          return;
        }
        if (config.local_resolution_max_anchor_bytes > 0 &&
            panel.selected_anchor_bytes + targetAnchorBytes >
                config.local_resolution_max_anchor_bytes) {
          ++panel.budget_skipped_targets;
          panel.budget_skipped_anchor_records += targetAnchorRecords;
          panel.budget_skipped_anchor_bytes += targetAnchorBytes;
          return;
        }
        panel.targets.push_back(ChimeraClassify::LocalResolutionTarget{
            row.genus, row.species, row.target_len, row.anchor_count,
            row.anchor_byte_offset, row.anchor_byte_size, row.target_name});
        ++panel.selected_targets;
        panel.selected_anchor_records += targetAnchorRecords;
        panel.selected_anchor_bytes += targetAnchorBytes;
        selectedTargetsByGroup[row.genus] = currentGroupTargets + 1;
      };

  if (config.local_resolution_max_anchor_bytes == 0) {
    for (uint32_t species : selectedSpecies) {
      const auto found = targetsBySpecies.find(species);
      if (found == targetsBySpecies.end()) {
        continue;
      }
      for (const auto &row : found->second) {
        admitTarget(row);
      }
    }
    finalize_local_resolution_panel_shadow(panel);
    return panel;
  }

  std::vector<uint32_t> targetAdmissionSpecies = selectedSpecies;
  std::sort(targetAdmissionSpecies.begin(), targetAdmissionSpecies.end(),
            [&](uint32_t lhs, uint32_t rhs) {
              const double lhsScore = selectedSpeciesScores[lhs];
              const double rhsScore = selectedSpeciesScores[rhs];
              if (lhsScore != rhsScore) {
                return lhsScore > rhsScore;
              }
              return lhs < rhs;
            });

  size_t maxTargetsPerSpecies = 0;
  for (uint32_t species : targetAdmissionSpecies) {
    const auto found = targetsBySpecies.find(species);
    if (found != targetsBySpecies.end()) {
      maxTargetsPerSpecies = std::max(maxTargetsPerSpecies,
                                      found->second.size());
    }
  }
  for (size_t targetRank = 0; targetRank < maxTargetsPerSpecies; ++targetRank) {
    for (uint32_t species : targetAdmissionSpecies) {
      const auto found = targetsBySpecies.find(species);
      if (found == targetsBySpecies.end() || targetRank >= found->second.size()) {
        continue;
      }
      admitTarget(found->second[targetRank]);
    }
  }
  finalize_local_resolution_panel_shadow(panel);
  return panel;
}

static LocalResolutionCallMap make_local_resolution_call_map(
    const ChimeraClassify::LocalResolutionResult &localResult) {
  LocalResolutionCallMap out;
  out.reserve(localResult.reads.size());
  for (const auto &call : localResult.reads) {
    if (!call.read_id.empty()) {
      out.emplace(call.read_id, call);
    }
  }
  return out;
}

static bool apply_local_resolution_result(
    ChimeraClassify::classifyResult &result,
    const LocalResolutionCallMap *localCalls, double sampleDivergence,
    double divergenceThreshold) {
  if (localCalls == nullptr || sampleDivergence < divergenceThreshold) {
    return false;
  }
  const auto found = localCalls->find(result.id);
  if (found == localCalls->end()) {
    return false;
  }
  const auto &call = found->second;
  if (call.candidates.empty()) {
    std::string hint;
    if (!result.taxidCount.empty()) {
      hint = result.taxidCount.front().first;
    }
    result.taxidCount.clear();
    result.taxidCount.emplace_back("unclassified", 1.0);
    result.posteriors.clear();
    result.reject_reason = "local_resolution_absent";
    if (!hint.empty() && hint != "unclassified") {
      result.best_taxid_hint = hint;
    }
    return true;
  }

  const auto &top = call.candidates.front();
  if (top.taxid.empty() || top.taxid == "0") {
    return false;
  }
  result.taxidCount.clear();
  result.taxidCount.emplace_back(top.taxid, static_cast<double>(top.score));
  result.posteriors.clear();
  const double denom = std::max<uint32_t>(1, top.score);
  const size_t topk = std::min<size_t>(16, call.candidates.size());
  result.posteriors.reserve(topk);
  for (size_t i = 0; i < topk; ++i) {
    result.posteriors.emplace_back(
        call.candidates[i].taxid,
        static_cast<double>(call.candidates[i].score) /
            static_cast<double>(denom));
  }
  result.reject_reason.clear();
  return true;
}

static void cleanup_part_paths(const std::vector<std::string> &partPaths) {
  for (const auto &path : partPaths) {
    std::error_code ec;
    std::filesystem::remove(path, ec);
  }
}

template <typename ChunkConsumer>
static void for_each_spool_postem_chunk(
    const std::string &spoolPath, const SpoolEMFit &fit,
    const SpoolSampleMixtureFit &sampleMixtureFit,
    const ChimeraClassify::EMOptions &options,
    const ChimeraClassify::DecisionConfig &decisionConfig,
    const ChimeraClassify::TaxDict &tax,
    const ChimeraClassify::PresenceDecision *presenceDecision,
    bool keepRecords, ChunkConsumer consume_chunk) {
  constexpr size_t kSpoolOutputChunkSize = 4096;
  std::vector<ChimeraClassify::SpoolReadRecord> chunkRecords;
  std::vector<ChimeraClassify::classifyResult> chunk;
  if (keepRecords) {
    chunkRecords.reserve(kSpoolOutputChunkSize);
  }
  chunk.reserve(kSpoolOutputChunkSize);

  auto flush_chunk = [&]() {
    if (chunk.empty()) {
      return;
    }
    ChimeraClassify::postEmDecision(chunk, decisionConfig, fit.class_weights,
                                    tax, presenceDecision);
    consume_chunk(keepRecords ? &chunkRecords : nullptr, chunk);
    if (keepRecords) {
      chunkRecords.clear();
    }
    chunk.clear();
  };

  for_each_spool_record(
      spoolPath, [&](const ChimeraClassify::SpoolReadRecord &record) {
        if (keepRecords) {
          chunkRecords.push_back(record);
        }
        chunk.push_back(materialize_spool_result(record, fit, sampleMixtureFit,
                                                tax, options));
        if (chunk.size() >= kSpoolOutputChunkSize) {
          flush_chunk();
        }
      });
  flush_chunk();
}

static LocalResolutionPostTopkScores collect_local_resolution_post_topk_scores(
    const std::vector<std::string> &spoolPaths, const SpoolEMFit &fit,
    const SpoolSampleMixtureFit &sampleMixtureFit,
    const ChimeraClassify::EMOptions &options,
    const ChimeraClassify::DecisionConfig &decisionConfig,
    const ChimeraClassify::TaxDict &tax,
    const ChimeraClassify::PresenceDecision *presenceDecision,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  LocalResolutionPostTopkScores scores;
  for (const auto &spoolPath : spoolPaths) {
    for_each_spool_postem_chunk(
        spoolPath, fit, sampleMixtureFit, options, decisionConfig, tax,
        presenceDecision, false,
        [&](const std::vector<ChimeraClassify::SpoolReadRecord> *,
            std::vector<ChimeraClassify::classifyResult> &chunk) {
          for (auto &result : chunk) {
            ++scores.rows;
            ChimeraClassify::readout::apply_selective_readout(result,
                                                              ncbiTaxdump);
            if (result.posteriors.empty()) {
              continue;
            }
            ++scores.rows_with_post_topk;
            uint64_t positive_items = 0;
            double top1 = 0.0;
            std::unordered_set<uint32_t> seenSpecies;
            seenSpecies.reserve(result.posteriors.size());
            for (const auto &[taxidText, weight] : result.posteriors) {
              ++scores.post_topk_items;
              if (!(weight > 0.0)) {
                continue;
              }
              ++positive_items;
              if (weight > top1) {
                top1 = weight;
              }
              uint32_t taxid = 0;
              if (!chimera::utils::try_parse_u32(taxidText, taxid) ||
                  taxid == 0) {
                continue;
              }
              uint32_t species = taxid;
              if (ncbiTaxdump && ncbiTaxdump->enabled()) {
                species = ncbiTaxdump->to_species(taxid);
              }
              if (species == 0) {
                continue;
              }
              ++scores.usable_post_topk_items;
              scores.species_scores[species] += weight;
              seenSpecies.insert(species);
            }
            scores.top1_sum += top1;
            if (positive_items == 1) {
              ++scores.singleton_rows;
            }
            if (positive_items > 1 || top1 < 0.95) {
              ++scores.ambiguous_rows;
            }
            for (uint32_t species : seenSpecies) {
              ++scores.species_read_support[species];
            }
          }
        });
  }
  return scores;
}

static LocalResolutionEligibility derive_local_resolution_eligibility(
    const LocalResolutionPostTopkScores &scores) {
  LocalResolutionEligibility eligibility;
  if (scores.rows_with_post_topk == 0) {
    return eligibility;
  }
  const double denom = static_cast<double>(scores.rows_with_post_topk);
  eligibility.mean_post_items =
      static_cast<double>(scores.usable_post_topk_items) / denom;
  eligibility.singleton_rate =
      static_cast<double>(scores.singleton_rows) / denom;
  eligibility.mean_top1 = scores.top1_sum / denom;
  eligibility.ambiguous_rate =
      static_cast<double>(scores.ambiguous_rows) / denom;

  eligibility.hard_ambiguity =
      eligibility.mean_post_items >= 3.0 &&
      eligibility.singleton_rate <= 0.5 &&
      eligibility.mean_top1 <= 0.85 &&
      eligibility.ambiguous_rate >= 0.5;
  return eligibility;
}

static std::vector<std::string> make_output_part_paths(
    const std::string &outputFile, size_t partCount) {
  std::vector<std::string> partPaths(partCount);
  for (size_t i = 0; i < partCount; ++i) {
    partPaths[i] = outputFile + ".part." + std::to_string(i) + ".tmp";
  }
  return partPaths;
}

template <typename PartWriter>
static void write_spool_parts_parallel(
    const std::vector<std::string> &spoolPaths,
    const std::vector<std::string> &partPaths, PartWriter write_part) {
  std::vector<std::exception_ptr> partErrors(spoolPaths.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::int64_t part_i = 0;
       part_i < static_cast<std::int64_t>(spoolPaths.size()); ++part_i) {
    const size_t part_idx = static_cast<size_t>(part_i);
    try {
      write_part(part_idx);
    } catch (...) {
      partErrors[part_idx] = std::current_exception();
    }
  }

  for (const auto &error : partErrors) {
    if (error) {
      cleanup_part_paths(partPaths);
      std::rethrow_exception(error);
    }
  }
}

static void write_spool_output_part(
    const std::string &spoolPath, const std::string &partPath,
    const SpoolEMFit &fit,
    const SpoolSampleMixtureFit &sampleMixtureFit,
    const ChimeraClassify::EMOptions &options,
    const ChimeraClassify::DecisionConfig &decisionConfig,
    const ChimeraClassify::TaxDict &tax,
    const ChimeraClassify::PresenceDecision *presenceDecision,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump,
    const LocalResolutionCallMap *localCalls, double sampleDivergence,
    double divergenceThreshold,
    SpoolOutputPartStats &partStats,
    EvidenceAggregateMap &postemPrimaryEvidenceAggregate,
    EvidenceAggregateMap &postemDecisionEvidenceAggregate) {
  std::ofstream partOs(partPath, std::ios::out | std::ios::binary);
  if (!partOs.is_open()) {
    throw std::runtime_error("Failed to open classify output part: " +
                             partPath);
  }
  std::vector<char> partBuffer(1 << 20, '\0');
  partOs.rdbuf()->pubsetbuf(
      partBuffer.data(), static_cast<std::streamsize>(partBuffer.size()));
  std::ostringstream postTopkOss;

  for_each_spool_postem_chunk(
      spoolPath, fit, sampleMixtureFit, options, decisionConfig, tax,
      presenceDecision, false,
      [&](const std::vector<ChimeraClassify::SpoolReadRecord> *,
          std::vector<ChimeraClassify::classifyResult> &chunk) {
        for (auto &result : chunk) {
          for (const auto &[taxid, count] : result.taxidCount) {
            if (!taxid.empty() && taxid != "unclassified" && count > 0.0) {
              postemPrimaryEvidenceAggregate[taxid] += count;
            }
          }
          if (!result.taxidCount.empty()) {
            const std::string &decisionTaxid = result.taxidCount.front().first;
            if (!decisionTaxid.empty() && decisionTaxid != "unclassified") {
              postemDecisionEvidenceAggregate[decisionTaxid] += 1.0;
            }
          }
          ChimeraClassify::readout::apply_selective_readout(result,
                                                            ncbiTaxdump);
          apply_local_resolution_result(result, localCalls, sampleDivergence,
                                        divergenceThreshold);
          if (!result.taxidCount.empty() &&
              result.taxidCount.front().first == "unclassified") {
            ++partStats.unclassified;
          } else {
            ++partStats.classified;
          }
          ChimeraClassify::writeResultRecord(partOs, result, postTopkOss);
        }
      });
  partOs.close();
  if (!partOs.good()) {
    throw std::runtime_error("Failed to close classify output part: " +
                             partPath);
  }
}

static void merge_classify_output_parts(
    const std::string &outputFile, const std::vector<std::string> &partPaths,
    const std::vector<SpoolOutputPartStats> &partStats,
    ChimeraClassify::FileInfo &fileInfo) {
  std::ofstream os(outputFile, std::ios::out | std::ios::binary);
  if (!os.is_open()) {
    throw std::runtime_error("Failed to open file: " + outputFile);
  }
  std::vector<char> outputBuffer(1 << 20, '\0');
  os.rdbuf()->pubsetbuf(outputBuffer.data(),
                        static_cast<std::streamsize>(outputBuffer.size()));
  fileInfo.classifiedNum = 0;
  fileInfo.unclassifiedNum = 0;
  for (size_t i = 0; i < partPaths.size(); ++i) {
    std::ifstream partIs(partPaths[i], std::ios::in | std::ios::binary);
    if (!partIs.is_open()) {
      throw std::runtime_error("Failed to open classify output part: " +
                               partPaths[i]);
    }
    std::error_code partSizeEc;
    if (std::filesystem::file_size(partPaths[i], partSizeEc) == 0 &&
        !partSizeEc) {
      continue;
    }
    os << partIs.rdbuf();
    if (!os.good()) {
      throw std::runtime_error("Failed to append classify output part: " +
                               partPaths[i]);
    }
    fileInfo.classifiedNum += partStats[i].classified;
    fileInfo.unclassifiedNum += partStats[i].unclassified;
  }
  os.close();
}

static const char *presence_state(double logPosterior, double threshold) {
  if (logPosterior >= threshold) {
    return "accepted";
  }
  if (logPosterior <= -threshold) {
    return "rejected";
  }
  return "unknown";
}

static void write_presence_evidence(
    const std::string &presencePath,
    const ChimeraClassify::PresenceEvidenceTable &table,
    const ChimeraClassify::TaxDict &tax) {
  std::filesystem::path parent =
      std::filesystem::path(presencePath).parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::vector<const ChimeraClassify::PresenceEvidenceRow *> rows;
  rows.reserve(table.rows.size());
  for (const auto &row : table.rows) {
    rows.push_back(&row);
  }
  std::sort(rows.begin(), rows.end(), [&tax](const auto *lhs, const auto *rhs) {
    const std::string lt = spool_taxid_to_string(lhs->tid, tax);
    const std::string rt = spool_taxid_to_string(rhs->tid, tax);
    return lt < rt;
  });

  std::ofstream os(presencePath, std::ios::out | std::ios::binary);
  if (!os.is_open()) {
    throw std::runtime_error("Failed to open presence evidence output: " +
                             presencePath);
  }
  os << "taxid\tpresence_state\tpresence_posterior\tpresence_log_odds"
     << "\tpresence_log_bf\tscore\tunique_score\thits\tunique_hits"
     << "\tread_hits\tunique_reads\tunique_obs\tbreadth_obs"
     << "\tbreadth_ratio\tunique_effective\tlocal_factor\texposure"
     << "\tunique_reference\ttotal_signatures\tgenome_length"
     << "\tunique_density\texpected_unique_per_ref_read\n";
  os << std::setprecision(17);
  for (const auto *row : rows) {
    const std::string taxid = spool_taxid_to_string(row->tid, tax);
    if (taxid.empty() || taxid == "unclassified") {
      continue;
    }
    os << taxid << '\t' << presence_state(row->log_posterior, table.threshold)
       << '\t' << row->posterior << '\t' << row->log_posterior << '\t'
       << row->log_bf << '\t' << row->score << '\t' << row->unique_score
       << '\t' << row->hits << '\t' << row->unique_hits << '\t'
       << row->read_hits << '\t' << row->unique_reads << '\t'
       << row->unique_obs << '\t' << row->breadth_obs << '\t'
       << row->breadth_ratio << '\t' << row->unique_effective << '\t'
       << row->local_factor << '\t' << row->exposure << '\t'
       << row->unique_reference << '\t' << row->total_signatures << '\t'
       << row->genome_length << '\t' << row->unique_density << '\t'
       << row->expected_unique_per_ref_read << '\n';
  }
  os.close();
  if (!os.good()) {
    throw std::runtime_error("Failed to close presence evidence output: " +
                             presencePath);
  }
}

static ChimeraClassify::PresenceDecision presence_decision_from_table(
    const ChimeraClassify::PresenceEvidenceTable &table) {
  ChimeraClassify::PresenceDecision decision;
  decision.threshold = table.threshold;
  decision.posteriors.reserve(table.rows.size());
  decision.logPosteriors.reserve(table.rows.size());
  for (const auto &row : table.rows) {
    decision.posteriors[row.tid] = row.posterior;
    decision.logPosteriors[row.tid] = row.log_posterior;
  }
  return decision;
}

static EvidenceAggregateMap merge_evidence_aggregates(
    const std::vector<EvidenceAggregateMap> &parts) {
  EvidenceAggregateMap merged;
  for (const auto &partMap : parts) {
    for (const auto &[taxid, count] : partMap) {
      merged[taxid] += count;
    }
  }
  return merged;
}

static void write_spool_em_results(
    const std::vector<std::string> &spoolPaths, const SpoolEMFit &fit,
    const SpoolSampleMixtureFit &sampleMixtureFit,
    const ChimeraClassify::EMOptions &options,
    const ChimeraClassify::DecisionConfig &decisionConfig,
    const ChimeraClassify::TaxDict &tax,
    const ChimeraClassify::PresenceDecision *presenceDecision,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump,
    const ChimeraClassify::ClassifyConfig &config,
    const LocalResolutionCallMap *localCalls, double sampleDivergence,
    ChimeraClassify::FileInfo &fileInfo) {
  const std::string outputFile = resolve_tsv_output_path(config.outputFile);

  const size_t part_count = spoolPaths.size();
  std::vector<std::string> partPaths =
      make_output_part_paths(outputFile, part_count);
  const std::string evidencePath = resolve_evidence_output_path(outputFile);
  {
    std::filesystem::path evidenceParent =
        std::filesystem::path(evidencePath).parent_path();
    if (!evidenceParent.empty()) {
      std::filesystem::create_directories(evidenceParent);
    }
  }
  const ChimeraClassify::AbundanceAutoPolicy abundancePolicy =
      ChimeraClassify::derive_abundance_auto_policy(config, options);
  std::cout << "[classify][auto] abundance-em"
            << " prune_ratio=" << abundancePolicy.prune_ratio
            << " iterations=" << abundancePolicy.iterations
            << " dispersion_s=" << config.community_dispersion_s << "\n";
  SpoolAbundanceFit abundanceFit = fit_spool_abundance_em(
      spoolPaths, tax, abundancePolicy.iterations, options,
      abundancePolicy.prune_ratio);
  std::vector<EvidenceAggregateMap> postemPrimaryEvidenceAggregates(part_count);
  std::vector<EvidenceAggregateMap> postemDecisionEvidenceAggregates(part_count);
  std::vector<SpoolOutputPartStats> partStats(part_count);
  write_spool_parts_parallel(
      spoolPaths, partPaths, [&](size_t part_idx) {
        write_spool_output_part(
            spoolPaths[part_idx], partPaths[part_idx], fit, sampleMixtureFit,
            options, decisionConfig, tax, presenceDecision, ncbiTaxdump,
            localCalls, sampleDivergence,
            config.local_resolution_divergence_threshold,
            partStats[part_idx],
            postemPrimaryEvidenceAggregates[part_idx],
            postemDecisionEvidenceAggregates[part_idx]);
      });

  try {
    merge_classify_output_parts(outputFile, partPaths, partStats, fileInfo);
  } catch (...) {
    cleanup_part_paths(partPaths);
    throw;
  }
  EvidenceAggregateMap postemPrimaryEvidenceMass =
      merge_evidence_aggregates(postemPrimaryEvidenceAggregates);
  EvidenceAggregateMap postemDecisionEvidenceMass =
      merge_evidence_aggregates(postemDecisionEvidenceAggregates);
  // ChimeraEvidence.tsv is the sample-level profiling evidence head.  It is
  // intentionally aggregated before final selective readout rewrites per-read
  // labels, so it must not be interpreted as a histogram of ChimeraClassify.tsv.
  write_abundance_em_evidence(evidencePath, abundanceFit,
                              postemPrimaryEvidenceMass,
                              postemDecisionEvidenceMass, tax,
                              abundancePolicy.abundance_weight);
  cleanup_part_paths(partPaths);
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

static ChimeraClassify::CommunityDispersionStats
run_community_dispersion_calibration(
    ChimeraBuild::IMCFConfig &imcfConfig,
    const ChimeraClassify::ClassifyConfig &config,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    const ChimeraClassify::TaxDict &tax,
    const chimera::feature::Params &feature_params, size_t feature_min_len,
    const ChimeraClassify::WeightingContext &weightCtx) {
  constexpr uint32_t kDispersionCalibrationReads = 200000;
  ChimeraClassify::FileInfo calibrationInfo;
  std::vector<moodycamel::ConcurrentQueue<ChimeraClassify::batchReads>>
      calibrationQueues(
      static_cast<size_t>(std::max<uint16_t>(1, config.threads)));
  std::vector<ChimeraClassify::QueueThrottle> calibrationThrottles(
      calibrationQueues.size());
  std::vector<ChimeraClassify::classifyResult> calibrationResults;
  std::atomic<bool> calibration_done{false};
  std::thread calibrationProducer([&]() {
    ChimeraClassify::parseReads(calibrationQueues, config, calibrationInfo,
                                kDispersionCalibrationReads,
                                &calibrationThrottles);
    calibration_done.store(true, std::memory_order_release);
  });

  ChimeraClassify::ClassifyConfig calibrationConfig = config;
  calibrationConfig.sample_state_calibration = true;
  ChimeraClassify::classify_streaming(
      imcfConfig, calibrationQueues, calibrationConfig, imcf, tax,
      calibrationResults, calibrationInfo, calibration_done, feature_params,
      feature_min_len, weightCtx, nullptr, &calibrationThrottles);
  calibrationProducer.join();

  if (!calibrationResults.empty()) {
    std::sort(calibrationResults.begin(), calibrationResults.end(),
              [](const ChimeraClassify::classifyResult &a,
                 const ChimeraClassify::classifyResult &b) {
                return a.id < b.id;
              });
  }

  bool calibrationEm = false;
  if (!calibrationResults.empty()) {
    ChimeraClassify::EMOptions options;
    options.temp = 1.05;
    options.prune_ratio = config.em_prune_ratio;
    options.conf_power = config.em_conf_power;
    auto [posterior, weights] =
        ChimeraClassify::EMAlgorithm(std::move(calibrationResults),
                                     config.emIter, 0.0, options, nullptr);
    calibrationResults = std::move(posterior);
    calibrationEm = true;
  }

  std::unordered_map<std::string, double> species_counts;
  species_counts.reserve(calibrationResults.size() / 4 + 8);
  std::size_t unclassified_reads = 0;

  for (const auto &res : calibrationResults) {
    std::string top;
    if (calibrationEm && !res.posteriors.empty()) {
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
  return ChimeraClassify::compute_community_dispersion_stats(
      counts, unclassified_reads);
}

static std::unordered_map<std::string, double> build_em_prior_scale(
    const chimera::presence::CoverageMeta &coverageMeta,
    const ChimeraClassify::PresenceDecision &presenceDecision,
    const ChimeraClassify::TaxDict &tax,
    const ChimeraClassify::FileInfo &fileInfo) {
  std::unordered_map<std::string, double> emPriorScale;
  if (coverageMeta.entries.empty()) {
    return emPriorScale;
  }
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
    } else if (entry.expected_unique_per_ref_read > 0.0 && window_ref > 0.0) {
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
  return emPriorScale;
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
  const std::optional<LocalResolutionArtifacts> resolvedLocalArtifacts =
      config.local_resolution_enabled
          ? resolve_local_resolution_artifacts(config.dbFile)
          : std::nullopt;
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
    rebuild_bin_slot_rep_lookup(tax, weightCtx.tid2speciesRep);
  }
  std::vector<std::vector<std::string>>().swap(indexToTaxid);
  PresenceSummary presenceSummary(config.presence_breadth_bits);
  PresenceSummary *presencePtr = &presenceSummary;

  config.community_dispersion_u = 1.0;
  config.community_dispersion_s = 1.0;

  CommunityDispersionStats dispersionStats = run_community_dispersion_calibration(
      imcfConfig, config, imcf, tax, feature_params, feature_min_len,
      weightCtx);

  const CommunityAutoPolicy communityPolicy =
      derive_community_auto_policy(config, dispersionStats);
  config.community_dispersion_u = communityPolicy.community_dispersion_u;
  config.community_dispersion_s = communityPolicy.community_dispersion_s;
  config.firstFilterBeta = communityPolicy.first_filter_beta;
  config.em_conf_power = communityPolicy.em_conf_power;
  std::cout << "[classify][auto] community-dispersion"
            << " u=" << std::fixed << std::setprecision(4)
            << config.community_dispersion_u
            << " s=" << config.community_dispersion_s
            << " top_mass=" << dispersionStats.top_mass
            << " eff_species=" << dispersionStats.eff_species
            << " simpson_species=" << dispersionStats.simpson_species
            << " uncls=" << dispersionStats.unclassified << "\n";
  std::cout << "[classify][auto] firstFilterBeta=" << config.firstFilterBeta
            << " em_conf_power=" << config.em_conf_power
            << " dispersion_gate=" << config.community_dispersion_s << "\n";
  std::cout << std::defaultfloat;

  FileInfo fileInfo;
  seqan3::contrib::bgzf_thread_count = config.threads;
  std::vector<moodycamel::ConcurrentQueue<batchReads>> readQueues(
      static_cast<size_t>(std::max<uint16_t>(1, config.threads)));
  std::vector<QueueThrottle> queueThrottles(readQueues.size());
  std::unordered_map<std::string, double> classWeights;
  bool posteriorModelUsed = false;
  const std::filesystem::path spoolDir = make_spool_dir(config.outputFile);
  const bool keepClassifySpool =
      env_flag_enabled("CHIMERA_KEEP_CLASSIFY_SPOOL");
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
                           feature_min_len, weightCtx, presencePtr,
                           &queueThrottles);
  producer.join();
  std::vector<moodycamel::ConcurrentQueue<batchReads>>().swap(readQueues);
  if (fileInfo.sequenceNum > 0) {
    fileInfo.avgLen = fileInfo.bpLength / fileInfo.sequenceNum;
  }
  size_t presenceTotalReads = fileInfo.sequenceNum;
  size_t presenceMeanReadLen = fileInfo.avgLen;

  const PresenceEvidenceTable presenceEvidence =
      build_presence_evidence_table(presenceSummary, tax, config, coverageMeta,
                                    presenceTotalReads, presenceMeanReadLen);
  PresenceDecision presenceDecision =
      presence_decision_from_table(presenceEvidence);
  write_presence_evidence(
      resolve_presence_output_path(resolve_tsv_output_path(config.outputFile)),
      presenceEvidence, tax);
  PresenceSummary().stats.swap(presenceSummary.stats);
  presenceSummary.sketchBits = 0;
  presenceSummary.sketchWords = 0;
  presencePtr = nullptr;

  std::unordered_map<std::string, double> emPriorScale =
      build_em_prior_scale(coverageMeta, presenceDecision, tax, fileInfo);
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
    const PostPiAutoPolicy postPiPolicy = derive_post_pi_auto_policy(
        config, fileInfo, coverageMeta, dispersionStats, classWeights);
    decisionConfig.min_class_weight = postPiPolicy.tuned;
    std::cout << "[classify][auto] post-pi"
              << " s=" << std::fixed << std::setprecision(4)
              << config.community_dispersion_s
              << " dominance=" << postPiPolicy.dominance
              << " relax=" << postPiPolicy.relax
              << " base=" << config.post_pi_min
              << " low_weight_mass=" << postPiPolicy.low_weight_mass
              << " low_weight_target=" << postPiPolicy.low_weight_target
              << " tuned=" << postPiPolicy.tuned
              << " evidence_t=" << postPiPolicy.evidence_t
              << " avg_len=" << fileInfo.avgLen
              << " evidence_strength="
              << postPiPolicy.sample_evidence_strength << "\n";
    std::cout << std::defaultfloat;
    SpoolSampleMixtureFit sampleMixtureFit =
        fit_spool_sample_mixture_em(spoolPaths, config.emIter);
    std::cout << "[classify][auto] sample-mixture"
              << " reads=" << sampleMixtureFit.read_count
              << " candidates=" << sampleMixtureFit.candidate_count
              << " taxa=" << sampleMixtureFit.taxid_list.size() << "\n";

    LocalResolutionCallMap localResolutionCalls;
    const LocalResolutionCallMap *localResolutionCallPtr = nullptr;
    double localResolutionDivergence = config.community_dispersion_s;
    const std::string localProfileOutput =
        resolve_tsv_output_path(config.outputFile);
    if (config.local_resolution_enabled && config.pairedFiles.empty()) {
      const auto postTopkScores = collect_local_resolution_post_topk_scores(
          spoolPaths, speciesFit, sampleMixtureFit, options, decisionConfig,
          tax, &presenceDecision, weightCtx.ncbiTaxdump);
      const LocalResolutionEligibility localEligibility =
          derive_local_resolution_eligibility(postTopkScores);
      const bool localMayChangeOutput =
          localResolutionDivergence >=
              config.local_resolution_divergence_threshold &&
          localEligibility.hard_ambiguity;
      if (!localMayChangeOutput) {
        write_local_resolution_profile_json(
            localProfileOutput,
            localResolutionDivergence <
                    config.local_resolution_divergence_threshold
                ? "skipped_low_sample_divergence"
                : "skipped_not_hard_ambiguity",
            localResolutionDivergence,
            config.local_resolution_divergence_threshold, nullptr, nullptr,
            0.0, 0.0, 0.0);
        std::cout << "[classify][local-resolution]"
                  << " skipped="
                  << (localResolutionDivergence <
                              config.local_resolution_divergence_threshold
                          ? "low_sample_divergence"
                          : "not_hard_ambiguity")
                  << " sample_divergence=" << std::fixed
                  << std::setprecision(4) << localResolutionDivergence
                  << " divergence_threshold="
                  << config.local_resolution_divergence_threshold
                  << " mean_post_items="
                  << localEligibility.mean_post_items
                  << " singleton_rate="
                  << localEligibility.singleton_rate
                  << " mean_top1=" << localEligibility.mean_top1
                  << " ambiguous_rate="
                  << localEligibility.ambiguous_rate << "\n";
        std::cout << std::defaultfloat;
      } else {
        if (resolvedLocalArtifacts.has_value()) {
          const std::filesystem::path localIndexPath =
              resolvedLocalArtifacts->index_path;
          const std::filesystem::path repMetadataPath =
              resolvedLocalArtifacts->rep_metadata_path;
          if (!std::filesystem::exists(repMetadataPath)) {
            throw std::runtime_error(
                "Local resolution metadata is missing next to database: " +
                repMetadataPath.string());
          }
          const auto metadataStarted = std::chrono::steady_clock::now();
          const auto repMetadata =
              chimera::local_resolution::RepMetadata::open(repMetadataPath);
          const auto metadataLoaded = std::chrono::steady_clock::now();
          const LocalResolutionPanel panel = build_local_resolution_panel(
              postTopkScores.species_scores,
              postTopkScores.species_read_support, repMetadata, config,
              resolvedLocalArtifacts->k);
          const auto panelBuilt = std::chrono::steady_clock::now();
          const double metadataSeconds =
              std::chrono::duration<double>(metadataLoaded - metadataStarted)
                  .count();
          const double panelSeconds =
              std::chrono::duration<double>(panelBuilt - metadataLoaded)
                  .count();
          if (panel.selected_targets > 0) {
            ChimeraClassify::LocalResolutionRequest localRequest;
            localRequest.read_files = config.singleFiles;
            localRequest.index_file = localIndexPath.string();
            localRequest.targets = panel.targets;
            localRequest.diag_bin = config.lpc_diag_bin;
            localRequest.max_occ = config.lpc_max_occ;
            localRequest.min_chain = config.lpc_min_chain;
            localRequest.threads = config.threads;
            const auto started = std::chrono::steady_clock::now();
            ChimeraClassify::LocalResolutionResult localResult =
                run_local_resolution_engine(localRequest);
            localResolutionCalls = make_local_resolution_call_map(localResult);
            localResolutionCallPtr = &localResolutionCalls;
            const auto finished = std::chrono::steady_clock::now();
            const double seconds =
                std::chrono::duration<double>(finished - started).count();
            std::cout << "[classify][local-resolution]"
                      << " sample_divergence=" << std::fixed
                      << std::setprecision(4) << localResolutionDivergence
                      << " divergence_threshold="
                      << config.local_resolution_divergence_threshold
                      << " hard_ambiguity="
                      << (localEligibility.hard_ambiguity ? 1 : 0)
                      << " mean_post_items="
                      << localEligibility.mean_post_items
                      << " singleton_rate="
                      << localEligibility.singleton_rate
                      << " mean_top1=" << localEligibility.mean_top1
                      << " ambiguous_rate="
                      << localEligibility.ambiguous_rate
                      << " groups=" << panel.selected_groups
                      << " species=" << panel.selected_species
                      << " targets=" << panel.selected_targets
                      << " post_topk_rows="
                      << postTopkScores.rows_with_post_topk
                      << " post_topk_items="
                      << postTopkScores.usable_post_topk_items
                      << " panel_anchor_records="
                      << panel.selected_anchor_records
                      << " panel_anchor_bytes="
                      << panel.selected_anchor_bytes
                      << " species_source_pairs="
                      << panel.selected_species_source_pairs
                      << " source_extra_targets="
                      << panel.selected_source_extra_targets
                      << " source_cap1_targets="
                      << panel.source_cap1_targets
                      << " source_cap1_anchor_bytes="
                      << panel.source_cap1_anchor_bytes
                      << " source_cap2_targets="
                      << panel.source_cap2_targets
                      << " source_cap2_anchor_bytes="
                      << panel.source_cap2_anchor_bytes
                      << " targets_per_species="
                      << panel.targets_per_species
                      << " max_targets_per_group="
                      << panel.max_targets_per_group
                      << " anchor_byte_budget="
                      << panel.anchor_byte_budget
                      << " group_cap_skipped_targets="
                      << panel.group_cap_skipped_targets
                      << " group_cap_skipped_anchor_bytes="
                      << panel.group_cap_skipped_anchor_bytes
                      << " budget_skipped_targets="
                      << panel.budget_skipped_targets
                      << " budget_skipped_anchor_bytes="
                      << panel.budget_skipped_anchor_bytes
                      << " reads=" << localResult.stats.reads
                      << " query_hashes=" << localResult.stats.query_hashes
                      << " target_filter=" << localResult.stats.target_filter
                      << " selected_targets="
                      << localResult.stats.selected_targets
                      << " direct_targets=" << localResult.stats.direct_targets
                      << " target_anchor_records_scanned="
                      << localResult.stats.target_anchor_records_scanned
                      << " target_anchor_bytes_read="
                      << localResult.stats.target_anchor_bytes_read
                      << " target_anchor_records_matched="
                      << localResult.stats.target_anchor_records_matched
                      << " target_hash_prefilter_rejects="
                      << localResult.stats.target_hash_prefilter_rejects
                      << " direct_load_batches="
                      << localResult.stats.direct_load_batches
                      << " pread_calls=" << localResult.stats.pread_calls
                      << " pread_bytes=" << localResult.stats.pread_bytes
                      << " raw_chain_records="
                      << localResult.stats.raw_chain_records
                      << " kept_chain_records="
                      << localResult.stats.kept_chain_records
                      << " index_hash_keys="
                      << localResult.stats.index_hash_keys
                      << " overflow_hash_keys="
                      << localResult.stats.overflow_hash_keys
                      << " dropped_broad_keys="
                      << localResult.stats.dropped_broad_keys
                      << " dropped_broad_records="
                      << localResult.stats.dropped_broad_records
                      << " local_hits=" << localResult.stats.local_hits
                      << " local_absent=" << localResult.stats.local_absent
                      << " metadata_seconds=" << metadataSeconds
                      << " panel_seconds=" << panelSeconds
                      << " read_seconds=" << localResult.stats.read_seconds
                      << " ref_seconds=" << localResult.stats.ref_seconds
                      << " target_io_seconds="
                      << localResult.stats.target_io_seconds
                      << " target_read_seconds="
                      << localResult.stats.target_read_seconds
                      << " target_filter_seconds="
                      << localResult.stats.target_filter_seconds
                      << " target_collect_seconds="
                      << localResult.stats.target_collect_seconds
                      << " posting_merge_seconds="
                      << localResult.stats.posting_merge_seconds
                      << " index_finalize_seconds="
                      << localResult.stats.index_finalize_seconds
                      << " chain_seconds=" << localResult.stats.chain_seconds
                      << " seconds=" << seconds << "\n";
            std::cout << std::defaultfloat;
            write_local_resolution_profile_json(
                localProfileOutput, "ran", localResolutionDivergence,
                config.local_resolution_divergence_threshold, &panel,
                &localResult.stats, metadataSeconds, panelSeconds, seconds);
          } else {
            write_local_resolution_profile_json(
                localProfileOutput, "skipped_no_sample_targets",
                localResolutionDivergence,
                config.local_resolution_divergence_threshold, &panel, nullptr,
                metadataSeconds, panelSeconds, 0.0);
            std::cout << "[classify][local-resolution]"
                      << " skipped=no_sample_targets"
                      << " metadata_seconds=" << metadataSeconds
                      << " panel_seconds=" << panelSeconds << "\n";
          }
        } else {
          write_local_resolution_profile_json(
              localProfileOutput, "skipped_no_local_resolution_data",
              localResolutionDivergence,
              config.local_resolution_divergence_threshold, nullptr, nullptr,
              0.0, 0.0, 0.0);
          std::cout << "[classify][local-resolution]"
                    << " skipped=no_local_resolution_data\n";
        }
      }
    } else if (config.local_resolution_enabled && !config.pairedFiles.empty()) {
      write_local_resolution_profile_json(
          localProfileOutput, "skipped_paired_input", localResolutionDivergence,
          config.local_resolution_divergence_threshold, nullptr, nullptr, 0.0,
          0.0, 0.0);
      std::cout << "[classify][local-resolution]"
                << " skipped=paired_input\n";
    }

    write_spool_em_results(spoolPaths, speciesFit, sampleMixtureFit, options,
                           decisionConfig, tax, &presenceDecision,
                           weightCtx.ncbiTaxdump, config, localResolutionCallPtr,
                           localResolutionDivergence, fileInfo);
  }
  if (keepClassifySpool) {
    std::cout << "[classify][debug] keeping spool=" << spoolDir.string()
              << "\n";
  } else {
    std::filesystem::remove_all(spoolDir);
  }
}

} // namespace ChimeraClassify
