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
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
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
maybe_load_ncbi_taxdump_for_profile();
static std::unique_ptr<ChimeraClassify::NcbiTaxdump>
maybe_load_profiledb_taxonomy_for_profile(const std::filesystem::path &dbPath);
static std::unique_ptr<ChimeraClassify::NcbiTaxdump>
maybe_load_tax_tsv_for_profile(const std::filesystem::path &dbPath);

static std::unique_ptr<ChimeraClassify::NcbiTaxdump>
maybe_load_ncbi_taxdump(const std::string &taxonomyKind) {
  if (taxonomyKind != "ncbi") {
    return nullptr;
  }
  return maybe_load_ncbi_taxdump_for_profile();
}

static std::unique_ptr<ChimeraClassify::NcbiTaxdump>
maybe_load_ncbi_taxdump_for_profile() {
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
      tax->scientific_name.resize(static_cast<size_t>(tid) + 1);
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
  const std::filesystem::path names = base / "names.dmp";
  if (std::filesystem::exists(names)) {
    std::ifstream names_is(names);
    if (names_is.is_open()) {
      while (std::getline(names_is, line)) {
        if (line.empty()) {
          continue;
        }
        std::vector<std::string> fields;
        fields.reserve(4);
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, '|')) {
          tok = trim_copy(tok);
          if (tok.empty()) {
            continue;
          }
          fields.push_back(tok);
          if (fields.size() >= 4) {
            break;
          }
        }
        if (fields.size() < 3 || fields.back() != "scientific name") {
          continue;
        }
        uint32_t tid = 0;
        if (!chimera::utils::try_parse_u32(fields[0], tid)) {
          continue;
        }
        if (tid >= tax->scientific_name.size()) {
          tax->scientific_name.resize(static_cast<size_t>(tid) + 1);
        }
        tax->scientific_name[tid] = fields[1];
      }
    }
  }
  return tax;
}

template <typename T>
bool read_profiledb_taxonomy_pod(std::istream &in, T &value) {
  return static_cast<bool>(in.read(reinterpret_cast<char *>(&value),
                                   static_cast<std::streamsize>(sizeof(T))));
}

bool read_profiledb_taxonomy_string(std::istream &in, std::string &value) {
  uint32_t size = 0;
  if (!read_profiledb_taxonomy_pod(in, size)) {
    return false;
  }
  if (size > (1u << 30)) {
    return false;
  }
  value.assign(static_cast<size_t>(size), '\0');
  if (size == 0) {
    return true;
  }
  return static_cast<bool>(
      in.read(value.data(), static_cast<std::streamsize>(size)));
}

std::filesystem::path profiledb_root_for_db(
    const std::filesystem::path &dbPath) {
  if (std::filesystem::is_directory(dbPath)) {
    return dbPath;
  }
  const auto corePath = chimera::local_resolution::core_archive_path_for(dbPath);
  const auto sidecarRoot =
      corePath.parent_path() / (corePath.stem().string() + ".profiledb");
  if (std::filesystem::is_directory(sidecarRoot)) {
    return sidecarRoot;
  }
  return corePath.parent_path();
}

static std::unique_ptr<ChimeraClassify::NcbiTaxdump>
maybe_load_tax_tsv_for_profile(const std::filesystem::path &dbPath) {
  const auto corePath = chimera::local_resolution::core_archive_path_for(dbPath);
  std::filesystem::path taxPath = corePath;
  taxPath.replace_extension(".tax");
  if (!std::filesystem::exists(taxPath)) {
    return nullptr;
  }

  std::ifstream in(taxPath);
  if (!in.is_open()) {
    return nullptr;
  }

  auto tax = std::make_unique<ChimeraClassify::NcbiTaxdump>();
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    std::vector<std::string> fields;
    size_t begin = 0;
    while (begin <= line.size() && fields.size() < 4) {
      const size_t end = line.find('\t', begin);
      if (end == std::string::npos) {
        fields.push_back(line.substr(begin));
        break;
      }
      fields.push_back(line.substr(begin, end - begin));
      begin = end + 1;
    }
    if (fields.size() < 4) {
      continue;
    }
    uint32_t tid = 0;
    uint32_t parent = 0;
    if (!chimera::utils::try_parse_u32(fields[0], tid) ||
        !chimera::utils::try_parse_u32(fields[1], parent)) {
      continue;
    }
    if (tid >= tax->parent.size()) {
      tax->parent.resize(static_cast<size_t>(tid) + 1, 0);
      tax->is_species.resize(static_cast<size_t>(tid) + 1, 0);
      tax->is_genus.resize(static_cast<size_t>(tid) + 1, 0);
      tax->scientific_name.resize(static_cast<size_t>(tid) + 1);
    }
    tax->parent[tid] = parent;
    tax->is_species[tid] = (fields[2] == "species") ? 1 : 0;
    tax->is_genus[tid] = (fields[2] == "genus") ? 1 : 0;
    tax->scientific_name[tid] = fields[3];
  }

  return tax->enabled() ? std::move(tax) : nullptr;
}

static std::unique_ptr<ChimeraClassify::NcbiTaxdump>
maybe_load_profiledb_taxonomy_for_profile(const std::filesystem::path &dbPath) {
  const std::filesystem::path taxonomyDir =
      profiledb_root_for_db(dbPath) / "taxonomy";
  const std::filesystem::path nodesPath = taxonomyDir / "nodes.bin";
  const std::filesystem::path namesPath = taxonomyDir / "names.bin";
  if (!std::filesystem::exists(nodesPath)) {
    return nullptr;
  }

  std::ifstream nodes(nodesPath, std::ios::binary);
  if (!nodes.is_open()) {
    return nullptr;
  }
  uint32_t version = 0;
  uint64_t count = 0;
  if (!read_profiledb_taxonomy_pod(nodes, version) ||
      !read_profiledb_taxonomy_pod(nodes, count) || version != 1) {
    return nullptr;
  }

  auto tax = std::make_unique<ChimeraClassify::NcbiTaxdump>();
  for (uint64_t i = 0; i < count; ++i) {
    uint32_t tid = 0;
    uint32_t parent = 0;
    std::string rank;
    if (!read_profiledb_taxonomy_pod(nodes, tid) ||
        !read_profiledb_taxonomy_pod(nodes, parent) ||
        !read_profiledb_taxonomy_string(nodes, rank)) {
      return nullptr;
    }
    if (tid >= tax->parent.size()) {
      tax->parent.resize(static_cast<size_t>(tid) + 1, 0);
      tax->is_species.resize(static_cast<size_t>(tid) + 1, 0);
      tax->is_genus.resize(static_cast<size_t>(tid) + 1, 0);
      tax->scientific_name.resize(static_cast<size_t>(tid) + 1);
    }
    tax->parent[tid] = parent;
    tax->is_species[tid] = (rank == "species") ? 1 : 0;
    tax->is_genus[tid] = (rank == "genus") ? 1 : 0;
  }

  if (std::filesystem::exists(namesPath)) {
    std::ifstream names(namesPath, std::ios::binary);
    if (names.is_open()) {
      uint32_t nameVersion = 0;
      uint64_t nameCount = 0;
      if (read_profiledb_taxonomy_pod(names, nameVersion) &&
          read_profiledb_taxonomy_pod(names, nameCount) && nameVersion == 1) {
        for (uint64_t i = 0; i < nameCount; ++i) {
          uint32_t tid = 0;
          std::string name;
          if (!read_profiledb_taxonomy_pod(names, tid) ||
              !read_profiledb_taxonomy_string(names, name)) {
            return nullptr;
          }
          if (tid >= tax->scientific_name.size()) {
            tax->scientific_name.resize(static_cast<size_t>(tid) + 1);
          }
          tax->scientific_name[tid] = std::move(name);
        }
      }
    }
  }

  return tax->enabled() ? std::move(tax) : nullptr;
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

struct SeqProfileRow {
  uint32_t species_taxid{0};
  uint32_t genus_taxid{0};
  double assigned_reads{0.0};
  double unique_reads{0.0};
  double assigned_features{0.0};
  double unique_features{0.0};
  double sequence_abundance{0.0};
  double feature_abundance{0.0};
  double effective_length{0.0};
  double effective_callable_signatures{0.0};
  double effective_unique_signatures{0.0};
  double length_normalized_abundance{0.0};
  double sqrt_length_normalized_abundance{0.0};
  double callable_normalized_abundance{0.0};
  double unique_callable_normalized_abundance{0.0};
  double feature_length_normalized_abundance{0.0};
  double feature_sqrt_length_normalized_abundance{0.0};
  uint32_t source_count{0};
};

struct SeqProfileFit {
  std::vector<SeqProfileRow> rows;
  uint64_t input_reads{0};
  uint64_t reads_with_candidates{0};
  uint64_t multi_candidate_reads{0};
  uint64_t candidate_edges{0};
  double assigned_reads{0.0};
  double assigned_features{0.0};
};

enum class PrimaryProfileScale {
  LengthNormalized,
  SqrtLengthNormalized,
  CallableNormalized,
  UniqueCallableNormalized,
  SequenceAbundance,
};

struct SpeciesProfileMasses {
  std::unordered_map<uint32_t, double> assigned_reads;
  std::unordered_map<uint32_t, double> unique_reads;
  std::unordered_map<uint32_t, double> assigned_features;
  std::unordered_map<uint32_t, double> unique_features;
  uint64_t input_reads{0};
  uint64_t reads_with_candidates{0};
  uint64_t multi_candidate_reads{0};
  uint64_t candidate_edges{0};
};

static std::string spool_taxid_to_string(
    uint32_t tid, const ChimeraClassify::TaxDict &tax) {
  if (tid == ChimeraClassify::kSpoolUnclassifiedTid ||
      tid >= tax.id2str.size()) {
    return "unclassified";
  }
  return tax.id2str[tid];
}

static uint32_t spool_tid_to_taxid(uint32_t tid,
                                   const ChimeraClassify::TaxDict &tax) {
  if (tid == ChimeraClassify::kSpoolUnclassifiedTid ||
      tid >= tax.id2str.size()) {
    return 0;
  }
  uint32_t taxid = 0;
  if (!chimera::utils::try_parse_u32(tax.id2str[tid], taxid)) {
    return 0;
  }
  return taxid;
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

static const std::vector<ChimeraClassify::SpoolCandidate> &
select_spool_sequence_profile_candidates(
    const ChimeraClassify::SpoolReadRecord &record) {
  return select_spool_abundance_candidates(record);
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

static uint32_t taxid_to_species(const ChimeraClassify::TaxDict &tax,
                                 const ChimeraClassify::NcbiTaxdump *ncbiTaxdump,
                                 uint32_t tid_id) {
  if (tid_id == ChimeraClassify::kSpoolUnclassifiedTid ||
      tid_id >= tax.id2str.size()) {
    return 0;
  }
  uint32_t taxid = 0;
  if (!chimera::utils::try_parse_u32(tax.id2str[tid_id], taxid) ||
      taxid == 0) {
    return 0;
  }
  if (ncbiTaxdump && ncbiTaxdump->enabled()) {
    taxid = ncbiTaxdump->to_species(taxid);
  }
  return taxid;
}

static bool is_strict_species_taxid(
    uint32_t taxid, const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  if (taxid == 0) {
    return false;
  }
  if (!(ncbiTaxdump && ncbiTaxdump->enabled())) {
    return true;
  }
  return taxid < ncbiTaxdump->is_species.size() &&
         ncbiTaxdump->is_species[taxid] != 0;
}

static uint32_t taxid_to_strict_species(
    const ChimeraClassify::TaxDict &tax,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump, uint32_t tid_id) {
  const uint32_t species = taxid_to_species(tax, ncbiTaxdump, tid_id);
  return is_strict_species_taxid(species, ncbiTaxdump) ? species : 0;
}

static std::unordered_map<uint32_t, std::vector<double>>
build_species_length_samples(
    const chimera::presence::CoverageMeta &coverageMeta,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  std::unordered_map<uint32_t, std::vector<double>> lengths;
  lengths.reserve(coverageMeta.entries.size());
  for (const auto &entry : coverageMeta.entries) {
    uint32_t taxid = 0;
    if (!chimera::utils::try_parse_u32(entry.taxid, taxid) || taxid == 0 ||
        entry.genome_length == 0) {
      continue;
    }
    uint32_t species = taxid;
    if (ncbiTaxdump && ncbiTaxdump->enabled()) {
      species = ncbiTaxdump->to_species(taxid);
    }
    if (species == 0) {
      continue;
    }
    lengths[species].push_back(static_cast<double>(entry.genome_length));
  }
  return lengths;
}

struct SpeciesCallableScaleSamples {
  std::vector<double> genome_lengths;
  std::vector<double> callable_signatures;
  std::vector<double> unique_signatures;
};

static std::unordered_map<uint32_t, SpeciesCallableScaleSamples>
build_species_callable_scale_samples(
    const chimera::presence::CoverageMeta &coverageMeta,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  std::unordered_map<uint32_t, SpeciesCallableScaleSamples> samples;
  samples.reserve(coverageMeta.entries.size());
  for (const auto &entry : coverageMeta.entries) {
    uint32_t taxid = 0;
    if (!chimera::utils::try_parse_u32(entry.taxid, taxid) || taxid == 0) {
      continue;
    }
    uint32_t species = taxid;
    if (ncbiTaxdump && ncbiTaxdump->enabled()) {
      species = ncbiTaxdump->to_species(taxid);
    }
    if (species == 0) {
      continue;
    }
    auto &row = samples[species];
    if (entry.genome_length > 0) {
      row.genome_lengths.push_back(static_cast<double>(entry.genome_length));
    }
    if (entry.total_signatures > 0) {
      row.callable_signatures.push_back(
          static_cast<double>(entry.total_signatures));
    }
    if (entry.unique_signatures > 0) {
      row.unique_signatures.push_back(
          static_cast<double>(entry.unique_signatures));
    }
  }
  return samples;
}

static double median_length(std::vector<double> values) {
  if (values.empty()) {
    return 0.0;
  }
  const size_t mid = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid),
                   values.end());
  double median = values[mid];
  if ((values.size() % 2) == 0 && mid > 0) {
    const auto max_lower =
        std::max_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid));
    if (max_lower != values.begin() + static_cast<std::ptrdiff_t>(mid)) {
      median = 0.5 * (median + *max_lower);
    }
  }
  return median;
}

static SeqProfileFit fit_spool_sequence_profile(
    const std::vector<std::string> &spoolPaths,
    const ChimeraClassify::TaxDict &tax,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump,
    const chimera::presence::CoverageMeta &coverageMeta) {
  constexpr double kProfileScoreRatio = 0.95;
  SeqProfileFit fit;
  std::unordered_map<uint32_t, double> assigned;
  std::unordered_map<uint32_t, double> unique;
  std::unordered_map<uint32_t, double> feature_assigned;
  std::unordered_map<uint32_t, double> feature_unique;

  for (const auto &path : spoolPaths) {
    for_each_spool_record(
        path, [&](const ChimeraClassify::SpoolReadRecord &record) {
          ++fit.input_reads;
          const auto &candidates =
              select_spool_sequence_profile_candidates(record);
          if (candidates.empty()) {
            return;
          }
          std::unordered_map<uint32_t, double> by_species;
          by_species.reserve(candidates.size());
          for (const auto &cand : candidates) {
            if (!(cand.score > 0.0)) {
              continue;
            }
            const uint32_t species =
                taxid_to_species(tax, ncbiTaxdump, cand.tid);
            if (species == 0) {
              continue;
            }
            auto it = by_species.find(species);
            if (it == by_species.end() || cand.score > it->second) {
              by_species[species] = cand.score;
            }
          }
          if (by_species.empty()) {
            return;
          }
          double best = 0.0;
          for (const auto &[species, score] : by_species) {
            (void)species;
            if (score > best) {
              best = score;
            }
          }
          if (!(best > 0.0)) {
            return;
          }
          std::vector<std::pair<uint32_t, double>> kept;
          kept.reserve(by_species.size());
          const double cutoff = best * kProfileScoreRatio;
          double kept_score_sum = 0.0;
          for (const auto &[species, score] : by_species) {
            if (score >= cutoff) {
              kept.emplace_back(species, score);
              kept_score_sum += score;
            }
          }
          if (kept.empty() || !(kept_score_sum > 0.0)) {
            return;
          }
          ++fit.reads_with_candidates;
          fit.candidate_edges += kept.size();
          if (kept.size() > 1) {
            ++fit.multi_candidate_reads;
          }
          const double feature_weight =
              record.evaluated > 0.0 ? record.evaluated : 0.0;
          for (const auto &[species, score] : kept) {
            const double share = score / kept_score_sum;
            assigned[species] += share;
            if (feature_weight > 0.0) {
              const double feature_share = feature_weight * share;
              feature_assigned[species] += feature_share;
            }
          }
          if (kept.size() == 1) {
            unique[kept.front().first] += 1.0;
            if (feature_weight > 0.0) {
              feature_unique[kept.front().first] += feature_weight;
            }
          }
        });
  }

  fit.assigned_reads = 0.0;
  for (const auto &[species, value] : assigned) {
    (void)species;
    fit.assigned_reads += value;
  }
  fit.assigned_features = 0.0;
  for (const auto &[species, value] : feature_assigned) {
    (void)species;
    fit.assigned_features += value;
  }
  const auto scale_samples =
      build_species_callable_scale_samples(coverageMeta, ncbiTaxdump);
  double length_scaled_sum = 0.0;
  double sqrt_length_scaled_sum = 0.0;
  double callable_scaled_sum = 0.0;
  double unique_callable_scaled_sum = 0.0;
  double feature_length_scaled_sum = 0.0;
  double feature_sqrt_length_scaled_sum = 0.0;
  fit.rows.reserve(assigned.size());
  for (const auto &[species, value] : assigned) {
    if (!(value > 0.0)) {
      continue;
    }
    SeqProfileRow row;
    row.species_taxid = species;
    row.assigned_reads = value;
    auto unique_it = unique.find(species);
    if (unique_it != unique.end()) {
      row.unique_reads = unique_it->second;
    }
    auto feature_it = feature_assigned.find(species);
    if (feature_it != feature_assigned.end()) {
      row.assigned_features = feature_it->second;
    }
    auto feature_unique_it = feature_unique.find(species);
    if (feature_unique_it != feature_unique.end()) {
      row.unique_features = feature_unique_it->second;
    }
    if (fit.assigned_reads > 0.0) {
      row.sequence_abundance = value / fit.assigned_reads;
    }
    if (fit.assigned_features > 0.0) {
      row.feature_abundance = row.assigned_features / fit.assigned_features;
    }
    auto scale_it = scale_samples.find(species);
    if (scale_it != scale_samples.end()) {
      row.source_count =
          static_cast<uint32_t>(scale_it->second.genome_lengths.size());
      row.effective_length = median_length(scale_it->second.genome_lengths);
      if (row.effective_length > 0.0) {
        length_scaled_sum += value / row.effective_length;
        sqrt_length_scaled_sum += value / std::sqrt(row.effective_length);
        if (row.assigned_features > 0.0) {
          feature_length_scaled_sum +=
              row.assigned_features / row.effective_length;
          feature_sqrt_length_scaled_sum +=
              row.assigned_features / std::sqrt(row.effective_length);
        }
      }
      row.effective_callable_signatures =
          median_length(scale_it->second.callable_signatures);
      if (row.effective_callable_signatures > 0.0) {
        callable_scaled_sum += value / row.effective_callable_signatures;
      }
      row.effective_unique_signatures =
          median_length(scale_it->second.unique_signatures);
      if (row.effective_unique_signatures > 0.0) {
        unique_callable_scaled_sum += value / row.effective_unique_signatures;
      }
    }
    if (ncbiTaxdump && ncbiTaxdump->enabled()) {
      row.genus_taxid = ncbiTaxdump->to_genus(species);
    }
    fit.rows.push_back(row);
  }
  if (length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0) {
        row.length_normalized_abundance =
            (row.assigned_reads / row.effective_length) / length_scaled_sum;
      }
    }
  }
  if (sqrt_length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0) {
        row.sqrt_length_normalized_abundance =
            (row.assigned_reads / std::sqrt(row.effective_length)) /
            sqrt_length_scaled_sum;
      }
    }
  }
  if (callable_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_callable_signatures > 0.0) {
        row.callable_normalized_abundance =
            (row.assigned_reads / row.effective_callable_signatures) /
            callable_scaled_sum;
      }
    }
  }
  if (unique_callable_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_unique_signatures > 0.0) {
        row.unique_callable_normalized_abundance =
            (row.assigned_reads / row.effective_unique_signatures) /
            unique_callable_scaled_sum;
      }
    }
  }
  if (feature_length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0 && row.assigned_features > 0.0) {
        row.feature_length_normalized_abundance =
            (row.assigned_features / row.effective_length) /
            feature_length_scaled_sum;
      }
    }
  }
  if (feature_sqrt_length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0 && row.assigned_features > 0.0) {
        row.feature_sqrt_length_normalized_abundance =
            (row.assigned_features / std::sqrt(row.effective_length)) /
            feature_sqrt_length_scaled_sum;
      }
    }
  }
  std::sort(fit.rows.begin(), fit.rows.end(),
            [](const SeqProfileRow &lhs, const SeqProfileRow &rhs) {
              const double la = lhs.length_normalized_abundance > 0.0
                                    ? lhs.length_normalized_abundance
                                    : lhs.sequence_abundance;
              const double ra = rhs.length_normalized_abundance > 0.0
                                    ? rhs.length_normalized_abundance
                                    : rhs.sequence_abundance;
              if (la != ra) {
                return la > ra;
              }
              return lhs.species_taxid < rhs.species_taxid;
            });
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
  result.query_length = record.query_length;
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

static std::string resolve_profile_output_path(const std::string &outputFile) {
  std::filesystem::path path(outputFile);
  if (path.filename() == "ChimeraClassify.tsv") {
    return (path.parent_path() / "ChimeraProfile.tsv").string();
  }
  path.replace_extension(".profile.tsv");
  return path.string();
}

static std::string
resolve_profile_cami_output_path(const std::string &outputFile) {
  std::filesystem::path path(outputFile);
  if (path.filename() == "ChimeraClassify.tsv") {
    return (path.parent_path() / "ChimeraProfile.cami.tsv").string();
  }
  path.replace_extension(".profile.cami.tsv");
  return path.string();
}

static std::string resolve_profile_variant_output_path(
    const std::string &outputFile, const std::string &standardFilename,
    const std::string &replacementExtension) {
  std::filesystem::path path(outputFile);
  if (path.filename() == "ChimeraClassify.tsv") {
    return (path.parent_path() / standardFilename).string();
  }
  path.replace_extension(replacementExtension);
  return path.string();
}

static std::string
resolve_native_profile_trace_output_path(const std::string &outputFile) {
  std::filesystem::path path(outputFile);
  if (path.filename() == "ChimeraClassify.tsv") {
    return (path.parent_path() / "ChimeraProfile.trace.tsv").string();
  }
  path.replace_extension(".profile.trace.tsv");
  return path.string();
}

static std::string
resolve_native_profile_debug_output_path(const std::string &outputFile) {
  std::filesystem::path path(outputFile);
  if (path.filename() == "ChimeraClassify.tsv") {
    return (path.parent_path() / "ChimeraProfile.debug.tsv").string();
  }
  path.replace_extension(".profile.debug.tsv");
  return path.string();
}

static void copy_profile_file_or_throw(const std::filesystem::path &src,
                                       const std::filesystem::path &dst) {
  if (!std::filesystem::exists(src)) {
    throw std::runtime_error("Profile response-calibrated source is missing: " +
                             src.string());
  }
  if (!dst.parent_path().empty()) {
    std::filesystem::create_directories(dst.parent_path());
  }
  std::error_code ec;
  std::filesystem::copy_file(src, dst,
                             std::filesystem::copy_options::overwrite_existing,
                             ec);
  if (ec) {
    throw std::runtime_error("Failed to write profile response-calibrated "
                             "output from " +
                             src.string() + " to " + dst.string() + ": " +
                             ec.message());
  }
}

struct SpoolOutputPartStats {
  uint64_t classified{0};
  uint64_t unclassified{0};
  std::unordered_map<uint32_t, double> decision_taxid_counts;
  std::unordered_map<uint32_t, double> profile_response_taxid_counts;
  std::unordered_map<uint32_t, double> decision_species_counts;
  SpeciesProfileMasses localmix_masses;
};

static double decision_support_fraction(
    const ChimeraClassify::classifyResult &result) {
  if (!(result.evaluated > 0.0)) {
    return 1.0;
  }
  if (result.query_length == 0) {
    return 1.0;
  }
  const double density =
      result.evaluated / static_cast<double>(result.query_length);
  return std::clamp(density, 0.0, 1.0);
}

static uint32_t taxid_text_to_u32_or_zero(const std::string &taxidText) {
  uint32_t taxid = 0;
  if (!chimera::utils::try_parse_u32(taxidText, taxid)) {
    return 0;
  }
  return taxid;
}

using EvidenceAggregateMap = std::unordered_map<std::string, double>;
using LocalResolutionCallMap =
    std::unordered_map<std::string, ChimeraClassify::LocalResolutionReadCall>;

static uint32_t taxid_text_to_species(
    const std::string &taxidText,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump);

static const ChimeraClassify::LocalResolutionReadCall *
find_local_resolution_call(const LocalResolutionCallMap *localCalls,
                           const std::string &readId) {
  if (localCalls == nullptr) {
    return nullptr;
  }
  const auto found = localCalls->find(readId);
  if (found == localCalls->end()) {
    size_t firstTokenEnd = 0;
    while (firstTokenEnd < readId.size() &&
           std::isspace(static_cast<unsigned char>(readId[firstTokenEnd])) ==
               0) {
      ++firstTokenEnd;
    }
    if (firstTokenEnd == readId.size()) {
      return nullptr;
    }
    const auto tokenFound = localCalls->find(readId.substr(0, firstTokenEnd));
    if (tokenFound == localCalls->end()) {
      return nullptr;
    }
    return &tokenFound->second;
  }
  return &found->second;
}

static void accumulate_localmix_profile_candidates(
    const ChimeraClassify::classifyResult &result,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump,
    SpeciesProfileMasses &masses) {
  ++masses.input_reads;
  if (result.taxidCount.empty() ||
      result.taxidCount.front().first == "unclassified") {
    return;
  }

  const uint32_t decisionSpecies =
      taxid_text_to_species(result.taxidCount.front().first, ncbiTaxdump);
  if (decisionSpecies == 0) {
    return;
  }

  uint32_t decisionGenus = 0;
  if (ncbiTaxdump && ncbiTaxdump->enabled()) {
    decisionGenus = ncbiTaxdump->to_genus(decisionSpecies);
  }

  std::unordered_map<uint32_t, double> speciesScore;
  speciesScore.reserve(result.posteriors.size() + 1);
  for (const auto &[taxidText, weight] : result.posteriors) {
    if (!(weight > 0.0) || taxidText.empty() ||
        taxidText == "unclassified") {
      continue;
    }
    const uint32_t species = taxid_text_to_species(taxidText, ncbiTaxdump);
    if (species == 0) {
      continue;
    }
    if (decisionGenus != 0 && ncbiTaxdump && ncbiTaxdump->enabled() &&
        ncbiTaxdump->to_genus(species) != decisionGenus) {
      continue;
    }
    auto &slot = speciesScore[species];
    slot = std::max(slot, weight);
  }

  if (speciesScore.empty()) {
    speciesScore.emplace(decisionSpecies, 1.0);
  }

  double total = 0.0;
  for (auto &[species, score] : speciesScore) {
    (void)species;
    score = std::sqrt(std::max(0.0, score));
    total += score;
  }
  if (!(total > 0.0)) {
    speciesScore.clear();
    speciesScore.emplace(decisionSpecies, 1.0);
    total = 1.0;
  }

  ++masses.reads_with_candidates;
  masses.candidate_edges += speciesScore.size();
  if (speciesScore.size() > 1) {
    ++masses.multi_candidate_reads;
  }

  const double feature_weight =
      result.evaluated > 0.0 ? result.evaluated : 0.0;
  for (const auto &[species, score] : speciesScore) {
    const double share = score / total;
    masses.assigned_reads[species] += share;
    if (feature_weight > 0.0) {
      masses.assigned_features[species] += feature_weight * share;
    }
  }
  if (speciesScore.size() == 1) {
    const uint32_t species = speciesScore.begin()->first;
    masses.unique_reads[species] += 1.0;
    if (feature_weight > 0.0) {
      masses.unique_features[species] += feature_weight;
    }
  }
}

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
  std::filesystem::path shard_manifest_path;
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
      chimera::local_resolution::materialize_manifest_path(
          corePath, manifest->shard_manifest),
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
  if (localCalls == nullptr) {
    return false;
  }
  if (sampleDivergence < divergenceThreshold) {
    return false;
  }
  const auto *callPtr = find_local_resolution_call(localCalls, result.id);
  if (callPtr == nullptr) {
    return false;
  }
  const auto &call = *callPtr;
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
      presenceDecision, true,
      [&](const std::vector<ChimeraClassify::SpoolReadRecord> *records,
          std::vector<ChimeraClassify::classifyResult> &chunk) {
        for (size_t resultIndex = 0; resultIndex < chunk.size();
             ++resultIndex) {
          auto &result = chunk[resultIndex];
          if (records != nullptr && resultIndex < records->size()) {
            const uint32_t observed =
                (*records)[resultIndex].profile_response_taxid;
            if (observed != 0) {
              partStats.profile_response_taxid_counts[observed] += 1.0;
            }
          }
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
          accumulate_localmix_profile_candidates(
              result, ncbiTaxdump, partStats.localmix_masses);
          if (!result.taxidCount.empty() &&
              result.taxidCount.front().first == "unclassified") {
            ++partStats.unclassified;
          } else {
            ++partStats.classified;
            if (!result.taxidCount.empty()) {
              const uint32_t decisionTaxid =
                  taxid_text_to_u32_or_zero(result.taxidCount.front().first);
              if (decisionTaxid != 0) {
                partStats.decision_taxid_counts[decisionTaxid] += 1.0;
              }
              const uint32_t decisionSpecies =
                  taxid_text_to_species(result.taxidCount.front().first,
                                        ncbiTaxdump);
              if (decisionSpecies != 0) {
                partStats.decision_species_counts[decisionSpecies] += 1.0;
              }
            }
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

static std::string profile_taxon_name(
    uint32_t taxid, const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  if (ncbiTaxdump && ncbiTaxdump->enabled()) {
    const std::string &name = ncbiTaxdump->name(taxid);
    if (!name.empty()) {
      return name;
    }
  }
  return std::to_string(taxid);
}

static uint32_t taxid_text_to_species(
    const std::string &taxidText,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  uint32_t taxid = 0;
  if (!chimera::utils::try_parse_u32(taxidText, taxid) || taxid == 0) {
    return 0;
  }
  if (ncbiTaxdump && ncbiTaxdump->enabled()) {
    taxid = ncbiTaxdump->to_species(taxid);
  }
  return taxid;
}

static SeqProfileFit make_profile_from_species_counts(
    const std::unordered_map<uint32_t, double> &speciesCounts,
    uint64_t inputReads, uint64_t classifiedReads,
    const chimera::presence::CoverageMeta &coverageMeta,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump);

static uint32_t strict_species_for_taxid(
    uint32_t taxid, const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  if (taxid == 0) {
    return 0;
  }
  if (!(ncbiTaxdump && ncbiTaxdump->enabled())) {
    return taxid;
  }
  const uint32_t species = ncbiTaxdump->to_species(taxid);
  return is_strict_species_taxid(species, ncbiTaxdump) ? species : 0;
}

static bool taxid_is_descendant_of(
    uint32_t taxid, uint32_t ancestor,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  if (taxid == 0 || ancestor == 0) {
    return false;
  }
  if (taxid == ancestor) {
    return true;
  }
  if (!(ncbiTaxdump && ncbiTaxdump->enabled()) ||
      taxid >= ncbiTaxdump->parent.size()) {
    return false;
  }
  uint32_t cur = taxid;
  for (int steps = 0; steps < 128; ++steps) {
    if (cur == 0 || cur >= ncbiTaxdump->parent.size()) {
      break;
    }
    const uint32_t parent = ncbiTaxdump->parent[cur];
    if (parent == 0 || parent == cur) {
      break;
    }
    if (parent == ancestor) {
      return true;
    }
    cur = parent;
  }
  return false;
}

static SeqProfileFit make_profile_from_species_masses(
    const SpeciesProfileMasses &masses,
    const chimera::presence::CoverageMeta &coverageMeta,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  SeqProfileFit fit;
  fit.input_reads = masses.input_reads;
  fit.reads_with_candidates = masses.reads_with_candidates;
  fit.multi_candidate_reads = masses.multi_candidate_reads;
  fit.candidate_edges = masses.candidate_edges;

  for (const auto &[species, value] : masses.assigned_reads) {
    if (species != 0 && value > 0.0) {
      fit.assigned_reads += value;
    }
  }
  for (const auto &[species, value] : masses.assigned_features) {
    if (species != 0 && value > 0.0) {
      fit.assigned_features += value;
    }
  }

  const auto scale_samples =
      build_species_callable_scale_samples(coverageMeta, ncbiTaxdump);
  double length_scaled_sum = 0.0;
  double sqrt_length_scaled_sum = 0.0;
  double callable_scaled_sum = 0.0;
  double unique_callable_scaled_sum = 0.0;
  double feature_length_scaled_sum = 0.0;
  double feature_sqrt_length_scaled_sum = 0.0;

  fit.rows.reserve(masses.assigned_reads.size());
  for (const auto &[species, value] : masses.assigned_reads) {
    if (species == 0 || !(value > 0.0)) {
      continue;
    }
    SeqProfileRow row;
    row.species_taxid = species;
    row.assigned_reads = value;
    if (fit.assigned_reads > 0.0) {
      row.sequence_abundance = value / fit.assigned_reads;
    }
    auto unique_it = masses.unique_reads.find(species);
    if (unique_it != masses.unique_reads.end()) {
      row.unique_reads = unique_it->second;
    }
    auto feature_it = masses.assigned_features.find(species);
    if (feature_it != masses.assigned_features.end()) {
      row.assigned_features = feature_it->second;
    }
    auto feature_unique_it = masses.unique_features.find(species);
    if (feature_unique_it != masses.unique_features.end()) {
      row.unique_features = feature_unique_it->second;
    }
    if (fit.assigned_features > 0.0) {
      row.feature_abundance = row.assigned_features / fit.assigned_features;
    }
    auto scale_it = scale_samples.find(species);
    if (scale_it != scale_samples.end()) {
      row.source_count =
          static_cast<uint32_t>(scale_it->second.genome_lengths.size());
      row.effective_length = median_length(scale_it->second.genome_lengths);
      if (row.effective_length > 0.0) {
        length_scaled_sum += value / row.effective_length;
        sqrt_length_scaled_sum += value / std::sqrt(row.effective_length);
        if (row.assigned_features > 0.0) {
          feature_length_scaled_sum +=
              row.assigned_features / row.effective_length;
          feature_sqrt_length_scaled_sum +=
              row.assigned_features / std::sqrt(row.effective_length);
        }
      }
      row.effective_callable_signatures =
          median_length(scale_it->second.callable_signatures);
      if (row.effective_callable_signatures > 0.0) {
        callable_scaled_sum += value / row.effective_callable_signatures;
      }
      row.effective_unique_signatures =
          median_length(scale_it->second.unique_signatures);
      if (row.effective_unique_signatures > 0.0) {
        unique_callable_scaled_sum += value / row.effective_unique_signatures;
      }
    }
    if (ncbiTaxdump && ncbiTaxdump->enabled()) {
      row.genus_taxid = ncbiTaxdump->to_genus(species);
    }
    fit.rows.push_back(row);
  }

  if (length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0) {
        row.length_normalized_abundance =
            (row.assigned_reads / row.effective_length) / length_scaled_sum;
      }
    }
  }
  if (sqrt_length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0) {
        row.sqrt_length_normalized_abundance =
            (row.assigned_reads / std::sqrt(row.effective_length)) /
            sqrt_length_scaled_sum;
      }
    }
  }
  if (callable_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_callable_signatures > 0.0) {
        row.callable_normalized_abundance =
            (row.assigned_reads / row.effective_callable_signatures) /
            callable_scaled_sum;
      }
    }
  }
  if (unique_callable_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_unique_signatures > 0.0) {
        row.unique_callable_normalized_abundance =
            (row.assigned_reads / row.effective_unique_signatures) /
            unique_callable_scaled_sum;
      }
    }
  }
  if (feature_length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0 && row.assigned_features > 0.0) {
        row.feature_length_normalized_abundance =
            (row.assigned_features / row.effective_length) /
            feature_length_scaled_sum;
      }
    }
  }
  if (feature_sqrt_length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0 && row.assigned_features > 0.0) {
        row.feature_sqrt_length_normalized_abundance =
            (row.assigned_features / std::sqrt(row.effective_length)) /
            feature_sqrt_length_scaled_sum;
      }
    }
  }

  std::sort(fit.rows.begin(), fit.rows.end(),
            [](const SeqProfileRow &lhs, const SeqProfileRow &rhs) {
              const double la = lhs.length_normalized_abundance > 0.0
                                    ? lhs.length_normalized_abundance
                                    : lhs.sequence_abundance;
              const double ra = rhs.length_normalized_abundance > 0.0
                                    ? rhs.length_normalized_abundance
                                    : rhs.sequence_abundance;
              if (la != ra) {
                return la > ra;
              }
              return lhs.species_taxid < rhs.species_taxid;
            });
  return fit;
}

static SeqProfileFit make_profile_from_species_counts(
    const std::unordered_map<uint32_t, double> &speciesCounts,
    uint64_t inputReads, uint64_t classifiedReads,
    const chimera::presence::CoverageMeta &coverageMeta,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  SeqProfileFit fit;
  fit.input_reads = inputReads;
  fit.reads_with_candidates = classifiedReads;
  fit.candidate_edges = classifiedReads;
  fit.rows.reserve(speciesCounts.size());
  double assigned_sum = 0.0;
  for (const auto &[species, value] : speciesCounts) {
    if (species != 0 && value > 0.0) {
      assigned_sum += value;
    }
  }
  fit.assigned_reads = assigned_sum;

  const auto scale_samples =
      build_species_callable_scale_samples(coverageMeta, ncbiTaxdump);
  double length_scaled_sum = 0.0;
  double sqrt_length_scaled_sum = 0.0;
  double callable_scaled_sum = 0.0;
  double unique_callable_scaled_sum = 0.0;
  for (const auto &[species, value] : speciesCounts) {
    if (species == 0 || !(value > 0.0)) {
      continue;
    }
    SeqProfileRow row;
    row.species_taxid = species;
    row.assigned_reads = value;
    row.unique_reads = value;
    if (assigned_sum > 0.0) {
      row.sequence_abundance = value / assigned_sum;
    }
    auto scale_it = scale_samples.find(species);
    if (scale_it != scale_samples.end()) {
      row.source_count =
          static_cast<uint32_t>(scale_it->second.genome_lengths.size());
      row.effective_length = median_length(scale_it->second.genome_lengths);
      if (row.effective_length > 0.0) {
        length_scaled_sum += value / row.effective_length;
        sqrt_length_scaled_sum += value / std::sqrt(row.effective_length);
      }
      row.effective_callable_signatures =
          median_length(scale_it->second.callable_signatures);
      if (row.effective_callable_signatures > 0.0) {
        callable_scaled_sum += value / row.effective_callable_signatures;
      }
      row.effective_unique_signatures =
          median_length(scale_it->second.unique_signatures);
      if (row.effective_unique_signatures > 0.0) {
        unique_callable_scaled_sum += value / row.effective_unique_signatures;
      }
    }
    if (ncbiTaxdump && ncbiTaxdump->enabled()) {
      row.genus_taxid = ncbiTaxdump->to_genus(species);
    }
    fit.rows.push_back(row);
  }
  if (length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0) {
        row.length_normalized_abundance =
            (row.assigned_reads / row.effective_length) / length_scaled_sum;
      }
    }
  }
  if (sqrt_length_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_length > 0.0) {
        row.sqrt_length_normalized_abundance =
            (row.assigned_reads / std::sqrt(row.effective_length)) /
            sqrt_length_scaled_sum;
      }
    }
  }
  if (callable_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_callable_signatures > 0.0) {
        row.callable_normalized_abundance =
            (row.assigned_reads / row.effective_callable_signatures) /
            callable_scaled_sum;
      }
    }
  }
  if (unique_callable_scaled_sum > 0.0) {
    for (auto &row : fit.rows) {
      if (row.effective_unique_signatures > 0.0) {
        row.unique_callable_normalized_abundance =
            (row.assigned_reads / row.effective_unique_signatures) /
            unique_callable_scaled_sum;
      }
    }
  }
  return fit;
}

static double primary_profile_abundance(const SeqProfileRow &row,
                                        PrimaryProfileScale primaryScale) {
  if (primaryScale == PrimaryProfileScale::LengthNormalized) {
    return row.length_normalized_abundance > 0.0
               ? row.length_normalized_abundance
               : row.sequence_abundance;
  }
  if (primaryScale == PrimaryProfileScale::SqrtLengthNormalized) {
    return row.sqrt_length_normalized_abundance > 0.0
               ? row.sqrt_length_normalized_abundance
               : row.sequence_abundance;
  }
  if (primaryScale == PrimaryProfileScale::CallableNormalized) {
    return row.callable_normalized_abundance > 0.0
               ? row.callable_normalized_abundance
               : row.sequence_abundance;
  }
  if (primaryScale == PrimaryProfileScale::UniqueCallableNormalized) {
    return row.unique_callable_normalized_abundance > 0.0
               ? row.unique_callable_normalized_abundance
               : row.sequence_abundance;
  }
  return row.sequence_abundance;
}

static const char *primary_profile_scale_label(
    PrimaryProfileScale primaryScale) {
  switch (primaryScale) {
  case PrimaryProfileScale::LengthNormalized:
    return "length_normalized_abundance";
  case PrimaryProfileScale::SqrtLengthNormalized:
    return "sqrt_length_normalized_abundance";
  case PrimaryProfileScale::CallableNormalized:
    return "callable_signature_normalized_abundance";
  case PrimaryProfileScale::UniqueCallableNormalized:
    return "unique_signature_normalized_abundance";
  case PrimaryProfileScale::SequenceAbundance:
    return "sequence_abundance";
  }
  return "sequence_abundance";
}

static std::string clean_profile_tsv_field(std::string value) {
  for (char &ch : value) {
    if (ch == '\t' || ch == '\n' || ch == '\r') {
      ch = ' ';
    }
  }
  return value;
}

static uint32_t
profile_parent_taxid(uint32_t taxid,
                     const ChimeraClassify::NcbiTaxdump *taxdump) {
  if (taxdump == nullptr || !taxdump->enabled() || taxid == 0) {
    return 0;
  }
  const uint32_t genus = taxdump->to_genus(taxid);
  if (genus != 0) {
    return genus;
  }
  return taxid < taxdump->parent.size() ? taxdump->parent[taxid] : 0;
}

struct ClassifierResponseReportabilityMeta {
  std::unordered_map<uint32_t, uint32_t> genus_species_count;
};

static ClassifierResponseReportabilityMeta build_classifier_response_meta(
    const chimera::presence::CoverageMeta &coverageMeta,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump) {
  ClassifierResponseReportabilityMeta meta;
  if (ncbiTaxdump == nullptr || !ncbiTaxdump->enabled()) {
    return meta;
  }
  std::unordered_map<uint32_t, std::unordered_set<uint32_t>> by_genus;
  by_genus.reserve(coverageMeta.entries.size());
  for (const auto &entry : coverageMeta.entries) {
    uint32_t taxid = 0;
    if (!chimera::utils::try_parse_u32(entry.taxid, taxid) || taxid == 0) {
      continue;
    }
    const uint32_t species = ncbiTaxdump->to_species(taxid);
    if (species == 0) {
      continue;
    }
    const uint32_t genus = ncbiTaxdump->to_genus(species);
    if (genus == 0) {
      continue;
    }
    by_genus[genus].insert(species);
  }
  for (const auto &[genus, species] : by_genus) {
    const uint32_t count = static_cast<uint32_t>(species.size());
    for (uint32_t taxid : species) {
      meta.genus_species_count[taxid] = count;
    }
  }
  return meta;
}

static std::string classifier_response_column_type(const SeqProfileRow &row) {
  if (!(row.effective_callable_signatures > 0.0)) {
    return "unusable_no_callable_exposure";
  }
  if (row.source_count <= 1) {
    return "single_source_limited";
  }
  return "species_resolvable";
}

static double classifier_response_effective_taxa(
    const std::vector<double> &normalizedAbundances) {
  double sum_sq = 0.0;
  for (double value : normalizedAbundances) {
    if (value > 0.0 && std::isfinite(value)) {
      sum_sq += value * value;
    }
  }
  if (!(sum_sq > 0.0)) {
    return 0.0;
  }
  return 1.0 / sum_sq;
}

static double classifier_response_dynamic_report_floor_percent(
    const std::vector<double> &normalizedAbundances) {
  const double effective_taxa =
      classifier_response_effective_taxa(normalizedAbundances);
  if (!(effective_taxa > 0.0) || !std::isfinite(effective_taxa)) {
    return 0.05;
  }
  return std::clamp(0.5 / effective_taxa, 0.01, 0.05);
}

static uint64_t classifier_response_min_read_support(const SeqProfileFit &fit) {
  const uint64_t reads = fit.reads_with_candidates > 0 ? fit.reads_with_candidates
                                                       : fit.input_reads;
  const double support = 0.01 * std::sqrt(static_cast<double>(reads));
  return static_cast<uint64_t>(std::max<double>(3.0, std::ceil(support)));
}

static std::unordered_map<uint32_t, double>
profile_abundance_by_species(const SeqProfileFit &fit,
                             PrimaryProfileScale scale) {
  std::unordered_map<uint32_t, double> mass;
  mass.reserve(fit.rows.size());
  double total = 0.0;
  for (const SeqProfileRow &row : fit.rows) {
    if (row.species_taxid == 0) {
      continue;
    }
    const double abundance = primary_profile_abundance(row, scale);
    if (!(abundance > 0.0) || !std::isfinite(abundance)) {
      continue;
    }
    mass[row.species_taxid] += abundance;
    total += abundance;
  }
  if (total > 0.0) {
    const double inv_total = 1.0 / total;
    for (auto &[taxid, value] : mass) {
      (void)taxid;
      value *= inv_total;
    }
  }
  return mass;
}

static void write_classifier_response_reportable_profile_outputs(
    const std::string &nativePath, const std::string &camiPath,
    const std::string &debugPath, const std::string &tracePath,
    const SeqProfileFit &fit, const ChimeraClassify::NcbiTaxdump *ncbiTaxdump,
    const chimera::presence::CoverageMeta &coverageMeta,
    const SeqProfileFit *localmixFit = nullptr,
    const char *responseSourceOverride = nullptr,
    PrimaryProfileScale primaryScale = PrimaryProfileScale::CallableNormalized) {
  struct CandidateRow {
    const SeqProfileRow *row{nullptr};
    double raw_abundance{0.0};
    double raw_percent{0.0};
    double response_abundance{0.0};
    double response_percent{0.0};
    bool kept{false};
    std::string column_type;
    std::string reportability;
    uint32_t genus_species_count{0};
  };

  constexpr const char *kProfileInput =
      "classifier_response_sparse_reportability_v1";
  const PrimaryProfileScale kScale = primaryScale;

  std::vector<CandidateRow> rows;
  rows.reserve(fit.rows.size());
  double raw_total = 0.0;
  for (const SeqProfileRow &row : fit.rows) {
    if (row.species_taxid == 0) {
      continue;
    }
    const double abundance = primary_profile_abundance(row, kScale);
    if (!(abundance > 0.0) || !std::isfinite(abundance)) {
      continue;
    }
    raw_total += abundance;
    rows.push_back(
        CandidateRow{&row, abundance, 0.0, 0.0, 0.0, false, {}, {}, 0});
  }

  std::vector<double> normalized;
  normalized.reserve(rows.size());
  if (raw_total > 0.0) {
    for (CandidateRow &entry : rows) {
      entry.raw_abundance /= raw_total;
      entry.raw_percent = 100.0 * entry.raw_abundance;
      normalized.push_back(entry.raw_abundance);
    }
  }
  const double effective_taxa =
      classifier_response_effective_taxa(normalized);
  const double report_floor_percent =
      classifier_response_dynamic_report_floor_percent(normalized);
  const uint64_t min_read_support =
      classifier_response_min_read_support(fit);
  const auto meta = build_classifier_response_meta(coverageMeta, ncbiTaxdump);
  const bool use_localmix_response = responseSourceOverride == nullptr &&
                                     localmixFit != nullptr &&
                                     effective_taxa >= 32.0;
  const auto localmix_mass =
      use_localmix_response
          ? profile_abundance_by_species(*localmixFit,
                                         PrimaryProfileScale::SequenceAbundance)
          : std::unordered_map<uint32_t, double>{};
  const char *abundance_response_source =
      responseSourceOverride != nullptr
          ? responseSourceOverride
          : (use_localmix_response ? "localmix_restricted_response"
                                   : "callable_diagonal_response");

  uint64_t no_exposure_rows = 0;
  uint64_t low_support_rows = 0;
  uint64_t below_floor_rows = 0;
  uint64_t single_source_rows = 0;
  uint64_t reportable_without_response_rows = 0;
  double kept_mass = 0.0;
  double kept_observation_mass = 0.0;
  double nuisance_mass = 0.0;
  double single_source_kept_mass = 0.0;

  for (CandidateRow &entry : rows) {
    const SeqProfileRow &row = *entry.row;
    entry.column_type = classifier_response_column_type(row);
    auto genus_it = meta.genus_species_count.find(row.species_taxid);
    if (genus_it != meta.genus_species_count.end()) {
      entry.genus_species_count = genus_it->second;
    }

    if (!(row.effective_callable_signatures > 0.0)) {
      entry.reportability = "nuisance_no_callable_exposure";
      ++no_exposure_rows;
    } else if (row.assigned_reads <
               static_cast<double>(min_read_support)) {
      entry.reportability = "nuisance_low_runtime_support";
      ++low_support_rows;
    } else if (entry.raw_percent < report_floor_percent) {
      entry.reportability = "nuisance_below_dynamic_detection_floor";
      ++below_floor_rows;
    } else if (row.source_count <= 1) {
      entry.kept = true;
      entry.reportability = "single_source_limited_reportable";
      ++single_source_rows;
    } else {
      entry.kept = true;
      entry.reportability = "classifier_response_reportable";
    }

    if (entry.kept) {
      entry.response_abundance = entry.raw_abundance;
      if (use_localmix_response) {
        auto response_it = localmix_mass.find(row.species_taxid);
        entry.response_abundance =
            response_it == localmix_mass.end() ? 0.0 : response_it->second;
      }
      if (!(entry.response_abundance > 0.0) ||
          !std::isfinite(entry.response_abundance)) {
        entry.kept = false;
        entry.reportability = "nuisance_no_abundance_response";
        ++reportable_without_response_rows;
        nuisance_mass += entry.raw_abundance;
        continue;
      }
      kept_mass += entry.response_abundance;
      kept_observation_mass += entry.raw_abundance;
      if (row.source_count <= 1) {
        single_source_kept_mass += entry.response_abundance;
      }
    } else {
      nuisance_mass += entry.raw_abundance;
    }
    entry.response_percent = 100.0 * entry.response_abundance;
  }

  std::sort(rows.begin(), rows.end(), [](const CandidateRow &lhs,
                                         const CandidateRow &rhs) {
    const double lhs_response = lhs.kept ? lhs.response_abundance : 0.0;
    const double rhs_response = rhs.kept ? rhs.response_abundance : 0.0;
    if (lhs_response != rhs_response) {
      return lhs_response > rhs_response;
    }
    if (lhs.kept != rhs.kept) {
      return lhs.kept;
    }
    if (lhs.raw_abundance != rhs.raw_abundance) {
      return lhs.raw_abundance > rhs.raw_abundance;
    }
    const uint32_t lt = lhs.row == nullptr ? 0 : lhs.row->species_taxid;
    const uint32_t rt = rhs.row == nullptr ? 0 : rhs.row->species_taxid;
    return lt < rt;
  });

  auto ensure_parent = [](const std::string &path) {
    const auto parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) {
      std::filesystem::create_directories(parent);
    }
  };

  ensure_parent(nativePath);
  {
    std::ofstream out(nativePath, std::ios::out | std::ios::binary);
    if (!out.is_open()) {
      throw std::runtime_error(
          "Failed to open classifier-response profile: " + nativePath);
    }
    out << "# profile_input=" << kProfileInput << "\n";
    out << "# count_unit=" << primary_profile_scale_label(kScale) << "\n";
    out << "# abundance_basis=reportability_constrained_response_normalized\n";
    out << "# abundance_response_source=" << abundance_response_source << "\n";
    out << "# classifier_response_effective_taxa=" << effective_taxa << "\n";
    out << "# report_floor_policy=clamp(0.5/effective_taxa,0.01,0.05)\n";
    out << "# dynamic_report_min_percent=" << report_floor_percent << "\n";
    out << "# min_runtime_read_support=" << min_read_support << "\n";
    out << "# source_profile_rows=" << rows.size() << "\n";
    out << "# reportable_profile_rows="
        << std::count_if(rows.begin(), rows.end(),
                         [](const CandidateRow &row) { return row.kept; })
        << "\n";
    out << "# raw_biological_mass=" << raw_total << "\n";
    out << "# kept_response_mass=" << kept_mass << "\n";
    out << "# kept_observation_mass=" << kept_observation_mass << "\n";
    out << "# nuisance_raw_mass=" << nuisance_mass << "\n";
    out << "# single_source_limited_kept_mass=" << single_source_kept_mass
        << "\n";
    out << "# no_callable_exposure_rows=" << no_exposure_rows << "\n";
    out << "# low_runtime_support_rows=" << low_support_rows << "\n";
    out << "# below_dynamic_floor_rows=" << below_floor_rows << "\n";
    out << "# reportable_without_response_rows="
        << reportable_without_response_rows << "\n";
    out << "rank\ttaxid\tname\tparent_taxid\trelative_abundance"
        << "\tassigned_mass\teffective_yield\tread_support"
        << "\tmeasurement_unit_support\tsingle_source_limited\treportability\n";
    out << std::setprecision(17);
    if (kept_mass > 0.0) {
      for (const CandidateRow &entry : rows) {
        if (!entry.kept || entry.row == nullptr) {
          continue;
        }
        const SeqProfileRow &row = *entry.row;
        out << "species\t" << row.species_taxid << '\t'
            << clean_profile_tsv_field(
                   profile_taxon_name(row.species_taxid, ncbiTaxdump))
            << '\t' << profile_parent_taxid(row.species_taxid, ncbiTaxdump)
            << '\t' << (entry.response_abundance / kept_mass) << '\t'
            << entry.response_abundance << '\t'
            << row.effective_callable_signatures << '\t'
            << row.assigned_reads << '\t' << row.source_count << '\t'
            << (row.source_count <= 1 ? 1 : 0) << '\t'
            << entry.reportability << '\n';
      }
    }
    out.close();
    if (!out.good()) {
      throw std::runtime_error(
          "Failed to close classifier-response profile: " + nativePath);
    }
  }

  ensure_parent(camiPath);
  {
    std::ofstream out(camiPath, std::ios::out | std::ios::binary);
    if (!out.is_open()) {
      throw std::runtime_error(
          "Failed to open classifier-response CAMI profile: " + camiPath);
    }
    out << "# profile_input=" << kProfileInput << "\n";
    out << "# count_unit=" << primary_profile_scale_label(kScale) << "\n";
    out << "# abundance_basis=reportability_constrained_response_normalized\n";
    out << "# abundance_response_source=" << abundance_response_source << "\n";
    out << "# classifier_response_effective_taxa=" << effective_taxa << "\n";
    out << "# dynamic_report_min_percent=" << report_floor_percent << "\n";
    out << "# min_runtime_read_support=" << min_read_support << "\n";
    out << "# kept_response_mass=" << kept_mass << "\n";
    out << "# kept_observation_mass=" << kept_observation_mass << "\n";
    out << "# nuisance_raw_mass=" << nuisance_mass << "\n";
    out << "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n";
    out << std::setprecision(12);
    if (kept_mass > 0.0) {
      for (const CandidateRow &entry : rows) {
        if (!entry.kept || entry.row == nullptr) {
          continue;
        }
        const SeqProfileRow &row = *entry.row;
        const std::string taxid = std::to_string(row.species_taxid);
        const std::string name =
            profile_taxon_name(row.species_taxid, ncbiTaxdump);
        out << taxid << "\tspecies\t" << taxid << '\t' << name << '\t'
            << (100.0 * entry.response_abundance / kept_mass) << '\n';
      }
    }
    out.close();
    if (!out.good()) {
      throw std::runtime_error(
          "Failed to close classifier-response CAMI profile: " + camiPath);
    }
  }

  auto write_debug_like = [&](const std::string &path, const char *kind) {
    ensure_parent(path);
    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out.is_open()) {
      throw std::runtime_error(
          "Failed to open classifier-response debug profile: " + path);
    }
    out << "# profile_input=" << kProfileInput << "\n";
    out << "# debug_kind=" << kind << "\n";
    out << "key\tvalue\n";
    out << "source_profile_rows\t" << rows.size() << "\n";
    out << "reportable_profile_rows\t"
        << std::count_if(rows.begin(), rows.end(),
                         [](const CandidateRow &row) { return row.kept; })
        << "\n";
    out << "dynamic_report_min_percent\t" << report_floor_percent << "\n";
    out << "abundance_response_source\t" << abundance_response_source << "\n";
    out << "classifier_response_effective_taxa\t" << effective_taxa << "\n";
    out << "min_runtime_read_support\t" << min_read_support << "\n";
    out << "raw_biological_mass\t" << raw_total << "\n";
    out << "kept_response_mass\t" << kept_mass << "\n";
    out << "kept_observation_mass\t" << kept_observation_mass << "\n";
    out << "nuisance_raw_mass\t" << nuisance_mass << "\n";
    out << "single_source_limited_kept_mass\t" << single_source_kept_mass
        << "\n";
    out << "no_callable_exposure_rows\t" << no_exposure_rows << "\n";
    out << "low_runtime_support_rows\t" << low_support_rows << "\n";
    out << "below_dynamic_floor_rows\t" << below_floor_rows << "\n";
    out << "reportable_without_response_rows\t"
        << reportable_without_response_rows << "\n";
    out << "\n";
    out << "taxid\tname\traw_percent\tresponse_percent\trelative_if_reported"
        << "\tassigned_reads\teffective_callable_exposure\tsource_count"
        << "\tgenus_species_count\tcolumn_type\tkept\treportability\n";
    out << std::setprecision(12);
    for (const CandidateRow &entry : rows) {
      if (entry.row == nullptr) {
        continue;
      }
      const SeqProfileRow &row = *entry.row;
      out << row.species_taxid << '\t'
          << clean_profile_tsv_field(
                 profile_taxon_name(row.species_taxid, ncbiTaxdump))
          << '\t' << entry.raw_percent << '\t'
          << entry.response_percent << '\t'
          << (entry.kept && kept_mass > 0.0
                  ? 100.0 * entry.response_abundance / kept_mass
                  : 0.0)
          << '\t' << row.assigned_reads << '\t'
          << row.effective_callable_signatures << '\t' << row.source_count
          << '\t' << entry.genus_species_count << '\t'
          << entry.column_type << '\t' << (entry.kept ? 1 : 0) << '\t'
          << entry.reportability << '\n';
    }
    out.close();
    if (!out.good()) {
      throw std::runtime_error(
          "Failed to close classifier-response debug profile: " + path);
    }
  };
  write_debug_like(debugPath, "debug");
  write_debug_like(tracePath, "trace");
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

static SpeciesProfileMasses
merge_localmix_profile_masses(const std::vector<SpoolOutputPartStats> &parts) {
  SpeciesProfileMasses merged;
  auto merge_map = [](std::unordered_map<uint32_t, double> &dst,
                      const std::unordered_map<uint32_t, double> &src) {
    for (const auto &[species, value] : src) {
      dst[species] += value;
    }
  };

  for (const auto &part : parts) {
    const SpeciesProfileMasses &src = part.localmix_masses;
    merged.input_reads += src.input_reads;
    merged.reads_with_candidates += src.reads_with_candidates;
    merged.multi_candidate_reads += src.multi_candidate_reads;
    merged.candidate_edges += src.candidate_edges;
    merge_map(merged.assigned_reads, src.assigned_reads);
    merge_map(merged.unique_reads, src.unique_reads);
    merge_map(merged.assigned_features, src.assigned_features);
    merge_map(merged.unique_features, src.unique_features);
  }
  return merged;
}

static std::string lowercase_copy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  return value;
}

static bool ends_with_suffix(const std::string &value,
                             const std::string &suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) ==
             0;
}

static bool path_looks_like_fasta(const std::string &inputPath) {
  std::string name =
      lowercase_copy(std::filesystem::path(inputPath).filename().string());
  for (const std::string &compressed : {".gz", ".zst", ".bz2", ".xz"}) {
    if (ends_with_suffix(name, compressed)) {
      name.resize(name.size() - compressed.size());
      break;
    }
  }
  return ends_with_suffix(name, ".fa") || ends_with_suffix(name, ".fasta") ||
         ends_with_suffix(name, ".fna") || ends_with_suffix(name, ".ffn");
}

static bool classify_inputs_are_unweighted_fasta(
    const ChimeraClassify::ClassifyConfig &config) {
  if (!config.pairedFiles.empty() || config.singleFiles.empty()) {
    return false;
  }
  for (const auto &path : config.singleFiles) {
    if (!path_looks_like_fasta(path)) {
      return false;
    }
  }
  return true;
}

static void write_spool_em_results(
    const std::vector<std::string> &spoolPaths, const SpoolEMFit &fit,
    const SpoolSampleMixtureFit &sampleMixtureFit,
    const ChimeraClassify::EMOptions &options,
    const ChimeraClassify::DecisionConfig &decisionConfig,
    const ChimeraClassify::TaxDict &tax,
    const ChimeraClassify::PresenceDecision *presenceDecision,
    const ChimeraClassify::NcbiTaxdump *ncbiTaxdump,
    const chimera::presence::CoverageMeta &coverageMeta,
    const ChimeraClassify::ClassifyConfig &config,
    const LocalResolutionCallMap *localCalls, double sampleDivergence,
    ChimeraClassify::FileInfo &fileInfo,
    std::unordered_map<uint32_t, double> *profileClassTaxonPriors = nullptr) {
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
            config.local_resolution_divergence_threshold, partStats[part_idx],
            postemPrimaryEvidenceAggregates[part_idx],
            postemDecisionEvidenceAggregates[part_idx]);
      });

  try {
    merge_classify_output_parts(outputFile, partPaths, partStats, fileInfo);
  } catch (...) {
    cleanup_part_paths(partPaths);
    throw;
  }
  std::unordered_map<uint32_t, double> decisionSpeciesCounts;
  uint64_t decisionClassifiedReads = 0;
  for (const auto &stats : partStats) {
    decisionClassifiedReads += stats.classified;
    for (const auto &[species, value] : stats.decision_species_counts) {
      decisionSpeciesCounts[species] += value;
    }
  }
  SeqProfileFit classifyDecisionProfile = make_profile_from_species_counts(
      decisionSpeciesCounts, fileInfo.classifiedNum + fileInfo.unclassifiedNum,
      decisionClassifiedReads, coverageMeta, ncbiTaxdump);
  if (profileClassTaxonPriors != nullptr) {
    profileClassTaxonPriors->clear();
    if (ncbiTaxdump != nullptr && ncbiTaxdump->enabled()) {
      for (const SeqProfileRow &row : classifyDecisionProfile.rows) {
        const double value = row.callable_normalized_abundance > 0.0
                                 ? row.callable_normalized_abundance
                                 : row.sequence_abundance;
        if (row.species_taxid == 0 || !(value > 0.0)) {
          continue;
        }
        const uint32_t species = ncbiTaxdump->to_species(row.species_taxid);
        const uint32_t genus = ncbiTaxdump->to_genus(species);
        if (species != 0) {
          (*profileClassTaxonPriors)[species] += value;
        }
        if (genus != 0) {
          (*profileClassTaxonPriors)[genus] += value;
        }
      }
    }
  }
  SpeciesProfileMasses localmixMasses =
      merge_localmix_profile_masses(partStats);
  SeqProfileFit localmixProfile =
      make_profile_from_species_masses(localmixMasses, coverageMeta,
                                       ncbiTaxdump);
  const SeqProfileFit *profileOutputFit = &classifyDecisionProfile;
  const SeqProfileFit *profileOutputLocalmixFit = &localmixProfile;
  const char *profileResponseSourceOverride = nullptr;
  PrimaryProfileScale profileOutputScale =
      PrimaryProfileScale::CallableNormalized;
  if (localCalls != nullptr && sampleDivergence >=
                                  config.local_resolution_divergence_threshold) {
    profileOutputLocalmixFit = nullptr;
    profileResponseSourceOverride = "lpc_final_readcount";
    profileOutputScale = PrimaryProfileScale::SequenceAbundance;
  }
  write_classifier_response_reportable_profile_outputs(
      resolve_profile_variant_output_path(
          outputFile, "ChimeraProfile.classifier_response.tsv",
          ".profile.classifier_response.tsv"),
      resolve_profile_variant_output_path(
          outputFile, "ChimeraProfile.classifier_response.cami.tsv",
          ".profile.classifier_response.cami.tsv"),
      resolve_profile_variant_output_path(
          outputFile, "ChimeraProfile.classifier_response.debug.tsv",
          ".profile.classifier_response.debug.tsv"),
      resolve_profile_variant_output_path(
          outputFile, "ChimeraProfile.classifier_response.trace.tsv",
          ".profile.classifier_response.trace.tsv"),
      *profileOutputFit, ncbiTaxdump, coverageMeta, profileOutputLocalmixFit,
      profileResponseSourceOverride, profileOutputScale);
  std::cout << "[classify][profile] classifier-response reportable profile"
            << " path="
            << resolve_profile_variant_output_path(
                   outputFile, "ChimeraProfile.classifier_response.tsv",
                   ".profile.classifier_response.tsv")
            << "\n";
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
  std::unique_ptr<NcbiTaxdump> profileTaxdump;
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
  if (taxonomyKind == "ncbi") {
    ncbiTaxdump = maybe_load_profiledb_taxonomy_for_profile(config.dbFile);
    if (!ncbiTaxdump) {
      ncbiTaxdump = maybe_load_ncbi_taxdump(taxonomyKind);
    }
  }
  if (!ncbiTaxdump) {
    ncbiTaxdump = maybe_load_tax_tsv_for_profile(config.dbFile);
  }
  if (ncbiTaxdump && ncbiTaxdump->enabled()) {
    weightCtx.ncbiTaxdump = ncbiTaxdump.get();
  }
  profileTaxdump = maybe_load_profiledb_taxonomy_for_profile(config.dbFile);
  if (!profileTaxdump) {
    profileTaxdump = maybe_load_ncbi_taxdump_for_profile();
  }
  if (!profileTaxdump) {
    profileTaxdump = maybe_load_tax_tsv_for_profile(config.dbFile);
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
          (localResolutionDivergence >=
           config.local_resolution_divergence_threshold) &&
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
          const std::filesystem::path shardManifestPath =
              resolvedLocalArtifacts->shard_manifest_path;
          if (!std::filesystem::exists(repMetadataPath)) {
            throw std::runtime_error(
                "Local resolution metadata is missing next to database: " +
                repMetadataPath.string());
          }
          if (!std::filesystem::exists(shardManifestPath)) {
            throw std::runtime_error(
                "Local resolution shard manifest is missing next to database: " +
                shardManifestPath.string());
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
            localRequest.shard_manifest_file = shardManifestPath.string();
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
                           profileTaxdump && profileTaxdump->enabled()
                               ? profileTaxdump.get()
                               : weightCtx.ncbiTaxdump,
                           coverageMeta, config,
                           localResolutionCallPtr, localResolutionDivergence,
                           fileInfo, nullptr);
  }
  {
    const std::string outputFile = resolve_tsv_output_path(config.outputFile);
    copy_profile_file_or_throw(
        resolve_profile_variant_output_path(
            outputFile, "ChimeraProfile.classifier_response.tsv",
            ".profile.classifier_response.tsv"),
        resolve_profile_output_path(outputFile));
    copy_profile_file_or_throw(
        resolve_profile_variant_output_path(
            outputFile, "ChimeraProfile.classifier_response.cami.tsv",
            ".profile.classifier_response.cami.tsv"),
        resolve_profile_cami_output_path(outputFile));
    copy_profile_file_or_throw(
        resolve_profile_variant_output_path(
            outputFile, "ChimeraProfile.classifier_response.debug.tsv",
            ".profile.classifier_response.debug.tsv"),
        resolve_native_profile_debug_output_path(outputFile));
    copy_profile_file_or_throw(
        resolve_profile_variant_output_path(
            outputFile, "ChimeraProfile.classifier_response.trace.tsv",
            ".profile.classifier_response.trace.tsv"),
        resolve_native_profile_trace_output_path(outputFile));
    std::cout << "[classify][profile] classifier-response reportable default"
              << " path=" << resolve_profile_output_path(outputFile) << "\n";
  }
  std::vector<chimera::presence::CoverageEntry>().swap(coverageMeta.entries);
  if (keepClassifySpool) {
    std::cout << "[classify][debug] keeping spool=" << spoolDir.string()
              << "\n";
  } else {
    std::filesystem::remove_all(spoolDir);
  }
}

} // namespace ChimeraClassify
