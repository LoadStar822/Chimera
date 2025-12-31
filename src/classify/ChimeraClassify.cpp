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
#include <cctype>
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
maybe_load_ncbi_taxdump(const std::string &taxonomyKind, bool verbose) {
  if (taxonomyKind != "ncbi") {
    return nullptr;
  }
  const char *env_dir = std::getenv("CHIMERA_NCBI_TAXDUMP_DIR");
  std::filesystem::path base =
      env_dir && *env_dir ? env_dir
                          : "/mnt/sda/tianqinzhong/project/SimDataset/taxdump";
  std::filesystem::path nodes = base / "nodes.dmp";
  if (!std::filesystem::exists(nodes)) {
    if (verbose) {
      std::cerr << "Warning: nodes.dmp not found, skip strain->species collapse: "
                << nodes.string() << std::endl;
    }
    return nullptr;
  }

  std::ifstream is(nodes);
  if (!is.is_open()) {
    if (verbose) {
      std::cerr << "Warning: cannot open nodes.dmp, skip strain->species collapse: "
                << nodes.string() << std::endl;
    }
    return nullptr;
  }

  auto tax = std::make_unique<ChimeraClassify::NcbiTaxdump>();
  std::string line;
  uint32_t max_id = 0;
  size_t species_nodes = 0;
  size_t genus_nodes = 0;
  size_t parsed = 0;
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
    try {
      tid = static_cast<uint32_t>(std::stoul(fields[0]));
      parent = static_cast<uint32_t>(std::stoul(fields[1]));
    } catch (...) {
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
    if (is_sp) {
      ++species_nodes;
    }
    const bool is_g = (rank == "genus");
    tax->is_genus[tid] = is_g ? 1 : 0;
    if (is_g) {
      ++genus_nodes;
    }
    if (tid > max_id) {
      max_id = tid;
    }
    ++parsed;
  }

  if (!tax->enabled()) {
    return nullptr;
  }
  if (verbose) {
    std::cout << "Loaded NCBI nodes.dmp for strain->species collapse: "
              << "nodes=" << parsed << ", max_id=" << max_id
              << ", species_nodes=" << species_nodes
              << ", genus_nodes=" << genus_nodes << std::endl;
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

  if (config.verbose) {
    std::cout << config << std::endl;
  }
  omp_set_num_threads(config.threads);
  auto TotalclassifyStart = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<std::string>> indexToTaxid;
  chimera::imcf::InterleavedMergedCuckooFilter imcf;
  ChimeraBuild::IMCFConfig imcfConfig;
  chimera::presence::CoverageMeta coverageMeta;
  auto indexStatus =
      loadFilter(config.dbFile, imcf, imcfConfig, indexToTaxid, &coverageMeta);
  WeightingContext weightCtx;
  long long rebuildActiveMs = 0;
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
    } catch (const std::exception &ex) {
      std::cerr << "Warning: 无法加载频率 Sketch，回退到旧权重模型: " << ex.what()
                << std::endl;
    }
  }
  size_t uniq_nonzero = 0;
  uint64_t uniq_max = 0;
  size_t density_nonzero = 0;
  double density_max = 0.0;
  size_t genome_nonzero = 0;
  size_t total_nonzero = 0;
  if (!coverageMeta.entries.empty()) {
    for (const auto &entry : coverageMeta.entries) {
      if (entry.unique_signatures > 0) {
        ++uniq_nonzero;
        if (entry.unique_signatures > uniq_max) {
          uniq_max = entry.unique_signatures;
        }
      }
      if (entry.unique_density > 0.0) {
        ++density_nonzero;
        if (entry.unique_density > density_max) {
          density_max = entry.unique_density;
        }
      }
      if (entry.total_signatures > 0) {
        ++total_nonzero;
      }
      if (entry.genome_length > 0) {
        ++genome_nonzero;
      }
    }
  }

  bool freq_trusted = true;
  if (coverageMeta.freq_model.enabled() && !coverageMeta.entries.empty() &&
      uniq_nonzero == 0) {
    freq_trusted = false;
    std::cout << "[presence][auto] unique_gate=off reason=unique_nonzero0 "
              << "entries=" << coverageMeta.entries.size()
              << " total_nonzero=" << total_nonzero
              << " genome_nonzero=" << genome_nonzero << std::endl;
    if (total_nonzero > 0 && genome_nonzero > 0) {
      std::cout << "[presence][auto] exposure_fallback=total_density "
                << "reason=unique_nonzero0 total_nonzero=" << total_nonzero
                << " genome_nonzero=" << genome_nonzero << std::endl;
    }
  }
  weightCtx.freq_trusted = freq_trusted;

  if (config.verbose) {
    std::cout << "DB presence meta: entries=" << coverageMeta.entries.size()
              << ", ref_read_len=" << coverageMeta.ref_read_length
              << ", span=" << coverageMeta.effective_span
              << ", unique_deg=" << coverageMeta.unique_deg_threshold
              << ", freq_model.enabled=" << coverageMeta.freq_model.enabled()
              << ", freq_model.depth=" << coverageMeta.freq_model.depth
              << ", freq_model.width=" << coverageMeta.freq_model.width
              << ", counters=" << coverageMeta.freq_model.counters.size()
              << std::endl;
    std::cout << "DB presence meta stats: unique_nonzero=" << uniq_nonzero
              << ", unique_max=" << uniq_max
              << ", density_nonzero=" << density_nonzero
              << ", density_max=" << std::scientific << density_max
              << std::defaultfloat << ", genome_nonzero=" << genome_nonzero
              << std::endl;
  }

  std::unordered_map<std::string, double> sampleWeights;
  if (!config.weight_map_file.empty()) {
    std::cout << "Loading sample weight map: " << config.weight_map_file
              << std::endl;
    sampleWeights = load_weight_map_file(config.weight_map_file);
    if (sampleWeights.empty()) {
      std::cerr << "Warning: weight map parsed empty; fallback to default weights."
                << std::endl;
    } else {
      std::cout << "Loaded " << sampleWeights.size()
                << " sample weights." << std::endl;
      weightCtx.sampleWeights = &sampleWeights;
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
  // Optional NCBI strain/subspecies -> species collapse.
  // This helps avoid candidate saturation by many strain taxids, which can
  // prevent sister species from entering pre-EM/posterior lists.
  ncbiTaxdump = maybe_load_ncbi_taxdump(config.taxonomyKind, config.verbose);
  if (ncbiTaxdump && ncbiTaxdump->enabled()) {
    weightCtx.ncbiTaxdump = ncbiTaxdump.get();
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

  auto stage_t0 = std::chrono::steady_clock::now();
  auto stage_last = stage_t0;
  auto log_stage = [&](const char *stage) {
    auto now = std::chrono::steady_clock::now();
    auto delta_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - stage_last)
            .count();
    auto total_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - stage_t0)
            .count();
    uint64_t rss_kb = current_rss_kb();
    std::cout << "[classify] stage=" << stage
              << " rss_gib=" << std::fixed << std::setprecision(2)
              << kb_to_gib(rss_kb) << std::defaultfloat
              << " delta_ms=" << delta_ms << " total_ms=" << total_ms
              << std::endl;
    stage_last = now;
  };

  log_stage("before_build_tax_dict");
  const TaxDict tax = build_tax_dict(indexToTaxid);
  std::cout << "[classify] tax_dict tid_count=" << tax.id2str.size()
            << " bins=" << tax.idx2id.size() << std::endl;
  log_stage("after_build_tax_dict");
  std::vector<uint32_t> tid2speciesGroup;
  if (config.deg_by_species && weightCtx.ncbiTaxdump &&
      weightCtx.ncbiTaxdump->enabled()) {
    tid2speciesGroup.resize(tax.id2str.size(), 0);
    auto parse_u32 = [](const std::string &s, uint32_t &out) -> bool {
      if (s.empty()) {
        return false;
      }
      for (unsigned char c : s) {
        if (!std::isdigit(c)) {
          return false;
        }
      }
      try {
        unsigned long v = std::stoul(s);
        if (v > std::numeric_limits<uint32_t>::max()) {
          return false;
        }
        out = static_cast<uint32_t>(v);
        return true;
      } catch (...) {
        return false;
      }
    };

    size_t numeric = 0;
    size_t mapped = 0;
    for (uint32_t tid_id = 0; tid_id < tax.id2str.size(); ++tid_id) {
      const std::string &taxid = tax.id2str[tid_id];
      uint32_t tid = 0;
      if (parse_u32(taxid, tid)) {
        ++numeric;
        uint32_t sid = weightCtx.ncbiTaxdump->to_species(tid);
        tid2speciesGroup[tid_id] = sid;
        ++mapped;
      } else {
        // Stable synthetic id to avoid collisions with real numeric taxids.
        tid2speciesGroup[tid_id] = 0x80000000u | tid_id;
      }
    }
    weightCtx.tid2speciesGroup = &tid2speciesGroup;
    if (config.verbose) {
      std::cout << "NCBI tid->species group mapping for deg: tids="
                << tax.id2str.size() << ", numeric=" << numeric
                << ", mapped=" << mapped << std::endl;
    }
  }
  std::vector<uint32_t> tid2speciesRep;
  if ((config.collapse_strain_hits || config.collapse_strain_candidates) &&
      weightCtx.ncbiTaxdump &&
      weightCtx.ncbiTaxdump->enabled()) {
    tid2speciesRep.resize(tax.id2str.size());
    for (uint32_t i = 0; i < tid2speciesRep.size(); ++i) {
      tid2speciesRep[i] = i;
    }

    auto parse_u32 = [](const std::string &s, uint32_t &out) -> bool {
      if (s.empty()) {
        return false;
      }
      for (unsigned char c : s) {
        if (!std::isdigit(c)) {
          return false;
        }
      }
      try {
        unsigned long v = std::stoul(s);
        if (v > std::numeric_limits<uint32_t>::max()) {
          return false;
        }
        out = static_cast<uint32_t>(v);
        return true;
      } catch (...) {
        return false;
      }
    };

    robin_hood::unordered_flat_map<uint32_t, uint32_t> species2rep;
    species2rep.reserve(tax.id2str.size() / 2 + 1);
    size_t numeric = 0;
    size_t collapsed = 0;
    for (uint32_t tid_id = 0; tid_id < tax.id2str.size(); ++tid_id) {
      const std::string &taxid = tax.id2str[tid_id];
      uint32_t tid = 0;
      if (!parse_u32(taxid, tid)) {
        continue;
      }
      ++numeric;
      uint32_t sid = weightCtx.ncbiTaxdump->to_species(tid);
      auto it = species2rep.find(sid);
      if (it == species2rep.end()) {
        species2rep.emplace(sid, tid_id);
        tid2speciesRep[tid_id] = tid_id;
      } else {
        tid2speciesRep[tid_id] = it->second;
        if (it->second != tid_id) {
          ++collapsed;
        }
      }
    }
    weightCtx.tid2speciesRep = &tid2speciesRep;
    if (config.verbose) {
      std::cout << "NCBI tid->species representative mapping: tids="
                << tax.id2str.size() << ", numeric=" << numeric
                << ", species=" << species2rep.size()
                << ", collapsed=" << collapsed << std::endl;
    }
  }
  if (config.decoy_reps == 0 && !(config.presence_noise > 0.0) &&
      config.verbose) {
    std::cerr << "Warning: decoy_reps=0，覆盖模型的噪声将退回默认 μ=1e-4；建议 decoy_reps>=1"
              << std::endl;
  }
  PresenceSummary presenceSummary(static_cast<size_t>(config.decoy_reps),
                                  config.presence_breadth_bits);
  PresenceSummary *presencePtr = &presenceSummary;
  uint64_t presenceSeed = 0;
  std::hash<std::string> hasher;
  // Stable seed: avoid run_dir/output path affecting decoy randomness so metrics
  // are reproducible across repeated runs on the same dataset+DB.
  if (!config.singleFiles.empty()) {
    presenceSeed ^= hasher(config.singleFiles.front());
  }
  if (!config.pairedFiles.empty()) {
    presenceSeed ^= (hasher(config.pairedFiles.front()) << 1);
  }
  presenceSeed ^= (hasher(config.dbFile) << 2);
  presenceSeed ^= static_cast<uint64_t>(imcfConfig.fpSalt);

  if (config.low_div_auto && config.low_div_probe_reads > 0) {
    ClassifyConfig probe_config = config;
    probe_config.max_reads = config.low_div_probe_reads;
    probe_config.verbose = false;

    FileInfo probeInfo;
    std::vector<moodycamel::ConcurrentQueue<batchReads>> probeQueues(
        static_cast<size_t>(std::max<uint16_t>(1, config.threads)));
    std::vector<classifyResult> probeResults;
    std::atomic<bool> probe_done{false};
    std::thread probeProducer([&]() {
      parseReads(probeQueues, probe_config, probeInfo);
      probe_done.store(true, std::memory_order_release);
    });

    classify_streaming(imcfConfig, probeQueues, probe_config, imcf, indexToTaxid,
                       tax, probeResults, probeInfo, probe_done, feature_params,
                       feature_min_len, weightCtx, nullptr, presenceSeed);
    probeProducer.join();

    if (!probeResults.empty()) {
      std::sort(probeResults.begin(), probeResults.end(),
                [](const classifyResult &a, const classifyResult &b) {
                  return a.id < b.id;
                });
    }

    bool probe_em = false;
    if (config.em && !probeResults.empty()) {
      EMOptions options;
      options.temp = 1.05;
      options.prior_strength = config.em_prior_strength;
      options.coexist_penalty = config.em_coexist_penalty;
      options.prune_ratio = config.em_prune_ratio;
      options.conf_power = config.em_conf_power;
      auto [posterior, weights] = EMAlgorithm(
          probeResults, config.emIter, config.emThreshold, options, nullptr);
      probeResults = std::move(posterior);
      probe_em = true;
    }

    std::unordered_map<std::string, double> species_counts;
    species_counts.reserve(probeResults.size() / 4 + 8);
    double unclassified_weight = 0.0;

    auto parse_u32 = [](const std::string &s, uint32_t &out) -> bool {
      if (s.empty()) {
        return false;
      }
      for (unsigned char c : s) {
        if (!std::isdigit(c)) {
          return false;
        }
      }
      try {
        unsigned long v = std::stoul(s);
        if (v > std::numeric_limits<uint32_t>::max()) {
          return false;
        }
        out = static_cast<uint32_t>(v);
        return true;
      } catch (...) {
        return false;
      }
    };

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
        if (parse_u32(top, tid)) {
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
    constexpr std::size_t kLowDivTopK = 10;
    constexpr double kLowDivTopMass = 0.75;
    constexpr double kLowDivEffSpeciesMax = 32.0;
    std::size_t unclassified_reads =
        static_cast<std::size_t>(std::llround(unclassified_weight));
    LowDivStats stats =
        compute_low_div_stats(counts, unclassified_reads, kLowDivTopK);
    bool low_div = is_low_diversity(stats, kLowDivTopMass, kLowDivEffSpeciesMax);

    std::cout << "[classify][auto] lowdiv probe: reads=" << probeInfo.sequenceNum
              << " eff_species=" << std::fixed << std::setprecision(2)
              << stats.eff_species << " top_mass=" << (stats.top_mass * 100.0)
              << " low_div=" << (low_div ? 1 : 0)
              << " thresh_top_mass=" << (kLowDivTopMass * 100.0)
              << " eff_species_max=" << kLowDivEffSpeciesMax << std::defaultfloat
              << std::endl;

    if (low_div) {
      apply_low_div_overrides(config);
      std::cout << "[classify][auto] lowdiv=1 overrides: preEmTopK="
                << config.preEmTopK
                << " exclusive_gamma=" << config.exclusive_gamma
                << " em_conf_power=" << config.em_conf_power
                << " em_coexist_penalty=" << config.em_coexist_penalty
                << " dump_post_topk=" << config.dump_post_topk << std::endl;
    } else {
      std::cout << "[classify][auto] lowdiv=0 use default classify config"
                << std::endl;
    }
  }

  FileInfo fileInfo;
  seqan3::contrib::bgzf_thread_count = config.threads;
  std::vector<moodycamel::ConcurrentQueue<batchReads>> readQueues(
      static_cast<size_t>(std::max<uint16_t>(1, config.threads)));
  std::vector<classifyResult> classifyResults;
  std::unordered_map<std::string, double> classWeights;
  bool posteriorModelUsed = false;

  auto readStart = std::chrono::high_resolution_clock::now();
  std::cout << "Reading input files..." << std::endl;
  std::atomic<bool> producer_done{false};
  auto readEnd = readStart;
  std::thread producer([&]() {
    parseReads(readQueues, config, fileInfo);
    readEnd = std::chrono::high_resolution_clock::now();
    if (config.verbose) {
      auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
          readEnd - readStart);
      std::cout << "\nRead time: ";
      print_classify_time(readDuration.count());
    }
    producer_done.store(true, std::memory_order_release);
  });

  auto classifyStart = std::chrono::high_resolution_clock::now();
  std::cout << "Classifying sequences by imcf (feature=" << config.feature
            << ")..." << std::endl;
  classify_streaming(imcfConfig, readQueues, config, imcf, indexToTaxid, tax,
                     classifyResults, fileInfo, producer_done, feature_params,
                     feature_min_len, weightCtx, presencePtr, presenceSeed);
  auto classifyEnd = std::chrono::high_resolution_clock::now();
  auto classifyDuration =
      std::chrono::duration_cast<std::chrono::milliseconds>(classifyEnd -
                                                            classifyStart);
  producer.join();
  // Determinism: classify_streaming merges thread-local batches in the order
  // threads finish, which is non-deterministic. Sort by read/contig id before
  // downstream EM/post-processing so repeated runs are stable.
  std::sort(classifyResults.begin(), classifyResults.end(),
            [](const classifyResult &a, const classifyResult &b) {
              return a.id < b.id;
            });
  const char *dump_pre_em_path = std::getenv("CHIMERA_DUMP_PRE_EM");
  if (dump_pre_em_path && *dump_pre_em_path) {
    std::ofstream dump_pre(dump_pre_em_path);
    if (dump_pre.good()) {
      for (const auto &result : classifyResults) {
        dump_pre << result.id;
        if (result.taxidCount.empty()) {
          dump_pre << "\tunclassified:1";
        } else {
          for (const auto &kv : result.taxidCount) {
            dump_pre << '\t' << kv.first << ':' << kv.second;
          }
        }
        dump_pre << '\n';
      }
    }
  }
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

  size_t presenceTotalReads = fileInfo.sequenceNum;
  size_t presenceMeanReadLen = fileInfo.avgLen;
  if (weightCtx.has_sample_weights()) {
    presenceMeanReadLen = (coverageMeta.ref_read_length > 0)
                              ? static_cast<size_t>(coverageMeta.ref_read_length)
                              : static_cast<size_t>(150);
  }

  PresenceDecision presenceDecision = evaluate_presence_coverage(
      presenceSummary, tax, config, coverageMeta, presenceTotalReads,
      presenceMeanReadLen);
  if (config.verbose) {
    auto oldFlags = std::cout.flags();
    auto oldPrecision = std::cout.precision();
    std::cout << "Presence caller (coverage): tests=" << presenceDecision.tested
              << ", accepted=" << presenceDecision.acceptedCount
              << ", mu=" << std::scientific << std::setprecision(6)
              << presenceDecision.noiseMu << ", pi=" << config.presence_pi
              << ", tau=" << config.presence_tau << std::defaultfloat
              << std::endl;

    if (!presenceDecision.logPosteriors.empty()) {
      std::vector<double> vals;
      vals.reserve(presenceDecision.logPosteriors.size());
      for (const auto &kv : presenceDecision.logPosteriors) {
        vals.push_back(kv.second);
      }
      std::sort(vals.begin(), vals.end());
      auto pick = [&](double q) -> double {
        if (vals.empty()) {
          return 0.0;
        }
        q = std::clamp(q, 0.0, 1.0);
        size_t idx = static_cast<size_t>(
            std::floor(q * static_cast<double>(vals.size() - 1)));
        if (idx >= vals.size()) {
          idx = vals.size() - 1;
        }
        return vals[idx];
      };
      std::cout << "Presence logPosterior: p50=" << std::fixed
                << std::setprecision(3) << pick(0.50)
                << ", p90=" << pick(0.90) << ", p99=" << pick(0.99)
                << ", max=" << vals.back() << std::defaultfloat << std::endl;
    }
    std::cout.flags(oldFlags);
    std::cout.precision(oldPrecision);
  }

  const char *dump_presence_path = std::getenv("CHIMERA_DUMP_PRESENCE");
  if (dump_presence_path && *dump_presence_path &&
      !presenceDecision.logPosteriors.empty()) {
    std::ofstream dump_presence(dump_presence_path);
    if (dump_presence.good()) {
      double pi = std::clamp(config.presence_pi, 1e-9, 1.0 - 1e-6);
      double logPriorOdds = std::log(pi) - std::log1p(-pi);
      dump_presence << "# presence_model=coverage\n";
      dump_presence << "# tested=" << presenceDecision.tested
                    << "\n# accepted=" << presenceDecision.acceptedCount
                    << "\n# mu=" << std::scientific << presenceDecision.noiseMu
                    << "\n# pi=" << std::defaultfloat << pi
                    << "\n# tau=" << presenceDecision.threshold << "\n";
      dump_presence << "taxid\tlogPosterior\tlogBF\tposterior\tlambda_hat\taccepted\n";
      std::vector<std::pair<uint32_t, double>> rows;
      rows.reserve(presenceDecision.logPosteriors.size());
      for (const auto &kv : presenceDecision.logPosteriors) {
        rows.emplace_back(kv.first, kv.second);
      }
      std::sort(rows.begin(), rows.end(),
                [](const auto &a, const auto &b) { return a.second > b.second; });
      dump_presence.setf(std::ios::fixed);
      dump_presence << std::setprecision(6);
      for (const auto &[tid, logPosterior] : rows) {
        std::string taxid = std::to_string(tid);
        if (tid < tax.id2str.size()) {
          taxid = tax.id2str[tid];
        }
        double posteriorProb = 0.0;
        if (auto itp = presenceDecision.posteriors.find(tid);
            itp != presenceDecision.posteriors.end()) {
          posteriorProb = itp->second;
        }
        double lambda_hat = 0.0;
        if (auto itl = presenceDecision.lambdaHats.find(tid);
            itl != presenceDecision.lambdaHats.end()) {
          lambda_hat = itl->second;
        }
        double logBF = logPosterior - logPriorOdds;
        bool accepted =
            presenceDecision.accepted.find(tid) != presenceDecision.accepted.end();
        dump_presence << taxid << '\t' << logPosterior << '\t' << logBF << '\t'
                      << posteriorProb << '\t' << lambda_hat << '\t'
                      << (accepted ? 1 : 0) << '\n';
      }
    }
  }
  PresenceFilterStats filterStats{};
  if (!config.em || config.presence_pre_filter) {
    // Presence filter is a strong precision gate for short-read datasets.
    filterStats = apply_presence_filter(presenceDecision, tax, classifyResults,
                                        fileInfo);
    if (config.verbose &&
        (filterStats.trimmedAssignments > 0 ||
         filterStats.forcedUnclassified > 0)) {
      std::cout << "Presence filter: trimmed " << filterStats.trimmedAssignments
                << " assignments, forced " << filterStats.forcedUnclassified
                << " reads to unclassified" << std::endl;
    }
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
    options.prior_strength = config.em_prior_strength; // 可选先验强度，默认 0
    options.coexist_penalty = config.em_coexist_penalty;
    options.prune_ratio = config.em_prune_ratio;
    options.conf_power = config.em_conf_power;
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
	    // Auto-tune post_pi_min for high-diversity samples:
	    // a fixed global pi cut (e.g., 5e-4) can hard-kill low-abundance true taxa,
	    // causing massive unclassified reads on datasets like CAMI-long.
	    // We only relax the threshold when the EM weight distribution shows a
	    // heavy tail AND the head is not dominating (to avoid low-diversity mocks).
	    double tuned_post_pi_min = config.post_pi_min;
	    if (tuned_post_pi_min > 0.0 && !classWeights.empty()) {
	      std::vector<double> weights;
	      weights.reserve(classWeights.size());
	      double total_mass = 0.0;
	      for (const auto &kv : classWeights) {
	        double w = kv.second;
	        if (w > 0.0) {
	          weights.push_back(w);
	          total_mass += w;
	        }
	      }
	      if (total_mass > 0.0 && weights.size() > 1) {
	        std::sort(weights.begin(), weights.end(),
	                  [](double a, double b) { return a > b; });
	        double top10_mass = 0.0;
	        const size_t topK = std::min<size_t>(10, weights.size());
	        for (size_t i = 0; i < topK; ++i) {
	          top10_mass += weights[i];
	        }
	        top10_mass /= total_mass;
	        double top50_mass = 0.0;
	        const size_t topK50 = std::min<size_t>(50, weights.size());
	        for (size_t i = 0; i < topK50; ++i) {
	          top50_mass += weights[i];
	        }
	        top50_mass /= total_mass;
	
	        auto tail_mass = [&](double pi_min) -> double {
	          if (!(pi_min > 0.0)) {
	            return 0.0;
	          }
	          double tail = 0.0;
	          // weights are sorted desc, so scan from the end.
	          for (size_t i = weights.size(); i-- > 0;) {
	            double w = weights[i];
	            if (w >= pi_min) {
	              break;
	            }
	            tail += w;
	          }
	          return tail / total_mass;
	        };
	
	        const double base = tuned_post_pi_min;
	        const double tail_base = tail_mass(base);
	        const double tail_trigger = 0.02;     // only relax when base kills >2% mass
	        // Low-diversity guard: avoid relaxing global pi in head-heavy samples,
	        // otherwise many ambiguous reads get forced into random taxa, exploding FPs.
	        const double head_guard10 = 0.65;
	        const double head_guard50 = 0.90;
	        const bool head_heavy =
	            (top10_mass >= head_guard10) || (top50_mass >= head_guard50);
	        if (tail_base > tail_trigger && !head_heavy) {
	          const double target_tail = 0.01;    // aim to kill <=1% mass by global pi cut
	          const std::array<double, 4> candidates = {5e-4, 2e-4, 1e-4, 5e-5};
	          double tuned = base;
	          for (double c : candidates) {
	            if (c > base) {
	              continue;
	            }
	            tuned = c;
	            if (tail_mass(c) <= target_tail) {
	              break;
	            }
	          }
	          tuned_post_pi_min = tuned;
	        }
	        if (config.verbose) {
	          std::cout << "Auto post_pi_min: base=" << config.post_pi_min
	                    << " tuned=" << tuned_post_pi_min
	                    << " top10_mass=" << top10_mass
	                    << " top50_mass=" << top50_mass
	                    << " tail_mass=" << tail_base << std::endl;
	        }
	      }
	    }
	    decisionConfig.min_class_weight = tuned_post_pi_min;
	    // postEmDecision tail controls:
	    // - They govern how much posterior mass is allowed to fan out into multiple
	    //   taxids for *taxidCount* (abundance/presence is very sensitive).
	    // - POST_TOPK dump uses the pruned+renormalized posterior view (same as decision).
	    decisionConfig.posterior_min_fraction = config.post_min_fraction;
	    decisionConfig.posterior_power = config.post_power;

	    // Default/auto behavior: when we relax global pi to recover recall in
	    // high-diversity samples, tighten per-read tail allocation to avoid
	    // exploding the number of tiny non-zero taxa.
	    double auto_head_mass = 0.95;
	    uint32_t auto_max_taxa = 8;
	    if (tuned_post_pi_min > 0.0 && tuned_post_pi_min < config.post_pi_min) {
	      auto_head_mass = 0.90;
	      auto_max_taxa = 5;
	    }
		    decisionConfig.posterior_head_mass =
		        (config.post_head_mass > 0.0) ? config.post_head_mass : auto_head_mass;
		    decisionConfig.posterior_max_taxa =
		        (config.post_max_taxa > 0) ? config.post_max_taxa : auto_max_taxa;
		    decisionConfig.allow_fallback_on_reject = !config.low_div_active;

		    if (config.verbose) {
		      std::cout << "PostEM tail: min_fraction=" << decisionConfig.posterior_min_fraction
		                << " power=" << decisionConfig.posterior_power
		                << " head_mass=" << decisionConfig.posterior_head_mass
		                << (config.post_head_mass > 0.0 ? " (cli)" : " (auto)")
		                << " max_taxa=" << decisionConfig.posterior_max_taxa
		                << (config.post_max_taxa > 0 ? " (cli)" : " (auto)")
		                << " fallback_on_reject=" << (decisionConfig.allow_fallback_on_reject ? 1 : 0)
		                << " fallback_gap_min=" << decisionConfig.fallback_gap_min
		                << std::endl;
		    }

    postEmDecision(classifyResults, decisionConfig, classWeights, tax,
                   &presenceDecision, weightCtx.ncbiTaxdump);
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
  auto print_preem_stats = [&]() {
    if (fileInfo.preem_route_reads == 0 && fileInfo.preem_cap_checks == 0 &&
        fileInfo.preem_floor_checks == 0) {
      return;
    }
    double rare_frac = 0.0;
    if (fileInfo.preem_route_reads > 0) {
      rare_frac = static_cast<double>(fileInfo.preem_route_rare_reads) /
                  static_cast<double>(fileInfo.preem_route_reads);
    }
    double cap_frac = 0.0;
    if (fileInfo.preem_cap_checks > 0) {
      cap_frac = static_cast<double>(fileInfo.preem_cap_expanded) /
                 static_cast<double>(fileInfo.preem_cap_checks);
    }
    double floor_frac = 0.0;
    if (fileInfo.preem_floor_checks > 0) {
      floor_frac = static_cast<double>(fileInfo.preem_floor_applied) /
                   static_cast<double>(fileInfo.preem_floor_checks);
    }
    double fill_frac = 0.0;
    if (fileInfo.preem_underfull_fill_checks > 0) {
      fill_frac =
          static_cast<double>(fileInfo.preem_underfull_fill_applied) /
          static_cast<double>(fileInfo.preem_underfull_fill_checks);
    }
    double relax_frac = 0.0;
    if (config.preem_beta_relax && fileInfo.preem_beta_relax_checks > 0) {
      relax_frac = static_cast<double>(fileInfo.preem_beta_relax_applied) /
                   static_cast<double>(fileInfo.preem_beta_relax_checks);
    }
    double relax_thr_drop_avg = 0.0;
    if (config.preem_beta_relax && fileInfo.preem_beta_relax_applied > 0) {
      relax_thr_drop_avg =
          static_cast<double>(fileInfo.preem_beta_relax_thr_drop_sum) /
          static_cast<double>(fileInfo.preem_beta_relax_applied);
    }
    std::cout << "[classify][auto] preem_route_rare="
              << fileInfo.preem_route_rare_reads << "/"
              << fileInfo.preem_route_reads << " (" << std::fixed
              << std::setprecision(3) << rare_frac << ")"
              << " cap_expand=" << fileInfo.preem_cap_expanded << "/"
              << fileInfo.preem_cap_checks << " (" << cap_frac << ")"
              << " floor=" << fileInfo.preem_floor_applied << "/"
              << fileInfo.preem_floor_checks << " (" << floor_frac << ")"
              << " floor_added=" << fileInfo.preem_floor_added
              << " floor_skip_dominant=" << fileInfo.preem_floor_skipped_dominant
              << " floor_filtered_weak=" << fileInfo.preem_floor_filtered_weak
              << " floor_target=" << config.preem_floor_target
              << " floor_min_ratio=" << config.preem_floor_min_ratio
              << " floor_add_min_ratio=" << config.preem_floor_add_min_ratio
              << " keepalive=" << fileInfo.preem_keepalive_applied << "/"
              << fileInfo.preem_keepalive_attempt
              << " keepalive_blk_ratio="
              << fileInfo.preem_keepalive_blocked_low_ratio
              << " keepalive_blk_gain="
              << fileInfo.preem_keepalive_blocked_low_gain
              << " keepalive_blk_abs=" << fileInfo.preem_keepalive_blocked_low_abs
              << " keepalive_min_ratio=" << config.preem_keepalive_min_ratio
              << " keepalive_replace_ratio=" << config.preem_keepalive_replace_ratio
              << " keepalive_abs_min=" << config.preem_keepalive_abs_min
              << " fill="
              << (config.preem_underfull_fill ?
                      (std::to_string(fileInfo.preem_underfull_fill_applied) +
                       "/" +
                       std::to_string(fileInfo.preem_underfull_fill_checks)) :
                      std::string("off"))
              << (config.preem_underfull_fill ?
                      (" (" + std::to_string(fill_frac) + ")") :
                      std::string(""))
              << (config.preem_underfull_fill ?
                      (" fill_added=" +
                       std::to_string(fileInfo.preem_underfull_fill_added)) :
                      std::string(""))
              << (config.preem_underfull_fill ?
                      (" fill_stage2_added=" +
                       std::to_string(fileInfo.preem_underfull_fill_stage2_added)) :
                      std::string(""))
              << " beta_relax="
              << (config.preem_beta_relax ?
                      (std::to_string(fileInfo.preem_beta_relax_applied) + "/" +
                       std::to_string(fileInfo.preem_beta_relax_checks)) :
                      std::string("off"))
              << (config.preem_beta_relax ? (" (" + std::to_string(relax_frac) +
                                            ")") :
                                           std::string(""))
              << (config.preem_beta_relax ? (" beta_relax_thr_drop_avg=" +
                                            std::to_string(relax_thr_drop_avg)) :
                                           std::string(""))
              << (config.preem_beta_relax ? (" beta_relax_suppress_overflow=" +
                                            std::to_string(fileInfo
                                                               .preem_beta_relax_suppressed_overflow)) :
                                           std::string(""))
              << (config.preem_beta_relax ? (" beta_seen=" +
                                            std::to_string(fileInfo
                                                               .preem_beta_relax_seen)) :
                                           std::string(""))
              << (config.preem_beta_relax ? (" beta_dom=" +
                                            std::to_string(fileInfo
                                                               .preem_beta_relax_dom_beta)) :
                                           std::string(""))
              << (config.preem_beta_relax ? (" beta_le_halfk=" +
                                            std::to_string(fileInfo
                                                               .preem_beta_relax_strict_le_halfk)) :
                                           std::string(""))
              << (config.preem_beta_relax ? (" beta_efflt=" +
                                            std::to_string(fileInfo
                                                               .preem_beta_relax_eff_lt_min)) :
                                           std::string(""))
              << (config.preem_beta_relax ? (" beta_user=" +
                                            std::to_string(fileInfo
                                                               .preem_beta_relax_beta_user)) :
                                           std::string(""))
              << " dyn96=" << fileInfo.preem_dynamic_topk_96
              << " near_gap=" << fileInfo.preem_neartie_gap
              << " near_ratio=" << fileInfo.preem_neartie_ratio
              << " ov_bins=" << fileInfo.preem_overflow_topbins
              << " ov_size=" << fileInfo.preem_overflow_size
              << " finalK<=16=" << fileInfo.preem_finalk_le16
              << " finalK17_32=" << fileInfo.preem_finalk_17_32
              << " finalK33_64=" << fileInfo.preem_finalk_33_64
              << " finalK65_96=" << fileInfo.preem_finalk_65_96
              << " finalKeq96=" << fileInfo.preem_finalk_eq96
              << " head_mass_thresh=0.5 base_cap=256 max_cap=512"
              << std::defaultfloat << std::endl;
  };
  if (config.verbose) {
    print_rejects();
    print_preem_stats();
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
