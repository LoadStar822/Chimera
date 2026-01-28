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
		              << ", freq_model.quantile=" << coverageMeta.freq_model.quantile
		              << ", freq.df_high_threshold=" << weightCtx.freqStats.df_high_threshold
		              << ", freq.df_max_observed=" << weightCtx.freqStats.df_max_observed
		              << ", freq.nonzero_counters=" << weightCtx.freqStats.nonzero_counters
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
  PresenceSummary presenceSummary(config.presence_breadth_bits);
  PresenceSummary *presencePtr = &presenceSummary;

  size_t avg_len_hint = 0;
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
                       feature_min_len, weightCtx, nullptr);
    probeProducer.join();

    if (probeInfo.sequenceNum > 0) {
      avg_len_hint = probeInfo.bpLength / probeInfo.sequenceNum;
    }

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
      std::cout << "[classify][auto] lowdiv=1 overrides: em_conf_power="
                << config.em_conf_power
                << " dump_post_topk=" << config.dump_post_topk << std::endl;
    } else {
      std::cout << "[classify][auto] lowdiv=0 use default classify config"
                << std::endl;
    }
  }

  FileInfo fileInfo;
  if (avg_len_hint > 0) {
    fileInfo.avgLen = avg_len_hint;
  }
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
                     feature_min_len, weightCtx, presencePtr);
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

  PresenceFilterStats filterStats{};
  if (!config.em) {
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

  // Optional: saturate per-read sample weights to avoid extreme long tails
  // dominating EM in high-diversity datasets (e.g., CAMI-long).
  //
  // Default: auto (enabled only when weights are present and highly skewed).
  {
    const bool low_div = config.low_div_active;
    bool applied = false;
    double scale = 1.0;
    std::string reason = "disabled";

    std::vector<double> raw;
    raw.reserve(classifyResults.size());
    for (const auto &res : classifyResults) {
      if (res.sample_weight > 0.0) {
        raw.push_back(res.sample_weight);
      }
    }
    if (low_div) {
      reason = "low_div";
    } else if (raw.empty()) {
      reason = "no_weights";
    } else {
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
      if (!skewed) {
        reason = "not_skewed";
      } else if (!(med > 0.0) || !(std::log1p(med) > 0.0)) {
        reason = "bad_median";
      } else {
        scale = med / std::log1p(med);
        for (auto &res : classifyResults) {
          if (res.sample_weight > 0.0) {
            res.sample_weight = scale * std::log1p(res.sample_weight);
          }
        }
        applied = true;
        reason = "ok";

        if (config.verbose) {
          const double eff_p50 = scale * std::log1p(med);
          const double eff_p90 = scale * std::log1p(p90);
          const double eff_p99 = scale * std::log1p(p99);
          const double eff_max = scale * std::log1p(mx);

          auto oldFlags = std::cout.flags();
          auto oldPrecision = std::cout.precision();
          std::cout.setf(std::ios::fixed);
          std::cout << std::setprecision(4);
          std::cout << "[em][wsat] requested=auto"
                    << " mode=log1p"
                    << " applied=" << (applied ? 1 : 0)
                    << " low_div=" << (low_div ? 1 : 0)
                    << " avgLen=" << fileInfo.avgLen
                    << " raw(p50/p90/p99/max)=" << med << '/' << p90 << '/'
                    << p99 << '/' << mx
                    << " raw_ratio(max/med)=" << (mx / denom)
                    << " raw_ratio(p99/med)=" << (p99 / denom)
                    << " eff(p50/p90/p99/max)=" << eff_p50 << '/' << eff_p90
                    << '/' << eff_p99 << '/' << eff_max
                    << " scale=" << scale << " reason=" << reason
                    << std::endl;
          std::cout.flags(oldFlags);
          std::cout.precision(oldPrecision);
        }
      }
    }
    (void)applied;
    (void)scale;
    (void)reason;
  }

  if (config.em) {
    auto EMstart = std::chrono::high_resolution_clock::now();
    std::cout << "Running EM algorithm..." << std::endl;
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
		    // Auto-tune post_pi_min for high-diversity samples:
		    // a fixed global pi cut (e.g., 5e-4) can hard-kill low-abundance true taxa,
		    // causing massive unclassified reads on datasets like CAMI-long/short.
		    // We only relax the threshold under conservative guards (low-div stays
		    // untouched), and further apply an avgLen-driven smooth curve for short
		    // samples to reduce posterior_weight/em_post rejects.
		    double tuned_post_pi_min = config.post_pi_min;
		    double tuned_post_pi_min_weights = tuned_post_pi_min;
		    double top10_mass = 0.0;
		    double top50_mass = 0.0;
		    double tail_base = 0.0;
		    bool weight_stats_ok = false;
		    bool head_heavy = false;
		    if (!config.low_div_active && tuned_post_pi_min > 0.0 &&
		        !classWeights.empty()) {
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
		        weight_stats_ok = true;
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
		        tail_base = tail_mass(base);
		        const double tail_trigger = 0.02; // relax only when base kills >2% mass
		        // Extra low-div guard: avoid relaxing global pi in head-heavy samples,
		        // otherwise ambiguous reads can get forced into random taxa, exploding FPs.
		        const double head_guard10 = 0.65;
		        const double head_guard50 = 0.90;
			        head_heavy =
			            (top10_mass >= head_guard10) || (top50_mass >= head_guard50);
			        if (tail_base > tail_trigger && !head_heavy) {
			          double target_tail = 0.01; // aim to kill <=1% mass by global pi cut
			          // For long-read, high-diversity samples we can be slightly more
			          // permissive to avoid pruning true low-abundance taxa.
			          if (fileInfo.avgLen >= 2000) {
			            target_tail = 0.005; // kill <=0.5% mass
			          }
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
		      }
		    }
		    tuned_post_pi_min_weights = tuned_post_pi_min;
		
		    constexpr double kAutoPostPiMinLoDefault = 1e-4;
		    constexpr size_t kAutoPostPiMinL0 = 2000;
		    constexpr size_t kAutoPostPiMinL1 = 800;
		    double auto_pi_lo = kAutoPostPiMinLoDefault;
		
		    auto len_tune = ChimeraClassify::tune_post_pi_min_by_avg_len(
		        tuned_post_pi_min_weights, auto_pi_lo, fileInfo.avgLen,
		        config.low_div_active || head_heavy, kAutoPostPiMinL0,
		        kAutoPostPiMinL1);
		    tuned_post_pi_min = len_tune.tuned;
		
		    std::string reason = "base";
		    if (!(config.post_pi_min > 0.0)) {
		      reason = "disabled";
		    } else if (config.low_div_active) {
		      reason = "lowdiv";
		    } else {
		      const bool weights_tuned =
		          (tuned_post_pi_min_weights > 0.0 &&
		           tuned_post_pi_min_weights < config.post_pi_min);
		      if (weights_tuned && len_tune.applied) {
		        reason = "weights+avgLen";
		      } else if (len_tune.applied) {
		        reason = "avgLen";
		      } else if (weights_tuned) {
		        reason = "weights";
		      }
		    }
		
		    {
		      auto oldFlags = std::cout.flags();
		      auto oldPrecision = std::cout.precision();
		      std::cout.setf(std::ios::scientific);
		      std::cout << std::setprecision(6);
		      std::cout << "[em][auto] post_pi_min=" << tuned_post_pi_min
		                << " reason=" << reason
		                << " low_div=" << (config.low_div_active ? 1 : 0)
		                << " avgLen=" << fileInfo.avgLen
		                << " base=" << config.post_pi_min
		                << " tuned_w=" << tuned_post_pi_min_weights
		                << " pi_lo=" << auto_pi_lo;
		      if (len_tune.applied) {
		        std::cout << " t=" << std::fixed << std::setprecision(3) << len_tune.t;
		      }
		      if (weight_stats_ok) {
		        std::cout << std::defaultfloat
		                  << " top10_mass=" << top10_mass
		                  << " top50_mass=" << top50_mass
		                  << " tail_mass=" << tail_base
		                  << " head_heavy=" << (head_heavy ? 1 : 0);
		      }
		      std::cout << std::endl;
		      std::cout.flags(oldFlags);
		      std::cout.precision(oldPrecision);
		    }
		
		    decisionConfig.min_class_weight = tuned_post_pi_min;
		    // postEmDecision tail controls:
		    // - They govern how much posterior mass is allowed to fan out into multiple
		    //   taxids for *taxidCount* (abundance/presence is very sensitive).
		    // - POST_TOPK dump uses the pruned+renormalized posterior view (same as decision).

	    // Default/auto behavior: when we relax global pi to recover recall in
	    // high-diversity samples, tighten per-read tail allocation to avoid
	    // exploding the number of tiny non-zero taxa.
	    double auto_head_mass = 0.95;
	    uint32_t auto_max_taxa = 8;
	    if (tuned_post_pi_min > 0.0 && tuned_post_pi_min < config.post_pi_min) {
	      auto_head_mass = 0.90;
	      auto_max_taxa = 5;
	    }
		    decisionConfig.posterior_head_mass = auto_head_mass;
		    decisionConfig.posterior_max_taxa = auto_max_taxa;
		    decisionConfig.allow_fallback_on_reject = !config.low_div_active;

		    if (config.verbose) {
		      std::cout << "PostEM tail: head_mass=" << decisionConfig.posterior_head_mass
		                << " max_taxa=" << decisionConfig.posterior_max_taxa
		                << " fallback_on_reject=" << (decisionConfig.allow_fallback_on_reject ? 1 : 0)
		                << " fallback_gap_min=" << decisionConfig.fallback_gap_min
		                << std::endl;
		    }

    postEmDecision(classifyResults, decisionConfig, classWeights, tax,
                   &presenceDecision, weightCtx.ncbiTaxdump, fileInfo.avgLen);
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
