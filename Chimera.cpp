/*
 * -----------------------------------------------------------------------------
 * Filename:      Chimera.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-07-09
 *
 * Last Modified: 2024-11-18
 *
 * Description:
 *  This is the main entry for Chimera
 *
 * Version:
 *  1.3
 * -----------------------------------------------------------------------------
 */
#include <CLI11.hpp>
#include <ChimeraBuild.hpp>
#include <ChimeraClassify.hpp>
#include <buildConfig.hpp>
#include <classifyConfig.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <cctype>
#include <thread>
#include <vector>

#if defined(__unix__) || defined(__APPLE__) || defined(__linux__)
#include <sys/utsname.h>
#endif

#ifdef CHIMERA_VERSION
#define VERSION_INFO CHIMERA_VERSION
#else
#define VERSION_INFO "unknown"
#endif

namespace {

std::string detect_os() {
#if defined(_WIN32)
  return "Windows";
#elif defined(__APPLE__)
  struct utsname info {};
  if (uname(&info) == 0) {
    return std::string(info.sysname) + " " + info.release;
  }
  return "macOS";
#elif defined(__linux__)
  std::ifstream osrelease("/etc/os-release");
  if (osrelease) {
    std::string line;
    while (std::getline(osrelease, line)) {
      constexpr char key[] = "PRETTY_NAME=";
      if (line.rfind(key, 0) == 0) {
        std::string value = line.substr(sizeof(key) - 1);
        if (!value.empty() && value.front() == '"' && value.back() == '"') {
          value = value.substr(1, value.size() - 2);
        }
        if (!value.empty()) {
          return value;
        }
      }
    }
  }
  struct utsname info {};
  if (uname(&info) == 0) {
    return std::string(info.sysname) + " " + info.release;
  }
  return "Linux";
#elif defined(__FreeBSD__)
  struct utsname info {};
  if (uname(&info) == 0) {
    return std::string(info.sysname) + " " + info.release;
  }
  return "FreeBSD";
#else
  return "Unknown";
#endif
}

} // namespace

int main(int argc, char **argv) {
  // Create the main application object
  CLI::App app{"Chimera - A versatile tool for metagenomic classification"};
  ChimeraBuild::BuildConfig buildConfig;
  ChimeraClassify::ClassifyConfig classifyConfig;

  bool show_version = false;
  app.add_flag("-v,--version", show_version, "Show version information");
  // Create subcommands
  auto build = app.add_subcommand("build", "Build a sequence database");
  auto classify = app.add_subcommand("classify", "Classify sequences");

  unsigned int hardware_threads = std::thread::hardware_concurrency();
  if (hardware_threads == 0) {
    hardware_threads = 1;
  }
  const auto max_threads =
      static_cast<unsigned int>(std::numeric_limits<uint16_t>::max());
  if (hardware_threads > max_threads) {
    hardware_threads = max_threads;
  }
  const auto default_threads = static_cast<uint16_t>(
      std::min<unsigned int>(hardware_threads, 192u));

  bool buildQuietRequested = false;
  bool classifyQuietRequested = false;

  // Build
  build
      ->add_option("-i,--input", buildConfig.input_file,
                   "Input file for building")
      ->required()
      ->check(CLI::ExistingFile);
  build
      ->add_option("-o,--output", buildConfig.output_file,
                   "Output file for building")
      ->default_val("ChimeraDB");
  build
      ->add_option("-k,--kmer", buildConfig.kmer_size, "Kmer size for building")
      ->default_val(31)
      ->check(CLI::Range(1, 50));
  build
      ->add_option("-s,--syncmer-s", buildConfig.smer_size,
                   "Syncmer s-mer size (must be < k)")
      ->default_val(16);
  build
      ->add_option("-P,--syncmer-pos", buildConfig.syncmer_position,
                   "Syncmer minimal s-mer offset (0-based)")
      ->default_val(7);
  auto *min_length_option = build
      ->add_option("-l,--min-length", buildConfig.min_length,
                   "Minimum length sequence for building (auto => k-mer size)");
  min_length_option->default_val(0);
  min_length_option->default_str("auto");
  CLI::Validator min_length_validator;
  min_length_validator.description("auto or non-negative integer");
  min_length_validator.operation([](std::string &input) -> std::string {
    std::string lowered = input;
    std::transform(lowered.begin(), lowered.end(), lowered.begin(), [](unsigned char ch) {
      return static_cast<char>(std::tolower(ch));
    });
    if (lowered == "auto") {
      input = "0";
      return {};
    }
    try {
      long long value = std::stoll(input);
      if (value < 0) {
        return "Minimum length must be >= 0 or 'auto'";
      }
    } catch (const std::exception &) {
      return "Minimum length must be an integer or 'auto'";
    }
    return {};
  });
  min_length_option->check(min_length_validator);
  build
      ->add_option("-t,--threads", buildConfig.threads,
                   "Number of threads for building")
      ->default_val(default_threads);
  build
      ->add_option("--load-factor", buildConfig.load_factor,
                   "IMCF 滤器负载因子")
      ->default_val(0.8);
  build
      ->add_option("-M,--max-hashes,--k-max", buildConfig.max_hashes_per_taxid,
                   "每 taxid 最多写入的 hash 数 (0=>auto=2e6)")
      ->default_val(0);
  build
      ->add_option("--k-base", buildConfig.k_base,
                   "每 taxid 的基准 hash 预算 (5Mb 时)")
      ->default_val(400000);
  build
      ->add_option("--k-min", buildConfig.k_min,
                   "每 taxid 最少写入的 hash 数")
      ->default_val(100000);
  build
      ->add_option("--core-alpha", buildConfig.core_alpha,
                   "Core-IDF 中的 core 平滑项 alpha")
      ->default_val(1.0);
  build
      ->add_option("--core-beta", buildConfig.core_beta,
                   "Core-IDF 中的 core 指数 beta")
      ->default_val(2.0);
  build
      ->add_option("--taxid-file-cap", buildConfig.taxid_file_cap,
                   "每个 taxid 参与 Core-IDF 的 genome 文件上限 (0=>不限制)")
      ->default_val(256);
  build
      ->add_option("--sig-oversample", buildConfig.sig_oversample,
                   "候选池倍率: 目标 K 的 oversample 倍")
      ->default_val(6.0);
  build
      ->add_option("--sig-s-min", buildConfig.sig_s_min,
                   "每文件签名最小大小")
      ->default_val(512);
  build
      ->add_option("--sig-s-max", buildConfig.sig_s_max,
                   "每文件签名最大大小")
      ->default_val(131072);
  build
      ->add_option("--presence-unique-deg", buildConfig.presence_unique_deg,
                   "Degree cutoff (<=) treated as unique signature for coverage meta")
      ->default_val(1)
      ->check(CLI::Range(1, 65535));
  build
      ->add_option("--feature", buildConfig.feature,
                   "Feature 提取方式 (syncmer|strobemer|auto)")
      ->default_val("strobemer");
  build
      ->add_option("--strobe-k", buildConfig.strobemer_k,
                   "Strobemer k-mer 长度")
      ->default_val(28);
  build
      ->add_option("--strobe-order", buildConfig.strobemer_order,
                   "Strobemer 阶数 (目前仅支持 2)")
      ->default_val(2);
  build
      ->add_option("--strobe-w-min", buildConfig.strobemer_w_min,
                   "Strobemer 最小窗口")
      ->default_val(12);
  build
      ->add_option("--strobe-w-max", buildConfig.strobemer_w_max,
                   "Strobemer 最大窗口")
      ->default_val(32);
  build
      ->add_option("--taxonomy-kind", buildConfig.taxonomy_kind,
                   "taxonomy 数据源标识 (auto|ncbi|gtdb)")
      ->default_val("auto");
  build
      ->add_option("--taxonomy-version", buildConfig.taxonomy_version,
                   "taxonomy 数据版本标识，例如 ncbi-taxdump-2025-09-15 或 gtdb-rs226")
      ->default_val("auto");
  build->add_flag("-q,--quiet", buildQuietRequested, "Quiet output");

  build->callback([&buildConfig, min_length_option, &buildQuietRequested]() {
    if (buildConfig.smer_size == 0) {
      throw CLI::ValidationError("--syncmer-s must be greater than 0");
    }
    if (buildConfig.smer_size >= buildConfig.kmer_size) {
      throw CLI::ValidationError("--syncmer-s must be smaller than k-mer size");
    }
    const uint16_t window_span = static_cast<uint16_t>(buildConfig.kmer_size - buildConfig.smer_size + 1);
    if (buildConfig.syncmer_position >= window_span) {
      throw CLI::ValidationError("--syncmer-pos must be < k - s + 1");
    }
    auto normalize_kind = [](std::string &value) {
      std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
      });
      if (value.empty()) {
        value = "auto";
      }
    };
    auto normalize_feature = [](std::string &value) {
      std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
      });
      if (value.empty()) {
        value = "auto";
      }
    };
    auto sanitize_version = [](std::string &value) {
      if (value.empty()) {
        value = "auto";
      }
    };
    normalize_kind(buildConfig.taxonomy_kind);
    normalize_feature(buildConfig.feature);
    sanitize_version(buildConfig.taxonomy_version);
    if (min_length_option->count() == 0) {
      buildConfig.min_length = buildConfig.kmer_size;
    }
    buildConfig.min_length = std::max<uint64_t>(buildConfig.min_length, buildConfig.kmer_size);
    if (!(buildConfig.feature == "auto" || buildConfig.feature == "syncmer" ||
          buildConfig.feature == "strobemer")) {
      throw CLI::ValidationError("--feature must be one of auto|syncmer|strobemer");
    }
    if (buildConfig.strobemer_w_min == 0 || buildConfig.strobemer_w_max == 0) {
      throw CLI::ValidationError("--strobe-w-min/w-max must be greater than 0");
    }
    if (buildConfig.strobemer_w_min > buildConfig.strobemer_w_max) {
      throw CLI::ValidationError("--strobe-w-min must be <= --strobe-w-max");
    }
    if (buildConfig.strobemer_k < 8) {
      throw CLI::ValidationError("--strobe-k must be >= 8");
    }
    if (buildConfig.strobemer_order != 2) {
      throw CLI::ValidationError("--strobe-order 当前仅支持取 2");
    }
    if (buildConfig.sig_s_min == 0) {
      throw CLI::ValidationError("--sig-s-min must be > 0");
    }
    if (buildConfig.sig_s_max < buildConfig.sig_s_min) {
      throw CLI::ValidationError("--sig-s-max must be >= --sig-s-min");
    }
    if (buildConfig.sig_oversample < 1.0) {
      throw CLI::ValidationError("--sig-oversample must be >= 1");
    }
    if (buildConfig.core_alpha < 0.0 || buildConfig.core_beta < 0.0) {
      throw CLI::ValidationError("--core-alpha/--core-beta must be >= 0");
    }
    buildConfig.verbose = !buildQuietRequested;
  });

  // Classify
  // Add --single option
  auto singleOpt = classify
                       ->add_option("-i,--single", classifyConfig.singleFiles,
                                    "Input file for classifying")
                       ->check(CLI::ExistingFile);

  // Add --paired option
  auto pairedOpt =
      classify
          ->add_option("-p,--paired", classifyConfig.pairedFiles,
                       "Paired input files for classifying")
          ->check(CLI::ExistingFile)
          ->excludes(
              singleOpt) // Ensure that single and paired are mutually exclusive
          ->each([](const std::string &) {
          }); // Use each to allow multiple inputs for paired

  classify->add_option("--weight-map", classifyConfig.weight_map_file,
                       "Optional contig/read weight map: id<TAB>weight or CAMI mapping.tsv")
      ->check(CLI::ExistingFile);

  // Custom validation function to ensure that the --paired option must have an
  // even number of files
  classify->callback([pairedOpt, &classifyConfig, &classifyQuietRequested]() {
    if (pairedOpt->count() > 0 && pairedOpt->count() % 2 != 0) {
      throw CLI::ValidationError(
          "--paired option must have an even number of input files");
    }
    auto normalize_kind = [](std::string &value) {
      std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
      });
      if (value.empty()) {
        value = "auto";
      }
    };
    auto normalize_feature = [](std::string &value) {
      std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
      });
      if (value.empty()) {
        value = "auto";
      }
    };
    auto normalize_presence = [](std::string &value, const char *fallback) {
      value = fallback; // presence 模式固定为 coverage，参数已移除
    };
    auto sanitize_version = [](std::string &value) {
      if (value.empty()) {
        value = "auto";
      }
    };
    normalize_kind(classifyConfig.taxonomyKind);
    sanitize_version(classifyConfig.taxonomyVersion);
    normalize_feature(classifyConfig.feature);
    normalize_presence(classifyConfig.decoy_mode, "imcf-edge-shuffle");
    if (!(classifyConfig.feature == "auto" || classifyConfig.feature == "syncmer" ||
          classifyConfig.feature == "strobemer")) {
      throw CLI::ValidationError("--feature must be one of auto|syncmer|strobemer");
    }
    if (classifyConfig.strobemer_w_min != 0 && classifyConfig.strobemer_w_max != 0 &&
        classifyConfig.strobemer_w_min > classifyConfig.strobemer_w_max) {
      throw CLI::ValidationError("--strobe-w-min must be <= --strobe-w-max");
    }
    if (classifyConfig.strobemer_order != 0 &&
        classifyConfig.strobemer_order != 2) {
      throw CLI::ValidationError("--strobe-order 当前仅支持取 2");
    }
    if (classifyConfig.strobemer_k != 0 && classifyConfig.strobemer_k < 8) {
      throw CLI::ValidationError("--strobe-k must be >= 8 when specified");
    }
    if (classifyConfig.decoy_mode != "imcf-edge-shuffle") {
      throw CLI::ValidationError("--decoy-mode 当前仅支持 imcf-edge-shuffle");
    }
    classifyConfig.verbose = !classifyQuietRequested;
  });

  classify
      ->add_option("-o,--output", classifyConfig.outputFile,
                   "Output file for classifying")
      ->default_val("ChimeraClassify");
  classify
      ->add_option("-d,--database", classifyConfig.dbFile,
                   "Database file for classifying")
      ->required()
      ->check(CLI::ExistingFile);
  classify
      ->add_option("--feature", classifyConfig.feature,
                   "Feature 提取方式 (syncmer|strobemer|auto；auto 表示跟随数据库)")
      ->default_val("auto");
  auto *strobeKOpt =
      classify
          ->add_option("--strobe-k", classifyConfig.strobemer_k,
                       "Strobemer k-mer 长度 (默认跟随数据库)")
          ->default_val(0);
  strobeKOpt->default_str("inherit");
  auto *strobeOrderOpt =
      classify
      ->add_option("--strobe-order", classifyConfig.strobemer_order,
                       "Strobemer 阶数 (默认跟随数据库，当前仅支持 2)")
          ->default_val(0);
  strobeOrderOpt->default_str("inherit");
  auto *strobeWminOpt =
      classify
          ->add_option("--strobe-w-min", classifyConfig.strobemer_w_min,
                       "Strobemer 最小窗口 (默认跟随数据库)")
          ->default_val(0);
  strobeWminOpt->default_str("inherit");
  auto *strobeWmaxOpt =
      classify
          ->add_option("--strobe-w-max", classifyConfig.strobemer_w_max,
                       "Strobemer 最大窗口 (默认跟随数据库)")
          ->default_val(0);
  strobeWmaxOpt->default_str("inherit");
  classify
      ->add_option("-s,--shot-threshold", classifyConfig.shotThreshold,
                   "Shot threshold for classifying")
      ->default_val(0.70);
  classify
      ->add_flag("--adaptive-shot,!--no-adaptive-shot",
                 classifyConfig.adaptive_shot,
                 "Scale thresholds by actually evaluated syncmers")
      ->default_val(true);
  auto *firstBetaOpt =
      classify
          ->add_option("--first-filter-beta", classifyConfig.firstFilterBeta,
                       "Keep bins >= beta * max count in first filter")
          ->default_val(0.80);
  classify
      ->add_option("--pre-em-topk", classifyConfig.preEmTopK,
                   "Keep top-K candidates per read before EM/VEM")
      ->default_val(16);
  classify->add_flag(
      "--collapse-strain-hits", classifyConfig.collapse_strain_hits,
      "[NCBI-only] Collapse per-hash hit lists to 1 representative taxid per species (reduces strain saturation; changes scoring).");
  classify->add_flag(
      "--collapse-strain-candidates", classifyConfig.collapse_strain_candidates,
      "[NCBI-only] Collapse pre-EM candidates to 1 representative per species before topK truncation (reduces strain saturation; lighter than --collapse-strain-hits).");
  classify->add_flag(
      "--deg-by-species", classifyConfig.deg_by_species,
      "[NCBI-only] Compute deg/exclusivity by species groups (mitigates strain saturation) while keeping output taxids unchanged.");
  classify
      ->add_option("--presence-pi", classifyConfig.presence_pi,
                   "Presence prior P(z=1) for coverage模型 (0-1)")
      ->check(CLI::Range(1e-6, 0.5))
      ->default_val(1e-3);
  classify
      ->add_option("--presence-tau", classifyConfig.presence_tau,
                   "Log posterior odds threshold for presence (coverage)")
      ->default_val(4.6);
  classify
      ->add_option("--presence-noise", classifyConfig.presence_noise,
                   "Override noise μ (per-signature) for coverage; <=0 auto")
      ->default_val(0.0);
  classify
      ->add_option("--presence-u-min", classifyConfig.presence_u_min,
                   "Minimum U_j used in coverage model to防止极小基因组过拟合 (1-16)")
      ->default_val(1)
      ->check(CLI::Range(1u, 16u));
  classify
      ->add_option("--presence-breadth-bits", classifyConfig.presence_breadth_bits,
                   "Sketch bits for presence breadth (power-of-two suggested)")
      ->check(CLI::Range(64u, 1048576u))
      ->default_val(2048);
  classify
      ->add_option("--presence-breadth-min-ratio",
                   classifyConfig.presence_breadth_min_ratio,
                   "Minimum breadth ratio (0 disables)")
      ->check(CLI::Range(0.0, 1.0))
      ->default_val(0.0);
  classify
      ->add_option("--presence-breadth-min-obs",
                   classifyConfig.presence_breadth_min_obs,
                   "Minimum observed unique signatures (0 disables)")
      ->check(CLI::Range(0u, 100000000u))
      ->default_val(0);
  classify
      ->add_option("--presence-breadth-penalty",
                   classifyConfig.presence_breadth_penalty,
                   "Penalty subtracted from logPosterior when breadth is low (0 disables)")
      ->check(CLI::Range(0.0, 1e6))
      ->default_val(0.0);
  classify
      ->add_option("--decoy-mode", classifyConfig.decoy_mode,
                   "Decoy generation mode (imcf-edge-shuffle)")
      ->default_val("imcf-edge-shuffle");
  classify
      ->add_option("--decoy-reps", classifyConfig.decoy_reps,
                   "Number of IMCF edge-shuffle decoy replicates (0-5)")
      ->check(CLI::Range(0, 5))
      ->default_val(3);
  classify
      ->add_option("--exclusive-gamma", classifyConfig.exclusive_gamma,
                   "Exclusive edge weighting gamma (0.0-2.0)")
      ->check(CLI::Range(0.0, 2.0))
      ->default_val(1.2);
  classify
      ->add_option("-t,--threads", classifyConfig.threads,
                   "Number of threads for classifying")
      ->default_val(default_threads);
  classify
      ->add_option("-b,--batch-size", classifyConfig.batchSize,
                   "Batch size for classifying")
      ->default_val(400);
  bool classifyNoEm = false;
  auto emFlag = classify->add_flag("-e,--EM", classifyConfig.em,
                                   "Enable EM mode (default; use --no-em to disable)");
  auto noEmFlag =
      classify->add_flag("--no-em", classifyNoEm, "Disable EM mode");
  noEmFlag->excludes(emFlag);
  emFlag->excludes(noEmFlag);
  classify
      ->add_option("--em-threshold", classifyConfig.emThreshold, "EM threshold")
      ->default_val(1e-4);
  classify->add_option("--em-iter", classifyConfig.emIter, "EM iteration")
      ->default_val(100);
  classify
      ->add_option("--em-prune-ratio", classifyConfig.em_prune_ratio,
                   "Relative sparsity prune ratio to max expected count in EM")
      ->check(CLI::Range(1e-6, 1.0))
      ->default_val(2e-4);
  classify
      ->add_option("--em-prior-strength", classifyConfig.em_prior_strength,
                   "Dirichlet prior strength (pseudo-count mass); 0 uses alpha only")
      ->check(CLI::Range(0.0, 10.0))
      ->default_val(1.0);
  classify
      ->add_option("--em-coexist-penalty", classifyConfig.em_coexist_penalty,
                   "Penalty when multiple taxa tie in EM (0 disables)")
      ->check(CLI::Range(0.0, 1.0))
      ->default_val(0.6);
  classify
      ->add_option("--em-conf-power", classifyConfig.em_conf_power,
                   "Confidence exponent for EM M-step weighting (0 disables)")
      ->check(CLI::Range(0.0, 3.0))
      ->default_val(2.0);
  classify
      ->add_option("--hash-sample-max", classifyConfig.hash_sample_max,
                   "Max sampled hashes for coarse candidate selection")
      ->check(CLI::Range(32.0, 4096.0))
      ->default_val(96);
  classify
      ->add_option("--hash-sample-min", classifyConfig.hash_sample_min,
                   "Min sampled hashes for coarse candidate selection")
      ->check(CLI::Range(8.0, 1024.0))
      ->default_val(16);
  classify
      ->add_option("--idf-max", classifyConfig.idf_max,
                   "Upper cap for IDF weight in evaluate_minimizer")
      ->check(CLI::Range(1.0, 50.0))
      ->default_val(8.0);
  classify
      ->add_option("--post-thres", classifyConfig.post_thres,
                   "Posterior acceptance threshold")
      ->default_val(0.56);
  classify
      ->add_option("--post-pi-min", classifyConfig.post_pi_min,
                   "Minimum global class weight")
      ->default_val(5e-4);
  classify
      ->add_option("--post-min-fraction", classifyConfig.post_min_fraction,
                   "postEmDecision: ignore taxa with posterior < fraction (0 disables)")
      ->check(CLI::Range(0.0, 1.0))
      ->default_val(0.01);
  classify
      ->add_option("--post-power", classifyConfig.post_power,
                   "postEmDecision: posterior^alpha sharpening (>=1)")
      ->check(CLI::Range(0.0, 10.0))
      ->default_val(1.5);
  auto *postHeadMassOpt =
      classify
          ->add_option("--post-head-mass", classifyConfig.post_head_mass,
                       "postEmDecision: keep top head_mass of sum(posterior^alpha) (0=>auto)")
          ->check(CLI::Range(0.0, 1.0))
          ->default_val(0.0);
  postHeadMassOpt->default_str("auto");
  auto *postMaxTaxaOpt =
      classify
          ->add_option("--post-max-taxa", classifyConfig.post_max_taxa,
                       "postEmDecision: max taxa per read in taxidCount (0=>auto)")
          ->check(CLI::Range(0u, 512u))
          ->default_val(0);
  postMaxTaxaOpt->default_str("auto");
  classify
      ->add_option("--dump-post-topk", classifyConfig.dump_post_topk,
                   "Dump top-K posterior candidates as POST_TOPK=... token in output TSV (0 disables)")
      ->check(CLI::Range(0u, 512u))
      ->default_val(256);
  classify->add_flag(
      "--presence-pre-filter", classifyConfig.presence_pre_filter,
      "Apply presence filter before EM/VEM (hard-prune candidates; may increase unclassified)");
  // TODO: 后处理相关参数暂时废弃，内部逻辑维持默认行为
  classify->add_flag("-q,--quiet", classifyQuietRequested, "Quiet output");

  if (argc == 1) {
    std::cout << app.help() << std::endl;
    return 0;
  }

  CLI11_PARSE(app, argc, argv);

  if (classifyNoEm) {
    classifyConfig.em = false;
  }
  if (firstBetaOpt && firstBetaOpt->count() > 0 &&
      classifyConfig.firstFilterBeta > 0.0) {
    classifyConfig.firstFilterBeta_user = true;
  }

  if (show_version) {
    std::cout << "======================================" << std::endl;
    std::cout << "        Chimera - Metagenomic Tool" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Version      : " << VERSION_INFO << std::endl;
    std::cout << "Build Date   : " << __DATE__ << " " << __TIME__ << std::endl;
#ifdef __clang__
    std::cout << "Compiled with: Clang " << __clang_major__ << '.'
              << __clang_minor__ << '.' << __clang_patchlevel__ << std::endl;
#elif defined(__GNUC__)
    std::cout << "Compiled with: GCC " << __GNUC__ << '.' << __GNUC_MINOR__ << '.'
              << __GNUC_PATCHLEVEL__ << std::endl;
#elif defined(__VERSION__)
    std::cout << "Compiled with: " << __VERSION__ << std::endl;
#else
    std::cout << "Compiled with: Unknown compiler" << std::endl;
#endif
    std::cout << "OS           : " << detect_os() << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Developed by : Qinzhong Tian" << std::endl;
    std::cout << "Team         : MalabZ" << std::endl;
    std::cout << "Homepage     : https://github.com/LoadStar822/Chimera"
              << std::endl;
    std::cout << "======================================" << std::endl;
    return 0;
  }

  if (*build) {
    ChimeraBuild::run(buildConfig);
  } else if (*classify) {
    ChimeraClassify::run(classifyConfig);
  }

  return 0;
}
