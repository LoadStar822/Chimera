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
  const auto default_threads = static_cast<uint16_t>(hardware_threads);

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
      ->default_val(0.85);
  build
      ->add_option("-M,--max-hashes", buildConfig.max_hashes_per_taxid,
                   "Maximum number of hashes per taxid")
      ->default_val(2000000);
  build->add_flag("--adaptive-cutoff", buildConfig.adaptive_cutoff,
                  "启用基于文件规模的自适应 cutoff");
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
      std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
      });
      if (value.empty()) {
        value = fallback;
      }
    };
    auto sanitize_version = [](std::string &value) {
      if (value.empty()) {
        value = "auto";
      }
    };
    normalize_kind(classifyConfig.taxonomyKind);
    sanitize_version(classifyConfig.taxonomyVersion);
    normalize_feature(classifyConfig.feature);
    normalize_presence(classifyConfig.presence_caller, "tdfdr");
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
    if (!(classifyConfig.presence_caller == "tdfdr" ||
          classifyConfig.presence_caller == "hard_cutoff")) {
      throw CLI::ValidationError("--presence-caller must be tdfdr or hard_cutoff");
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
      ->default_val(0.62);
  classify
      ->add_flag("--adaptive-shot,!--no-adaptive-shot",
                 classifyConfig.adaptive_shot,
                 "Scale thresholds by actually evaluated syncmers")
      ->default_val(true);
  classify
      ->add_option("--first-filter-beta", classifyConfig.firstFilterBeta,
                   "Keep bins >= beta * max count in first filter")
      ->default_val(0.80);
  classify
      ->add_option("--pre-em-topk", classifyConfig.preEmTopK,
                   "Keep top-K candidates per read before EM/VEM")
      ->default_val(32);
  classify
      ->add_option("--presence-caller", classifyConfig.presence_caller,
                   "Presence caller strategy (tdfdr|hard_cutoff)")
      ->check(CLI::IsMember({"tdfdr", "hard_cutoff"}))
      ->default_val("tdfdr");
  classify
      ->add_option("--presence-q", classifyConfig.presence_q,
                   "Target q-value cutoff for presence calling")
      ->check(CLI::Range(0.001, 0.2))
      ->default_val(0.05);
  classify
      ->add_flag("--presence-report-only,!--presence-hard-filter",
                 classifyConfig.presence_report_only,
                 "Presence caller gating mode (default hard filter; use --presence-report-only to skip trimming)")
      ->default_val(false);
  classify
      ->add_option("--presence-abundance-prior",
                   classifyConfig.presence_abundance_prior,
                   "Down-weight p-values based on abundance (0-1)")
      ->check(CLI::Range(0.0, 1.0))
      ->default_val(0.30);
  classify
      ->add_flag("--auto-q-tune,!--no-auto-q-tune",
                 classifyConfig.auto_q_tune,
                 "Auto-tune presence q from {0.05,0.02,0.01,0.005} to cap decoy positives")
      ->default_val(true);
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
                   "Exclusive edge weighting gamma (0.5-2.0)")
      ->check(CLI::Range(0.5, 2.0))
      ->default_val(1.2);
  classify
      ->add_option("--min-unique-evidence", classifyConfig.min_unique_evidence,
                   "Minimum IMCF unique-edge hits required before testing (1-10)")
      ->check(CLI::Range(1, 10))
      ->default_val(5);
  classify
      ->add_flag("--adaptive-fdr,!--no-adaptive-fdr",
                 classifyConfig.adaptive_fdr,
                 "Enable per-read background-calibrated threshold")
      ->default_val(true);
  classify
      ->add_option("--fdr-z", classifyConfig.fdr_z,
                   "Z for Poisson(mu)+Z*sqrt(mu) threshold")
      ->default_val(3.0);
  classify
      ->add_option("--min-eval-count", classifyConfig.min_eval_count,
                   "Minimum #evaluated syncmers to classify")
      ->default_val(0);
  classify
      ->add_option("-t,--threads", classifyConfig.threads,
                   "Number of threads for classifying")
      ->default_val(default_threads);
  classify
      ->add_option("-b,--batch-size", classifyConfig.batchSize,
                   "Batch size for classifying")
      ->default_val(400);
  auto emFlag = classify->add_flag("-e,--EM", classifyConfig.em,
                                   "Enable EM mode (default)");
  auto vemFlag =
      classify->add_flag("-V,--VEM", classifyConfig.vem, "Enable VEM mode")
          ->excludes(emFlag);
  classify
      ->add_option("--em-threshold", classifyConfig.emThreshold, "EM threshold")
      ->default_val(0.001);
  classify->add_option("--em-iter", classifyConfig.emIter, "EM iteration")
      ->default_val(80);
  classify
      ->add_option("--post-thres", classifyConfig.post_thres,
                   "Posterior acceptance threshold")
      ->default_val(0.56);
  classify
      ->add_option("--post-margin", classifyConfig.post_margin,
                   "Minimum gap between top posteriors")
      ->default_val(0.03);
  classify->add_option("--post-ratio", classifyConfig.post_ratio,
                       "Minimum ratio between top1 and top2 posteriors")
      ->default_val(1.30);
  classify
      ->add_option("--post-relax-abs", classifyConfig.post_relax_abs,
                   "Absolute posterior floor for relaxed gate")
      ->default_val(0.50);
  classify
      ->add_option("--post-relax-ratio", classifyConfig.post_relax_ratio,
                   "Ratio required when using relaxed gate")
      ->default_val(1.50);
  classify
      ->add_option("--post-relax-delta", classifyConfig.post_relax_delta,
                   "Gap required when using relaxed gate")
      ->default_val(0.08);
  classify
      ->add_option("--post-relax-delta-abs",
                   classifyConfig.post_relax_delta_abs,
                   "Absolute posterior floor for gap-based relaxed gate")
      ->default_val(0.52);
  classify
      ->add_option("--post-pi-min", classifyConfig.post_pi_min,
                   "Minimum global class weight")
      ->default_val(1e-4);
  classify
      ->add_flag("--evidence-override,!--no-evidence-override",
                 classifyConfig.evidence_override,
                 "Allow strong pre-EM evidence to bypass posterior filter")
      ->default_val(false);
  classify
      ->add_option("--em-temp", classifyConfig.em_temp,
                   "Softmax temperature inside EM/VEM")
      ->default_val(1.10);
  classify
      ->add_option("--em-prior-strength", classifyConfig.em_prior_strength,
                   "Dirichlet pseudo-count mass for abundance prior")
      ->default_val(0.25);
  classify
      ->add_option("--em-coexist-penalty",
                   classifyConfig.em_coexist_penalty,
                   "Penalty term when candidates are near-tied")
      ->default_val(0.20);
  classify
      ->add_flag("--output-posterior,!--no-output-posterior",
                 classifyConfig.output_posterior,
                 "Write posterior probabilities to the TSV output")
      ->default_val("true");
  // TODO: 后处理相关参数暂时废弃，内部逻辑维持默认行为
  classify->add_flag("-q,--quiet", classifyQuietRequested, "Quiet output");

  if (argc == 1) {
    std::cout << app.help() << std::endl;
    return 0;
  }

  CLI11_PARSE(app, argc, argv);

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
