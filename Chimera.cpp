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
  build->add_flag("--disable-strobemer",
                  [&buildConfig](size_t count) {
                    if (count > 0) {
                      buildConfig.enable_strobemers = false;
                    }
                  },
                  "禁用 strobemer 双特征流");
  build
      ->add_option("--strobe-wmin", buildConfig.strobemer_w_min,
                   "Strobemer 窗口下界 (0 表示自动)")
      ->default_val(0);
  build
      ->add_option("--strobe-wmax", buildConfig.strobemer_w_max,
                   "Strobemer 窗口上界 (0 表示自动)")
      ->default_val(0);
  build
      ->add_option("--strobe-q", buildConfig.strobemer_q,
                   "Strobemer popcount q 掩码 (0 表示自动)")
      ->default_val(0);
  build
      ->add_option("--strobe-maxdist", buildConfig.strobemer_max_dist,
                   "Strobemer 最大跨距 (0 表示自动)")
      ->default_val(0);
  build
      ->add_option("--strobe-aux-len", buildConfig.strobemer_aux_len,
                   "Strobemer aux_len (默认 15，用于主哈希掩码)")
      ->default_val(15);
  build
      ->add_option("--strobe-weight", buildConfig.strobemer_weight,
                   "Strobemer 事件权重 (0 表示自动)")
      ->default_val(0.0);
  build
      ->add_option("--strobe-ratio", buildConfig.strobemer_ratio,
                   "Strobemer 采样配额占比 [0,1]，默认 0.6")
      ->default_val(0.6)
      ->check(CLI::Range(0.0, 1.0));
  build
      ->add_option("--strobe-seed", buildConfig.strobemer_seed,
                   "Strobemer 额外哈希种子 (0 表示自动派生)")
      ->default_val(0ULL);
  build
      ->add_option("-l,--min-length", buildConfig.min_length,
                   "Minimum length sequence for building")
      ->default_val(0);
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
                   "每个 taxid bottom-k 采样数量")
      ->default_val(2000000);
  build
      ->add_option("--toxic-quantile", buildConfig.toxic_quantile,
                   "过滤高频 syncmer 的分位阈值 (0-1]")
      ->default_val(0.999)
      ->check(CLI::Range(0.0, 1.0));
  build
      ->add_option("--toxic-topN", buildConfig.toxic_top_n,
                   "优先剔除全库计数前 N 的 hash")
      ->default_val(0);
  build
      ->add_option("--toxic-min-fraction", buildConfig.toxic_min_fraction,
                   "毒性候选的最低全库占比")
      ->default_val(0.0001)
      ->check(CLI::Range(0.0, 1.0));
  build
      ->add_option("--toxic-safety-frac", buildConfig.toxic_safety_fraction,
                   "每个 taxid 最少保留比例")
      ->default_val(0.1)
      ->check(CLI::Range(0.0, 1.0));
  build
      ->add_option("--toxic-safety-min", buildConfig.toxic_safety_min,
                   "每个 taxid 最少保留数量")
      ->default_val(1024);
  build->add_flag("--adaptive-cutoff", buildConfig.adaptive_cutoff,
                  "启用基于文件规模的自适应 cutoff");
  build
      ->add_option("-f,--filter", buildConfig.filter,
                   "构建使用的滤器类型 (imcf)")
      ->check(CLI::IsMember({"imcf"}))
      ->default_val("imcf");
  build
      ->add_option("--taxonomy-kind", buildConfig.taxonomy_kind,
                   "taxonomy 数据源 (ncbi/gtdb 等)")
      ->default_val("auto");
  build
      ->add_option("--taxonomy-version", buildConfig.taxonomy_version,
                   "taxonomy 数据版本标识，例如 ncbi-taxdump-2025-09-15 或 gtdb-rs226")
      ->default_val("auto");
  build->add_flag("-q,--quiet", buildConfig.verbose, "Quiet output")
      ->default_val(true)
      ->disable_flag_override();

  build->callback([&buildConfig]() {
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
    if (buildConfig.toxic_top_n == 0 &&
        (buildConfig.toxic_quantile <= 0.0 || buildConfig.toxic_quantile > 1.0)) {
      throw CLI::ValidationError("当未指定 --toxic-topN 时，--toxic-quantile 必须落在 (0, 1] 区间");
    }
    if (buildConfig.toxic_min_fraction < 0.0 || buildConfig.toxic_min_fraction > 1.0) {
      throw CLI::ValidationError("--toxic-min-fraction must be between 0 and 1");
    }
    if (buildConfig.toxic_safety_fraction < 0.0 || buildConfig.toxic_safety_fraction > 1.0) {
      throw CLI::ValidationError("--toxic-safety-frac must be between 0 and 1");
    }
    if (buildConfig.toxic_top_n > 0 && buildConfig.toxic_quantile == 0.0) {
      buildConfig.toxic_quantile = 1.0;
    }
    if (buildConfig.strobemer_ratio < 0.0 || buildConfig.strobemer_ratio > 1.0) {
      throw CLI::ValidationError("--strobe-ratio must be between 0 and 1");
    }
    const bool manual_strobe =
        buildConfig.strobemer_w_min > 0 || buildConfig.strobemer_w_max > 0 ||
        buildConfig.strobemer_q > 0 || buildConfig.strobemer_max_dist > 0 ||
        buildConfig.strobemer_weight > 0.0 || buildConfig.strobemer_seed != 0;
    buildConfig.strobemer_auto = !manual_strobe;
    if (buildConfig.enable_strobemers && !buildConfig.strobemer_auto) {
      if (buildConfig.strobemer_w_min == 0 || buildConfig.strobemer_w_max == 0) {
        throw CLI::ValidationError("手动指定 strobemer 参数时必须提供 --strobe-wmin 与 --strobe-wmax");
      }
      if (buildConfig.strobemer_w_min > buildConfig.strobemer_w_max) {
        throw CLI::ValidationError("--strobe-wmin must be <= --strobe-wmax");
      }
      if (buildConfig.strobemer_q == 0) {
        throw CLI::ValidationError("手动模式下 --strobe-q 不能为 0");
      }
      if (buildConfig.strobemer_max_dist == 0) {
        throw CLI::ValidationError("手动模式下 --strobe-maxdist 不能为 0");
      }
      if (buildConfig.strobemer_weight <= 0.0) {
        throw CLI::ValidationError("手动模式下 --strobe-weight 必须大于 0");
      }
    }
    if (!buildConfig.enable_strobemers) {
      buildConfig.strobemer_ratio = 0.0;
    } else if (buildConfig.strobemer_ratio == 0.0 && buildConfig.strobemer_auto) {
      buildConfig.strobemer_ratio = 0.6;
    }
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
  classify->callback([pairedOpt, &classifyConfig]() {
    if (pairedOpt->count() > 0 && pairedOpt->count() % 2 != 0) {
      throw CLI::ValidationError(
          "--paired option must have an even number of input files");
    }
    const bool manual_strobe =
        classifyConfig.strobemer_w_min > 0 || classifyConfig.strobemer_w_max > 0 ||
        classifyConfig.strobemer_q > 0 || classifyConfig.strobemer_max_dist > 0 ||
        classifyConfig.strobemer_weight > 0.0 || classifyConfig.strobemer_seed != 0;
    classifyConfig.strobemer_auto = !manual_strobe;
    if (classifyConfig.enable_strobemers && !classifyConfig.strobemer_auto) {
      if (classifyConfig.strobemer_w_min == 0 || classifyConfig.strobemer_w_max == 0) {
        throw CLI::ValidationError("手动指定 strobemer 参数时必须提供 --strobe-wmin 与 --strobe-wmax");
      }
      if (classifyConfig.strobemer_w_min > classifyConfig.strobemer_w_max) {
        throw CLI::ValidationError("--strobe-wmin must be <= --strobe-wmax");
      }
      if (classifyConfig.strobemer_q == 0) {
        throw CLI::ValidationError("手动模式下 --strobe-q 不能为 0");
      }
      if (classifyConfig.strobemer_max_dist == 0) {
        throw CLI::ValidationError("手动模式下 --strobe-maxdist 不能为 0");
      }
      if (classifyConfig.strobemer_weight <= 0.0) {
        throw CLI::ValidationError("手动模式下 --strobe-weight 必须大于 0");
      }
    }
    if (!classifyConfig.enable_strobemers) {
      classifyConfig.strobemer_weight = 0.0;
    }
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
      ->add_option("-f,--filter", classifyConfig.filter,
                   "分类使用的滤器类型 (imcf)")
      ->default_val("imcf")
      ->check(CLI::IsMember({"imcf"}));
  classify
      ->add_option("--taxonomy-kind", classifyConfig.taxonomyKind,
                   "期望的 taxonomy 数据源 (ncbi/gtdb 等)")
      ->default_val("auto");
  classify
      ->add_option("--taxonomy-version", classifyConfig.taxonomyVersion,
                   "期望的 taxonomy 数据版本号，例如 ncbi-taxdump-2025-09-15")
      ->default_val("auto");
  classify->add_flag("--disable-strobemer",
                     [&classifyConfig](size_t count) {
                       if (count > 0) {
                         classifyConfig.enable_strobemers = false;
                       }
                     },
                     "禁用 strobemer 特征");
  classify
      ->add_option("--strobe-wmin", classifyConfig.strobemer_w_min,
                   "Strobemer 窗口下界 (0 表示按读长自适应)")
      ->default_val(0);
  classify
      ->add_option("--strobe-wmax", classifyConfig.strobemer_w_max,
                   "Strobemer 窗口上界 (0 表示按读长自适应)")
      ->default_val(0);
  classify
      ->add_option("--strobe-q", classifyConfig.strobemer_q,
                   "Strobemer popcount 掩码 (0 表示按读长自适应)")
      ->default_val(0);
  classify
      ->add_option("--strobe-maxdist", classifyConfig.strobemer_max_dist,
                   "Strobemer 最大跨距 (0 表示按读长自适应)")
      ->default_val(0);
  classify
      ->add_option("--strobe-aux-len", classifyConfig.strobemer_aux_len,
                   "Strobemer aux_len (0 表示使用数据库配置)")
      ->default_val(0);
  classify
      ->add_option("--strobe-weight", classifyConfig.strobemer_weight,
                   "Strobemer 事件权重 (0 表示按读长自适应)")
      ->default_val(0.0);
  classify
      ->add_option("--strobe-seed", classifyConfig.strobemer_seed,
                   "Strobemer 混合哈希种子 (0 表示使用数据库里的配置)")
      ->default_val(0ULL);
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
      ->add_option("--post-pi-min", classifyConfig.post_pi_min,
                   "Minimum global class weight")
      ->default_val(1e-4);
  classify
      ->add_flag("--evidence-override,!--no-evidence-override",
                 classifyConfig.evidence_override,
                 "Allow strong pre-EM evidence to bypass posterior filter")
      ->default_val(true);
  classify
      ->add_flag("--output-posterior,!--no-output-posterior",
                 classifyConfig.output_posterior,
                 "Write posterior probabilities to the TSV output")
      ->default_val("true");
  // TODO: 后处理相关参数暂时废弃，内部逻辑维持默认行为
  classify->add_flag("-q,--quiet", classifyConfig.verbose, "Quiet output")
      ->default_val(true)
      ->disable_flag_override();

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
