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
      ->default_val(19)
      ->check(CLI::Range(1, 50));
  build
      ->add_option("-w,--window", buildConfig.window_size,
                   "Window size for building")
      ->default_val(31);
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
      ->default_val(0.9);
  build
      ->add_option("-M,--max-hashes", buildConfig.max_hashes_per_taxid,
                   "Maximum number of hashes per taxid")
      ->default_val(2000000);
  build->add_option("-c,--fixed-cutoff", buildConfig.fixed_cutoff,
                    "Fixed cutoff for building (0 - 255)");
  build
      ->add_option("-f,--filter", buildConfig.filter,
                   "构建使用的滤器类型 (imcf)")
      ->check(CLI::IsMember({"imcf"}))
      ->default_val("imcf");
  build->add_flag("-q,--quiet", buildConfig.verbose, "Quiet output")
      ->default_val(true)
      ->disable_flag_override();

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
  classify->callback([pairedOpt]() {
    if (pairedOpt->count() > 0 && pairedOpt->count() % 2 != 0) {
      throw CLI::ValidationError(
          "--paired option must have an even number of input files");
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
      ->default_val(0.65);
  classify
      ->add_flag("--adaptive-shot,!--no-adaptive-shot",
                 classifyConfig.adaptive_shot,
                 "Scale thresholds by actually evaluated minimizers")
      ->default_val(true);
  classify
      ->add_option("--first-filter-beta", classifyConfig.firstFilterBeta,
                   "Keep bins >= beta * max count in first filter")
      ->default_val(0.75);
  classify
      ->add_option("--pre-em-topk", classifyConfig.preEmTopK,
                   "Keep top-K candidates per read before EM/VEM")
      ->default_val(24);
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
                   "Minimum #evaluated minimizers to classify")
      ->default_val(24);
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
      ->add_option("-b,--batch-size", classifyConfig.batchSize,
                   "Batch size for classifying")
      ->default_val(400);
  classify->add_flag("--lca", classifyConfig.lca, "Enable LCA mode");
  classify
      ->add_option("--tax-file", classifyConfig.taxFile,
                   "Taxonomy file for LCA mode")
      ->check(CLI::ExistingFile);
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
      ->default_val(0.9);
  classify
      ->add_option("--post-margin", classifyConfig.post_margin,
                   "Minimum gap between top posteriors")
      ->default_val(0.2);
  classify->add_option("--post-ratio", classifyConfig.post_ratio,
                       "Minimum ratio between top1 and top2 posteriors");
  classify
      ->add_option("--post-pi-min", classifyConfig.post_pi_min,
                   "Minimum global class weight")
      ->default_val(1e-4);
  classify->add_flag("--lca-fallback", classifyConfig.lca_fallback,
                     "Use LCA fallback when EM/VEM confidence is low");
  classify
      ->add_flag("--output-posterior,!--no-output-posterior",
                 classifyConfig.output_posterior,
                 "Write posterior probabilities to the TSV output")
      ->default_val("true");
  classify
      ->add_flag("--skip-postfilter,!--no-skip-postfilter",
                 classifyConfig.skip_post_filter,
                 "Skip legacy filtering steps when EM/VEM is enabled")
      ->default_val("true");
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
