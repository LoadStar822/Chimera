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

std::string normalize_auto_min_length_help() {
  return "Minimum length sequence for building (auto => strobemer minimum span)";
}

uint16_t default_cli_threads() {
  unsigned int hardware_threads = std::thread::hardware_concurrency();
  if (hardware_threads == 0) {
    hardware_threads = 1;
  }
  const auto max_threads =
      static_cast<unsigned int>(std::numeric_limits<uint16_t>::max());
  if (hardware_threads > max_threads) {
    hardware_threads = max_threads;
  }
  return static_cast<uint16_t>(std::min<unsigned int>(hardware_threads, 192u));
}

std::string validate_min_length_option(std::string &input) {
  std::string lowered = input;
  std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                 [](unsigned char ch) {
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
}

void normalize_auto_string(std::string &value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char ch) {
                   return static_cast<char>(std::tolower(ch));
                 });
  if (value.empty()) {
    value = "auto";
  }
}

void validate_build_config(ChimeraBuild::BuildConfig &buildConfig,
                           bool buildQuietRequested) {
  normalize_auto_string(buildConfig.taxonomy_kind);
  buildConfig.min_length =
      std::max<uint64_t>(buildConfig.min_length, buildConfig.strobemer_k);
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
    throw CLI::ValidationError("--strobe-order currently only supports value 2");
  }
  buildConfig.verbose = !buildQuietRequested;
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

  const uint16_t default_threads = default_cli_threads();

  bool buildQuietRequested = false;

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
      ->add_option("--strobe-k", buildConfig.strobemer_k,
                   "Strobemer k-mer length")
      ->default_val(28);
  build
      ->add_option("--strobe-order", buildConfig.strobemer_order,
                   "Strobemer order (currently only 2 is supported)")
      ->default_val(2);
  build
      ->add_option("--strobe-w-min", buildConfig.strobemer_w_min,
                   "Strobemer minimum window")
      ->default_val(12);
  build
      ->add_option("--strobe-w-max", buildConfig.strobemer_w_max,
                   "Strobemer maximum window")
      ->default_val(32);
  auto *min_length_option = build
      ->add_option("-l,--min-length", buildConfig.min_length,
                   normalize_auto_min_length_help());
  min_length_option->default_val(0);
  min_length_option->default_str("auto");
  CLI::Validator min_length_validator;
  min_length_validator.description("auto or non-negative integer");
  min_length_validator.operation(validate_min_length_option);
  min_length_option->check(min_length_validator);
  build
      ->add_option("-t,--threads", buildConfig.threads,
                   "Number of threads for building")
      ->default_val(default_threads);
  build
      ->add_option("--load-factor", buildConfig.load_factor,
                   "IMCF filter load factor")
      ->default_val(0.85);
  build
      ->add_option("--presence-unique-deg", buildConfig.presence_unique_deg,
                   "Degree cutoff (<=) treated as unique signature for coverage meta")
      ->default_val(1)
      ->check(CLI::Range(1, 65535));
  build
      ->add_option("--taxonomy-kind", buildConfig.taxonomy_kind,
                   "Taxonomy source identifier (auto|ncbi|gtdb)")
      ->default_val("auto");
  build
      ->add_option("--taxonomy-version", buildConfig.taxonomy_version,
                   "Taxonomy version label, for example ncbi-taxdump-2025-09-15 or gtdb-rs226")
      ->default_val("auto");
  build->add_flag("-q,--quiet", buildQuietRequested, "Quiet output");

  build->callback([&buildConfig, &buildQuietRequested]() {
    validate_build_config(buildConfig, buildQuietRequested);
  });

  // Classify
  // Add --single option
  auto singleOpt = classify
                       ->add_option("-i,--single", classifyConfig.singleFiles,
                                    "Input file for classifying")
                       ->check(CLI::ExistingFile);

  // Add --paired option
  classify
      ->add_option("-p,--paired", classifyConfig.pairedFiles,
                   "Paired input files for classifying")
      ->check(CLI::ExistingFile)
      ->excludes(singleOpt); // Ensure that single and paired are mutually exclusive

  // Custom validation function to ensure that the --paired option must have an
  // even number of files
  classify->callback([&classifyConfig]() {
    if (!classifyConfig.pairedFiles.empty() &&
        classifyConfig.pairedFiles.size() % 2 != 0) {
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
      ->default_val(0.70);
  auto *firstBetaOpt =
      classify
          ->add_option("--first-filter-beta", classifyConfig.firstFilterBeta,
                       "Keep bins >= beta * max count in first filter")
          ->default_val(0.80);
  classify
      ->add_option("--presence-pi", classifyConfig.presence_pi,
                   "Presence prior P(z=1) for the coverage model (0-1)")
      ->check(CLI::Range(1e-6, 0.5))
      ->default_val(1e-3);
  classify
      ->add_option("--presence-tau", classifyConfig.presence_tau,
                   "Log posterior odds threshold for presence (coverage)")
      ->default_val(4.6);
  classify
      ->add_option("--presence-breadth-bits", classifyConfig.presence_breadth_bits,
                   "Sketch bits for presence breadth (power-of-two suggested)")
      ->check(CLI::Range(64u, 1048576u))
      ->default_val(2048);
  classify
      ->add_option("-t,--threads", classifyConfig.threads,
                   "Number of threads for classifying")
      ->default_val(default_threads);
  classify
      ->add_option("-b,--batch-size", classifyConfig.batchSize,
                   "Batch size for classifying")
      ->default_val(400);
  classify->add_option("--em-iter", classifyConfig.emIter, "EM iteration")
      ->default_val(100);
  classify
      ->add_option("--em-prune-ratio", classifyConfig.em_prune_ratio,
                   "Relative sparsity prune ratio to max expected count in EM")
      ->check(CLI::Range(1e-6, 1.0))
      ->default_val(2e-4);
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
      ->add_option("--post-pi-min", classifyConfig.post_pi_min,
                   "Minimum global class weight")
      ->default_val(5e-4);
  classify->add_flag("--read-evidence", classifyConfig.write_read_evidence,
                     "Write ChimeraReadEvidence.tsv for read-resolved audits");
  // TODO: Deprecated post-processing knobs remain fixed to internal defaults.

  if (argc == 1) {
    std::cout << app.help() << std::endl;
    return 0;
  }

  CLI11_PARSE(app, argc, argv);

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
