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
 * Last Modified: 2024-10-05
 *
 * Description:
 *  This is the main entry for Chimera
 *
 * Version:
 *  1.3
 * -----------------------------------------------------------------------------
 */
#include <CLI11.hpp>
#include <buildConfig.hpp>
#include <iostream>
#include <ChimeraBuild.hpp>
#include <classifyConfig.hpp>
#include <ChimeraClassify.hpp>
#include <vector>
#include <string>
#include <stdexcept>

#ifdef CHIMERA_VERSION
#define VERSION_INFO CHIMERA_VERSION
#else
#define VERSION_INFO "unknown"
#endif

int main(int argc, char** argv)
{
	// Create the main application object
	CLI::App app{ "Chimera - A versatile tool for metagenomic classification" };
	ChimeraBuild::BuildConfig buildConfig;
	ChimeraClassify::ClassifyConfig classifyConfig;

	bool show_version = false;
	app.add_flag("-v,--version", show_version, "Show version information");
	// Create subcommands
	auto build = app.add_subcommand("build", "Build a sequence database");
	auto classify = app.add_subcommand("classify", "Classify sequences");

	// Build
	build->add_option("-i,--input", buildConfig.input_file, "Input file for building")
		->required()
		->check(CLI::ExistingFile);
	build->add_option("-o,--output", buildConfig.output_file, "Output file for building")
		->default_val("ChimeraDB");
	build->add_option("-m,--mode", buildConfig.mode, "Mode for building")
		->default_val("normal");
	build->add_option("-k,--kmer", buildConfig.kmer_size, "Kmer size for building")
		->default_val(19)
		->check(CLI::Range(1, 50));
	build->add_option("-w,--window", buildConfig.window_size, "Window size for building")
		->default_val(31);
	build->add_option("-l,--min-length", buildConfig.min_length, "Minimum length sequence for building")
		->default_val(0);
	build->add_option("-t,--threads", buildConfig.threads, "Number of threads for building")
		->default_val(32);
	build->add_option("--load-factor", buildConfig.load_factor, "Loading ratio of ICF")
		->default_val(0.95);
	build->add_option("-M,--max-hashes", buildConfig.maxHashesPerTaxid, "Maximum number of hashes per taxid")
		->default_val(2000000);
	build->add_flag("-q,--quiet", buildConfig.verbose, "Quiet output")->default_val(true)->disable_flag_override();

	// Classify
	// Add --single option
	auto singleOpt = classify->add_option("-i,--single", classifyConfig.singleFiles, "Input file for classifying")
		->check(CLI::ExistingFile);

	// Add --paired option
	auto pairedOpt = classify->add_option("-p,--paired", classifyConfig.pairedFiles, "Paired input files for classifying")
		->check(CLI::ExistingFile)
		->excludes(singleOpt)  // Ensure that single and paired are mutually exclusive
		->each([](const std::string&) {});  // Use each to allow multiple inputs for paired

	// Custom validation function to ensure that the --paired option must have an even number of files
	classify->callback([pairedOpt]() {
		if (pairedOpt->count() > 0 && pairedOpt->count() % 2 != 0) {
			throw CLI::ValidationError("--paired option must have an even number of input files");
		}
		});

	classify->add_option("-o,--output", classifyConfig.outputFile, "Output file for classifying")
		->default_val("ChimeraClassify");
	classify->add_option("-d,--database", classifyConfig.dbFile, "Database file for classifying")
		->required()
		->check(CLI::ExistingFile);
	classify->add_option("-s,--shot-threshold", classifyConfig.shotThreshold, "Shot threshold for classifying")
		->default_val(0.7);
	classify->add_option("-t,--threads", classifyConfig.threads, "Number of threads for classifying")
		->default_val(32);
	classify->add_option("-m,--mode", classifyConfig.mode, "Mode for classifying")
		->default_val("fast");
	classify->add_option("-b,--batch-size", classifyConfig.batchSize, "Batch size for classifying")
		->default_val(400);
	classify->add_flag("--lca", classifyConfig.lca, "Enable LCA mode");
	classify->add_option("--tax-file", classifyConfig.taxFile, "Taxonomy file for LCA mode")
		->check(CLI::ExistingFile);
	auto emFlag = classify->add_flag("-e,--EM", classifyConfig.em, "Enable EM mode");
	auto vemFlag = classify->add_flag("-V,--VEM", classifyConfig.vem, "Enable VEM mode")->excludes(emFlag);
	classify->add_option("--em-threshold", classifyConfig.emThreshold, "EM threshold")
		->default_val(0.001);
	classify->add_option("--em-iter", classifyConfig.emIter, "EM iteration")
		->default_val(100);
	classify->add_flag("-q,--quiet", classifyConfig.verbose, "Quiet output")->default_val(true)->disable_flag_override();

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
		std::cout << "Compiled with: " << "GCC " << __VERSION__ << std::endl;
		std::cout << "OS           : Ubuntu 20.04" << std::endl;
		std::cout << "======================================" << std::endl;
		std::cout << "Developed by : Qinzhong Tian" << std::endl;
		std::cout << "Team         : MalabZ" << std::endl;
		std::cout << "Homepage     : https://github.com/LoadStar822/Chimera" << std::endl;
		std::cout << "======================================" << std::endl;
		return 0;
	}

	if (*build) {
		ChimeraBuild::run(buildConfig);
	}
	else if (*classify) {
		ChimeraClassify::run(classifyConfig);
	}

	return 0;
}