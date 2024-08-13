﻿/*
 * -----------------------------------------------------------------------------
 * Filename:      Chimera.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-07-09
 *
 * Last Modified: 2024-07-09
 *
 * Description:
 *  This is a simple C++ program that outputs "Hello, World!".
 *
 * Version:
 *  1.0
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
	CLI::App app{ "Chimera - A tool for sequence analysis" };
	ChimeraBuild::BuildConfig buildConfig;
	ChimeraClassify::ClassifyConfig classifyConfig;

	bool show_version = false;
	app.add_flag("-v,--version", show_version, "Show version information");
	// Create subcommands
	auto build = app.add_subcommand("build", "Build a sequence database");
	auto classify = app.add_subcommand("classify", "Classify sequences");
	auto profile = app.add_subcommand("profile", "Generate sequence profile");

	// Build
	build->add_option("-i,--input", buildConfig.input_file, "Input file for building")
		->required()
		->check(CLI::ExistingFile);
	build->add_option("-o,--output", buildConfig.output_file, "Output file for building")
		->default_val("ChimeraDB");
	build->add_option("-m,--mode", buildConfig.mode, "Mode for building")
		->default_val("default");
	build->add_option("-k,--kmer", buildConfig.kmer_size, "Kmer size for building")
		->default_val(19)
		->check(CLI::Range(1, 31));
	build->add_option("-w,--window", buildConfig.window_size, "Window size for building")
		->default_val(31);
	build->add_option("-l,--min-length", buildConfig.min_length, "Minimum length sequence for building")
		->default_val(0);
	build->add_option("-t,--threads", buildConfig.threads, "Number of threads for building")
		->default_val(32);
	build->add_option("--load-factor", buildConfig.load_factor, "Loading ratio of ICF")
		->default_val(0.95);
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
		->default_val(10);
	classify->add_option("-n,--batch-num", classifyConfig.batchNum, "Batch number for classifying")
		->default_val(1000);
	classify->add_flag("-q,--quiet", classifyConfig.verbose, "Quiet output")->default_val(true)->disable_flag_override();




	// Profile
	std::string profile_input_file;
	profile->add_option("-i,--input", profile_input_file, "Input file for profiling")
		->required()
		->check(CLI::ExistingFile);



	if (argc == 1) {
		std::cout << app.help() << std::endl;
		return 0;
	}

	CLI11_PARSE(app, argc, argv);

	if (show_version) {
		std::cout << "Chimera version: " << VERSION_INFO << std::endl;
		return 0;
	}

	if (*build) {
		ChimeraBuild::run(buildConfig);
	}
	else if (*classify) {
		ChimeraClassify::run(classifyConfig);
	}
	else if (*profile) {
		std::cout << "Profiling sequences..." << std::endl;
		std::cout << "Input file: " << profile_input_file << std::endl;
	}

	return 0;
}