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

	build->add_flag("-V,--verbose", buildConfig.verbose, "Verbose output");

	// Classify
	std::string classify_input_file;
	std::string model_file;
	classify->add_option("-i,--input", classify_input_file, "Input file for classification")
		->required()
		->check(CLI::ExistingFile);
	classify->add_option("-m,--model", model_file, "Model file for classification")
		->required()
		->check(CLI::ExistingFile);

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
		std::cout << "Classifying sequences..." << std::endl;
		std::cout << "Input file: " << classify_input_file << std::endl;
		std::cout << "Model file: " << model_file << std::endl;
	}
	else if (*profile) {
		std::cout << "Profiling sequences..." << std::endl;
		std::cout << "Input file: " << profile_input_file << std::endl;
	}

	return 0;
}