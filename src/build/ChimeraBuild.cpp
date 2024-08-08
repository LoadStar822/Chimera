/*
 * -----------------------------------------------------------------------------
 * Filename:      ChimeraBuild.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-07-30
 *
 * Last Modified: 2024-08-06
 *
 * Description:
 *  This is a simple C++ program that outputs "Hello, World!".
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#include <ChimeraBuild.hpp>

namespace ChimeraBuild {
	/**
	* Print the build time in a human-readable format.
	*
	* @param milliseconds The build time in milliseconds.
	*/
	void print_build_time(long long milliseconds) {
		// Calculate seconds, minutes, and hours
		long long total_seconds = milliseconds / 1000;
		long long seconds = total_seconds % 60;
		long long total_minutes = total_seconds / 60;
		long long minutes = total_minutes % 60;
		long long hours = total_minutes / 60;

		// Output different formats based on the length of time
		if (hours > 0) {
			std::cout << hours << "h " << minutes << "min " << seconds << "s" << std::endl;
		}
		else {
			std::cout << minutes << "min " << seconds << "s" << std::endl;
		}
	}

	/**
	* Parse the input file and populate the inputFiles and hashCount maps.
	*
	* @param filePath The path to the input file.
	* @param inputFiles The map to store the input files.
	* @param hashCount The map to store the hash count.
	* @param fileInfo The struct to store file information.
	*/
	void parseInputFile(const std::string& filePath, robin_hood::unordered_map<std::string, std::vector<std::string>>& inputFiles, robin_hood::unordered_map<std::string, uint64_t>& hashCount, FileInfo& fileInfo) {
		// Open the input file
		std::ifstream inputFile(filePath);
		if (!inputFile.is_open()) {
			std::cerr << "Failed to open input file: " << filePath << std::endl;
			return;
		}
		std::string line;
		while (std::getline(inputFile, line)) {
			std::istringstream iss(line);
			std::string filePath;
			std::string taxidStr;
			if (!(iss >> filePath >> taxidStr)) {
				std::cerr << "Failed to parse line: " << line << std::endl;
				fileInfo.invalidNum++;
				continue;
			}
			hashCount[taxidStr] = 0;
			inputFiles[taxidStr].push_back(filePath);
			fileInfo.fileNum++;
		}

		// Close the input file
		inputFile.close();
	}

	/**
	* Create or reset the directory specified by the dir parameter.
	*
	* @param dir The directory path.
	* @param config The build configuration.
	*/
	void createOrResetDirectory(const std::string& dir, const BuildConfig& config) {
		try {
			std::filesystem::path directoryPath = dir;

			// Check if the directory exists
			if (std::filesystem::exists(directoryPath)) {
				if (std::filesystem::is_directory(directoryPath)) {
					// If it is a directory, remove it
					std::filesystem::remove_all(directoryPath);
					if (config.verbose) {
						std::cout << "Directory '" << dir << "' existed and was removed." << std::endl;
					}
				}
				else {
					// If it is a file instead of a directory, print an error message and return
					if (config.verbose) {
						std::cerr << "'" << dir << "' exists but is not a directory, can't be replaced." << std::endl;
					}
					return;
				}
			}

			// Create a new directory
			if (std::filesystem::create_directory(directoryPath)) {
				if (config.verbose) {
					std::cout << "Directory '" << dir << "' created successfully." << std::endl;
				}
			}
			else {
				if (config.verbose) {
					std::cerr << "Failed to create directory '" << dir << "'." << std::endl;
				}
			}
		}
		catch (const std::filesystem::filesystem_error& e) {
			if (config.verbose) {
				std::cerr << "Error: " << e.what() << std::endl;
			}
		}
	}

	/**
	 * Adjust the seed value based on the kmer size.
	 *
	 * @param kmer_size The size of the kmer.
	 * @param seed The seed value to adjust (default: 0x8F3F73B5CF1C9ADEULL).
	 * @return The adjusted seed value.
	 */
	inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept {
		// Right shift the seed by (64 - 2 * kmer_size) bits
		return seed >> (64u - 2u * kmer_size);
	}

	/**
	 * Count the minimisers for each taxid and its associated files.
	 *
	 * @param config The build configuration.
	 * @param inputFiles The map of input files.
	 * @param hashCount The map to store the hash count.
	 * @param fileInfo The struct to store file information.
	 */
	void minimiser_count(
		BuildConfig& config,
		robin_hood::unordered_map<std::string, std::vector<std::string>>& inputFiles,
		robin_hood::unordered_map<std::string, uint64_t>& hashCount,
		FileInfo& fileInfo)
	{
		// Create a view for generating minimiser hashes
		auto minimiser_view = seqan3::views::minimiser_hash(
			seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
			seqan3::window_size{ config.window_size },
			seqan3::seed{ adjust_seed(config.kmer_size) });

		// Parallelize the task using dynamic scheduling
#pragma omp parallel
		{
#pragma omp single
			{
				// Iterate over each taxid and its associated files
				for (auto& [taxid, files] : inputFiles)
				{
#pragma omp task firstprivate(taxid, files)
					{
						// Create a private FileInfo object for each thread
						FileInfo threadFileInfo;

						// Process each file
						for (auto& file : files)
						{
							// Open the sequence file for reading
							seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq >> fin{ file };
							robin_hood::unordered_set<uint64_t> hashes;

							// Iterate over each sequence in the file
							for (auto const& [header, seq] : fin)
							{
								// Skip sequences that are shorter than the minimum length
								if (seq.size() < config.min_length)
								{
									threadFileInfo.skippedNum++;
									continue;
								}
								threadFileInfo.sequenceNum++;
								threadFileInfo.bpLength += seq.size();

								// Generate minimiser hashes for the sequence and insert them into the set
								const auto minihash = seq | minimiser_view | std::views::common;
								hashes.insert(minihash.begin(), minihash.end());
							}

							// Update the hashCount map using a critical section
#pragma omp critical
							{
								hashCount[taxid] += hashes.size();
							}

							// Write the minimiser hashes to a file
							std::ofstream ofile("tmp/" + taxid + ".mini", std::ios::binary | std::ios::app);
							for (auto hash : hashes)
							{
								ofile.write(reinterpret_cast<char*>(&hash), sizeof(hash));
							}
							ofile.close();
							if (ofile.fail())
							{
								std::cerr << "Failed to write minimiser file: " << taxid << ".mini" << std::endl;
							}
						}

						// Merge the thread's private FileInfo into the global FileInfo
						#pragma omp critical
						{
							fileInfo.skippedNum += threadFileInfo.skippedNum;
							fileInfo.sequenceNum += threadFileInfo.sequenceNum;
							fileInfo.bpLength += threadFileInfo.bpLength;
						}
					}
				}
			}
		}
	}

	/**
	* Run the build process with the given build configuration.
	*
	* @param config The build configuration.
	*/
	void run(BuildConfig config) {
		if (config.verbose) {
			std::cout << config << std::endl;
		}
		omp_set_num_threads(config.threads);
		auto build_start = std::chrono::high_resolution_clock::now();

		auto read_start = std::chrono::high_resolution_clock::now();
		std::cout << "Reading input files..." << std::endl;
		ICFConfig icfConfig;
		FileInfo fileInfo;
		icfConfig.kmer_size = config.kmer_size;
		icfConfig.window_size = config.window_size;
		robin_hood::unordered_map<std::string, uint64_t> hashCount;
		robin_hood::unordered_map<std::string, std::vector<std::string>> inputFiles;
		parseInputFile(config.input_file, inputFiles, hashCount, fileInfo);
		auto read_end = std::chrono::high_resolution_clock::now();
		auto read_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(read_end - read_start).count();
		if (config.verbose) {
			std::cout << "Read time: ";
			print_build_time(read_total_time);
			std::cout << std::endl;
		}


		auto calculate_start = std::chrono::high_resolution_clock::now();
		std::cout << "Calculating minimizers..." << std::endl;
		std::filesystem::path dir = "tmp";
		createOrResetDirectory(dir, config);
		minimiser_count(config, inputFiles, hashCount, fileInfo);
		auto calculate_end = std::chrono::high_resolution_clock::now();
		auto calculate_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(calculate_end - calculate_start).count();
		if (config.verbose) {
			std::cout << "Calculate time: ";
			print_build_time(calculate_total_time);
			//Êä³öfileinfoµÄÄÚÈÝ
			std::cout << "File information:" << std::endl;
			std::cout << "Number of files: " << fileInfo.fileNum << std::endl;
			std::cout << "Number of invalid files: " << fileInfo.invalidNum << std::endl;
			std::cout << "Number of sequences: " << fileInfo.sequenceNum << std::endl;
			std::cout << "Number of skipped sequences: " << fileInfo.skippedNum << std::endl;
			std::cout << "Total base pairs: " << fileInfo.bpLength << std::endl << std::endl;
		}


		auto calculate_filter_size_start = std::chrono::high_resolution_clock::now();
		std::cout << "Calculating filter size..." << std::endl;


		auto calculate_filter_size_end = std::chrono::high_resolution_clock::now();
		auto calculate_filter_size_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(calculate_filter_size_end - calculate_filter_size_start).count();
		if (config.verbose) {
			std::cout << "Calculate filter size time: ";
			print_build_time(calculate_filter_size_total_time);
		}

		auto build_end = std::chrono::high_resolution_clock::now();

		// Calculate the total build time in milliseconds
		auto build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
	}
}