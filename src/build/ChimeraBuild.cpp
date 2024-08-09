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
            std::cout << hours << "h " << minutes << "min " << seconds << "s " << milliseconds % 1000 << "ms" << std::endl;
        }
        else if (minutes > 0) {
            std::cout << minutes << "min " << seconds << "s " << milliseconds % 1000 << "ms" << std::endl;
        }
        else {
            std::cout << seconds << "s " << milliseconds % 1000 << "ms" << std::endl;
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
    * Get the maximum value from the hashCount map.
    *
    * @param hashCount The map containing the values.
    * @return The maximum value.
    */
	uint64_t getMaxValue(const robin_hood::unordered_map<std::string, uint64_t>& hashCount) {
		uint64_t maxValue = 0;
		for (const auto& kv : hashCount) {
			if (kv.second > maxValue) {
				maxValue = kv.second;
			}
		}
		return maxValue;
	}

    /**
    * Calculate the total size of all values in the hashCount map.
    *
    * @param hashCount The map containing the values.
    * @return The total size.
    */
    uint64_t calculateTotalSize(const robin_hood::unordered_map<std::string, uint64_t>& hashCount) {
    uint64_t totalSize = 0;
    for (const auto& kv : hashCount) {
    totalSize += kv.second;
    }
    return totalSize;
    }

    /**
     * Calculate the filter size based on the hash count and configuration.
     *
     * @param hashCount The map containing the hash count for each taxid.
     * @param icfConfig The configuration for the Interleaved Cuckoo Filter.
     * @param load_factor The desired load factor for the filter.
     * @param mode The mode of the calculation.
     */
    void calculateFilterSize(robin_hood::unordered_map<std::string, uint64_t>& hashCount, ICFConfig& icfConfig, double load_factor, std::string mode) {
        // Get the maximum value from the hashCount map
        uint64_t maxValue = getMaxValue(hashCount);
        // Calculate the total size of all values in the hashCount map
        uint64_t totalSize = calculateTotalSize(hashCount);
        // Set the initial iteration value
        size_t it = 100;
        if (maxValue < it) {
            it = maxValue;
        }
        uint64_t binSize, binNum, oldbinSize, oldbinNum;
		bool ifTheFirstIsFull = true;
        double load;
        // Iterate from the maximum value to the iteration value
        for (size_t i = maxValue + 1; i > it; i -= it) {
            binSize = i - 1;
            binNum = 0;
            // Calculate the number of bins for each taxid
            for (auto const& [taxid, count] : hashCount) {
                binNum += std::ceil(count / static_cast<double>(binSize));
            }
            // Calculate the load factor
            load = totalSize / static_cast<double>(binNum * binSize);
            // Check if the load factor exceeds the specified load factor
            if (load > load_factor) {
				// If  find that the load exceeds the load factor for the first time, multiply i by 2 first
				if (i == maxValue + 1 && ifTheFirstIsFull)
				{
					i = maxValue * 2;
					ifTheFirstIsFull = false;
					continue;
				}
                // Set the configuration values to the previous values
                icfConfig.bins = oldbinNum;
                icfConfig.bin_size = oldbinSize;
                break;
            } else {
                // Update the previous values
                oldbinSize = binSize;
                oldbinNum = binNum;
                continue;
            }
        }
    }
	// 计算 taxid 映射到的 bins
	robin_hood::unordered_map<std::string, std::size_t> calculateTaxidMapBins(
		const ICFConfig& config,
		const robin_hood::unordered_map<std::string, uint64_t>& hashCount) {

		robin_hood::unordered_map<std::string, std::size_t> taxidBins;
		size_t totalBins = 0;

		for (const auto& [taxid, count] : hashCount) {
			size_t binNum = static_cast<size_t>(std::ceil(count / static_cast<double>(config.bin_size)));
			totalBins += binNum;
			taxidBins[taxid] = totalBins;
		}

		return taxidBins;
	}
	void build(robin_hood::unordered_map<std::string, std::size_t> taxidBins, 
		ICFConfig config, chimera::InterleavedCuckooFilter& icf, 
		const robin_hood::unordered_map<std::string, uint64_t>& hashCount,
		robin_hood::unordered_map<std::string, std::vector<std::string>> inputFiles){
		size_t old = 0;
		for (const auto& [taxid, count] : hashCount) {
			//打开文件
			std::ifstream ifile("tmp/" + taxid + ".mini", std::ios::binary);
			if (ifile.fail()) {
				std::cerr << "Failed to open minimiser file: " << taxid << ".mini" << std::endl;
				continue;
			}
			//读取文件
			uint64_t hash;
			size_t currentPos = old; // 用于追踪当前存储位置

			while (ifile.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
				// 插入到 filter 中，将 hash 存储在 currentPos 位置
				icf.insertTag(currentPos, hash);

				// 更新存储位置到下一个位置
				currentPos++;

				// 如果 currentPos 超过 taxidBins[taxid]，则重置为 old
				if (currentPos == taxidBins[taxid]) {
					currentPos = old;
				}
			}

			// 更新 old 的值为当前 taxid 的下一个位置，以备后续使用
			old = taxidBins[taxid];
			ifile.close();
			//删除minimiser文件
			std::filesystem::remove("tmp/" + taxid + ".mini");
			}
		
	}

	void saveFilter(const std::string& output_file,
		const chimera::InterleavedCuckooFilter& icf,
		ICFConfig& icfConfig,
		const robin_hood::unordered_map<std::string, uint64_t>& hashCount,
		robin_hood::unordered_map<std::string, std::size_t> taxidBins) {

		// 打开输出文件
		std::ofstream os(output_file, std::ios::binary);

		// 检查文件是否成功打开
		if (!os.is_open()) {
			throw std::runtime_error("Failed to open file: " + output_file);
		}

		// 创建cereal的二进制档案
		cereal::BinaryOutputArchive archive(os);

		// 序列化 InterleavedCuckooFilter
		archive(icf);

		// 序列化 ICFConfig
		archive(icfConfig);

		// 手动将 robin_hood::unordered_map 数据转换并提取为 vector
		std::vector<std::pair<std::string, uint64_t>> hashCountData;
		for (const auto& kv : hashCount) {
			hashCountData.emplace_back(kv.first, kv.second);
		}

		std::vector<std::pair<std::string, std::size_t>> taxidBinsData;
		for (const auto& kv : taxidBins) {
			taxidBinsData.emplace_back(kv.first, kv.second);
		}

		// 序列化 vector
		archive(hashCountData);
		archive(taxidBinsData);

		// 关闭文件
		os.close();
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
			//输出fileinfo的内容
			std::cout << "File information:" << std::endl;
			std::cout << "Number of files: " << fileInfo.fileNum << std::endl;
			std::cout << "Number of invalid files: " << fileInfo.invalidNum << std::endl;
			std::cout << "Number of sequences: " << fileInfo.sequenceNum << std::endl;
			std::cout << "Number of skipped sequences: " << fileInfo.skippedNum << std::endl;
			std::cout << "Total base pairs: " << fileInfo.bpLength << std::endl << std::endl;
		}


		auto calculate_filter_size_start = std::chrono::high_resolution_clock::now();
		std::cout << "Calculating filter size..." << std::endl;
		calculateFilterSize(hashCount, icfConfig, config.load_factor, config.mode);
		auto calculate_filter_size_end = std::chrono::high_resolution_clock::now();
		auto calculate_filter_size_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(calculate_filter_size_end - calculate_filter_size_start).count();
		if (config.verbose) {
			std::cout << "Calculate filter size time: ";
			print_build_time(calculate_filter_size_total_time);
			std::cout << std::endl;
		}

		auto create_filter_start = std::chrono::high_resolution_clock::now();
		std::cout << "Creating filter..." << std::endl;
		chimera::InterleavedCuckooFilter icf(icfConfig.bins, icfConfig.bin_size);
		robin_hood::unordered_map<std::string, std::size_t> taxidBins = calculateTaxidMapBins(icfConfig, hashCount);
		build(taxidBins, icfConfig, icf, hashCount, inputFiles);
		saveFilter(config.output_file, icf, icfConfig, hashCount, taxidBins);
		auto create_filter_end = std::chrono::high_resolution_clock::now();
		auto create_filter_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(create_filter_end - create_filter_start).count();
		if (config.verbose) {
			std::cout << "Create filter time: ";
			print_build_time(create_filter_total_time);
			std::cout << std::endl;
		}

		auto build_end = std::chrono::high_resolution_clock::now();

		// Calculate the total build time in milliseconds
		auto build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
		if (config.verbose) {
			std::cout << "Total build time: ";
			print_build_time(build_total_time);
		}
	}
}