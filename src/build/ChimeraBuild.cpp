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
 * Last Modified: 2024-09-19
 *
 * Description:
 *  The main program of ChimeraBuild.
 *
 * Version:
 *  1.2
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

	/*
	 * Check if the file is compressed based on the file extension.
	 *
	 * @param filepath The path to the file.
	 * @return True if the file is compressed, false otherwise.
	 */
	static inline bool file_is_compressed(const std::filesystem::path& filepath)
	{
		std::string extension = filepath.extension().string();
		return extension == ".gz" || extension == ".bgzf" || extension == ".bz2";
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
		// Define the minimiser view
		auto minimiser_view = seqan3::views::minimiser_hash(
			seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
			seqan3::window_size{ config.window_size },
			seqan3::seed{ adjust_seed(config.kmer_size) });

		// Expand the pairs of taxid and files
		std::vector<std::pair<std::string, std::string>> taxid_file_pairs;
		robin_hood::unordered_map<std::string, size_t> taxid_to_index;
		std::vector<std::string> index_to_taxid;
		size_t index = 0;
		for (auto& [taxid, files] : inputFiles)
		{
			for (auto& file : files)
			{
				taxid_file_pairs.emplace_back(taxid, file);
			}
			taxid_to_index[taxid] = index++;
			index_to_taxid.push_back(taxid);
		}
		std::vector<std::mutex> taxidMutexes(index);
		// Global file information statistics
		FileInfo globalFileInfo = {};

		// Mutex to protect global file information
		std::mutex fileInfo_mutex;

		// OpenMP parallel processing
#pragma omp parallel for schedule(dynamic)
		for (size_t idx = 0; idx < taxid_file_pairs.size(); ++idx)
		{
			const auto& [taxid, filename] = taxid_file_pairs[idx];

			// Thread-local file information
			FileInfo localFileInfo = {};

			// Thread-local hash count
			robin_hood::unordered_map<uint64_t, uint8_t> local_hash_counts;

			// Calculate the threshold
			uint8_t cutoff = 1;
			if (config.fixed_cutoff > 0)
			{
				cutoff = config.fixed_cutoff;
			}
			else
			{
				// Calculate the threshold based on file size and type
				std::filesystem::path filepath(filename);
				size_t filesize = std::filesystem::file_size(filepath);
				bool is_compressed = file_is_compressed(filepath);

				size_t adjusted_filesize = filesize * 2 / (is_compressed ? 1 : 3);

				if (adjusted_filesize <= 314'572'800ULL) // 300 MB
					cutoff = 1;
				else if (adjusted_filesize <= 524'288'000ULL) // 500 MB
					cutoff = 3;
				else if (adjusted_filesize <= 1'073'741'824ULL) // 1 GB
					cutoff = 10;
				else if (adjusted_filesize <= 3'221'225'472ULL) // 3 GB
					cutoff = 20;
				else
					cutoff = 50;
			}

			// Set the maximum hash limit
			size_t max_hashes = config.max_hashes_per_taxid;

			// Open the sequence file
			seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq >> fin{ filename };

			// Iterate over the records in the sequence file
			for (auto& record : fin)
			{
				auto& seq = record.sequence();

				if (seq.size() < config.min_length)
				{
					localFileInfo.skippedNum++;
					continue;
				}
				localFileInfo.sequenceNum++;
				localFileInfo.bpLength += seq.size();

				// Compute minimizers and count hash values
				for (uint64_t hash : seq | minimiser_view)
				{
					uint8_t& count = local_hash_counts[hash];
					if (count < 255)
						++count;
				}
			}

			// Compute minimizers and count hash values
			std::vector<uint64_t> filtered_hashes;
			filtered_hashes.reserve(local_hash_counts.size());
			for (const auto& [hash, count] : local_hash_counts)
			{
				if (count >= cutoff)
				{
					filtered_hashes.push_back(hash);
				}
			}

			// Apply maximum hash limit if set
			if (max_hashes > 0 && filtered_hashes.size() > max_hashes)
			{
				// Optionally, sort hashes by value or count before truncating
				// Here, simply truncate
				filtered_hashes.resize(max_hashes);
			}

			// Write to the minimiser file
			size_t taxid_index = taxid_to_index[taxid];
			{
				// Write to the minimiser file
				std::string output_filename = "tmp/" + taxid + ".mini";

				// Lock is required to prevent multiple threads from writing to the same taxid file
				std::lock_guard<std::mutex> lock(taxidMutexes[taxid_index]);


				// Lock is required to prevent multiple threads from writing to the same taxid file
				std::ofstream ofile(output_filename, std::ios::binary | std::ios::app);
				if (!ofile.is_open())
				{
					std::cerr << "Unable to open the minimiser file: " << output_filename << std::endl;
					continue;
				}

				for (uint64_t hash : filtered_hashes)
				{
					ofile.write(reinterpret_cast<const char*>(&hash), sizeof(hash));
				}
			}

			// Update global hash count
			{
				std::lock_guard<std::mutex> lock(fileInfo_mutex);
				hashCount[taxid] += filtered_hashes.size();
				fileInfo.skippedNum += localFileInfo.skippedNum;
				fileInfo.sequenceNum += localFileInfo.sequenceNum;
				fileInfo.bpLength += localFileInfo.bpLength;
			}
		}
	}

	/**
	 * Count and filter minimiser hashes for each taxonomic ID, store in HyperLogLog sketches, and
	 * output filtered hashes to temporary files. This function iterates over input files, extracts
	 * minimiser hashes based on k-mer and window size, and applies a cutoff to control the number of
	 * stored hashes.
	 *
	 * @param config The configuration object containing parameters like k-mer size, window size,
	 *               cutoff values, and maximum hashes per taxid.
	 * @param inputFiles A map where each key is a taxonomic ID (taxid) associated with a list of file paths
	 *                   for sequences associated with that taxid.
	 * @param hashCount A map to store the total number of hashes processed for each taxid.
	 * @param fileInfo A reference to a FileInfo structure, updated with statistics on skipped, processed
	 *                 sequences, and total base pairs.
	 * @param hllVec A vector of HyperLogLog objects used to estimate cardinality for each taxid, updated
	 *               with minimiser hashes.
	 */
	std::vector<std::string> minimiser_count(
		BuildConfig& config,
		robin_hood::unordered_map<std::string, std::vector<std::string>>& inputFiles,
		size_t& hashCounts,
		FileInfo& fileInfo,
		std::vector<HyperLogLog>& hllVec)
	{
		auto minimiser_view = seqan3::views::minimiser_hash(
			seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
			seqan3::window_size{ config.window_size },
			seqan3::seed{ adjust_seed(config.kmer_size) });

		std::vector<std::pair<std::string, std::string>> taxid_file_pairs;
		// build taxid to hllVec index map
		robin_hood::unordered_map<std::string, size_t> taxid_to_index;
		std::vector<std::string> index_to_taxid;
		size_t index = 0;
		for (auto& [taxid, files] : inputFiles)
		{
			for (auto& file : files)
			{
				taxid_file_pairs.emplace_back(taxid, file);
			}
			taxid_to_index[taxid] = index++;
			index_to_taxid.push_back(taxid);
		}

		std::vector<std::mutex> hllMutexes(hllVec.size());

		FileInfo globalFileInfo = {};

		std::mutex fileInfo_mutex;

#pragma omp parallel for schedule(dynamic)
		for (size_t idx = 0; idx < taxid_file_pairs.size(); ++idx)
		{
			const auto& [taxid, filename] = taxid_file_pairs[idx];

			FileInfo localFileInfo = {};

			robin_hood::unordered_map<uint64_t, uint8_t> local_hash_counts;

			uint8_t cutoff = 1;
			if (config.fixed_cutoff > 0)
			{
				cutoff = config.fixed_cutoff;
			}
			else
			{
				std::filesystem::path filepath(filename);
				size_t filesize = std::filesystem::file_size(filepath);
				bool is_compressed = file_is_compressed(filepath);

				size_t adjusted_filesize = filesize * 2 / (is_compressed ? 1 : 3);

				if (adjusted_filesize <= 314'572'800ULL)
					cutoff = 1;
				else if (adjusted_filesize <= 524'288'000ULL)
					cutoff = 3;
				else if (adjusted_filesize <= 1'073'741'824ULL)
					cutoff = 10;
				else if (adjusted_filesize <= 3'221'225'472ULL)
					cutoff = 20;
				else
					cutoff = 50;
			}

			size_t max_hashes = config.max_hashes_per_taxid;

			seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq >> fin{ filename };

			for (auto& record : fin)
			{
				auto& seq = record.sequence();

				if (seq.size() < config.min_length)
				{
					localFileInfo.skippedNum++;
					continue;
				}
				localFileInfo.sequenceNum++;
				localFileInfo.bpLength += seq.size();

				for (uint64_t hash : seq | minimiser_view)
				{
					uint8_t& count = local_hash_counts[hash];
					if (count < 255)
						++count;
				}
			}

			std::vector<uint64_t> filtered_hashes;
			filtered_hashes.reserve(local_hash_counts.size());
			for (const auto& [hash, count] : local_hash_counts)
			{
				if (count >= cutoff)
				{
					filtered_hashes.push_back(hash);
				}
			}

			if (max_hashes > 0 && filtered_hashes.size() > max_hashes)
			{
				filtered_hashes.resize(max_hashes);
			}

			size_t hll_index = taxid_to_index[taxid];
			{
				std::string output_filename = "tmp/" + taxid + ".mini";

				std::lock_guard<std::mutex> lock(hllMutexes[hll_index]);

				std::ofstream ofile(output_filename, std::ios::binary | std::ios::app);
				if (!ofile.is_open())
				{
					std::cerr << "Unable to open the minimiser file: " << output_filename << std::endl;
					continue;
				}

				for (uint64_t hash : filtered_hashes)
				{
					hllVec[hll_index].add(hash);
					ofile.write(reinterpret_cast<const char*>(&hash), sizeof(hash));
				}
			}

			{
				std::lock_guard<std::mutex> lock(fileInfo_mutex);
				hashCounts += filtered_hashes.size();
				fileInfo.skippedNum += localFileInfo.skippedNum;
				fileInfo.sequenceNum += localFileInfo.sequenceNum;
				fileInfo.bpLength += localFileInfo.bpLength;
			}
		}

		return index_to_taxid;
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
		uint64_t maxValue = getMaxValue(hashCount);
		uint64_t totalSize = calculateTotalSize(hashCount);

		uint64_t minBinSize = 1;
		uint64_t maxBinSize = maxValue * 2;
		uint64_t bestBinSize = maxBinSize;
		uint64_t bestBinNum = 0;
		double bestLoad = 0.0;
		double bestLoadDiff = std::numeric_limits<double>::max();

		std::vector<uint64_t> counts;
		counts.reserve(hashCount.size());
		for (const auto& kv : hashCount) {
			counts.push_back(kv.second);
		}

		while (minBinSize <= maxBinSize) {
			uint64_t binSize = (minBinSize + maxBinSize) / 2;
			uint64_t binNum = 0;

			// could choice use openmp to parallel
			//#pragma omp parallel for reduction(+:binNum)
			for (size_t idx = 0; idx < counts.size(); ++idx) {
				uint64_t count = counts[idx];
				binNum += (count + binSize - 1) / binSize;
			}

			double load = static_cast<double>(totalSize) / (binNum * binSize);
			double loadDiff = std::abs(load - load_factor);

			if (loadDiff < bestLoadDiff) {
				bestLoadDiff = loadDiff;
				bestBinSize = binSize;
				bestBinNum = binNum;
				bestLoad = load;
			}

			if (load < load_factor) {
				maxBinSize = binSize - 1;
			}
			else if (load > load_factor) {
				minBinSize = binSize + 1;
			}
			else {
				break;
			}
		}

		icfConfig.bins = bestBinNum;
		icfConfig.bin_size = bestBinSize;
	}

	/**
	* Calculate the bins mapped to the taxid.
	*
	* @param config The configuration for the Interleaved Cuckoo Filter.
	* @param hashCount The map containing the hash count for each taxid.
	* @return The map containing the bins mapped to the taxid.
	*/
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

	/**
	 * Process the taxid by inserting the hashes into the Interleaved Cuckoo Filter.
	 *
	 * @param taxid The taxid to process.
	 * @param start The starting position in the filter.
	 * @param end The ending position in the filter.
	 * @param icf The Interleaved Cuckoo Filter to insert the hashes into.
	 */
	void processTaxid(
		const std::string& taxid,
		size_t start,
		size_t end,
		chimera::InterleavedCuckooFilter& icf) {
		// Open the minimiser file for reading
		std::ifstream ifile("tmp/" + taxid + ".mini", std::ios::binary);
		if (ifile.fail()) {
			std::cerr << "Failed to open minimiser file: " << taxid << ".mini" << std::endl;
			return;
		}

		uint64_t hash;
		size_t currentPos = start;

		// Read the hashes from the file and insert them into the filter
		while (ifile.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
			// Insert the hash into the filter at the current position
			//std::cout << "Inserting hash: " << hash << " at position: " << currentPos << "\n";
			icf.insertTag(currentPos, hash);
			// Update the current position to the next position
			currentPos++;

			// If the current position reaches the end, reset it to the start
			if (currentPos == end) {
				currentPos = start;
			}
		}

		// Close the minimiser file
		ifile.close();
	}

	/**
	 * Build the Interleaved Cuckoo Filter.
	 *
	 * @param taxidBins The map containing the bins mapped to the taxid.
	 * @param config The configuration for the Interleaved Cuckoo Filter.
	 * @param icf The Interleaved Cuckoo Filter object to build.
	 * @param hashCount The map containing the hash count for each taxid.
	 * @param inputFiles The map containing the input files for each taxid.
	 * @param numThreads The number of threads to use for building.
	 */
	void build(
		const robin_hood::unordered_map<std::string, std::size_t>& taxidBins,
		ICFConfig config,
		chimera::InterleavedCuckooFilter& icf,
		const robin_hood::unordered_map<std::string, uint64_t>& hashCount,
		robin_hood::unordered_map<std::string, std::vector<std::string>> inputFiles) {
		std::vector<std::tuple<std::string, size_t, size_t, uintmax_t>> taxid_info_pairs;
		taxid_info_pairs.reserve(hashCount.size());

		size_t previousEnd = 0;
		for (const auto& [taxid, count] : hashCount) {
			size_t currentEnd = taxidBins.at(taxid);

			std::string miniFilePath = "tmp/" + taxid + ".mini";
			uintmax_t fileSize = 0;
			if (fs::exists(miniFilePath)) {
				fileSize = fs::file_size(miniFilePath);
			}
			else {
				std::cerr << "File does not exist: " << miniFilePath << std::endl;
				continue;
			}

			taxid_info_pairs.emplace_back(taxid, previousEnd, currentEnd, fileSize);
			previousEnd = currentEnd;
		}

		std::sort(taxid_info_pairs.begin(), taxid_info_pairs.end(),
			[](const auto& a, const auto& b) {
				return std::get<3>(a) > std::get<3>(b);
			});

		int numThreads = omp_get_max_threads();

		std::vector<std::vector<std::tuple<std::string, size_t, size_t>>> threadTasks(numThreads);
		std::vector<uintmax_t> threadLoad(numThreads, 0);

		for (const auto& task : taxid_info_pairs) {
			const auto& [taxid, start, end, fileSize] = task;
			int minThread = std::distance(threadLoad.begin(), std::min_element(threadLoad.begin(), threadLoad.end()));
			threadTasks[minThread].emplace_back(taxid, start, end);
			threadLoad[minThread] += fileSize;
		}

#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();
			for (const auto& task : threadTasks[threadId]) {
				const auto& [taxid, start, end] = task;
			processTaxid(taxid, start, end, icf);
		}
	}
	}

	/**
	 * Save the Interleaved Cuckoo Filter, ICFConfig, hashCount, and taxidBins to the output file.
	 *
	 * @param output_file The path to the output file.
	 * @param icf The Interleaved Cuckoo Filter to be saved.
	 * @param icfConfig The configuration for the Interleaved Cuckoo Filter.
	 * @param hashCount The map containing the hash count for each taxid.
	 * @param taxidBins The map containing the number of bins for each taxid.
	 */
	void saveFilter(const std::string& output_file,
		const chimera::InterleavedCuckooFilter& icf,
		ICFConfig& icfConfig,
		const robin_hood::unordered_map<std::string, uint64_t>& hashCount,
		robin_hood::unordered_map<std::string, std::size_t> taxidBins) {
		// Open the output file
		std::ofstream os(output_file, std::ios::binary);

		// Check if the file is successfully opened
		if (!os.is_open()) {
			throw std::runtime_error("Failed to open file: " + output_file);
		}

		// Create a cereal binary archive
		cereal::BinaryOutputArchive archive(os);

		// Serialize the Interleaved Cuckoo Filter
		archive(icf);

		// Serialize the ICFConfig
		archive(icfConfig);

		// Manually convert and extract the robin_hood::unordered_map data to vector
		std::vector<std::pair<std::string, uint64_t>> hashCountData;
		for (const auto& kv : hashCount) {
			hashCountData.emplace_back(kv.first, kv.second);
		}

		std::vector<std::pair<std::string, std::size_t>> taxidBinsData;
		for (const auto& kv : taxidBins) {
			taxidBinsData.emplace_back(kv.first, kv.second);
		}

		// Serialize the vector
		archive(hashCountData);
		archive(taxidBinsData);

		// Close the file
		os.close();

		// Get the file size
		std::uintmax_t fileSize = std::filesystem::file_size(output_file);

		// Output the file size, formatted as KB, MB, or GB
		if (fileSize >= 1024 * 1024 * 1024) {
			std::cout << "Filter file size: " << std::fixed << std::setprecision(2)
				<< static_cast<double>(fileSize) / (1024 * 1024 * 1024) << " GB" << std::endl;
		}
		else if (fileSize >= 1024 * 1024) {
			std::cout << "Filter file size: " << std::fixed << std::setprecision(2)
				<< static_cast<double>(fileSize) / (1024 * 1024) << " MB" << std::endl;
		}
		else if (fileSize >= 1024) {
			std::cout << "Filter file size: " << std::fixed << std::setprecision(2)
				<< static_cast<double>(fileSize) / 1024 << " KB" << std::endl;
		}
		else {
			std::cout << "Filter file size: " << fileSize << " bytes" << std::endl;
		}
	}


	/**
	 * Estimate the k-mer counts for each bin and store them in the provided data structure. This function
	 * iterates through the HyperLogLog (HLL) sketches stored in `data.hllVec` and calculates an estimated
	 * k-mer count for each bin, saving the results in `data.estimateKmerCounts`.
	 *
	 * @param data A DataStore object containing pointers to the HLL vector (`hllVec`) and the vector to store
	 *             estimated k-mer counts (`estimateKmerCounts`).
	 */
	void estimateKmerCount(chimera::hicf::DataStore& data)
	{
		{
#pragma omp parallel for
			for (size_t i = 0; i < data.hllVec->size(); ++i)
			{
				(*data.estimateKmerCounts)[i] = (*data.hllVec)[i].estimate();
			}
		}
	}

	/**
	 * Rearrange user-defined bins based on estimated k-mer counts and sort them to optimize data layout.
	 * This function initializes the location vector, sorts it in descending order of estimated k-mer
	 * counts, and optionally calls `rearrange_bins` if not in "fast" mode.
	 *
	 * @param data A reference to a DataStore object that holds vectors for the locations, estimated k-mer
	 *             counts, and HyperLogLog (HLL) sketches.
	 * @param config The BuildConfig object containing settings for the rearrangement, such as mode and
	 *               number of threads.
	 */
	void rearrangeUserbins(chimera::hicf::DataStore& data, BuildConfig& config)
	{
		// Sort the location based on the estimated kmer counts
		std::sort(data.location.begin(), data.location.end(),
			[&](size_t a, size_t b) -> bool {
				return (*data.estimateKmerCounts)[a] > (*data.estimateKmerCounts)[b];
			});
		if (config.mode != "fast")
		{
			auto start = std::chrono::high_resolution_clock::now();
			toolbox::rearrange_bins(*data.hllVec, *data.estimateKmerCounts, data.location, 0.5, config.threads);
			auto end = std::chrono::high_resolution_clock::now();
			data.totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}
	}


	/**
	 * @brief Calculate the maximum merge level required to allocate user-defined bins into technical bins.
	 *
	 * This function determines the appropriate merge level by calculating the logarithmic ratio between the
	 * number of user-defined bins and the adjusted number of technical bins. The adjusted technical bins number
	 * is rounded down to the nearest multiple of 64 to ensure alignment and optimal merging.
	 *
	 * @param userBinsNum The total number of user-defined bins that need to be allocated.
	 * @param technicalBinsMaxNum The maximum number of technical bins available for allocation.
	 *
	 * @return The maximum merge level as a `size_t`, representing the depth of merging required to fit all user bins
	 *         within the technical bins constraint.
	 *
	 * @details
	 * - Adjusted Technical Bins Number (`lowerTechnicalBinsNum`):
	 *   The function first adjusts the `technicalBinsMaxNum` by rounding it down to the nearest multiple of 64.
	 *   This is achieved by adding 63 to `userBinsNum`, performing a right bitwise shift by 6 (equivalent to dividing by 64),
	 *   and then left-shifting by 6 to multiply by 64. The `std::min` function ensures that the adjusted number does not exceed
	 *   `technicalBinsMaxNum`.
	 *
	 * - Merge Level Calculation (`level`):
	 *   The merge level is calculated as the logarithm of `userBinsNum` divided by the logarithm of `lowerTechnicalBinsNum`.
	 *   This represents the number of times bins need to be merged to fit all user bins within the available technical bins.
	 *
	 * - Ceiling of Merge Level:
	 *   The resulting `level` is rounded up using `std::ceil` to ensure that the calculated merge level is sufficient.
	 *   The final value is cast to `size_t` before being returned.
	 *
	 * @note
	 * - This function is marked as `inline` to suggest to the compiler that it should attempt to embed the function's
	 *   code at each point of call, potentially reducing the overhead of function calls for small, frequently called functions.
	 * - The alignment to multiples of 64 in `lowerTechnicalBinsNum` is likely chosen for performance optimization,
	 *   leveraging cache line sizes or other hardware-specific optimizations.
	 *
	 * @example
	 * ```cpp
	 * size_t userBins = 150;
	 * size_t maxTechBins = 256;
	 * size_t mergeLevel = maxMergeLevel(userBins, maxTechBins);
	 * std::cout << "Maximum Merge Level: " << mergeLevel << std::endl;
	 * ```
	 * This code calculates and prints the maximum merge level required when allocating 150 user bins into a maximum of 256 technical bins.
	 */
	inline size_t maxMergeLevel(size_t userBinsNum, size_t technicalBinsMaxNum)
	{
		size_t lowerTechnicalBinsNum = std::min(technicalBinsMaxNum, userBinsNum);
		double level = std::log(static_cast<double>(userBinsNum)) / std::log(static_cast<double>(lowerTechnicalBinsNum));
		return static_cast<size_t>(std::ceil(level));
	}


	/**
	 * Handles the merge operation during the backtracking phase of the Hierarchical Interleaved Cuckoo Filters (HICF) layout computation.
	 * This function merges multiple user-defined bins into a single technical bin, updates the layout with the merged bins,
	 * and tracks the maximum bin size and its corresponding identifier.
	 *
	 * @param nextJ      The index of the next user bin to consider for merging.
	 * @param traceJ     The current index of the user bin being processed.
	 * @param binNum     The number of user bins to merge into the current technical bin.
	 * @param binId      The identifier of the current technical bin where user bins are being merged.
	 * @param maxId      Reference to the variable tracking the identifier of the largest technical bin.
	 * @param maxSize    Reference to the variable tracking the size of the largest technical bin.
	 * @param data       Reference to the `DataStore` object containing layout information, estimated k-mer counts,
	 *                   HyperLogLog (HLL) vectors, union estimates, and timing information.
	 * @param config     Reference to the `BuildConfig` object that holds build settings such as mode and alpha value.
	 * @param hicfConfig Reference to the `HICFConfig` object containing HICF-specific configuration parameters like load factors.
	 * @param isFirstRow A boolean flag indicating whether the current operation is on the first row of the layout matrix.
	 *
	 * @details
	 * This function performs the following steps:
	 * 1. K-mer Count Retrieval and HLL Initialization:
	 *    - Retrieves the estimated k-mer count for the current user bin.
	 *    - Initializes a HyperLogLog (HLL) object. If the build mode is not "fast", it uses the existing HLL from `hllVec`; otherwise, it initializes a new HLL.
	 *
	 * 2. Preparing Lower-Level DataStore:
	 *    - Creates a `lowerData` DataStore to hold information for the lower-level bin assignments.
	 *    - Copies relevant pointers and initializes the location vector with the current user bin.
	 *    - Adjusts `traceJ` based on whether it is the first row.
	 *
	 * 3. Merging User Bins:
	 *    - Iterates backwards from `traceJ` to `nextJ`, merging user bins into the current technical bin.
	 *    - If not in "fast" mode, merges the HLLs and updates the total processing time.
	 *    - If in "fast" mode, updates the k-mer count directly.
	 *    - Adds the current location to `lowerData.location`.
	 *
	 * 4. Updating the Layout and Tracking Maximum Bin Size:
	 *    - Determines if the current operation is at the top level (i.e., no previous bins).
	 *    - Updates the `previous` bin indices and the number of bins.
	 *    - If the number of locations exceeds the maximum technical bins, it rearranges user bins and recomputes the layout.
	 *    - Otherwise, it performs a dynamic programming approach to allocate user bins to technical bins, updating the `matrix` and `trace` accordingly.
	 *    - Backtracks through the `trace` matrix to finalize bin assignments, handling both merge and split cases.
	 *    - Updates the layout's `maxBins` with the newly computed maximum bin information.
	 *    - Estimates the final k-mer count and updates the maximum bin size and identifier if necessary.
	 *    - Accumulates the total processing time from lower levels.
	 *
	 * @note
	 * - This function is typically called during the backtracking process to finalize the assignment of user bins to technical bins.
	 * - Proper synchronization is assumed to be handled externally if this function is called in a multithreaded context.
	 * - Auxiliary functions such as `rearrangeUserbins`, `computeLayoutHICF`, `backtrackMerge`, and `backtraceSplit` are used to manage bin assignments and layout computations.
	 */
	void backtrackMerge(size_t nextJ,
						size_t traceJ,
						size_t binId,
						size_t& maxId,
						size_t& maxSize,
						chimera::hicf::DataStore& data,
						BuildConfig& config,
						HICFConfig& hicfConfig,
						bool isFirstRow)
	{
		size_t currentCount = (*data.estimateKmerCounts)[data.location[traceJ]];
		HyperLogLog hll = (config.mode == "fast") ? HyperLogLog() : (*data.hllVec)[data.location[traceJ]];
		chimera::hicf::DataStore lowerData;
		lowerData.estimateKmerCounts = data.estimateKmerCounts;
		lowerData.hllVec = data.hllVec;
		lowerData.layout = data.layout;
		lowerData.location = { data.location[traceJ] };
		lowerData.loadFactor = data.loadFactor;
		traceJ -= !isFirstRow;
		while(traceJ != nextJ)
		{
			traceJ -= isFirstRow;
			if (config.mode != "fast")
			{
				auto Start = std::chrono::high_resolution_clock::now();
				hll.merge((*data.hllVec)[data.location[traceJ]]);
				auto End = std::chrono::high_resolution_clock::now();
				data.totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Start).count();
			}
			else
			{
				currentCount += (*data.estimateKmerCounts)[data.location[traceJ]];
			}
			lowerData.location.push_back(data.location[traceJ]);

			traceJ -= !isFirstRow;
		}

		bool isTopLevel = data.previous.empty();
		lowerData.previous = data.previous;
		lowerData.previous.binIndices.push_back(binId);
		lowerData.previous.binsNum += (isTopLevel ? "" : ";") + std::string{ "1" };
		size_t lowerMaxBin;
		if (lowerData.location.size() > hicfConfig.technicalBinsMaxNum)
		{
			rearrangeUserbins(lowerData, config);
			lowerMaxBin = computeLayoutHICF(hicfConfig, config, lowerData);
		}
		else
		{
			size_t userBinNum = lowerData.location.size();
			size_t technicalBinNum = (lowerData.previous.empty()) ? hicfConfig.technicalBinsMaxNum :
				std::min(userBinNum + 63, hicfConfig.technicalBinsMaxNum);

			std::vector<std::vector<size_t>> matrix(technicalBinNum);
			for (size_t i = 0; i < technicalBinNum; i++)
			{
				matrix[i].resize(userBinNum, std::numeric_limits<size_t>::max());
			}
			std::vector<std::vector<size_t>> trace(technicalBinNum);
			for (size_t i = 0; i < technicalBinNum; i++)
			{
				trace[i].resize(userBinNum, std::numeric_limits<size_t>::max());
			}
			size_t extraBins = technicalBinNum - userBinNum + 1;
			double userBin0 = static_cast<double>((*lowerData.estimateKmerCounts)[lowerData.location[0]]);
			for (size_t i = 0; i < extraBins; i++)
			{
				matrix[i][0] = static_cast<size_t>(std::ceil(userBin0 / (hicfConfig.loadFactor * (i + 1))));
			}
			for (size_t j = 1; j < userBinNum; j++)
			{
				double userBin = static_cast<double>((*lowerData.estimateKmerCounts)[lowerData.location[j]]);
				for (size_t i = j; i < j + extraBins; i++)
				{
					size_t minimum{ std::numeric_limits<size_t>::max() };
					for (size_t Split = j - 1; Split < j; Split++)
					{
						size_t correctedValue = static_cast<size_t>(std::ceil(userBin / (hicfConfig.relaxedLoadFactor * (i - Split))));
						size_t score = std::max(matrix[Split][j - 1], correctedValue);
						if (score < minimum)
						{
							minimum = score;
							trace[i][j] = Split;
						}
					}
					matrix[i][j] = minimum;
				}
			}
			size_t traceILow = technicalBinNum - 1;
			size_t traceJLow = userBinNum - 1;
			size_t maxIdLow;
			size_t maxSizeLow = 0;
			size_t binIdLow = 0;
			while (traceJLow > 0)
			{
				size_t nextILow = trace[traceILow][traceJLow];
				size_t binsNum = traceILow - nextILow;
				size_t BinSize = (*lowerData.estimateKmerCounts)[lowerData.location[traceJLow]];
				size_t coreectedSize = static_cast<size_t>(std::ceil(BinSize / (hicfConfig.loadFactor * binsNum)));
				lowerData.layout->userBins.emplace_back(lowerData.previous.binIndices, binIdLow, binsNum, lowerData.location[traceJLow]);
				if (coreectedSize > maxSizeLow)
				{
					maxSizeLow = coreectedSize;
					maxIdLow = binIdLow;
				}
				binIdLow += binsNum;
				traceILow = trace[traceILow][traceJLow];
				--traceJLow;
			}
			++traceILow;
			size_t BinSize = (*lowerData.estimateKmerCounts)[lowerData.location[0]];
			size_t coreectedSize = static_cast<size_t>(std::ceil(BinSize / (hicfConfig.loadFactor * traceILow)));
			lowerData.layout->userBins.emplace_back(lowerData.previous.binIndices, binIdLow, traceILow, lowerData.location[0]);
			if (coreectedSize > maxSizeLow)
			{
				maxSizeLow = coreectedSize;
				maxIdLow = binIdLow;
			}
			lowerMaxBin = maxIdLow;
		}
		data.layout->maxBins.emplace_back(lowerData.previous.binIndices, lowerMaxBin);
		if (config.mode != "fast")
		{
			auto Start = std::chrono::high_resolution_clock::now();
			currentCount = hll.estimate();
			auto End = std::chrono::high_resolution_clock::now();
			data.totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Start).count();
		}
		currentCount /= hicfConfig.relaxedLoadFactor;
		if (maxSize < currentCount)
		{
			maxSize = currentCount;
			maxId = binId;
		}
		data.totalTime += lowerData.totalTime;
	}


	/**
	 * Handles the split operation during the backtracking phase of the HICF layout computation.
	 * This function assigns a specified number of user-defined bins to a technical bin, updates the layout
	 * with the new bin assignment, and tracks the maximum bin size and its corresponding identifier.
	 *
	 * @param traceJ    The index of the current user bin in the `location` vector being processed.
	 * @param binNum    The number of user bins to assign to the current technical bin.
	 * @param binId     The identifier of the technical bin to which the user bins are being assigned.
	 * @param maxId     Reference to the variable tracking the identifier of the largest technical bin.
	 * @param maxSize   Reference to the variable tracking the size of the largest technical bin.
	 * @param data      Reference to the `DataStore` object containing layout information, estimated k-mer counts,
	 *                  load factors, and other relevant data structures.
	 *
	 * @details
	 * This function performs the following steps:
	 * 1. Bin Assignment: Adds a new entry to the `userBins` vector within the layout, specifying the indices of the
	 *    user bins being assigned, the technical bin identifier (`binId`), the number of user bins (`binNum`), and the
	 *    location index (`traceJ`).
	 *
	 * 2. K-mer Count Retrieval: Retrieves the estimated k-mer count for the current user bin using the `estimateKmerCounts`
	 *    vector and the `location` index.
	 *
	 * 3. Size Calculation: Calculates the corrected size of the technical bin by dividing the current k-mer count by the
	 *    product of the load factor and the number of user bins assigned. The result is rounded up using `std::ceil` to ensure
	 *    sufficient capacity.
	 *
	 * 4. Maximum Size Tracking: If the calculated `correctedSize` exceeds the current `maxSize`, the function updates `maxSize`
	 *    and `maxId` to reflect the new largest technical bin.
	 *
	 * @note
	 * - This function is typically called during the backtracking process to finalize the assignment of user bins to technical bins.
	 * - Proper synchronization is assumed to be handled externally if this function is called in a multithreaded context.
	 */
	void backtraceSplit(size_t traceJ, 
						size_t binNum,
						size_t binId, 
						size_t& maxId,
						size_t& maxSize,
						chimera::hicf::DataStore& data)
	{
		// Assign the specified number of user bins to the technical bin with identifier `binId`
		data.layout->userBins.emplace_back(data.previous.binIndices, binId, binNum, data.location[traceJ]);
		// Retrieve the estimated k-mer count for the current user bin
		size_t currentCount = (*data.estimateKmerCounts)[data.location[traceJ]];
		// Calculate the corrected size of the technical bin based on the load factor and number of user bins
		size_t corectedSize = static_cast<size_t>(std::ceil(currentCount / (data.loadFactor * binNum)));
		// Update the maximum bin size and its identifier if the current bin exceeds the existing maximum
		if (corectedSize > maxSize)
		{
			maxSize = corectedSize;
			maxId = binId;
		}
	}


	/**
	 * Compute the layout for Hierarchical Interleaved Cuckoo Filters (HICF) based on user-defined and technical bins.
	 * This function determines the optimal allocation of user bins into technical bins by minimizing the maximum
	 * technical bin size and the space consumption of lower-level Inverted Bloom Filters (IBF). It utilizes dynamic
	 * programming to evaluate different allocation strategies and performs backtracking to finalize the bin assignments.
	 *
	 * @param hicfConfig A reference to the HICFConfig object containing configuration parameters such as the number
	 *                   of user bins, maximum technical bins, load factors, and relaxed load factors.
	 * @param config A reference to the BuildConfig object that holds build settings including mode, alpha value,
	 *               and other relevant parameters.
	 * @param data A reference to the DataStore object that contains data structures like location indices,
	 *             estimated k-mer counts, HyperLogLog (HLL) vectors, union estimates, and timing information.
	 * @return The maximum bin ID (`maxId`) assigned during the layout computation, representing the highest
	 *         technical bin identifier used.
	 */
	size_t computeLayoutHICF(	HICFConfig& hicfConfig,
								BuildConfig& config,
								chimera::hicf::DataStore& data)
	{
		// Initialize the number of user and technical bins
		hicfConfig.userBinsNum = data.location.size();
		hicfConfig.technicalBinsNum = (data.previous.empty()) ? hicfConfig.technicalBinsMaxNum :
			std::min(hicfConfig.userBinsNum, hicfConfig.technicalBinsMaxNum);

		

		// matrix[i][j] represents the maximum technical bin size when allocating the first j+1 user bins into i+1 technical bins.
		// Formula: M[i][j] = min { max(M[k][j-1], ceil(Size(j) / (i - k)) ) | 0 ≤ k < i }
		std::vector<std::vector<size_t>> matrix(hicfConfig.technicalBinsNum, std::vector<size_t>(hicfConfig.userBinsNum, std::numeric_limits<size_t>::max()));
		// matrixL[i][j] represents the space consumption of the lower-level IBF required after allocating the first j+1 user bins into i+1 technical bins.
		// Formula: L[i][j] = min { L[k][j-1] + AdditionalSpace(j) | 0 ≤ k < i }
		std::vector<std::vector<size_t>> matrixL(hicfConfig.technicalBinsNum, std::vector<size_t>(hicfConfig.userBinsNum, 0u));
		// trace[i][j] is used to track the allocation path, recording each decision point (the indices of the previous technical bin and user bin).
		// Formula: T[i][j] = argmin { max(M[k][j-1], ceil(Size(j) / (i - k)) ) | 0 ≤ k < i }
		std::vector<std::vector<std::pair<size_t, size_t>>> trace(hicfConfig.technicalBinsNum, std::vector<std::pair<size_t, size_t>>(hicfConfig.userBinsNum,
			{ std::numeric_limits<size_t>::max(),std::numeric_limits<size_t>::max() }));

		// Initialization
		size_t userBin0 = (*data.estimateKmerCounts)[data.location[0]];
		for (size_t i = 0; i < hicfConfig.technicalBinsNum; i++)
		{
			double correctedValue = static_cast<double>(userBin0) / (hicfConfig.loadFactor * static_cast<double>(i + 1));
			matrix[i][0] = static_cast<size_t>(std::ceil(correctedValue));
			trace[i][0] = { 0u,0u };
		}

		size_t sum = userBin0;
		if (config.mode != "fast")
		{
			auto start = std::chrono::high_resolution_clock::now();
			toolbox::precompute_initial_union_estimates(data.estimates, *data.hllVec, *data.estimateKmerCounts, data.location);
			auto end = std::chrono::high_resolution_clock::now();
			data.totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			for (size_t j = 1; j < hicfConfig.userBinsNum; j++)
			{
				sum += (*data.estimateKmerCounts)[data.location[j]];
				double correctedValue = static_cast<double>(data.estimates[j]) / hicfConfig.relaxedLoadFactor;
				matrix[0][j] = static_cast<size_t>(std::ceil(correctedValue));
				matrixL[0][j] = maxMergeLevel(j + 1, hicfConfig.technicalBinsMaxNum) * sum * 16;
				trace[0][j] = { 0u,j - 1 };
			}
		}
		else
		{
			for (size_t j = 1; j < hicfConfig.userBinsNum; j++)
			{
				sum += (*data.estimateKmerCounts)[data.location[j]];
				double correctedValue = static_cast<double>(sum) / hicfConfig.relaxedLoadFactor;
				matrix[0][j] = static_cast<size_t>(std::ceil(correctedValue));
				matrixL[0][j] = maxMergeLevel(j + 1, hicfConfig.technicalBinsMaxNum) * sum * 16;
				trace[0][j] = { 0u,j - 1 };
			}
		}

		// Recursion
		for (size_t j = 1; j < hicfConfig.userBinsNum; j++)
		{
			size_t weight = (*data.estimateKmerCounts)[data.location[j]];

			if (config.mode != "fast")
			{
				auto start = std::chrono::high_resolution_clock::now();
				toolbox::precompute_union_estimates_for(data.estimates, *data.hllVec, *data.estimateKmerCounts, data.location, j);
				auto end = std::chrono::high_resolution_clock::now();
				data.totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			for (size_t i = 1; i < hicfConfig.technicalBinsNum; i++)
			{
				size_t minTechnicalBinSize = std::numeric_limits<size_t>::max();
				size_t minSpaceConsumption = std::numeric_limits<size_t>::max();

				for (size_t Split = 0; Split < i; Split++)
				{
					size_t correctedValue = static_cast<size_t>(std::ceil(static_cast<double>(weight) / hicfConfig.loadFactor));
					size_t score = std::max((correctedValue + i - Split) / (i - Split), matrix[Split][j - 1]);
					size_t fullScore = score * (i + 1) * 16 + config.alpha * matrixL[Split][j - 1];
					if (fullScore < minSpaceConsumption)
					{
						minSpaceConsumption = fullScore;
						minTechnicalBinSize = score;
						trace[i][j] = { Split,j - 1 };
						matrixL[i][j] = matrixL[Split][j - 1];
					}
				}

				size_t Split = j - 1;
				auto getWeight = [&]() -> size_t
					{
						size_t const uncorrected = (config.mode == "fast") ? weight : data.estimates[Split + 1];
						return  static_cast<double>(uncorrected) / hicfConfig.relaxedLoadFactor;
					};
				while (Split != 0 && ((i - trace[i][Split].first) < 2) && getWeight() < minTechnicalBinSize)
				{
					weight += (*data.estimateKmerCounts)[data.location[Split]];
					Split--;
					size_t score = std::max<size_t>(matrix[i - 1][Split], getWeight());
					size_t level = matrixL[i - 1][Split] + maxMergeLevel(j - Split, hicfConfig.technicalBinsMaxNum) * weight * 16;
					size_t fullScore = score * (i + 1) * 16 + config.alpha * level;
					if (fullScore < minSpaceConsumption)
					{
						minSpaceConsumption = fullScore;
						minTechnicalBinSize = score;
						trace[i][j] = { i - 1,Split };
						matrixL[i][j] = level;
					}
				}
				matrix[i][j] = minTechnicalBinSize;
			}
		}

		// Backtracking
		size_t traceI = hicfConfig.technicalBinsNum - 1;
		size_t traceJ = hicfConfig.userBinsNum - 1;
		size_t maxId = 0;
		size_t maxSize = 0;
		size_t binId = 0;
		while (traceJ > 0 && traceI > 0)
		{
			size_t nextI = trace[traceI][traceJ].first;
			size_t nextJ = trace[traceI][traceJ].second;
			size_t binNum = traceI - nextI;
			if (binNum == 1 && nextJ != traceJ - 1) // Merge
			{
				backtrackMerge(nextJ, traceJ, binId, maxId, maxSize, data, config, hicfConfig, false);
				traceI = nextI;
				traceJ = nextJ;
			}
			else // Split
			{
				backtraceSplit(traceJ, binNum, binId, maxId, maxSize, data);
				traceI = trace[traceI][traceJ].first;
				--traceJ;
			}
			binId += binNum;
		}

		// Check
		if (traceI == 0 && traceJ > 0)
		{
			backtrackMerge(0, traceJ, binId, maxId, maxSize, data, config, hicfConfig, true);
		}
		else if (traceJ == 0)
		{
			backtraceSplit(traceJ, traceI + 1, binId, maxId, maxSize, data);
		}

		return maxId;
	}


	void saveFilter(chimera::hicf::HierarchicalInterleavedCuckooFilter& hicf,
		const std::string& outputFile,
		std::vector<std::string>& indexToTaxid,
		HICFConfig& hicfConfig)
	{
		std::string outputDir = outputFile + ".hicf";
		std::ofstream os(outputDir, std::ios::binary);
		if (!os.is_open()) {
			throw std::runtime_error("Failed to open file: " + outputDir);
		}
		cereal::BinaryOutputArchive archive(os);
		archive(hicf);
		archive(indexToTaxid);
		archive(hicfConfig);
		os.close();
		// Get the file size
		std::uintmax_t fileSize = std::filesystem::file_size(outputDir);

		// Output the file size, formatted as KB, MB, or GB
		if (fileSize >= 1024 * 1024 * 1024) {
			std::cout << "Filter file size: " << std::fixed << std::setprecision(2)
				<< static_cast<double>(fileSize) / (1024 * 1024 * 1024) << " GB" << std::endl;
		}
		else if (fileSize >= 1024 * 1024) {
			std::cout << "Filter file size: " << std::fixed << std::setprecision(2)
				<< static_cast<double>(fileSize) / (1024 * 1024) << " MB" << std::endl;
		}
		else if (fileSize >= 1024) {
			std::cout << "Filter file size: " << std::fixed << std::setprecision(2)
				<< static_cast<double>(fileSize) / 1024 << " KB" << std::endl;
		}
		else {
			std::cout << "Filter file size: " << fileSize << " bytes" << std::endl;
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
		FileInfo fileInfo;
		if (config.window_size < config.kmer_size) {
			std::cerr << "Window size must be greater than or equal to kmer size." << std::endl;
			return;
		}
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

		std::filesystem::path dir = "tmp";
		createOrResetDirectory(dir, config);

		if (config.filter == "hicf")
		{
			std::cout << "Calculating minimizers and hyperloglog..." << std::endl;
			auto calculate_start = std::chrono::high_resolution_clock::now();

			std::vector<HyperLogLog> hllVec;
			hllVec.resize(hashCount.size());
			size_t hashCounts{};
			std::vector<std::string> indexToTaxid = minimiser_count(config, inputFiles, hashCounts, fileInfo, hllVec);

			auto calculate_end = std::chrono::high_resolution_clock::now();
			auto calculate_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(calculate_end - calculate_start).count();
			if (config.verbose) {
				std::cout << "Calculate time: ";
				print_build_time(calculate_total_time);
				std::cout << "File information:" << std::endl;
				std::cout << "Number of files: " << fileInfo.fileNum << std::endl;
				std::cout << "Number of invalid files: " << fileInfo.invalidNum << std::endl;
				std::cout << "Number of sequences: " << fileInfo.sequenceNum << std::endl;
				std::cout << "Number of skipped sequences: " << fileInfo.skippedNum << std::endl;
				std::cout << "Total base pairs: " << fileInfo.bpLength << std::endl << std::endl;
			}

			std::cout << "Rearranging user bins..." << std::endl;
			auto rearrange_start = std::chrono::high_resolution_clock::now();

			std::vector<size_t> estimateKmerCounts(hllVec.size());
			std::vector<size_t> location(hllVec.size());
			HICFConfig hicfConfig;
			hicfConfig.loadFactor = config.load_factor * 4.0;
			hicfConfig.relaxedLoadFactor = config.relaxedLoadFactor * 4.0;
			hicfConfig.userBinsNum = hllVec.size();
			hicfConfig.technicalBinsMaxNum = std::sqrt(hllVec.size());
			hicfConfig.windowSize = config.window_size;
			hicfConfig.kmerSize = config.kmer_size;
			hicfConfig.totalHashes = hashCounts;
			chimera::hicf::Layout layout;
			chimera::hicf::DataStore data;
			data.location = std::move(location);
			data.estimateKmerCounts = &estimateKmerCounts;
			data.hllVec = &hllVec;
			data.layout = &layout;
			data.loadFactor = hicfConfig.loadFactor;
			data.relaxedLoadFactor = hicfConfig.relaxedLoadFactor;
			estimateKmerCount(data);
			data.createLocation();
			rearrangeUserbins(data, config);

			auto rearrange_end = std::chrono::high_resolution_clock::now();
			auto rearrange_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(rearrange_end - rearrange_start).count();
			if (config.verbose) {
				std::cout << "Rearrange time: ";
				print_build_time(rearrange_total_time);
				std::cout << std::endl;
			}

			std::cout << "Start calculating HICF layout..." << std::endl;
			auto calculate_layout_start = std::chrono::high_resolution_clock::now();

			data.layout->topMaxBinId = computeLayoutHICF(hicfConfig, config, data);
			std::ranges::sort(data.layout->maxBins,
				[](auto const& l, auto const& r)
				{
					return l.previousTechnical.size() < r.previousTechnical.size();
				});
			hicfConfig.userBinsNum = hllVec.size();

			auto calculate_layout_end = std::chrono::high_resolution_clock::now();
			auto calculate_layout_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(calculate_layout_end - calculate_layout_start).count();
			if (config.verbose) {
				std::cout << "Union estimate time: ";
				print_build_time(data.totalTime);
				std::cout << std::endl;
				std::cout << "Calculate layout time: ";
				print_build_time(calculate_layout_total_time);
				std::cout << std::endl;
			}

			std::cout << "Building HICF..." << std::endl;
			auto build_start = std::chrono::high_resolution_clock::now();
			chimera::hicf::HierarchicalInterleavedCuckooFilter hicf(*data.layout, hicfConfig, indexToTaxid);
			auto build_end = std::chrono::high_resolution_clock::now();
			auto build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
			if (config.verbose) {
				std::cout << "Build time: ";
				print_build_time(build_total_time);
				std::cout << std::endl;
				std::cout << hicf << std::endl;
			}
			saveFilter(hicf, config.output_file, indexToTaxid, hicfConfig);

		}
		else if (config.filter == "icf")
		{
			auto calculate_start = std::chrono::high_resolution_clock::now();
			std::cout << "Calculating minimizers..." << std::endl;
			minimiser_count(config, inputFiles, hashCount, fileInfo);
			auto calculate_end = std::chrono::high_resolution_clock::now();
			auto calculate_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(calculate_end - calculate_start).count();
			if (config.verbose) {
				std::cout << "Calculate time: ";
				print_build_time(calculate_total_time);
				std::cout << "File information:" << std::endl;
				std::cout << "Number of files: " << fileInfo.fileNum << std::endl;
				std::cout << "Number of invalid files: " << fileInfo.invalidNum << std::endl;
				std::cout << "Number of sequences: " << fileInfo.sequenceNum << std::endl;
				std::cout << "Number of skipped sequences: " << fileInfo.skippedNum << std::endl;
				std::cout << "Total base pairs: " << fileInfo.bpLength << std::endl << std::endl;
			}
			ICFConfig icfConfig;
			icfConfig.kmer_size = config.kmer_size;
			icfConfig.window_size = config.window_size;
			if (config.mode == "normal")
			{
				icfConfig.bitNum = 16;
			}
			else if (config.mode == "fast")
			{
				icfConfig.bitNum = 8;
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
			chimera::InterleavedCuckooFilter icf(icfConfig.bins, icfConfig.bin_size, icfConfig.bitNum);
			auto calculate_bins_start = std::chrono::high_resolution_clock::now();
			robin_hood::unordered_map<std::string, std::size_t> taxidBins = calculateTaxidMapBins(icfConfig, hashCount);
			auto calculate_bins_end = std::chrono::high_resolution_clock::now();
			auto calculate_bins_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(calculate_bins_end - calculate_bins_start).count();
			if (config.verbose) {
				std::cout << "Calculated bins time: ";
				print_build_time(calculate_bins_total_time);
				std::cout << std::endl;
			}
			build(taxidBins, icfConfig, icf, hashCount, inputFiles);
			saveFilter(config.output_file, icf, icfConfig, hashCount, taxidBins);
			auto create_filter_end = std::chrono::high_resolution_clock::now();
			auto create_filter_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(create_filter_end - create_filter_start).count();
			if (config.verbose) {
				std::cout << "Create filter time: ";
				print_build_time(create_filter_total_time);
				std::cout << std::endl;
				std::cout << icf << std::endl;
			}
		}

		auto build_end = std::chrono::high_resolution_clock::now();
		// Calculate the total build time in milliseconds
		auto build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
		if (config.verbose) {
			std::cout << "Total build time: ";
			print_build_time(build_total_time);
		}
		std::filesystem::remove_all(dir);
	}
}