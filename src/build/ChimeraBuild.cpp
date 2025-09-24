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
 * Last Modified: 2024-11-18
 *
 * Description:
 *  The main program of ChimeraBuild.
 *
 * Version:
 *  1.2
 * -----------------------------------------------------------------------------
 */
#include <ChimeraBuild.hpp>
#include <limits>
#include <future>
#include <sstream>
#include <iomanip>
#include <string_view>
#include <chrono>

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
	void parseInputFile(const std::string& filePath, robin_hood::unordered_flat_map<std::string, std::vector<std::string>>& inputFiles, robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount, FileInfo& fileInfo) {
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
		robin_hood::unordered_flat_map<std::string, std::vector<std::string>>& inputFiles,
		robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount,
		FileInfo& fileInfo)
	{
		// Define the minimiser view
		auto minimiser_view = seqan3::views::minimiser_hash(
			seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
			seqan3::window_size{ config.window_size },
			seqan3::seed{ adjust_seed(config.kmer_size) });

		// Expand the pairs of taxid and files
		std::vector<std::pair<std::string, std::string>> taxid_file_pairs;
		robin_hood::unordered_flat_map<std::string, size_t> taxid_to_index;
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

			size_t taxid_index = taxid_to_index[taxid];

			bool skip_processing = false;
			{
				std::lock_guard<std::mutex> lock(taxidMutexes[taxid_index]);
				if (config.max_hashes_per_taxid > 0 && hashCount[taxid] >= config.max_hashes_per_taxid)
				{
					skip_processing = true;
				}
			}
			if (skip_processing)
			{
				{
					std::lock_guard<std::mutex> lock(fileInfo_mutex);
					fileInfo.skippedNum++;
				}
				continue;
			}

			// Thread-local file information
			FileInfo localFileInfo = {};

			// Thread-local hash count
			robin_hood::unordered_flat_map<uint64_t, uint8_t> local_hash_counts;

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

			{
				// Open the sequence file
				seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq >> fin{ filename };

				// Iterate over the records in the sequence file
				for (auto& record : fin)
				{
					auto& seq = record.sequence();

					if (seq.size() < config.min_length)
					{
						localFileInfo.skippedSeqNum++;
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
			}
			// Compute minimizers and count hash values
			std::vector<uint64_t> filtered_hashes;
			filtered_hashes.reserve(local_hash_counts.size());
			for (const auto& [hash, count] : local_hash_counts)
			{
				if (count >= cutoff)
				{
					filtered_hashes.emplace_back(hash);
				}
			}

			local_hash_counts.clear();

			if (max_hashes > 0 && filtered_hashes.size() > max_hashes)
			{
				filtered_hashes.resize(max_hashes);
			}

			// Write to the minimiser file
			{


				// Lock is required to prevent multiple threads from writing to the same taxid file
				std::lock_guard<std::mutex> lock(taxidMutexes[taxid_index]);


				if (max_hashes > 0)
				{
					size_t current_count = hashCount[taxid];
					if (current_count >= max_hashes)
					{
						filtered_hashes.clear();
					}
					else
					{
						size_t remaining = max_hashes - current_count;
						if (filtered_hashes.size() > remaining)
						{
							filtered_hashes.resize(remaining);
						}
						hashCount[taxid] += filtered_hashes.size();
					}
				}
				else
				{
					hashCount[taxid] += filtered_hashes.size();
				}
				if (!filtered_hashes.empty())
				{
					std::string output_filename = "tmp/" + taxid + ".mini";
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

			}


			// Update global hash count
			{
				std::lock_guard<std::mutex> lock(fileInfo_mutex);
				fileInfo.skippedNum += localFileInfo.skippedNum;
				fileInfo.skippedSeqNum += localFileInfo.skippedSeqNum;
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
		robin_hood::unordered_flat_map<std::string, std::vector<std::string>>& inputFiles,
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
		robin_hood::unordered_flat_map<std::string, size_t> taxid_to_index;
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

			robin_hood::unordered_flat_map<uint64_t, uint8_t> local_hash_counts;

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
					localFileInfo.skippedSeqNum++;
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
				fileInfo.skippedSeqNum += localFileInfo.skippedSeqNum;
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
	uint64_t getMaxValue(const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount) {
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
	uint64_t calculateTotalSize(const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount) {
		uint64_t totalSize = 0;
		for (const auto& kv : hashCount) {
			totalSize += kv.second;
		}
		return totalSize;
	}

	/**
	 * @brief Constructs the Interleaved Merged Cuckoo Filter (IMCF) using provided groups and returns a mapping of indices to taxids.
	 *
	 * This function builds the IMCF by iterating over groups of taxids, reading minimizer files, and inserting minimizer hashes
	 * into the IMCF. Each group represents a set of taxids that are processed together, with each minimizer hash being inserted
	 * into the corresponding index in the filter.
	 *
	 * @param imcf A reference to the Interleaved Merged Cuckoo Filter object to be built.
	 * @param groups A vector of Group objects, where each Group contains a collection of taxids.
	 *
	 * @return A vector of vectors, where each sub-vector maps the index in the IMCF to the corresponding taxid strings.
	 *
	 * @details
	 * The function first resizes `indexToTaxid` to match the number of groups and resizes each sub-vector to the number of taxids
	 * in the group. It then iterates over each group, reads the minimizer hashes from a file named in the format "tmp/{taxid}.mini",
	 * and inserts each hash into the IMCF using the `insertTag` method.
	 *
	 * The function is parallelized using OpenMP to process each group concurrently, speeding up the construction process.
	 *
	 * If the minimizer file for a taxid cannot be opened, an error message is printed to `std::cerr`.
	 *
	 * @note
	 * Ensure that the temporary files containing minimizer hashes ("tmp/{taxid}.mini") are pre-generated and accessible.
	 * The function assumes that `imcf.insertTag` is implemented and capable of inserting 64-bit minimizer hashes into the filter.
	 */
	std::vector<std::vector<std::string>> buildIMCF(
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		std::vector<chimera::imcf::Group>& groups)
	{
		// Initialize the index-to-taxid mapping structure with the number of groups
		std::vector<std::vector<std::string>> indexToTaxid;
		indexToTaxid.resize(groups.size());
		// Parallel loop to process each group concurrently
#pragma omp parallel for
		for (size_t i = 0; i < groups.size(); i++)
		{
			// Resize the sub-vector to hold the taxids for the current group
			indexToTaxid[i].resize(groups[i].taxids.size());
			// Iterate over each taxid in the current group
			for (size_t index = 0; index < groups[i].taxids.size(); index++)
			{
				indexToTaxid[i][index] = groups[i].taxids[index];
				auto taxid = groups[i].taxids[index];
				// Open the corresponding minimizer file for reading
				std::ifstream ifile("tmp/" + taxid + ".mini", std::ios::binary);
				if (ifile.fail()) {
					std::cerr << "Failed to open minimiser file: " << taxid << ".mini" << std::endl;
				}
				// Read minimizer hashes and insert them into the IMCF
				uint64_t hash;
				while (ifile.read(reinterpret_cast<char*>(&hash), sizeof(hash)))
				{
					imcf.insertTag(i, hash, index);
				}
				ifile.close();
			}
		}
		// Return the constructed index-to-taxid mapping
		return indexToTaxid;
	}

	void saveIMCF(chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		const std::string& output_file,
		std::vector<std::vector<std::string>>& indexToTaxid,
		IMCFConfig& imcfConfig)
	{
	// Open the output file
	std::ofstream os(output_file + ".imcf", std::ios::binary);

	// Check if the file is successfully opened
	if (!os.is_open()) {
		throw std::runtime_error("Failed to open file: " + output_file);
	}

	using Clock = std::chrono::steady_clock;
	auto formatDuration = [](long long ms) {
		if (ms >= 60000) {
			std::ostringstream oss;
			long long totalSeconds = ms / 1000;
			long long hours = totalSeconds / 3600;
			long long minutes = (totalSeconds % 3600) / 60;
			long long seconds = totalSeconds % 60;
			if (hours > 0) {
				oss << hours << "h ";
			}
			if (minutes > 0 || hours > 0) {
				oss << minutes << "min ";
			}
			oss << seconds << "s";
			return oss.str();
		}
		if (ms >= 1000) {
			std::ostringstream oss;
			oss << std::fixed << std::setprecision(2) << (ms / 1000.0) << "s";
			return oss.str();
		}
		return std::to_string(ms) + "ms";
	};

	auto logStep = [&](std::string_view title, auto&& work) {
		std::cout << "  - " << title << "... " << std::flush;
		auto start = Clock::now();
		bool ok = work();
		auto end = Clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << (ok ? "done" : "failed") << " (" << formatDuration(elapsed) << ")" << std::endl;
		return ok;
	};

	bool runAsyncIndexBuild = false;
#ifdef _OPENMP
	if (omp_get_max_threads() > 1) {
		runAsyncIndexBuild = true;
	}
#else
	if (std::thread::hardware_concurrency() > 1) {
		runAsyncIndexBuild = true;
	}
#endif

	std::future<void> routerFuture;
	Clock::time_point routerLaunchTime{};
	if (runAsyncIndexBuild) {
		routerLaunchTime = Clock::now();
		routerFuture = std::async(std::launch::async, [&imcf]() {
			imcf.buildRouterIndex();
		});
	}

	logStep("Writing filter archive (.imcf)", [&]() {
		cereal::BinaryOutputArchive archive(os);
		archive(imcf);
		archive(indexToTaxid);
		archive(imcfConfig);
		os.close();
		return true;
	});

	const char* indexTitle = runAsyncIndexBuild ? "Building index structures (.idx/.rtr) [async]" : "Building index structures (.idx/.rtr)";
	std::cout << "  - " << indexTitle << "... " << std::flush;
	long long indexElapsedMs = 0;
	if (routerFuture.valid()) {
		routerFuture.get();
		auto end = Clock::now();
		indexElapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(end - routerLaunchTime).count();
	}
	else {
		auto start = Clock::now();
		imcf.buildRouterIndex();
		auto end = Clock::now();
		indexElapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	}
	std::cout << "done (" << formatDuration(indexElapsedMs) << ")" << std::endl;

	std::string idxPath = output_file + ".imcf.idx";
	logStep("Writing active index (.imcf.idx)", [&]() {
		if (!imcf.saveActiveIndex(idxPath)) {
			std::cerr << "Warning: failed to write IMCF index file: " << idxPath << std::endl;
			return false;
		}
		return true;
	});

	std::string routerPath = output_file + ".imcf.rtr";
	logStep("Writing router index (.imcf.rtr)", [&]() {
		if (!imcf.saveRouterIndex(routerPath)) {
			std::cerr << "Warning: failed to write IMCF router file: " << routerPath << std::endl;
			return false;
		}
		return true;
	});

		// Get the file size
		std::uintmax_t fileSize = std::filesystem::file_size(output_file + ".imcf");

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
		if (config.threads == 0) {
			unsigned int hardwareThreads = std::thread::hardware_concurrency();
			if (hardwareThreads == 0) {
				hardwareThreads = 1;
			}
			const auto maxThreads = static_cast<unsigned int>(std::numeric_limits<uint16_t>::max());
			if (hardwareThreads > maxThreads) {
				hardwareThreads = maxThreads;
			}
			config.threads = static_cast<uint16_t>(hardwareThreads);
		}
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
		robin_hood::unordered_flat_map<std::string, uint64_t> hashCount;
		robin_hood::unordered_flat_map<std::string, std::vector<std::string>> inputFiles;
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
		if (config.filter == "imcf")
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
				std::cout << "Number of skipped files: " << fileInfo.skippedNum << std::endl;
				std::cout << "Number of sequences: " << fileInfo.sequenceNum << std::endl;
				std::cout << "Number of skipped sequences: " << fileInfo.skippedSeqNum << std::endl;
				std::cout << "Total base pairs: " << fileInfo.bpLength << std::endl << std::endl;
			}

			auto partition_start = std::chrono::high_resolution_clock::now();
			std::cout << "Partitioning hashcout..." << std::endl;
			std::vector<chimera::imcf::Group> groups = chimera::imcf::partitionHashCount(hashCount, 16);
			auto partition_end = std::chrono::high_resolution_clock::now();
			auto partition_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(partition_end - partition_start).count();
			if (config.verbose) {
				std::cout << "Partition time: ";
				print_build_time(partition_total_time);
				std::cout << std::endl;
			}

			auto build_start = std::chrono::high_resolution_clock::now();
			std::cout << "Building IMCF..." << std::endl;
			IMCFConfig imcfConfig;
			imcfConfig.loadFactor = config.load_factor;
			imcfConfig.kmerSize = config.kmer_size;
			imcfConfig.windowSize = config.window_size;
			chimera::imcf::InterleavedMergedCuckooFilter imcf(groups, imcfConfig);
			std::vector<std::vector<std::string>> indexToTaxid = buildIMCF(imcf, groups);
			auto build_end = std::chrono::high_resolution_clock::now();
			auto build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
			if (config.verbose) {
				std::cout << "Build time: ";
				print_build_time(build_total_time);
				std::cout << std::endl;
				std::cout << imcf << std::endl;
			}

			auto save_start = std::chrono::high_resolution_clock::now();
			std::cout << "Saving IMCF..." << std::endl;
			saveIMCF(imcf, config.output_file, indexToTaxid, imcfConfig);
			auto save_end = std::chrono::high_resolution_clock::now();
			auto save_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(save_end - save_start).count();
			if (config.verbose) {
				std::cout << "Save time: ";
				print_build_time(save_total_time);
				std::cout << std::endl;
			}
		}
		else {
			throw std::runtime_error("Unsupported filter type: " + config.filter +
			                             ". 当前仅支持 imcf");
		}



		auto build_end = std::chrono::high_resolution_clock::now();
		// Calculate the total build time in milliseconds
		auto build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
		if (config.verbose) {
			std::cout << "Total build time: ";
			print_build_time(build_total_time);
			std::cout << "Remove temporary files..." << std::endl;
		}
		std::filesystem::remove_all(dir);
	}
}
