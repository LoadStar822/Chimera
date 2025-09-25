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
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <string_view>
#include <chrono>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <mutex>
#include <span>

namespace ChimeraBuild {
	static inline uint64_t mix64(uint64_t value) {
		return XXH3_64bits(&value, sizeof(value));
	}

	struct TaxidShardPlan {
		size_t groupIndex;
		size_t slotIndex;
		uint64_t expectedCount;
	};

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
			seqan3::seed{ ChimeraBuild::adjust_seed(config.kmer_size) });

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
		const std::vector<chimera::imcf::Group>& groups,
		const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount)
	{
		std::vector<std::vector<std::string>> indexToTaxid(groups.size());
		robin_hood::unordered_flat_map<std::string, std::vector<TaxidShardPlan>> shardPlan;
		shardPlan.reserve(hashCount.size());

		for (size_t groupIdx = 0; groupIdx < groups.size(); ++groupIdx) {
			const auto& group = groups[groupIdx];
			indexToTaxid[groupIdx] = group.taxids;
			if (group.taxids.size() != group.assignedHashes.size()) {
				throw std::runtime_error("IMCF build: taxid list and assigned hash quotas length mismatch");
			}
			for (size_t slot = 0; slot < group.taxids.size(); ++slot) {
				uint64_t expected = group.assignedHashes[slot];
				if (expected == 0) {
					continue;
				}
				auto& plans = shardPlan[group.taxids[slot]];
				plans.push_back({ groupIdx, slot, expected });
			}
		}

		const double kShardToleranceRatio = 0.001; // 0.1%
		for (const auto& [taxid, plans] : shardPlan) {
			uint64_t shardTotal = 0;
			for (const auto& plan : plans) {
				shardTotal += plan.expectedCount;
			}
			auto it = hashCount.find(taxid);
			if (it == hashCount.end()) {
				std::cerr << "Warning: shard plan references unknown taxid " << taxid << std::endl;
				continue;
			}
			uint64_t expected = it->second;
			if (expected == 0 && shardTotal == 0) {
				continue;
			}
			int64_t delta = static_cast<int64_t>(shardTotal) - static_cast<int64_t>(expected);
			uint64_t absDelta = delta >= 0 ? static_cast<uint64_t>(delta) : static_cast<uint64_t>(-delta);
			double ratio = expected == 0 ? (delta == 0 ? 0.0 : std::numeric_limits<double>::infinity())
				: static_cast<double>(absDelta) / static_cast<double>(expected);
				if (ratio > kShardToleranceRatio) {
					std::cerr << "Warning: minimiser count mismatch for taxid " << taxid
						<< ": expected=" << expected << ", planned=" << shardTotal
						<< ", delta=" << delta << " (" << std::fixed << std::setprecision(4)
						<< ratio * 100.0 << "%)" << std::endl;
				}
		}

		struct ShardEntry {
			std::string taxid;
			const std::vector<TaxidShardPlan>* plans;
		};
		std::vector<ShardEntry> shardEntries;
		shardEntries.reserve(shardPlan.size());
		for (auto& kv : shardPlan) {
			shardEntries.push_back({ kv.first, &kv.second });
		}

		std::vector<std::mutex> groupLocks(groups.size());

		std::vector<uint64_t> groupSuccess(groups.size(), 0);
		std::vector<uint64_t> groupFailures(groups.size(), 0);
		std::vector<std::vector<uint64_t>> slotSuccess(groups.size());
		std::vector<std::vector<uint64_t>> slotFailures(groups.size());
		for (size_t g = 0; g < groups.size(); ++g) {
			size_t slotCount = groups[g].taxids.size();
			slotSuccess[g].assign(slotCount, 0);
			slotFailures[g].assign(slotCount, 0);
		}

		constexpr size_t kFlushThreshold = 16384;
		constexpr size_t kReadBlockBytes = 512 * 1024;

		const auto shardEntriesView = std::span<const ShardEntry>(shardEntries.data(), shardEntries.size());

#pragma omp parallel for schedule(dynamic)
		for (size_t entryIdx = 0; entryIdx < shardEntriesView.size(); ++entryIdx) {
			const auto& entry = shardEntriesView[entryIdx];
			const std::string& taxid = entry.taxid;
			const auto& plans = *entry.plans;
			std::ifstream ifile("tmp/" + taxid + ".mini", std::ios::binary);
			if (!ifile) {
#pragma omp critical(imcf_log)
				std::cerr << "Failed to open minimiser file: " << taxid << ".mini" << std::endl;
				continue;
			}

			std::vector<uint64_t> remainingCounts(plans.size());
			for (size_t i = 0; i < plans.size(); ++i) {
				remainingCounts[i] = plans[i].expectedCount;
			}
			std::vector<std::vector<uint64_t>> slotBuffers(plans.size());
			for (size_t i = 0; i < plans.size(); ++i) {
				size_t reserveSize = static_cast<size_t>(std::min<uint64_t>(remainingCounts[i], static_cast<uint64_t>(kFlushThreshold)));
				if (reserveSize > 0) {
					slotBuffers[i].reserve(reserveSize);
				}
			}

			auto flushSlot = [&](size_t planIdx) {
				auto& buffer = slotBuffers[planIdx];
				if (buffer.empty()) {
					return;
				}
				const auto& plan = plans[planIdx];
				std::lock_guard<std::mutex> guard(groupLocks[plan.groupIndex]);
				uint64_t success = 0;
				uint64_t failure = 0;
				for (uint64_t value : buffer) {
					if (imcf.insertTag(plan.groupIndex, value, plan.slotIndex)) {
						++success;
					}
					else {
						++failure;
					}
				}
				groupSuccess[plan.groupIndex] += success;
				groupFailures[plan.groupIndex] += failure;
				slotSuccess[plan.groupIndex][plan.slotIndex] += success;
				slotFailures[plan.groupIndex][plan.slotIndex] += failure;
				buffer.clear();
			};

			uint64_t totalRemaining = 0;
			for (uint64_t v : remainingCounts) {
				totalRemaining += v;
			}

			auto pickShard = [&](uint64_t hash) -> size_t {
				if (totalRemaining == 0) {
					return plans.size();
				}
				uint64_t r = mix64(hash) % totalRemaining;
				for (size_t idx = 0; idx < remainingCounts.size(); ++idx) {
					uint64_t quota = remainingCounts[idx];
					if (quota == 0) {
						continue;
					}
					if (r < quota) {
						return idx;
					}
					r -= quota;
				}
				return plans.size();
			};

			std::vector<uint64_t> readBlock(kReadBlockBytes / sizeof(uint64_t));
			while (ifile) {
				ifile.read(reinterpret_cast<char*>(readBlock.data()), static_cast<std::streamsize>(readBlock.size() * sizeof(uint64_t)));
				std::streamsize bytesRead = ifile.gcount();
				if (bytesRead <= 0) {
					break;
				}
				size_t hashesRead = static_cast<size_t>(bytesRead) / sizeof(uint64_t);
				bool overflow = false;
				for (size_t h = 0; h < hashesRead; ++h) {
					size_t shard = pickShard(readBlock[h]);
					if (shard >= plans.size()) {
						// Extra minimisers beyond planned quota
	#pragma omp critical(imcf_log)
						std::cerr << "Warning: taxid " << taxid
							<< " produced more minimisers than expected; extra entries ignored" << std::endl;
						overflow = true;
						break;
					}
					slotBuffers[shard].push_back(readBlock[h]);
					if (slotBuffers[shard].size() >= kFlushThreshold) {
						flushSlot(shard);
					}
					if (remainingCounts[shard] > 0) {
						--remainingCounts[shard];
						--totalRemaining;
					}
					if (remainingCounts[shard] == 0) {
						flushSlot(shard);
					}
				}
				if (overflow) {
					break;
				}
			}

			if (!ifile.eof() && ifile.fail()) {
	#pragma omp critical(imcf_log)
				std::cerr << "Warning: failed to read minimiser file completely for taxid " << taxid << std::endl;
			}
			ifile.close();

			for (size_t shardIdx = 0; shardIdx < plans.size(); ++shardIdx) {
				flushSlot(shardIdx);
			}

			uint64_t missing = totalRemaining;
			if (missing != 0) {
#pragma omp critical(imcf_log)
				std::cerr << "Warning: taxid " << taxid
					<< " ended with " << missing << " minimisers short of expected quota" << std::endl;
			}
		}

		uint64_t totalInserted = 0;
		uint64_t totalFailed = 0;
		for (size_t i = 0; i < groups.size(); ++i) {
			totalInserted += groupSuccess[i];
			totalFailed += groupFailures[i];
		}
		uint64_t totalAttempts = totalInserted + totalFailed;
		if (totalAttempts > 0) {
			auto formatRate = [](double rate) {
				std::ostringstream oss;
				oss << std::fixed << std::setprecision(4) << rate * 100.0 << '%';
				return oss.str();
			};
				double failureRate = static_cast<double>(totalFailed) / static_cast<double>(totalAttempts);
				std::cout << "  - IMCF insertion result: " << totalFailed << "/" << totalAttempts
					<< " failed (" << formatRate(failureRate) << ")" << std::endl;
				constexpr double kFailureWarningThreshold = 0.001; // 0.1%
				constexpr uint64_t kFailureMinPrint = 100;
				size_t groupsWithFailures = 0;
				size_t groupsAboveWarn = 0;
				double maxGroupRate = 0.0;
				size_t maxGroupIndex = std::numeric_limits<size_t>::max();
				struct GroupSummary {
					size_t index;
					uint64_t attempts;
					uint64_t failures;
					double rate;
					size_t worstSlot;
					uint64_t worstSlotAttempts;
					uint64_t worstSlotFailures;
					double worstSlotRate;
				};
				std::vector<GroupSummary> summaries;
				summaries.reserve(groups.size());

				for (size_t i = 0; i < groups.size(); ++i) {
					uint64_t attempts = groupSuccess[i] + groupFailures[i];
					if (attempts == 0) {
						continue;
					}
					double groupFailureRate = static_cast<double>(groupFailures[i]) / static_cast<double>(attempts);
					if (groupFailures[i] > 0) {
						++groupsWithFailures;
					}
					if (groupFailureRate > kFailureWarningThreshold) {
						++groupsAboveWarn;
					}
					if (groupFailureRate > maxGroupRate) {
						maxGroupRate = groupFailureRate;
						maxGroupIndex = i;
					}
					if (groupFailures[i] == 0) {
						continue;
					}
					if (groupFailureRate <= kFailureWarningThreshold && groupFailures[i] < kFailureMinPrint) {
						continue;
					}

					size_t worstSlot = std::numeric_limits<size_t>::max();
					double worstSlotRate = 0.0;
					uint64_t worstSlotAttempts = 0;
					uint64_t worstSlotFailures = 0;
					for (size_t slot = 0; slot < slotSuccess[i].size(); ++slot) {
						uint64_t slotAttempts = slotSuccess[i][slot] + slotFailures[i][slot];
						if (slotAttempts == 0) {
							continue;
						}
						double slotRate = static_cast<double>(slotFailures[i][slot]) / static_cast<double>(slotAttempts);
						if (slotRate > worstSlotRate) {
							worstSlotRate = slotRate;
							worstSlot = slot;
							worstSlotAttempts = slotAttempts;
							worstSlotFailures = slotFailures[i][slot];
						}
					}

					GroupSummary summary{ i, attempts, groupFailures[i], groupFailureRate, worstSlot, worstSlotAttempts, worstSlotFailures, worstSlotRate };
					summaries.push_back(summary);
				}

				std::cout << "    groups with failures: " << groupsWithFailures << "/" << groups.size();
				if (maxGroupIndex != std::numeric_limits<size_t>::max()) {
					std::cout << ", worst group=" << maxGroupIndex << " (" << formatRate(maxGroupRate) << ")";
				}
				std::cout << std::endl;
				if (groupsAboveWarn > 0) {
					std::cout << "    groups above warning threshold (" << formatRate(kFailureWarningThreshold) << "): "
						<< groupsAboveWarn << std::endl;
				}

				std::sort(summaries.begin(), summaries.end(), [](const GroupSummary& a, const GroupSummary& b) {
					return a.rate > b.rate;
				});

				const size_t maxDetails = 5;
				for (size_t idx = 0; idx < summaries.size() && idx < maxDetails; ++idx) {
					const auto& summary = summaries[idx];
					std::ostringstream line;
					line << "    group " << summary.index << ": " << summary.failures << "/" << summary.attempts
						 << " failed (" << formatRate(summary.rate) << ")";
					if (summary.rate > kFailureWarningThreshold) {
						line << " — 建议降低 load_factor 或提高 MaxCuckooCount/扩大 binSize";
					}
					std::cout << line.str() << std::endl;
					if (summary.worstSlot != std::numeric_limits<size_t>::max() && summary.worstSlot < indexToTaxid[summary.index].size() && summary.worstSlotFailures > 0 && summary.worstSlotRate > kFailureWarningThreshold) {
						std::cout << "      worst slot " << summary.worstSlot << " (taxid="
							<< indexToTaxid[summary.index][summary.worstSlot] << ") : "
							<< summary.worstSlotFailures << "/" << summary.worstSlotAttempts
							<< " failed (" << formatRate(summary.worstSlotRate) << ")" << std::endl;
					}
					auto bucketStats = imcf.computeBucketStats(summary.index, 256);
					std::ostringstream statsLine;
					statsLine << std::fixed << std::setprecision(4)
						<< "      bucket load mean=" << bucketStats.meanLoad
						<< ", stdev=" << bucketStats.stddevLoad
						<< ", full=" << bucketStats.occupancy[4] << "/" << bucketStats.bucketCount
						<< " (" << formatRate(bucketStats.percentFull) << ")";
					statsLine << std::setprecision(3)
						<< ", low-bit load ratio range=" << bucketStats.lowBitMinRatio << " - " << bucketStats.lowBitMaxRatio;
					if (bucketStats.stashEntries > 0) {
						statsLine << std::setprecision(0)
							<< ", stash_species=" << bucketStats.stashEntries
							<< " (entries=" << bucketStats.stashEntryCount
							<< ", max_per_entry=" << bucketStats.stashMaxPerEntry << ")";
					}
					std::cout << statsLine.str() << std::endl;
				}
			uint64_t imcfFailureTotal = imcf.getInsertFailureTotal();
			uint64_t imcfFailureSaturated = imcf.getInsertFailureSaturatedFingerprint();
			if (imcfFailureTotal > 0) {
				double saturatedRatio = static_cast<double>(imcfFailureSaturated) / static_cast<double>(imcfFailureTotal);
				std::cout << "  - IMCF failure counters: total=" << imcfFailureTotal
					<< ", saturated_fp=" << imcfFailureSaturated << " ("
					<< formatRate(saturatedRatio) << ")" << std::endl;
			}
		}

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
			imcfConfig.seed64 = ChimeraBuild::adjust_seed(config.kmer_size);
			imcfConfig.fpSalt = IMCFConfig::DefaultFingerprintSalt;
			imcfConfig.hashVersion = IMCFConfig::CurrentHashVersion;
			chimera::imcf::InterleavedMergedCuckooFilter imcf(groups, imcfConfig);
			std::vector<std::vector<std::string>> indexToTaxid = buildIMCF(imcf, groups, hashCount);
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
