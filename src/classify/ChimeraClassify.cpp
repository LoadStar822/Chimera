/*
 * -----------------------------------------------------------------------------
 * Filename:      ChimeraClassify.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-08-09
 *
 * Last Modified: 2024-11-18
 *
 * Description:
 *  Classify functions for Chimera
 *
 * Version:
 *  1.4
 * -----------------------------------------------------------------------------
 */
#include <ChimeraClassify.hpp>

namespace ChimeraClassify {
	/**
	 * @brief Print the time taken for classification in a human-readable format.
	 *
	 * This function prints the time taken for classification in a human-readable format.
	 * It calculates the seconds, minutes, and hours from the milliseconds.
	 * The function outputs different formats based on the length of time.
	 *
	 * @param milliseconds The time taken for classification in milliseconds.
	 */
	void print_classify_time(long long milliseconds) {
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
	 * @brief Parse the reads from input files and store them in a queue.
	 *
	 * This function reads the sequences from input files and creates batchReads objects.
	 * It processes single files in parallel or paired files in parallel based on the configuration.
	 * The sequences are read in chunks and stored in batchReads objects.
	 * The batchReads objects are then pushed to the global queue for further processing.
	 * The number of sequences for each thread is counted and accumulated to the fileInfo.sequenceNum.
	 *
	 * @param readQueue The queue to store the batchReads objects.
	 * @param config The configuration for parsing the reads.
	 * @param fileInfo The information about the files and sequences.
	 */
	void parseReads(moodycamel::ConcurrentQueue<batchReads>& readQueue, ClassifyConfig config, FileInfo& fileInfo)
	{
		// Count the number of sequences for each thread and accumulate it to fileInfo.sequenceNum
		size_t totalSequenceNum = 0;

		if (!config.singleFiles.empty()) {
#pragma omp parallel for reduction(+:totalSequenceNum)
			for (size_t i = 0; i < config.singleFiles.size(); ++i) {
				const auto& file = config.singleFiles[i];

				fileInfo.fileNum++;

				seqan3::sequence_file_input<raptor::dna4_traits,
					seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{ file };
				std::vector<batchReads> localBatches;

				for (auto&& rec : fin | seqan3::views::chunk(config.batchSize)) {
					batchReads batch;
					for (auto&& r : rec) {
						batch.ids.emplace_back(std::move(r.id()));
						batch.seqs.emplace_back(std::move(r.sequence()));
					}

					totalSequenceNum += batch.ids.size();
					localBatches.emplace_back(std::move(batch));
				}

				for (auto& b : localBatches) {
					readQueue.enqueue(std::move(b));
				}
				fileInfo.sequenceNum += totalSequenceNum;
			}
		}
		else if (config.pairedFiles.size() > 0)
		{
#pragma omp parallel for reduction(+:totalSequenceNum)
			for (size_t i = 0; i < config.pairedFiles.size(); i += 2) {
				fileInfo.fileNum += 2;

				seqan3::sequence_file_input<raptor::dna4_traits,
					seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin1{ config.pairedFiles[i] };
				seqan3::sequence_file_input<raptor::dna4_traits,
					seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin2{ config.pairedFiles[i + 1] };
				std::vector<batchReads> localBatches;

				auto it1 = fin1 | seqan3::views::chunk(config.batchSize);
				auto it2 = fin2 | seqan3::views::chunk(config.batchSize);

				while (true) {
					auto chunk1 = it1.begin();
					auto end1 = it1.end();
					if (chunk1 == end1) break;

					auto chunk2 = it2.begin();
					auto end2 = it2.end();
					if (chunk2 == end2) break;

					batchReads batch;
					for (auto&& r : *chunk1) {
						batch.ids.emplace_back(std::move(r.id()));
						batch.seqs.emplace_back(std::move(r.sequence()));
					}
					for (auto&& r : *chunk2) {
						batch.seqs2.emplace_back(std::move(r.sequence()));
					}

					totalSequenceNum += batch.ids.size();
					localBatches.emplace_back(std::move(batch));
				}

				for (auto& b : localBatches) {
					readQueue.enqueue(std::move(b));
				}
			}

			fileInfo.sequenceNum += totalSequenceNum;
		}
		else
		{
			throw std::runtime_error("No input files specified");
		}

		if (config.verbose) {
			std::cout << "Number of files: " << fileInfo.fileNum << std::endl;
			std::cout << "Number of sequences: " << fileInfo.sequenceNum << std::endl << std::endl;
		}
	}

	/**
	 * @brief Load the filter from the input file.
	 *
	 * This function opens the input file and deserializes the filter and related data.
	 * It checks if the file is successfully opened and throws an exception if it fails.
	 * The function uses a cereal binary input archive to deserialize the filter and data.
	 * After deserialization, the function clears the hashCount and taxidBins containers.
	 * It then reinserts the data from the deserialized vectors into the hashCount and taxidBins containers.
	 * Finally, the function closes the file.
	 *
	 * @param input_file The path to the input file.
	 * @param icf The InterleavedCuckooFilter object to store the deserialized filter.
	 * @param icfConfig The ICFConfig object to store the deserialized configuration.
	 * @param hashCount The unordered_map to store the deserialized hash count data.
	 * @param taxidBins The unordered_map to store the deserialized taxid bins data.
	 */
	void loadFilter(const std::string& input_file,
		chimera::InterleavedCuckooFilter& icf,
		ChimeraBuild::ICFConfig& icfConfig,
		robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount,
		std::vector<std::pair<std::string, std::size_t>>& taxidBinsData
	) {
		// Open the input file
		std::ifstream is(input_file, std::ios::binary);

		// Check if the file is successfully opened
		if (!is.is_open()) {
			throw std::runtime_error("Failed to open file: " + input_file);
		}

		// Create a cereal binary input archive
		cereal::BinaryInputArchive archive(is);

		// Deserialize InterleavedCuckooFilter
		archive(icf);

		// Deserialize ICFConfig
		archive(icfConfig);

		// Deserialize vectors
		std::vector<std::pair<std::string, uint64_t>> hashCountData;
		archive(hashCountData);
		archive(taxidBinsData);

		// Clear hashCount and taxidBins
		hashCount.clear();

		// Reinsert data from vectors into robin_hood::unordered_flat_map
		for (const auto& kv : hashCountData) {
			hashCount[kv.first] = kv.second;
		}

		// Close the file
		is.close();
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
	 * @brief Process a sequence for classification.
	 *
	 * This function processes a sequence for classification using the Interleaved Cuckoo Filter (ICF).
	 * It takes the sequence's hash values, the ICF configuration, the taxid bins, the classification configuration,
	 * the ICF object, the sequence's ID, and the vector to store the classification results.
	 *
	 * The function calculates the number of hash values and the threshold for classification based on the configuration.
	 * It then performs a bulk count operation on the hash values using the ICF object.
	 * The function creates a classifyResult object to store the classification result for the sequence.
	 *
	 * The function iterates over the taxid bins and calculates the count for each bin based on the bulk count result.
	 * If the count for a bin is above the threshold, the bin and count are added to the classifyResult object.
	 * The function also keeps track of the bin with the maximum count.
	 *
	 * If the maximum count is greater than 0, the function checks the classification mode.
	 * If the mode is "fast" or there is only one bin with a count above the threshold, the classifyResult object is updated
	 * to only contain the bin with the maximum count.
	 *
	 * Finally, the function locks the resultMutex to ensure thread safety and adds the classifyResult object to the
	 * classifyResults vector.
	 *
	 * @param hashs1 The hash values of the sequence.
	 * @param icfConfig The configuration for the Interleaved Cuckoo Filter.
	 * @param taxidBins The taxid bins.
	 * @param config The configuration for the classification.
	 * @param icf The Interleaved Cuckoo Filter object.
	 * @param id The ID of the sequence.
	 * @param classifyResults The vector to store the classification results.
	 * @param fileInfo The information about the files and sequences.
	 */
	inline void processSequence(const std::vector<size_t>& hashs1,
		ChimeraBuild::ICFConfig& icfConfig,
		const std::vector<std::pair<std::string, std::size_t>>& taxidBins,
		ClassifyConfig& config,
		chimera::InterleavedCuckooFilter& icf,
		const std::string& id,
		std::vector<classifyResult>& classifyResults,
		FileInfo& fileInfo,
		LCA& lca)
	{
		// Calculate the number of hash values and the threshold for classification
		size_t hashNum = hashs1.size();
		size_t threshold = std::ceil(hashNum * config.shotThreshold);
		if (threshold == 0)
		{
			threshold = 1;
		}

		// Perform bulk count operation on the hash values using the Interleaved Cuckoo Filter
		auto count = icf.bulk_count(hashs1);
		classifyResult result;
		result.id = id;

		size_t oldIndex = 0;
		size_t maxBinCount = 0;
		std::pair<std::string, std::size_t> maxCount;
		robin_hood::unordered_flat_map<std::string, size_t> taxidToCount;

		// Iterate over the taxid bins and calculate the count for each bin based on the bulk count result
		for (auto const& [taxid, bins] : taxidBins)
		{
			size_t binCount = 0;
			for (size_t i = oldIndex; i < bins; i++)
			{
				binCount += kv_A(count, i);
			}
			oldIndex = bins;
			if (binCount > hashNum)
			{
				binCount = hashNum;
			}
			taxidToCount[taxid] = binCount;
			if (binCount > maxBinCount)
			{
				maxBinCount = binCount;
			}
		}

		// First filter step: remove bins with count less than 80% of the maximum count
		size_t thresholdCount = static_cast<size_t>(0.8 * static_cast<double>(maxBinCount));
		for (const auto& [taxid, binCount] : taxidToCount)
		{
			if (binCount >= thresholdCount && binCount >= threshold)
			{
				result.taxidCount.emplace_back(taxid, binCount);

				if (binCount == maxBinCount)
				{
					maxCount = std::make_pair(taxid, binCount);
				}
				fileInfo.taxidTotalMatches[taxid] += 1;
			}
		}

		// If the maximum count is greater than 0, update the classifyResult object based on the classification mode
		if (!result.taxidCount.empty())
		{

			bool isUniqueMapping = (result.taxidCount.size() == 1);
			if (isUniqueMapping)
			{
				fileInfo.uniqueTaxids.insert(result.taxidCount.front().first);
				const std::string& taxid = result.taxidCount.front().first;
				fileInfo.taxidUniqueMatches[taxid] += 1;
			}

			fileInfo.classifiedNum++;
			if (config.lca && result.taxidCount.size() > 1)
			{
				std::vector<std::string> taxids;
				for (auto& [taxid, count] : result.taxidCount)
				{
					taxids.push_back(taxid);
				}
				std::string lcaTaxid = lca.getLCA(taxids);
				result.taxidCount.clear();
				result.taxidCount.emplace_back(lcaTaxid, 0);
			}
			else if (config.mode == "fast" || result.taxidCount.size() == 1)
			{
				result.taxidCount.clear();
				result.taxidCount.emplace_back(maxCount);
			}
			else if (config.mode == "normal")
			{
				//// Sort the taxidCount vector based on the count in descending order
				//std::sort(result.taxidCount.begin(), result.taxidCount.end(), [](const auto& a, const auto& b) {
				//	return a.second > b.second;  // Sort in descending order
			}
		}
		else
		{
			fileInfo.unclassifiedNum++;
			result.taxidCount.emplace_back("unclassified", 1);
		}

		// Lock the resultMutex to ensure thread safety and add the classifyResult object to the classifyResults vector
		kv_destroy(count);
		classifyResults.emplace_back(std::move(result));
	}

	/**
	 * @brief Process a batch of reads for classification.
	 *
	 * This function processes a batch of reads for classification using the Interleaved Cuckoo Filter (ICF).
	 * It takes the batch of reads, the ICF configuration, the taxid bins, the classification configuration,
	 * the ICF object, the vector to store the classification results, and the minimiser view for hashing.
	 *
	 * If the batch contains paired-end reads, the function generates minimizer hash values for each sequence.
	 * If the sequence lengths meet the window size requirement, the function generates the hash values using the minimiser view.
	 * For paired-end reads, the function combines the hash values from both sequences.
	 *
	 * If the batch contains single-end reads, the function generates minimizer hash values for each sequence.
	 * If the sequence length meets the window size requirement, the function generates the hash values using the minimiser view.
	 *
	 * The function then calls the processSequence function to perform classification on the generated hash values.
	 *
	 * @param batch The batch of reads to process.
	 * @param icfConfig The configuration for the Interleaved Cuckoo Filter.
	 * @param taxidBins The taxid bins.
	 * @param config The configuration for the classification.
	 * @param icf The Interleaved Cuckoo Filter object.
	 * @param classifyResults The vector to store the classification results.
	 * @param minimiser_view The minimiser view for hashing.
	 * @param fileInfo The information about the files and sequences.
	 */
	inline void processBatch(batchReads batch,
		ChimeraBuild::ICFConfig& icfConfig,
		const std::vector<std::pair<std::string, std::size_t>>& taxidBins,
		ClassifyConfig& config,
		chimera::InterleavedCuckooFilter& icf,
		std::vector<classifyResult>& classifyResults,
		const auto& minimiser_view,
		FileInfo& fileInfo,
		LCA& lca)
	{
		// Process batch of reads
		std::vector<size_t> hashs1;
		if (!batch.seqs2.empty())
		{
			// Process paired-end reads
			for (size_t i = 0; i < batch.seqs2.size(); i += 2)
			{
				hashs1.clear();
				if (batch.seqs2[i].size() >= icfConfig.window_size)
				{
					// Generate minimizer hash values for the first sequence
					hashs1 = batch.seqs2[i] | minimiser_view | seqan3::ranges::to<std::vector>();
					if (batch.seqs2[i + 1].size() >= icfConfig.window_size)
					{
						// Generate minimizer hash values for the second sequence
						std::vector<size_t> hashs2 = batch.seqs2[i + 1] | minimiser_view | seqan3::ranges::to<std::vector>();
						// Combine the hash values from both sequences
						hashs1.insert(hashs1.end(), hashs2.begin(), hashs2.end());
					}
				}
				// Process the combined hash values for classification
				processSequence(hashs1, icfConfig, taxidBins, config, icf, batch.ids[i >> 1], classifyResults, fileInfo, lca);
			}
		}
		else
		{
			// Process single-end reads
			for (size_t i = 0; i < batch.seqs.size(); i++)
			{
				hashs1.clear();
				if (batch.seqs[i].size() >= icfConfig.window_size)
				{
					// Generate minimizer hash values for the sequence
					hashs1 = batch.seqs[i] | minimiser_view | seqan3::ranges::to<std::vector>();
				}
				// Process the hash values for classification
				processSequence(hashs1, icfConfig, taxidBins, config, icf, batch.ids[i], classifyResults, fileInfo, lca);
			}
		}
	}

	/*
	 * @brief Build the LCA object from the taxonomy file.
	 *
	 * This function builds the LCA object from the taxonomy file.
	 * It reads the taxonomy file line by line and splits each line into taxid and parent taxid.
	 * The function then adds the edge to the LCA object.
	 * After reading the entire file, the function performs Euler walk on the LCA object.
	 *
	 * @param lca The LCA object to build.
	 * @param taxFile The path to the taxonomy file.
	*/
	void buildLCA(LCA& lca, const std::string& taxFile)
	{
		// Open the taxonomy file
		std::ifstream is(taxFile);

		// Check if the file is successfully opened
		if (!is.is_open())
		{
			throw std::runtime_error("Failed to open file: " + taxFile);
		}

		// Read the taxonomy file line by line
		std::string line;
		while (std::getline(is, line))
		{
			// Split the line into taxid and parent taxid
			std::istringstream iss(line);
			std::string childID, parentID, rank, name;
			// Read the taxid, parent taxid, rank, and name from the line
			if (!(iss >> childID >> parentID >> rank >> name)) {
				continue;
			}

			// Add the edge to the LCA object
			lca.addEdge(parentID, childID);
		}

		// Close the file
		is.close();

		// Perform Euler walk on the LCA object
		lca.doEulerWalk("1");
	}


	/**
	 * @brief Applies a second filtering step to the classification results based on unique taxids.
	 *
	 * This function filters out taxids from the classification results that are not present in the provided set of unique taxids.
	 * If all taxids in a result are removed during this filtering, the result is marked as "unclassified".
	 *
	 * @param classifyResults A vector containing the classification results for each sequence.
	 * @param uniqueTaxids An unordered set of unique taxids that should be retained in the classification results.
	 *
	 * @details
	 * The function iterates over each `classifyResult` in `classifyResults` and applies the following logic:
	 * - It removes any taxid from `result.taxidCount` if it is not found in `uniqueTaxids`.
	 * - If the `taxidCount` vector becomes empty after filtering, the result is marked as "unclassified" by adding
	 *   a placeholder entry ("unclassified", 1).
	 */
	void secondFilteringStep(std::vector<classifyResult>& classifyResults, const std::unordered_set<std::string>& uniqueTaxids)
	{
		for (auto& result : classifyResults)
		{
			// Remove taxids that are not in the uniqueTaxids set
			result.taxidCount.erase(
				std::remove_if(
					result.taxidCount.begin(),
					result.taxidCount.end(),
					[&uniqueTaxids](const auto& pair)
					{
						return uniqueTaxids.find(pair.first) == uniqueTaxids.end();
					}),
				result.taxidCount.end()
			);
			// If all taxids were removed, mark the result as "unclassified"
			if (result.taxidCount.empty())
			{
				result.taxidCount.emplace_back("unclassified", 1);
			}
		}
	}

	/**
	 * @brief Applies a third filtering step to the classification results, focusing on removing low-uniqueness taxids and reassigning them.
	 *
	 * This function identifies taxids with low unique matches relative to their total matches and checks if their associated reads
	 * are shared with other taxids. If a taxid has significant overlap (95% or more) of reads with another taxid, it is removed from
	 * the classification results, and the reads are reassigned to the other taxid.
	 *
	 * @param classifyResults A vector containing the classification results for each sequence.
	 * @param fileInfo A structure holding statistical data about the classification process, including taxid match counts.
	 *
	 * @details
	 * Step 1: Identify taxids with low unique matches.
	 * - The function calculates the unique match rate for each taxid using the data in `fileInfo`.
	 * - Taxids with unique matches less than 5% of their total matches are marked as low-uniqueness taxids.
	 *
	 * Step 2: Map reads to taxids.
	 * - The function creates a map (`taxidToReads`) that associates each taxid with the set of read IDs it appears in.
	 *
	 * Step 3: Compare read sharing between taxids.
	 * - For each low-uniqueness taxid, the function checks how many of its reads are shared with other taxids.
	 * - If 95% or more of its reads are found in another taxid's read set, the original taxid is removed from the results,
	 *   and its reads are reassigned to the other taxid.
	 *
	 * Step 4: Update classification results.
	 * - The classification results are modified by removing occurrences of the low-uniqueness taxid and replacing it with the new taxid.
	 */
	void thirdFilteringStep(std::vector<classifyResult>& classifyResults, FileInfo& fileInfo)
	{
		// Identify taxids with low unique matches (less than 5% of their total matches)
		std::unordered_set<std::string> lowUniqueTaxids;
		for (const auto& [taxid, totalMatches] : fileInfo.taxidTotalMatches)
		{
			size_t uniqueMatches = fileInfo.taxidUniqueMatches[taxid];
			if (uniqueMatches < 0.05 * totalMatches)
			{
				lowUniqueTaxids.insert(taxid);
			}
		}

		// Map reads to the taxids they are associated with
		robin_hood::unordered_flat_map<std::string, std::unordered_set<std::string>> taxidToReads;
		for (const auto& result : classifyResults)
		{
			for (const auto& [taxid, count] : result.taxidCount)
			{
				taxidToReads[taxid].insert(result.id);
			}
		}

		// Compare read sharing between low-uniqueness taxids and other taxids
		for (const auto& taxid : lowUniqueTaxids)
		{
			const auto& reads = taxidToReads[taxid];
			if (reads.empty()) continue; 
			for (const auto& [otherTaxid, otherReads] : taxidToReads)
			{
				if (taxid == otherTaxid) continue;
				// Count shared reads between the current taxid and another taxid
				size_t sharedReads = 0;
				for (const auto& readId : reads)
				{
					if (otherReads.find(readId) != otherReads.end())
					{
						sharedReads++;
					}
				}
				// Calculate the percentage of shared reads
				double sharedPercentage = static_cast<double>(sharedReads) / reads.size();
				if (sharedPercentage >= 0.95)
				{
					// Update classification results: replace occurrences of the low-uniqueness taxid
					for (auto& result : classifyResults)
					{
						if (reads.find(result.id) != reads.end())
						{
							result.taxidCount.erase(
								std::remove_if(
									result.taxidCount.begin(),
									result.taxidCount.end(),
									[&taxid](const auto& pair)
									{
										return pair.first == taxid;
									}),
								result.taxidCount.end()
							);
							result.taxidCount.emplace_back(otherTaxid, 0);
						}
					}
					break;
				}
			}
		}
	}

	/**
	 * @brief Classify the reads using the Interleaved Cuckoo Filter.
	 *
	 * This function performs the classification of reads using the Interleaved Cuckoo Filter (ICF).
	 * It takes the ICF configuration, a queue of batchReads objects, the classification configuration,
	 * the ICF object, the taxid bins, and a vector to store the classification results.
	 * The function creates a minimiser view for hashing and processes batches from the read queue.
	 * Each batch is processed asynchronously using multiple threads, and the results are stored in the classifyResults vector.
	 * The function waits for all the tasks to complete before returning.
	 *
	 * @param icfConfig The configuration for the Interleaved Cuckoo Filter.
	 * @param readQueue The queue of batchReads objects containing the reads to classify.
	 * @param config The configuration for the classification.
	 * @param icf The Interleaved Cuckoo Filter object.
	 * @param taxidBins The taxid bins.
	 * @param classifyResults The vector to store the classification results.
	 * @param fileInfo The information about the files and sequences.
	 */
	void classify(ChimeraBuild::ICFConfig& icfConfig,
		moodycamel::ConcurrentQueue<batchReads>& readQueue,
		ClassifyConfig& config,
		chimera::InterleavedCuckooFilter& icf,
		std::vector<std::pair<std::string, std::size_t>>& taxidBins,
		std::vector<classifyResult>& classifyResults,
		FileInfo& fileInfo)
	{
		auto minimiser_view = seqan3::views::minimiser_hash(
			seqan3::shape{ seqan3::ungapped{ icfConfig.kmer_size } },
			seqan3::window_size{ icfConfig.window_size },
			seqan3::seed{ adjust_seed(icfConfig.kmer_size) });

		LCA lca;
		if (config.lca)
		{
			buildLCA(lca, config.taxFile);
		}

		using minimiser_view_t = decltype(minimiser_view);

		struct classify_thread_data {
			ChimeraBuild::ICFConfig* icfConfig;
			moodycamel::ConcurrentQueue<batchReads>* readQueue;
			ClassifyConfig* config;
			chimera::InterleavedCuckooFilter* icf;
			std::vector<std::pair<std::string, std::size_t>>* taxidBins;
			FileInfo localFileInfo;
			LCA* lca;
			minimiser_view_t minimiser_view;
			std::vector<classifyResult> localClassifyResults;

			classify_thread_data(
				ChimeraBuild::ICFConfig* icfConfig,
				moodycamel::ConcurrentQueue<batchReads>* readQueue,
				ClassifyConfig* config,
				chimera::InterleavedCuckooFilter* icf,
				std::vector<std::pair<std::string, std::size_t>>* taxidBins,
				LCA* lca,
				minimiser_view_t minimiser_view)
				: icfConfig(icfConfig),
				readQueue(readQueue),
				config(config),
				icf(icf),
				taxidBins(taxidBins),
				lca(lca),
				minimiser_view(std::move(minimiser_view))
			{
			}
		};

		int num_threads = config.threads;
		std::vector<classify_thread_data> thread_data;
		thread_data.reserve(num_threads);

		for (int i = 0; i < num_threads; ++i) {
			thread_data.emplace_back(
				&icfConfig,
				&readQueue,
				&config,
				&icf,
				&taxidBins,
				&lca,
				minimiser_view
			);
		}

		auto classify_worker = [](void* _data, long tid, int nthr) {
			classify_thread_data* thread_data = static_cast<classify_thread_data*>(_data);
			classify_thread_data& data = thread_data[tid];
			batchReads batch;
			while (data.readQueue->try_dequeue(batch)) {
				processBatch(batch,
					*data.icfConfig,
					*data.taxidBins,
					*data.config,
					*data.icf,
					data.localClassifyResults,
					data.minimiser_view,
					data.localFileInfo,
					*data.lca);
			}
			};

		void* thread_pool = kt_forpool_init(num_threads);

		kt_forpool(thread_pool, classify_worker, thread_data.data(), num_threads);

		kt_forpool_destroy(thread_pool);

		for (int i = 0; i < num_threads; ++i) {
			classifyResults.insert(classifyResults.end(),
				thread_data[i].localClassifyResults.begin(),
				thread_data[i].localClassifyResults.end());

			fileInfo.uniqueTaxids.insert(thread_data[i].localFileInfo.uniqueTaxids.begin(),
				thread_data[i].localFileInfo.uniqueTaxids.end());


			for (const auto& [taxid, count] : thread_data[i].localFileInfo.taxidTotalMatches) {
				fileInfo.taxidTotalMatches[taxid] += count;
			}
			for (const auto& [taxid, count] : thread_data[i].localFileInfo.taxidUniqueMatches) {
				fileInfo.taxidUniqueMatches[taxid] += count;
			}

			fileInfo.classifiedNum += thread_data[i].localFileInfo.classifiedNum;
			fileInfo.unclassifiedNum += thread_data[i].localFileInfo.unclassifiedNum;
		}

		secondFilteringStep(classifyResults, fileInfo.uniqueTaxids);

		thirdFilteringStep(classifyResults, fileInfo);

	}

	/**
	 * @brief Save the classification results to an output file.
	 *
	 * This function saves the classification results to an output file in TSV format.
	 * It takes the classification results and the configuration for saving the results.
	 * The function opens the output file and writes the classification results to the file.
	 * Each line in the file contains the sequence ID followed by the taxid and count pairs.
	 * The function closes the file after writing the results.
	 *
	 * @param classifyResults The classification results to save.
	 * @param config The configuration for saving the results.
	 */
	void saveResult(std::vector<classifyResult> classifyResults,
		ClassifyConfig config)
	{
		// Ensure the output file has a .tsv extension
		std::string outputFile = config.outputFile;
		if (std::filesystem::path(outputFile).extension() != ".tsv") {
			outputFile += ".tsv";
		}

		// Open the output file
		std::ofstream os(outputFile, std::ios::out);

		// Check if the file is successfully opened
		if (!os.is_open())
		{
			throw std::runtime_error("Failed to open file: " + config.outputFile);
		}

		// Write the classification results to the output file
		for (const auto& result : classifyResults)
		{
			os << result.id << '\t';
			for (const auto& [taxid, count] : result.taxidCount)
			{
				if (taxid == "unclassified")
				{
					os << taxid;
					continue;
				}
				os << taxid << ':' << count << '\t';
			}
			os << '\n';
		}

		// Close the file
		os.close();
	}


	void loadFilter(const std::string& inputFile,
		chimera::hicf::HierarchicalInterleavedCuckooFilter& hicf,
		ChimeraBuild::HICFConfig& hicfConfig,
		std::vector<std::string>& indexToTaxid)
	{
		// Open the input file
		std::ifstream is(inputFile, std::ios::binary);

		// Check if the file is successfully opened
		if (!is.is_open()) {
			throw std::runtime_error("Failed to open file: " + inputFile);
		}

		// Create a cereal binary input archive
		cereal::BinaryInputArchive archive(is);

		archive(hicf);
		archive(indexToTaxid);
		archive(hicfConfig);
		
		is.close();
	}


	inline void processSequence(const std::vector<size_t>& hashs1,
		ChimeraBuild::HICFConfig& hicfConfig,
		const std::vector<std::string>& indexToTaxid,
		ClassifyConfig& config,
		chimera::hicf::HierarchicalInterleavedCuckooFilter& hicf,
		const std::string& id,
		std::vector<classifyResult>& classifyResults,
		FileInfo& fileInfo,
		LCA& lca)
	{
		// Calculate the number of hash values and the threshold for classification
		size_t hashNum = hashs1.size();
		size_t threshold = std::ceil(hashNum * config.shotThreshold);
		if (threshold == 0)
		{
			threshold = 1;
		}

		// Perform bulk count operation on the hash values using the Interleaved Cuckoo Filter
		auto count = hicf.bulkCount(hashs1, threshold);
		classifyResult result;
		result.id = id;

		size_t oldIndex = 0;
		size_t maxBinCount = 0;
		std::pair<std::string, std::size_t> maxCount;
		robin_hood::unordered_flat_map<std::string, size_t> taxidToCount;
		size_t i = 0;
		for (const auto& taxid : indexToTaxid)
		{
			taxidToCount[taxid] = kv_A(count, i);
			if (kv_A(count, i) > maxBinCount)
			{
				maxBinCount = kv_A(count, i);
			}
			++i;
		}

		// First filter step: remove bins with count less than 80% of the maximum count
		size_t thresholdCount = static_cast<size_t>(0.8 * static_cast<double>(maxBinCount));
		for (const auto& [taxid, binCount] : taxidToCount)
		{
			if (binCount >= thresholdCount && binCount >= threshold)
			{
				result.taxidCount.emplace_back(taxid, binCount);

				if (binCount == maxBinCount)
				{
					maxCount = std::make_pair(taxid, binCount);
				}
				fileInfo.taxidTotalMatches[taxid] += 1;
			}
		}

		// If the maximum count is greater than 0, update the classifyResult object based on the classification mode
		if (!result.taxidCount.empty())
		{

			bool isUniqueMapping = (result.taxidCount.size() == 1);
			if (isUniqueMapping)
			{
				fileInfo.uniqueTaxids.insert(result.taxidCount.front().first);
				const std::string& taxid = result.taxidCount.front().first;
				fileInfo.taxidUniqueMatches[taxid] += 1;
			}

			fileInfo.classifiedNum++;
			if (config.lca && result.taxidCount.size() > 1)
			{
				std::vector<std::string> taxids;
				for (auto& [taxid, count] : result.taxidCount)
				{
					taxids.push_back(taxid);
				}
				std::string lcaTaxid = lca.getLCA(taxids);
				result.taxidCount.clear();
				result.taxidCount.emplace_back(lcaTaxid, 0);
			}
			else if (config.mode == "fast" || result.taxidCount.size() == 1)
			{
				result.taxidCount.clear();
				result.taxidCount.emplace_back(maxCount);
			}
			else if (config.mode == "normal")
			{
				//// Sort the taxidCount vector based on the count in descending order
				//std::sort(result.taxidCount.begin(), result.taxidCount.end(), [](const auto& a, const auto& b) {
				//	return a.second > b.second;  // Sort in descending order
			}
		}
		else
		{
			fileInfo.unclassifiedNum++;
			result.taxidCount.emplace_back("unclassified", 1);
		}

		// Lock the resultMutex to ensure thread safety and add the classifyResult object to the classifyResults vector
		kv_destroy(count);
		classifyResults.emplace_back(std::move(result));
	}


	inline void processBatch(batchReads batch,
		ChimeraBuild::HICFConfig& hicfConfig,
		const std::vector<std::string>& indexToTaxid,
		ClassifyConfig& config,
		chimera::hicf::HierarchicalInterleavedCuckooFilter& hicf,
		std::vector<classifyResult>& classifyResults,
		const auto& minimiser_view,
		FileInfo& fileInfo,
		LCA& lca)
	{
		// Process batch of reads
		std::vector<size_t> hashs1;
		if (!batch.seqs2.empty())
		{
			// Process paired-end reads
			for (size_t i = 0; i < batch.seqs2.size(); i += 2)
			{
				hashs1.clear();
				if (batch.seqs2[i].size() >= hicfConfig.windowSize)
				{
					// Generate minimizer hash values for the first sequence
					hashs1 = batch.seqs2[i] | minimiser_view | seqan3::ranges::to<std::vector>();
					if (batch.seqs2[i + 1].size() >= hicfConfig.windowSize)
					{
						// Generate minimizer hash values for the second sequence
						std::vector<size_t> hashs2 = batch.seqs2[i + 1] | minimiser_view | seqan3::ranges::to<std::vector>();
						// Combine the hash values from both sequences
						hashs1.insert(hashs1.end(), hashs2.begin(), hashs2.end());
					}
				}
				// Process the combined hash values for classification
				processSequence(hashs1, hicfConfig, indexToTaxid, config, hicf, batch.ids[i >> 1], classifyResults, fileInfo, lca);
			}
		}
		else
		{
			// Process single-end reads
			for (size_t i = 0; i < batch.seqs.size(); i++)
			{
				hashs1.clear();
				if (batch.seqs[i].size() >= hicfConfig.windowSize)
				{
					// Generate minimizer hash values for the sequence
					hashs1 = batch.seqs[i] | minimiser_view | seqan3::ranges::to<std::vector>();
				}
				// Process the hash values for classification
				processSequence(hashs1, hicfConfig, indexToTaxid, config, hicf, batch.ids[i], classifyResults, fileInfo, lca);
			}
		}
	}

	void classify(ChimeraBuild::HICFConfig& hicfConfig,
		moodycamel::ConcurrentQueue<batchReads>& readQueue,
		ClassifyConfig& config,
		chimera::hicf::HierarchicalInterleavedCuckooFilter& hicf,
		std::vector<std::string>& indexToTaxid,
		std::vector<classifyResult>& classifyResults,
		FileInfo& fileInfo)
	{
		auto minimiser_view = seqan3::views::minimiser_hash(
			seqan3::shape{ seqan3::ungapped{ hicfConfig.kmerSize } },
			seqan3::window_size{ hicfConfig.windowSize },
			seqan3::seed{ adjust_seed(hicfConfig.kmerSize) });

		LCA lca;
		if (config.lca)
		{
			buildLCA(lca, config.taxFile);
		}

		using minimiser_view_t = decltype(minimiser_view);

		struct classify_thread_data {
			ChimeraBuild::HICFConfig* hicfConfig;
			moodycamel::ConcurrentQueue<batchReads>* readQueue;
			ClassifyConfig* config;
			chimera::hicf::HierarchicalInterleavedCuckooFilter* hicf;
			std::vector<std::string>* indexToTaxid;
			FileInfo localFileInfo;
			LCA* lca;
			minimiser_view_t minimiser_view;
			std::vector<classifyResult> localClassifyResults;

			classify_thread_data(
				ChimeraBuild::HICFConfig* hicfConfig,
				moodycamel::ConcurrentQueue<batchReads>* readQueue,
				ClassifyConfig* config,
				chimera::hicf::HierarchicalInterleavedCuckooFilter* hicf,
				std::vector<std::string>* indexToTaxid,
				LCA* lca,
				minimiser_view_t minimiser_view)
				: hicfConfig(hicfConfig),
				readQueue(readQueue),
				config(config),
				hicf(hicf),
				indexToTaxid(indexToTaxid),
				lca(lca),
				minimiser_view(std::move(minimiser_view))
			{
			}
		};

		int num_threads = config.threads;
		std::vector<classify_thread_data> thread_data;
		thread_data.reserve(num_threads);

		for (int i = 0; i < num_threads; ++i) {
			thread_data.emplace_back(
				&hicfConfig,
				&readQueue,
				&config,
				&hicf,
				&indexToTaxid,
				&lca,
				minimiser_view
			);
		}

		auto classify_worker = [](void* _data, long tid, int nthr) {
			classify_thread_data* thread_data = static_cast<classify_thread_data*>(_data);
			classify_thread_data& data = thread_data[tid];
			batchReads batch;
			while (data.readQueue->try_dequeue(batch)) {
				processBatch(batch,
					*data.hicfConfig,
					*data.indexToTaxid,
					*data.config,
					*data.hicf,
					data.localClassifyResults,
					data.minimiser_view,
					data.localFileInfo,
					*data.lca);
			}
			};

		void* thread_pool = kt_forpool_init(num_threads);

		kt_forpool(thread_pool, classify_worker, thread_data.data(), num_threads);

		kt_forpool_destroy(thread_pool);

		for (int i = 0; i < num_threads; ++i) {
			classifyResults.insert(classifyResults.end(),
				thread_data[i].localClassifyResults.begin(),
				thread_data[i].localClassifyResults.end());

			fileInfo.uniqueTaxids.insert(thread_data[i].localFileInfo.uniqueTaxids.begin(),
				thread_data[i].localFileInfo.uniqueTaxids.end());


			for (const auto& [taxid, count] : thread_data[i].localFileInfo.taxidTotalMatches) {
				fileInfo.taxidTotalMatches[taxid] += count;
			}
			for (const auto& [taxid, count] : thread_data[i].localFileInfo.taxidUniqueMatches) {
				fileInfo.taxidUniqueMatches[taxid] += count;
			}

			fileInfo.classifiedNum += thread_data[i].localFileInfo.classifiedNum;
			fileInfo.unclassifiedNum += thread_data[i].localFileInfo.unclassifiedNum;
		}

		secondFilteringStep(classifyResults, fileInfo.uniqueTaxids);

		thirdFilteringStep(classifyResults, fileInfo);

	}


	void loadFilter(const std::string& inputFile,
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		ChimeraBuild::IMCFConfig& imcfConfig,
		std::vector<std::vector<std::string>>& indexToTaxid)
	{
		// Open the input file
		std::ifstream is(inputFile, std::ios::binary);

		// Check if the file is successfully opened
		if (!is.is_open()) {
			throw std::runtime_error("Failed to open file: " + inputFile);
		}

		// Create a cereal binary input archive
		cereal::BinaryInputArchive archive(is);

		// Deserialize
		archive(imcf);
		archive(indexToTaxid);
		archive(imcfConfig);

		// Close the file
		is.close();
	}



	/**
	 * @brief Processes a sequence of minimizer hashes for classification and applies filtering to determine taxonomic identity.
	 *
	 * This function handles the classification of a sequence by counting the number of matching hash values in the
	 * Interleaved Merged Cuckoo Filter (IMCF), performing filtering steps, and updating the classification results accordingly.
	 *
	 * @param hashs1 A vector containing the minimizer hash values generated from the sequence.
	 * @param imcfConfig Configuration settings for the IMCF, such as k-mer size and window size.
	 * @param indexToTaxid A mapping from IMCF indices to taxids used for classification.
	 * @param config Configuration settings for classification, including mode and filtering thresholds.
	 * @param imcf A reference to the IMCF used for counting hash matches.
	 * @param id The identifier for the current sequence being processed.
	 * @param classifyResults A vector to store the classification results for each processed sequence.
	 * @param fileInfo A structure to store information about processed reads, including counts of classified and unclassified reads.
	 * @param lca A reference to an LCA (Lowest Common Ancestor) structure for taxonomic resolution if LCA mode is enabled.
	 *
	 * @details
	 * Step 1: The function calculates the number of hash values (hashNum) and the minimum threshold for classification.
	 * If the threshold is calculated to be 0, it is set to 1.
	 *
	 * Step 2: The function performs a bulk count operation using the IMCF to count how many times each taxid appears for the given hash values.
	 *
	 * Step 3: The function aggregates the counts for each taxid and identifies the maximum count observed (maxBinCount).
	 *
	 * Step 4: The first filtering step removes bins with counts less than 80% of maxBinCount or below the minimum threshold.
	 * Taxids that pass the filter are added to the result.taxidCount vector.
	 *
	 * Step 5: If any taxids remain after filtering:
	 * - If only one taxid remains, it is marked as a unique mapping.
	 * - If LCA mode is enabled and multiple taxids remain, the LCA is determined and stored as the result.
	 * - If the classification mode is "fast" or only one taxid remains, the taxid with the maximum count is chosen.
	 * - If the classification mode is "normal," the results are left for further processing with the EM algorithm.
	 *
	 * Step 6: If no taxids pass the filters, the sequence is marked as unclassified.
	 *
	 * The function updates shared fileInfo fields for classified and unclassified reads and ensures thread safety when adding
	 * results to the classifyResults vector.
	 */
	inline void processSequence(const std::vector<size_t>& hashs1,
		ChimeraBuild::IMCFConfig& imcfConfig,
		std::vector<std::vector<std::string>>& indexToTaxid,
		ClassifyConfig& config,
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		const std::string& id,
		std::vector<classifyResult>& classifyResults,
		FileInfo& fileInfo,
		LCA& lca)
	{
		// Calculate the number of hash values and determine the threshold for classification
		size_t hashNum = hashs1.size();
		size_t threshold = std::ceil(hashNum * config.shotThreshold);
		if (threshold == 0)
		{
			threshold = 1;
		}

		// Perform bulk count operation using the IMCF
		std::vector<std::vector<size_t>> count;
		count.resize(indexToTaxid.size());
		for (int i = 0; i < count.size(); i++)
		{
			count[i].resize(indexToTaxid[i].size(), 0);
		}
		imcf.bulkCount(hashs1, count);
		classifyResult result;
		result.id = id;

		size_t oldIndex = 0;
		size_t maxBinCount = 0;
		std::pair<std::string, std::size_t> maxCount;
		robin_hood::unordered_flat_map<std::string, size_t> taxidToCount;

		// Aggregate counts for each taxid
		for (int i = 0; i < indexToTaxid.size(); i++)
		{
			for (int j = 0; j < indexToTaxid[i].size(); j++)
			{
				taxidToCount[indexToTaxid[i][j]] += count[i][j];
			}
		}

		// Determine the maximum count observed
		for (auto& [taxid, binCount] : taxidToCount)
		{
			if (binCount > hashNum)
			{
				binCount = hashNum;
			}
			if (binCount > maxBinCount)
			{
				maxBinCount = binCount;
			}
		}


		// First filtering step: keep bins with counts >= 80% of maxBinCount and above the threshold
		size_t thresholdCount = static_cast<size_t>(0.8 * static_cast<double>(maxBinCount));
		for (const auto& [taxid, binCount] : taxidToCount)
		{
			if (binCount >= thresholdCount && binCount >= threshold)
			{
				result.taxidCount.emplace_back(taxid, binCount);

				if (binCount == maxBinCount)
				{
					maxCount = std::make_pair(taxid, binCount);
				}
				fileInfo.taxidTotalMatches[taxid] += 1;
			}
		}

		// Update classifyResult based on the number of remaining taxids and classification mode
		if (!result.taxidCount.empty())
		{

			bool isUniqueMapping = (result.taxidCount.size() == 1);
			if (isUniqueMapping)
			{
				fileInfo.uniqueTaxids.insert(result.taxidCount.front().first);
				const std::string& taxid = result.taxidCount.front().first;
				fileInfo.taxidUniqueMatches[taxid] += 1;
			}

			fileInfo.classifiedNum++;
			if (config.lca && result.taxidCount.size() > 1)
			{
				std::vector<std::string> taxids;
				for (auto& [taxid, count] : result.taxidCount)
				{
					taxids.push_back(taxid);
				}
				std::string lcaTaxid = lca.getLCA(taxids);
				result.taxidCount.clear();
				result.taxidCount.emplace_back(lcaTaxid, 0);
			}
			else if (config.mode == "fast" || result.taxidCount.size() == 1)
			{
				result.taxidCount.clear();
				result.taxidCount.emplace_back(maxCount);
			}
			else if (config.mode == "normal")
			{
				//// Sort the taxidCount vector based on the count in descending order
				//std::sort(result.taxidCount.begin(), result.taxidCount.end(), [](const auto& a, const auto& b) {
				//	return a.second > b.second;  // Sort in descending order
				// Leave result for further EM algorithm processing
			}
		}
		else
		{
			fileInfo.unclassifiedNum++;
			result.taxidCount.emplace_back("unclassified", 1);
		}

		// Add the classifyResult to the shared classifyResults vector
		classifyResults.emplace_back(std::move(result));
	}

	
	/**
	 * @brief Processes a batch of reads for classification using the Interleaved Merged Cuckoo Filter (IMCF).
	 *
	 * This function handles both paired-end and single-end read processing. It generates minimizer hash values
	 * from the reads, combines them as needed, and passes the hash values to the `processSequence` function for
	 * classification. The results are stored in the `classifyResults` vector and file statistics are updated in `fileInfo`.
	 *
	 * @param batch A structure containing the batch of reads to be processed, including read sequences and their IDs.
	 * @param imcfConfig Configuration settings for the IMCF, including k-mer size and window size.
	 * @param indexToTaxid A mapping of indices in the IMCF to the corresponding taxids for classification.
	 * @param config Configuration settings for classification, including filtering options.
	 * @param imcf The Interleaved Merged Cuckoo Filter used for classification.
	 * @param classifyResults A vector to store the classification results for each read.
	 * @param minimiser_view A view used to generate minimizer hash values from the read sequences.
	 * @param fileInfo A structure to store information about the processed reads, such as the number of classified and unclassified reads.
	 * @param lca A reference to an LCA (Lowest Common Ancestor) structure, used if LCA classification is enabled.
	 *
	 * @details
	 * The function first checks if paired-end reads (`batch.seqs2`) are present. If so, it processes the reads in pairs,
	 * generating minimizer hash values for both sequences in each pair and combining the results. For single-end reads,
	 * it processes each sequence individually.
	 *
	 * The generated hash values are passed to the `processSequence` function for classification, which updates the
	 * classification results and the file information.
	 *
	 * - If the read length is smaller than the window size specified in `imcfConfig`, the read is skipped.
	 * - The function ensures that hash values are generated only for reads that meet the minimum length requirement.
	 */
	inline void processBatch(batchReads batch,
		ChimeraBuild::IMCFConfig& imcfConfig,
		std::vector<std::vector<std::string>>& indexToTaxid,
		ClassifyConfig& config,
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		std::vector<classifyResult>& classifyResults,
		const auto& minimiser_view,
		FileInfo& fileInfo,
		LCA& lca)
	{
		// Process batch of reads
		std::vector<size_t> hashs1;
		if (!batch.seqs2.empty())
		{
			// Process paired-end reads
			for (size_t i = 0; i < batch.seqs2.size(); i += 2)
			{
				hashs1.clear();
				if (batch.seqs2[i].size() >= imcfConfig.windowSize)
				{
					// Generate minimizer hash values for the first sequence
					hashs1 = batch.seqs2[i] | minimiser_view | seqan3::ranges::to<std::vector>();
					if (batch.seqs2[i + 1].size() >= imcfConfig.windowSize)
					{
						// Generate minimizer hash values for the second sequence
						std::vector<size_t> hashs2 = batch.seqs2[i + 1] | minimiser_view | seqan3::ranges::to<std::vector>();
						// Combine the hash values from both sequences
						hashs1.insert(hashs1.end(), hashs2.begin(), hashs2.end());
					}
				}
				// Process the combined hash values for classification
				processSequence(hashs1, imcfConfig, indexToTaxid, config, imcf, batch.ids[i >> 1], classifyResults, fileInfo, lca);
			}
		}
		else
		{
			// Process single-end reads
			for (size_t i = 0; i < batch.seqs.size(); i++)
			{
				hashs1.clear();
				if (batch.seqs[i].size() >= imcfConfig.windowSize)
				{
					// Generate minimizer hash values for the sequence
					hashs1 = batch.seqs[i] | minimiser_view | seqan3::ranges::to<std::vector>();
				}
				// Process the hash values for classification
				processSequence(hashs1, imcfConfig, indexToTaxid, config, imcf, batch.ids[i], classifyResults, fileInfo, lca);
			}
		}
	}

	/**
	 * @brief Classifies reads using the Interleaved Merged Cuckoo Filter (IMCF) and performs post-classification filtering.
	 *
	 * This function processes batches of reads using the IMCF, classifies them, and then applies three filtering steps
	 * (the first step is performed within the `processBatch` function). The classified results are collected, and
	 * information about the classified and unclassified reads is updated.
	 *
	 * @param imcfConfig Configuration parameters for the IMCF, including k-mer size and window size.
	 * @param readQueue A concurrent queue holding batches of reads to be classified.
	 * @param config Configuration settings for the classification, including LCA settings and taxonomic file paths.
	 * @param imcf A reference to the Interleaved Merged Cuckoo Filter used for classification.
	 * @param indexToTaxid A mapping of indices in the IMCF to taxids for resolving classification results.
	 * @param classifyResults A vector to store the classification results for each read.
	 * @param fileInfo A structure to store information about the processed reads, including counts of classified and unclassified reads.
	 *
	 * @details
	 * The function processes the input reads in parallel using OpenMP. Each thread dequeues a batch of reads from
	 * `readQueue`, processes it using `processBatch`, and stores the results in local variables. These local results
	 * are merged into the shared `classifyResults` vector and `fileInfo` structure in a critical section to ensure
	 * thread safety.
	 *
	 * The `processBatch` function performs the actual classification and the first filtering step. After all batches
	 * are processed, the function applies the `secondFilteringStep` and `thirdFilteringStep` functions to refine
	 * the classification results.
	 *
	 * If LCA (Lowest Common Ancestor) classification is enabled in `config`, an LCA structure is built using the
	 * `buildLCA` function with the provided taxonomic file.
	 *
	 * @note
	 * - The function assumes that the `processBatch`, `secondFilteringStep`, and `thirdFilteringStep` functions
	 *   are defined and correctly implemented.
	 * - OpenMP is used for parallel processing. Ensure proper thread synchronization and safety when using shared resources.
	 */
	void classify(ChimeraBuild::IMCFConfig& imcfConfig,
		moodycamel::ConcurrentQueue<batchReads>& readQueue,
		ClassifyConfig& config,
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		std::vector<std::vector<std::string>>& indexToTaxid,
		std::vector<classifyResult>& classifyResults,
		FileInfo& fileInfo)
	{
		// Create a minimiser hash view based on IMCF configuration
		auto minimiser_view = seqan3::views::minimiser_hash(
			seqan3::shape{ seqan3::ungapped{ imcfConfig.kmerSize } },
			seqan3::window_size{ imcfConfig.windowSize },
			seqan3::seed{ adjust_seed(imcfConfig.kmerSize) });

		// Initialize LCA structure if LCA mode is enabled in the configuration
		LCA lca;
		if (config.lca)
		{
			buildLCA(lca, config.taxFile);
		}

		// Parallel processing of batches using OpenMP
#pragma omp parallel
		{
			batchReads batch;
			std::vector<classifyResult> localClassifyResults;
			FileInfo localFileInfo;
			// Dequeue and process batches until the queue is empty
			while (readQueue.try_dequeue(batch))
			{
				processBatch(batch, imcfConfig, indexToTaxid, config, imcf, localClassifyResults, minimiser_view, localFileInfo, lca);
			}
			// Critical section to merge local results into shared resources
#pragma omp critical
			{
				classifyResults.insert(classifyResults.end(),
					localClassifyResults.begin(),
					localClassifyResults.end());

				fileInfo.classifiedNum += localFileInfo.classifiedNum;
				fileInfo.unclassifiedNum += localFileInfo.unclassifiedNum;

				// Merge unique taxids from local results
				fileInfo.uniqueTaxids.insert(localFileInfo.uniqueTaxids.begin(),
					localFileInfo.uniqueTaxids.end());

				// Update taxid match counts
				for (const auto& [taxid, count] : localFileInfo.taxidTotalMatches)
				{
					fileInfo.taxidTotalMatches[taxid] += count;
				}

				for (const auto& [taxid, count] : localFileInfo.taxidUniqueMatches)
				{
					fileInfo.taxidUniqueMatches[taxid] += count;
				}
			}
		}

		// Apply the second and third filtering steps
		secondFilteringStep(classifyResults, fileInfo.uniqueTaxids);
		thirdFilteringStep(classifyResults, fileInfo);

	}

	/**
	 * @brief Run the classification process.
	 *
	 * This function runs the classification process using the provided configuration.
	 * It prints the configuration if the verbose flag is set.
	 * The function sets the number of threads for OpenMP and starts the classification process.
	 * It measures the time taken for reading, classifying, and saving the results.
	 * The function prints the time taken for each step if the verbose flag is set.
	 *
	 * @param config The configuration for the classification process.
	 */
	void run(ClassifyConfig config) {
		if ((config.em || config.lca || config.vem) && config.mode == "fast")
		{
			config.mode = "normal";
			std::cout << "Warning: The mode is changed to 'normal' for EM algorithm or LCA classification\n";
		}
		if (config.verbose) {
			std::cout << config << std::endl;
		}
		omp_set_num_threads(config.threads);
		auto TotalclassifyStart = std::chrono::high_resolution_clock::now();
		auto readStart = std::chrono::high_resolution_clock::now();
		std::cout << "Reading input files..." << std::endl;

		FileInfo fileInfo;
		seqan3::contrib::bgzf_thread_count = config.threads;
		moodycamel::ConcurrentQueue<batchReads> readQueue;
		parseReads(readQueue, config, fileInfo);
		std::vector<classifyResult> classifyResults;
		if (config.filter == "imcf")
		{
			std::vector<std::vector<std::string>> indexToTaxid;
			chimera::imcf::InterleavedMergedCuckooFilter imcf;
			ChimeraBuild::IMCFConfig imcfConfig;
			loadFilter(config.dbFile, imcf, imcfConfig, indexToTaxid);
			auto readEnd = std::chrono::high_resolution_clock::now();
			auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(readEnd - readStart);
			if (config.verbose) {
				std::cout << "Read time: ";
				print_classify_time(readDuration.count());
				std::cout << std::endl;
			}

			auto classifyStart = std::chrono::high_resolution_clock::now();
			std::cout << "Classifying sequences by imcf..." << std::endl;
			classify(imcfConfig, readQueue, config, imcf, indexToTaxid, classifyResults, fileInfo);
			auto classifyEnd = std::chrono::high_resolution_clock::now();
			auto classifyDuration = std::chrono::duration_cast<std::chrono::milliseconds>(classifyEnd - classifyStart);
			if (config.verbose) {
				std::cout << "Classify time: ";
				print_classify_time(classifyDuration.count());
			}
		}
		else if (config.filter == "icf")
		{
			std::vector<std::pair<std::string, std::size_t>> taxidBins;
			chimera::InterleavedCuckooFilter icf;
			ChimeraBuild::ICFConfig icfConfig;
			robin_hood::unordered_flat_map<std::string, uint64_t> hashCount;
			loadFilter(config.dbFile, icf, icfConfig, hashCount, taxidBins);
			auto readEnd = std::chrono::high_resolution_clock::now();
			auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(readEnd - readStart);
			if (config.verbose) {
				std::cout << "Read time: ";
				print_classify_time(readDuration.count());
				std::cout << "Number of taxids: " << taxidBins.size() << std::endl << std::endl;
			}
			auto classifyStart = std::chrono::high_resolution_clock::now();
			std::cout << "Classifying sequences by icf..." << std::endl;
			classify(icfConfig, readQueue, config, icf, taxidBins, classifyResults, fileInfo);
			auto classifyEnd = std::chrono::high_resolution_clock::now();
			auto classifyDuration = std::chrono::duration_cast<std::chrono::milliseconds>(classifyEnd - classifyStart);
			if (config.verbose) {
				std::cout << "Classify time: ";
				print_classify_time(classifyDuration.count());
			}
		}
		else if (config.filter == "hicf")
		{
			std::vector<std::string> indexToTaxid;
			chimera::hicf::HierarchicalInterleavedCuckooFilter hicf;
			ChimeraBuild::HICFConfig hicfConfig;
			loadFilter(config.dbFile, hicf, hicfConfig, indexToTaxid);
			auto readEnd = std::chrono::high_resolution_clock::now();
			auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(readEnd - readStart);
			if (config.verbose) {
				std::cout << "Read time: ";
				print_classify_time(readDuration.count());
				std::cout << "Number of taxids: " << indexToTaxid.size() << std::endl << std::endl;
			}

			auto classifyStart = std::chrono::high_resolution_clock::now();
			std::cout << "Classifying sequences by hicf..." << std::endl;
			classify(hicfConfig, readQueue, config, hicf, indexToTaxid, classifyResults, fileInfo);
			auto classifyEnd = std::chrono::high_resolution_clock::now();
			auto classifyDuration = std::chrono::duration_cast<std::chrono::milliseconds>(classifyEnd - classifyStart);
			if (config.verbose) {
				std::cout << "Classify time: ";
				print_classify_time(classifyDuration.count());
			}
		}
		else
		{
			throw std::runtime_error("Invalid filter type: " + config.filter);
		}



		if (config.em)
		{
			auto EMstart = std::chrono::high_resolution_clock::now();
			std::cout << "Running EM algorithm..." << std::endl;
			classifyResults = EMAlgorithm(classifyResults, config.emIter, config.emThreshold);
			auto EMend = std::chrono::high_resolution_clock::now();
			auto EMduration = std::chrono::duration_cast<std::chrono::milliseconds>(EMend - EMstart);
			if (config.verbose) {
				std::cout << "EM time: ";
				print_classify_time(EMduration.count());
			}
		}
		if (config.vem)
		{
			auto VEMstart = std::chrono::high_resolution_clock::now();
			std::cout << "Running VEM algorithm..." << std::endl;
			classifyResults = VEMAlgorithm(classifyResults, config.emIter, config.emThreshold);
			auto VEMend = std::chrono::high_resolution_clock::now();
			auto VEMduration = std::chrono::duration_cast<std::chrono::milliseconds>(VEMend - VEMstart);
			if (config.verbose) {
				std::cout << "VEM time: ";
				print_classify_time(VEMduration.count());
			}
		}

		auto saveStart = std::chrono::high_resolution_clock::now();
		std::cout << "Saving classification results..." << std::endl;
		saveResult(classifyResults, config);
		auto saveEnd = std::chrono::high_resolution_clock::now();
		auto saveDuration = std::chrono::duration_cast<std::chrono::milliseconds>(saveEnd - saveStart);
		if (config.verbose) {
			std::cout << "\nSave time: ";
			print_classify_time(saveDuration.count());
			std::cout << "Total sequences: " << fileInfo.sequenceNum << std::endl;
			std::cout << "Classified sequences: " << fileInfo.classifiedNum << " (" << (static_cast<double>(fileInfo.classifiedNum) / fileInfo.sequenceNum) * 100 << "%)" << std::endl;
			std::cout << "Unclassified sequences: " << fileInfo.unclassifiedNum << " (" << (static_cast<double>(fileInfo.unclassifiedNum) / fileInfo.sequenceNum) * 100 << "%)" << std::endl;
			if (config.lca)
			{
				std::cout << "Total LCA classification: " << fileInfo.lcaNum << " (" << (static_cast<double>(fileInfo.lcaNum) / fileInfo.classifiedNum) * 100 << "%)" << std::endl;
			}
		}

		auto TotalclassifyEnd = std::chrono::high_resolution_clock::now();
		auto TotalclassifyDuration = std::chrono::duration_cast<std::chrono::milliseconds>(TotalclassifyEnd - TotalclassifyStart);

		if (config.verbose) {
			std::cout << "\nTotal classify time: ";
			print_classify_time(TotalclassifyDuration.count());
		}
	}
}