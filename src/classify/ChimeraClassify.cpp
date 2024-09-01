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
 * Last Modified: 2024-08-06
 *
 * Description:
 *  This is a simple C++ program that outputs "Hello, World!".
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#include <ChimeraClassify.hpp>


namespace ChimeraClassify {

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
    void parseReads(std::queue<batchReads>& readQueue, ClassifyConfig config, FileInfo& fileInfo)
    {
        // Count the number of sequences for each thread and accumulate it to fileInfo.sequenceNum
        std::atomic<size_t> localSequenceNum = 0;

        if (config.singleFiles.size() > 0)
        {
            // Process single files in parallel
    #pragma omp parallel for
            for (int i = 0; i < config.singleFiles.size(); ++i)
            {
                const auto& file = config.singleFiles[i];

                // Increase fileInfo.fileNum using atomic operation
    #pragma omp atomic
                fileInfo.fileNum++;

                // Open the input file using sequence_file_input
                seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{ file };
                std::vector<batchReads> localBatches;  // Local batch storage to reduce access to the global queue

                // Read sequences in chunks and create batchReads objects
                for (auto&& rec : fin | seqan3::views::chunk(config.batchSize))
                {
                    batchReads batch;
                    for (auto&& r : rec)
                    {
                        batch.ids.push_back(std::move(r.id()));
                        batch.seqs.push_back(std::move(r.sequence()));
                    }

                    // Update localSequenceNum
                    localSequenceNum += batch.ids.size();
                    // Store the batch in the localBatches vector
                    localBatches.push_back(std::move(batch));
                }

                // Push the local batches to the global queue
    #pragma omp critical
                {
                    for (auto& b : localBatches)
                    {
                        readQueue.push(std::move(b));
                    }
                }
            }

            // Update fileInfo.sequenceNum using atomic operation
    #pragma omp atomic
            fileInfo.sequenceNum += localSequenceNum.load();
        }
        else if (config.pairedFiles.size() > 0)
        {
            // Process paired files in parallel
    #pragma omp parallel for
            for (int i = 0; i < config.pairedFiles.size(); i += 2)
            {
                // Increase fileInfo.fileNum by 2 using atomic operation
    #pragma omp atomic
                fileInfo.fileNum += 2;

                // Open the input files using sequence_file_input
                seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin1{ config.pairedFiles[i] };
                seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin2{ config.pairedFiles[i + 1] };
                std::vector<batchReads> localBatches;  // Local batch storage

                // Read sequences in chunks and create batchReads objects
                for (auto&& rec1 : fin1 | seqan3::views::chunk(config.batchSize))
                {
                    batchReads batch;
                    for (auto&& r : rec1)
                    {
                        batch.ids.push_back(std::move(r.id()));
                        batch.seqs.push_back(std::move(r.sequence()));
                    }
                    for (auto&& rec2 : fin2 | seqan3::views::chunk(config.batchSize))
                    {
                        for (auto&& r : rec2)
                        {
                            batch.seqs2.push_back(std::move(r.sequence()));
                        }
                    }

                    // Update localSequenceNum
                    localSequenceNum += batch.ids.size();
                    // Store the batch in the localBatches vector
                    localBatches.push_back(std::move(batch));
                }

                // Push the local batches to the global queue
    #pragma omp critical
                {
                    for (auto& b : localBatches)
                    {
                        readQueue.push(std::move(b));
                    }
                }
            }

            // Update fileInfo.sequenceNum using atomic operation
    #pragma omp atomic
            fileInfo.sequenceNum += localSequenceNum.load();
        }
        else
        {
            throw std::runtime_error("No input files specified");
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
        robin_hood::unordered_map<std::string, uint64_t>& hashCount,
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

        // Reinsert data from vectors into robin_hood::unordered_map
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
	std::mutex resultMutex;
	inline void processSequence(const std::vector<size_t>& hashs1,
		ChimeraBuild::ICFConfig& icfConfig,
		const std::vector<std::pair<std::string, std::size_t>>& taxidBins,
		ClassifyConfig& config,
		chimera::InterleavedCuckooFilter& icf,
		const std::string& id,
		std::vector<classifyResult>& classifyResults,
		FileInfo& fileInfo)
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

        // Iterate over the taxid bins and calculate the count for each bin based on the bulk count result
        for (auto const& [taxid, bins] : taxidBins)
        {
            size_t binCount = 0;
            for (size_t i = oldIndex; i < bins; i++)
            {
                binCount += kv_A(count, i);
            }
            oldIndex = bins;
            if (binCount >= threshold)
            {
                // Add the bin and count to the classifyResult object if the count is above the threshold
                result.taxidCount.emplace_back(taxid, binCount);
                if (binCount > maxBinCount)
                {
                    maxBinCount = binCount;
                    maxCount = std::make_pair(taxid, binCount);
                }
            }
        }

        // If the maximum count is greater than 0, update the classifyResult object based on the classification mode
        if (maxBinCount > 0)
        {
            fileInfo.classifiedNum++;
            if (config.mode == "fast" || result.taxidCount.size() == 1)
            {
                result.taxidCount.clear();
                result.taxidCount.emplace_back(maxCount);
            }
            else if (config.mode == "normal")
            {
				// Sort the taxidCount vector based on the count in descending order
				std::sort(result.taxidCount.begin(), result.taxidCount.end(), [](const auto& a, const auto& b) {
					return a.second > b.second;  // Sort in descending order
					});
            }
        }
        else
        {
            fileInfo.unclassifiedNum++;
            result.taxidCount.emplace_back("unclassified", 1);
        }

        // Lock the resultMutex to ensure thread safety and add the classifyResult object to the classifyResults vector
        std::lock_guard<std::mutex> lock(resultMutex);
        classifyResults.emplace_back(result);
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
    void processBatch(batchReads batch,
        ChimeraBuild::ICFConfig& icfConfig,
        const std::vector<std::pair<std::string, std::size_t>>& taxidBins,
        ClassifyConfig& config,
        chimera::InterleavedCuckooFilter& icf,
        std::vector<classifyResult>& classifyResults,
        const auto& minimiser_view,
        FileInfo& fileInfo)
    {
        // Process batch of reads
        std::vector<size_t> hashs1;
        if (!batch.seqs2.empty())
        {
            // Process paired-end reads
            for (size_t i = 0; i < batch.seqs2.size(); i += 2)
            {

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
                processSequence(hashs1, icfConfig, taxidBins, config, icf, batch.ids[i >> 1], classifyResults, fileInfo);
            }
        }
        else
        {
            // Process single-end reads
            for (size_t i = 0; i < batch.seqs.size(); i++)
            {
                if (batch.seqs[i].size() >= icfConfig.window_size)
                {
                    // Generate minimizer hash values for the sequence
                    hashs1 = batch.seqs[i] | minimiser_view | seqan3::ranges::to<std::vector>();
                }
                // Process the hash values for classification
                processSequence(hashs1, icfConfig, taxidBins, config, icf, batch.ids[i], classifyResults, fileInfo);
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
        std::queue<batchReads>& readQueue,
        ClassifyConfig& config,
        chimera::InterleavedCuckooFilter& icf,
        std::vector<std::pair<std::string, std::size_t>>& taxidBins,
        std::vector<classifyResult>& classifyResults,
        FileInfo& fileInfo)
    {
        // Create a minimiser view for hashing
        auto minimiser_view = seqan3::views::minimiser_hash(
            seqan3::shape{ seqan3::ungapped{ icfConfig.kmer_size } },
            seqan3::window_size{ icfConfig.window_size },
            seqan3::seed{ adjust_seed(icfConfig.kmer_size) });

        BS::thread_pool pool(config.threads);
		//ctpl::thread_pool pool(config.threads);

        std::vector<std::future<void>> futures;

        // Process batches from the read queue
        while (!readQueue.empty())
        {
            batchReads batch = readQueue.front();
            readQueue.pop();


			futures.emplace_back(pool.submit_task([batch, &icfConfig, &taxidBins, &config, &icf, &classifyResults, &minimiser_view, &fileInfo]() {
                //auto timeStart = std::chrono::high_resolution_clock::now();
				processBatch(std::ref(batch), icfConfig, taxidBins, config, icf, classifyResults, minimiser_view, fileInfo);
    //            auto timeEnd = std::chrono::high_resolution_clock::now();
    //            auto timeDuration = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart);
    //            if (config.verbose) {
				//	std::cout << "Processed batch in: ";
				//	print_classify_time(timeDuration.count());
				//}
				}));
        }

        // Wait for the remaining tasks to complete
        for (auto& f : futures)
        {
            f.get();
        }
    }


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



	void run(ClassifyConfig config) {
		if (config.verbose) {
			std::cout << config << std::endl;
		}
		omp_set_num_threads(config.threads);
		auto TotalclassifyStart = std::chrono::high_resolution_clock::now();
		auto readStart = std::chrono::high_resolution_clock::now();
		std::cout << "Reading input files..." << std::endl;
		FileInfo fileInfo;
		seqan3::contrib::bgzf_thread_count = config.threads;
		std::queue<batchReads> readQueue;
		parseReads(readQueue, config, fileInfo);
        chimera::InterleavedCuckooFilter icf;
        ChimeraBuild::ICFConfig icfConfig;
        robin_hood::unordered_map<std::string, uint64_t> hashCount;
        std::vector<std::pair<std::string, std::size_t>> taxidBins;
        loadFilter(config.dbFile, icf, icfConfig, hashCount, taxidBins);
		auto readEnd = std::chrono::high_resolution_clock::now();
		auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(readEnd - readStart);
		if (config.verbose) {
			std::cout << "Read time: ";
			print_classify_time(readDuration.count());
            std::cout << "Number of taxids: " << taxidBins.size() << std::endl << std::endl;

		}
		auto classifyStart = std::chrono::high_resolution_clock::now();
        std::cout << "Classifying sequences..." << std::endl;
        std::vector<classifyResult> classifyResults;
        classify(icfConfig, readQueue, config, icf, taxidBins, classifyResults, fileInfo);

        auto classifyEnd = std::chrono::high_resolution_clock::now();
        auto classifyDuration = std::chrono::duration_cast<std::chrono::milliseconds>(classifyEnd - classifyStart);
        if (config.verbose) {
			std::cout << "Classify time: ";
			print_classify_time(classifyDuration.count());
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
        }
        


		auto TotalclassifyEnd = std::chrono::high_resolution_clock::now();
		auto TotalclassifyDuration = std::chrono::duration_cast<std::chrono::milliseconds>(TotalclassifyEnd - TotalclassifyStart);

        if(config.verbose) {
			std::cout << "\nTotal classify time: ";
			print_classify_time(TotalclassifyDuration.count());
		}
	}
}
	
