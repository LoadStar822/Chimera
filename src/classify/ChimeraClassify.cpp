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

    void classify(ChimeraBuild::ICFConfig& icfConfig,
        std::queue<batchReads>& readQueue,
        ClassifyConfig& config,
        chimera::InterleavedCuckooFilter& icf,
        std::vector<std::pair<std::string, std::size_t>>& taxidBins,
        std::vector<classifyResult>& classifyResults
         )
    {
        auto minimiser_view = seqan3::views::minimiser_hash(
            seqan3::shape{ seqan3::ungapped{ icfConfig.kmer_size } },
            seqan3::window_size{ icfConfig.window_size },
            seqan3::seed{ adjust_seed(icfConfig.kmer_size) });
        while (true)
        {
            batchReads batch;
            if (!readQueue.empty())
            {
                batch = readQueue.front();
                readQueue.pop();
            }
            else
            {
                break;
            }
            std::vector<size_t> hashs1;
            if (batch.seqs2.size() > 0)
            {
                for (size_t i = 0; i < batch.seqs2.size(); i += 2)
                {
                    classifyResult result;
                    result.id = batch.ids[i >> 1];
					std::vector<std::pair<std::string, std::size_t>> resultBuffer;
					std::pair<std::string, std::size_t> maxCount;
                    size_t seqLen1 = batch.seqs2[i].size();
                    size_t seqLen2 = batch.seqs2[i + 1].size();
                    if (seqLen1 >= icfConfig.window_size)
                    {
                        hashs1 = batch.seqs2[i] | minimiser_view | seqan3::ranges::to<std::vector>();
                        if(seqLen2 >= icfConfig.window_size)
						{
							std::vector<size_t> hashs2 = batch.seqs2[i + 1] | minimiser_view | seqan3::ranges::to<std::vector>();
                            hashs1.insert(hashs1.end(), hashs2.begin(), hashs2.end());
						}
                    }
					size_t hashNum = hashs1.size();
					size_t threshold = std::ceil(hashNum * config.shotThreshold);
					if (threshold == 0)
					{
						threshold = 1;
					}
					auto count = icf.bulk_count(hashs1);
					size_t oldIndex = 0;
					size_t maxBinCount = 0;
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
							resultBuffer.push_back(std::make_pair(taxid, binCount));
							if (binCount > maxBinCount)
							{
								maxBinCount = binCount;
								maxCount = std::make_pair(taxid, binCount);
							}
						}
					}
					if (maxBinCount > 0)
					{
                        if (config.mode == "fast" || resultBuffer.size() == 1)
                        {
                            result.taxidCount.push_back(maxCount);
                        }
					}
                }
            }
            else
            {
                for(size_t i = 0; i < batch.seqs.size(); i++)
				{
					size_t seqLen = batch.seqs[i].size();
					if (seqLen >= icfConfig.window_size)
					{
						hashs1 = batch.seqs[i] | minimiser_view | seqan3::ranges::to<std::vector>();
					}
                    size_t hashNum = hashs1.size();
                    size_t threshold = std::ceil(hashNum * config.shotThreshold);
                    if (threshold == 0)
					{
						threshold = 1;
					}
                    auto count = icf.bulk_count(hashs1);
                    classifyResult result;
                    result.id = batch.ids[i];
                    size_t oldIndex = 0;
                    size_t maxBinCount = 0;
                    std::pair<std::string, std::size_t> maxCount;
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
							result.taxidCount.push_back(std::make_pair(taxid, binCount));
							if (binCount > maxBinCount)
							{
								maxBinCount = binCount;
								maxCount = std::make_pair(taxid, binCount);
							}
						}
					}
                    if (maxBinCount > 0)
					{
						if (config.mode == "fast" || result.taxidCount.size() == 1)
						{
							result.taxidCount.push_back(maxCount);
						}
					}
					classifyResults.push_back(result);
				}
            }

        }
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
        classify(icfConfig, readQueue, config, icf, taxidBins, classifyResults);

        auto classifyEnd = std::chrono::high_resolution_clock::now();
        auto classifyDuration = std::chrono::duration_cast<std::chrono::milliseconds>(classifyEnd - classifyStart);
        if (config.verbose) {
			std::cout << "Classify time: ";
			print_classify_time(classifyDuration.count());
		}


		auto TotalclassifyEnd = std::chrono::high_resolution_clock::now();
		auto TotalclassifyDuration = std::chrono::duration_cast<std::chrono::milliseconds>(TotalclassifyEnd - TotalclassifyStart);

	}
}
