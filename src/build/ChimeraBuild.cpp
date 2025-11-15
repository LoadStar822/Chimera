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
#include <atomic>
#include <future>
#include <memory>
#include <cmath>
#include <algorithm>
#include <array>
#include <sstream>
#include <iomanip>
#include <string_view>
#include <chrono>
#include <cctype>
#include <stdexcept>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <mutex>
#include <span>
#include <unordered_map>

#include <utils/FeatureHasher.hpp>
#include <utils/CountMinSketch.hpp>

namespace ChimeraBuild {
	static inline uint64_t mix64(uint64_t value) {
		return XXH3_64bits(&value, sizeof(value));
	}

	struct TaxidShardPlan {
		size_t groupIndex;
		size_t slotIndex;
		uint64_t expectedCount;
	};

	struct HashFreqStats {
		uint32_t df_high_threshold{ std::numeric_limits<uint32_t>::max() };
		uint32_t df_max_observed{ 0 };
		uint64_t nonzero_counters{ 0 };
	};

	struct HashFrequencyContext {
		HashFreqMode mode{ HashFreqMode::Off };
		double quantile{ 0.999 };
		std::unique_ptr<CountMinSketch> sketch;
		HashFreqStats stats{};
		std::atomic<uint64_t> passA_total_hashes{ 0 };
		std::atomic<uint64_t> passB_total_hashes{ 0 };
		std::atomic<uint64_t> passB_filtered_hashes{ 0 };

		bool enabled() const {
			return mode != HashFreqMode::Off && static_cast<bool>(sketch);
		}
	};

	void print_build_time(long long milliseconds);

	static inline chimera::feature::Method resolve_feature_method(const BuildConfig &config)
	{
		if (config.feature == "syncmer")
			return chimera::feature::Method::Syncmer;
		if (config.feature == "strobemer")
			return chimera::feature::Method::Strobemer;
		// auto: 默认采用 strobemer，后续可根据配置扩展
		return chimera::feature::Method::Strobemer;
	}

	static inline chimera::feature::Params make_feature_params(const BuildConfig &config,
	                                                          chimera::feature::Method &selected,
	                                                          uint64_t &seed_out)
	{
		selected = resolve_feature_method(config);
#ifndef CHIMERA_HAS_STROBEMERS
		if (selected == chimera::feature::Method::Strobemer)
		{
			selected = chimera::feature::Method::Syncmer;
		}
#endif
		chimera::feature::Params params{};
		if (selected == chimera::feature::Method::Strobemer)
		{
			params.method = chimera::feature::Method::Strobemer;
			params.strobe.k = config.strobemer_k;
			params.strobe.order = config.strobemer_order;
			params.strobe.w_min = config.strobemer_w_min;
			params.strobe.w_max = config.strobemer_w_max;
			seed_out = ChimeraBuild::adjust_seed(params.strobe.k);
			params.strobe.seed = seed_out;
			params.strobe.canonical = true;
		}
		else
		{
			params.method = chimera::feature::Method::Syncmer;
			params.sync.k = config.kmer_size;
			params.sync.s = config.smer_size;
			params.sync.pos = config.syncmer_position;
			seed_out = ChimeraBuild::adjust_seed(config.kmer_size);
			params.sync.seed = seed_out;
			params.sync.canonical = true;
		}
		return params;
	}

	void print_build_time(long long milliseconds) {
		long long total_seconds = milliseconds / 1000;
		long long seconds = total_seconds % 60;
		long long total_minutes = total_seconds / 60;
		long long minutes = total_minutes % 60;
		long long hours = total_minutes / 60;

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

	static HashFreqStats compute_hash_freq_stats(const CountMinSketch &cms, double quantile)
	{
		HashFreqStats stats{};
		std::array<uint64_t, 32> histogram{};
		cms.forEachCounter([&](uint32_t value) {
			if (value == 0) {
				return;
			}
			++stats.nonzero_counters;
			if (value > stats.df_max_observed) {
				stats.df_max_observed = value;
			}
			int bucket = 0;
			if (value > 1) {
				double dv = static_cast<double>(value);
				bucket = static_cast<int>(std::log2(dv));
				if (bucket > 31) {
					bucket = 31;
				}
			}
			histogram[static_cast<size_t>(bucket)]++;
		});
		if (stats.nonzero_counters == 0) {
			stats.df_high_threshold = std::numeric_limits<uint32_t>::max();
			return stats;
		}
		const double clamped = std::clamp(quantile, 0.0, 0.999999);
		const double tail_fraction = 1.0 - clamped;
		uint64_t target = static_cast<uint64_t>(std::ceil(tail_fraction * static_cast<double>(stats.nonzero_counters)));
		if (target == 0) {
			target = 1;
		}
		uint64_t accumulated = 0;
		for (int bucket = static_cast<int>(histogram.size()) - 1; bucket >= 0; --bucket) {
			accumulated += histogram[static_cast<size_t>(bucket)];
			if (accumulated >= target) {
				uint32_t bucket_threshold = bucket == 0 ? 1u : (1u << bucket);
				stats.df_high_threshold = std::max<uint32_t>(1u, bucket_threshold);
				return stats;
			}
		}
		stats.df_high_threshold = std::max<uint32_t>(1u, stats.df_max_observed);
		return stats;
	}

	static void build_hash_frequency_sketch(const BuildConfig &config,
	                                       const robin_hood::unordered_flat_map<std::string, std::vector<std::string>> &inputFiles,
	                                       HashFrequencyContext &context)
	{
		context.passA_total_hashes.store(0, std::memory_order_relaxed);
		context.passB_total_hashes.store(0, std::memory_order_relaxed);
		context.passB_filtered_hashes.store(0, std::memory_order_relaxed);
		context.stats = {};
		if (context.mode == HashFreqMode::Off) {
			context.sketch.reset();
			return;
		}
		if (config.hash_sketch_depth == 0 || config.hash_sketch_width == 0) {
			throw std::invalid_argument("Hash sketch depth and width must be positive");
		}
		std::cout << "Building hash frequency sketch (mode="
		          << (context.mode == HashFreqMode::BasicFilter ? "basic" : "off")
		          << ", depth=" << config.hash_sketch_depth
		          << ", width=" << config.hash_sketch_width << ")..." << std::endl;
		context.sketch = std::make_unique<CountMinSketch>(config.hash_sketch_depth, config.hash_sketch_width);
		std::vector<std::string> all_files;
		for (const auto &entry : inputFiles) {
			for (const auto &filename : entry.second) {
				all_files.push_back(filename);
			}
		}
		if (all_files.empty()) {
			context.stats = {};
			return;
		}
		chimera::feature::Method feature_method{};
		uint64_t feature_seed = 0;
		auto feature_params = make_feature_params(config, feature_method, feature_seed);
		const size_t feature_min_length = chimera::feature::min_required_length(feature_params);
		const size_t min_required = std::max<size_t>(config.min_length, feature_min_length);
		auto sketch_start = std::chrono::high_resolution_clock::now();

	#pragma omp parallel
		{
			std::vector<uint64_t> hashes;
			hashes.reserve(4096);
			uint64_t local_hash_total = 0;

		#pragma omp for schedule(dynamic)
			for (size_t idx = 0; idx < all_files.size(); ++idx) {
				const std::string &filename = all_files[idx];
				try {
					seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq >> fin{ filename };
					for (auto &record : fin) {
						auto &seq = record.sequence();
						if (seq.size() < min_required) {
							continue;
						}
						hashes.clear();
						chimera::feature::compute_hashes_append(seq, feature_params, hashes);
						local_hash_total += hashes.size();
						for (uint64_t hash : hashes) {
							context.sketch->add(hash);
						}
					}
				}
				catch (const std::exception &ex) {
				#pragma omp critical(sketch_log)
					{
						std::cerr << "Failed to read sequence file for sketch: " << filename
						          << " (" << ex.what() << ")" << std::endl;
					}
				}
			}

			if (local_hash_total > 0) {
				context.passA_total_hashes.fetch_add(local_hash_total, std::memory_order_relaxed);
			}
		}
		auto sketch_end = std::chrono::high_resolution_clock::now();
		auto sketch_time = std::chrono::duration_cast<std::chrono::milliseconds>(sketch_end - sketch_start).count();
		context.stats = compute_hash_freq_stats(*context.sketch, context.quantile);
		std::cout << "Hash frequency sketch time: "
		          << sketch_time / 1000 << "s " << sketch_time % 1000 << "ms" << std::endl;
		const auto hashed = context.passA_total_hashes.load(std::memory_order_relaxed);
		std::cout << "Sketch summary:" << std::endl;
		std::cout << "  Streamed hashes: " << hashed << std::endl;
		std::cout << "  Active counters: " << context.stats.nonzero_counters << std::endl;
		std::cout << "  Max df estimate: " << context.stats.df_max_observed << std::endl;
		if (context.stats.df_high_threshold == std::numeric_limits<uint32_t>::max()) {
			std::cout << "  BasicFilter threshold: disabled (insufficient data)" << std::endl;
		} else {
			std::cout << "  BasicFilter threshold: df >= " << context.stats.df_high_threshold
			          << " (quantile=" << context.quantile << ")" << std::endl;
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

	static inline std::string tmp_hash_path(const std::string &taxid, std::string_view suffix)
	{
		std::string path = "tmp/";
		path.append(taxid);
		path.append(suffix);
		return path;
	}

	static inline std::string tmp_hash_thread_path(const std::string &taxid,
	                                               std::string_view suffix,
	                                               int thread_id)
	{
		std::string path = tmp_hash_path(taxid, suffix);
		path.push_back('.');
		path.append(std::to_string(thread_id));
		return path;
	}

	/**
	 * Count syncmers for each taxid and its associated files.
	 *
	 * @param config The build configuration.
	 * @param inputFiles The map of input files.
	 * @param hashCount The map to store the hash count.
	 * @param fileInfo The struct to store file information.
	 */
	void syncmer_count(
		BuildConfig& config,
		robin_hood::unordered_flat_map<std::string, std::vector<std::string>>& inputFiles,
		robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount,
		FileInfo& fileInfo,
		HashFrequencyContext* hashFreqContext)
	{
	chimera::feature::Method feature_method{};
	uint64_t feature_seed = 0;
	chimera::feature::Params feature_params = make_feature_params(config, feature_method, feature_seed);
	const std::string feature_suffix = (feature_method == chimera::feature::Method::Strobemer) ? ".strb" : ".sync";
	if (config.verbose) {
		if (feature_method == chimera::feature::Method::Strobemer) {
			std::cout << "Feature method: strobemer (k="
			          << static_cast<int>(feature_params.strobe.k)
			          << ", order=" << static_cast<int>(feature_params.strobe.order)
			          << ", w=[" << feature_params.strobe.w_min << ',' << feature_params.strobe.w_max
			          << "], seed=" << static_cast<unsigned long long>(feature_params.strobe.seed)
			          << ")" << std::endl;
		} else {
			std::cout << "Feature method: syncmer (k=" << static_cast<int>(feature_params.sync.k)
			          << ", s=" << feature_params.sync.s
			          << ", pos=" << feature_params.sync.pos
			          << ", seed=" << static_cast<unsigned long long>(feature_params.sync.seed)
			          << ")" << std::endl;
		}
	}

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
		const size_t max_hashes = config.max_hashes_per_taxid;
		constexpr size_t hash_buffer_flush_bytes = 4 * 1024 * 1024; // 4 MiB

		std::vector<std::mutex> taxidMutexes(index);
		// Track per-taxid buffered hashes so quota checks include pending data.
		std::vector<std::atomic<size_t>> pendingHashCounts(index);
		for (auto &pending : pendingHashCounts) {
			pending.store(0, std::memory_order_relaxed);
		}
		// Global file information statistics
		FileInfo globalFileInfo = {};

		// Mutex to protect global file information
		std::mutex fileInfo_mutex;

	const size_t feature_min_length = chimera::feature::min_required_length(feature_params);
	const size_t min_required = std::max<size_t>(config.min_length, feature_min_length);
	if (config.verbose) {
		std::cout << "Feature hash minimum read length: " << min_required << std::endl;
	}
	const bool freq_filter_enabled = (hashFreqContext != nullptr) && hashFreqContext->enabled();
	const uint32_t freq_threshold = freq_filter_enabled ? hashFreqContext->stats.df_high_threshold
		: std::numeric_limits<uint32_t>::max();
	auto keep_hash = [&](uint64_t hash) {
		if (!freq_filter_enabled) {
			return true;
		}
		hashFreqContext->passB_total_hashes.fetch_add(1, std::memory_order_relaxed);
		if (freq_threshold == std::numeric_limits<uint32_t>::max()) {
			return true;
		}
		const uint32_t df_est = hashFreqContext->sketch->estimate(hash);
		if (df_est >= freq_threshold) {
			hashFreqContext->passB_filtered_hashes.fetch_add(1, std::memory_order_relaxed);
			return false;
		}
		return true;
	};

		int used_threads = 1;

	#pragma omp parallel
		{
			std::vector<uint64_t> thread_buffer;
			thread_buffer.reserve(8192);
			size_t last_taxid_index = std::numeric_limits<size_t>::max();
			const int tid = omp_get_thread_num();
#ifdef _OPENMP
#pragma omp single
			{
				used_threads = omp_get_num_threads();
			}
#endif

			auto flush_buffer = [&](size_t taxid_index_local, std::vector<uint64_t> &buffer, int thread_id) {
				if (buffer.empty()) {
					return;
				}

				const size_t buffer_size = buffer.size();
				const std::string &taxid_ref = index_to_taxid[taxid_index_local];
				size_t write_count = buffer_size;
				bool should_write = true;

			size_t reserved = 0;
			{
				std::lock_guard<std::mutex> lock(taxidMutexes[taxid_index_local]);
				size_t &current_count_ref = hashCount[taxid_ref];
				if (max_hashes > 0)
				{
					if (current_count_ref >= max_hashes)
					{
						should_write = false;
					}
					else
					{
						size_t remaining = max_hashes - current_count_ref;
						if (write_count > remaining)
							write_count = remaining;
						should_write = (write_count > 0);
					}
				}
				else
				{
					should_write = true;
				}
				if (should_write)
				{
					current_count_ref += write_count;
					reserved = write_count;
				}
			}

			if (!should_write)
			{
				pendingHashCounts[taxid_index_local].fetch_sub(buffer_size, std::memory_order_relaxed);
				buffer.clear();
					return;
				}

				const std::string part_path = tmp_hash_thread_path(taxid_ref, feature_suffix, thread_id);
				std::ofstream stream(part_path, std::ios::binary | std::ios::app);
				stream.tie(nullptr);
				if (!stream.is_open())
				{
					std::cerr << "Unable to open the feature hash file: "
					          << part_path << std::endl;
					pendingHashCounts[taxid_index_local].fetch_sub(buffer_size, std::memory_order_relaxed);
					buffer.clear();
					return;
				}

				const auto bytes_to_write = static_cast<std::streamsize>(write_count * sizeof(uint64_t));
				stream.write(reinterpret_cast<const char*>(buffer.data()), bytes_to_write);
				if (!stream.good())
				{
					std::cerr << "Failed to write feature hashes for taxid " << taxid_ref << std::endl;
					if (reserved > 0)
					{
						std::lock_guard<std::mutex> lock(taxidMutexes[taxid_index_local]);
						hashCount[taxid_ref] -= reserved;
					}
				}

				pendingHashCounts[taxid_index_local].fetch_sub(buffer_size, std::memory_order_relaxed);
				buffer.clear();
			};

#pragma omp for schedule(dynamic)
			for (size_t idx = 0; idx < taxid_file_pairs.size(); ++idx)
			{
				const auto& [taxid, filename] = taxid_file_pairs[idx];

				size_t taxid_index = taxid_to_index[taxid];
				if (last_taxid_index != taxid_index && !thread_buffer.empty())
				{
					flush_buffer(last_taxid_index, thread_buffer, tid);
				}
				last_taxid_index = taxid_index;

				bool skip_processing = false;
				{
					std::lock_guard<std::mutex> lock(taxidMutexes[taxid_index]);
					size_t observed = hashCount[taxid] + pendingHashCounts[taxid_index].load(std::memory_order_relaxed);
					if (max_hashes > 0 && observed >= max_hashes)
					{
						skip_processing = true;
					}
				}
				if (skip_processing)
				{
					std::lock_guard<std::mutex> lock(fileInfo_mutex);
					fileInfo.skippedNum++;
					continue;
				}

				// Thread-local file information
				FileInfo localFileInfo = {};

			// Calculate the threshold
			uint8_t cutoff = 0;
			if (config.adaptive_cutoff)
			{
				// 基于文件大小估算 cutoff，仅在显式开启时触发
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

			// 未启用 adaptive_cutoff 时强制采用流式路径，max_hashes 仍由 flush_buffer 控制
			const bool streaming = (!config.adaptive_cutoff);

			// Open the sequence file
			seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq >> fin{ filename };
			std::vector<uint64_t> hashes;
			hashes.reserve(4096);

			if (streaming)
			{
				for (auto &record : fin)
				{
					auto &seq = record.sequence();

					if (seq.size() < min_required)
					{
						localFileInfo.skippedSeqNum++;
						continue;
					}
					localFileInfo.sequenceNum++;
					localFileInfo.bpLength += seq.size();

					hashes.clear();
					chimera::feature::compute_hashes_append(seq, feature_params, hashes);
					if (hashes.empty())
						continue;

					size_t appended = 0;
					if (freq_filter_enabled)
					{
						for (uint64_t hash : hashes)
						{
							if (keep_hash(hash))
							{
								thread_buffer.push_back(hash);
								++appended;
							}
						}
					}
					else
					{
						thread_buffer.insert(thread_buffer.end(), hashes.begin(), hashes.end());
						appended = hashes.size();
					}
					if (appended == 0)
					{
						continue;
					}
					pendingHashCounts[taxid_index].fetch_add(appended, std::memory_order_relaxed);
					if (thread_buffer.size() * sizeof(uint64_t) >= hash_buffer_flush_bytes)
					{
						flush_buffer(taxid_index, thread_buffer, tid);
					}
				}
			}
			else
			{
				robin_hood::unordered_flat_map<uint64_t, uint8_t> local_hash_counts;
				for (auto &record : fin)
				{
					auto &seq = record.sequence();

					if (seq.size() < min_required)
					{
						localFileInfo.skippedSeqNum++;
						continue;
					}
					localFileInfo.sequenceNum++;
					localFileInfo.bpLength += seq.size();

					hashes.clear();
					chimera::feature::compute_hashes_append(seq, feature_params, hashes);
					for (uint64_t hash : hashes)
					{
						uint8_t &count = local_hash_counts[hash];
						if (count < 255)
							++count;
					}
				}

				std::vector<uint64_t> filtered_hashes;
				filtered_hashes.reserve(local_hash_counts.size());
				for (const auto &[hash, count] : local_hash_counts)
				{
					if (count >= cutoff)
					{
						filtered_hashes.emplace_back(hash);
					}
				}

				local_hash_counts.clear();

				if (freq_filter_enabled && !filtered_hashes.empty())
				{
					size_t write_idx = 0;
					for (uint64_t hash : filtered_hashes)
					{
						if (keep_hash(hash))
						{
							filtered_hashes[write_idx++] = hash;
						}
					}
					filtered_hashes.resize(write_idx);
				}

				if (max_hashes > 0 && filtered_hashes.size() > max_hashes)
				{
					filtered_hashes.resize(max_hashes);
				}

				if (!filtered_hashes.empty())
				{
					thread_buffer.insert(thread_buffer.end(), filtered_hashes.begin(), filtered_hashes.end());
					pendingHashCounts[taxid_index].fetch_add(filtered_hashes.size(), std::memory_order_relaxed);
					if (thread_buffer.size() * sizeof(uint64_t) >= hash_buffer_flush_bytes)
					{
						flush_buffer(taxid_index, thread_buffer, tid);
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

			if (!thread_buffer.empty() && last_taxid_index != std::numeric_limits<size_t>::max())
			{
				flush_buffer(last_taxid_index, thread_buffer, tid);
			}
		}

			for (size_t idx = 0; idx < index_to_taxid.size(); ++idx)
			{
				const std::string &taxid = index_to_taxid[idx];
				std::string base_path = tmp_hash_path(taxid, feature_suffix);
				{
				std::ofstream base(base_path, std::ios::binary | std::ios::trunc);
				base.tie(nullptr);
				if (!base.is_open())
				{
					std::cerr << "Failed to open merge target: " << base_path << std::endl;
					continue;
				}
				std::array<char, 1 << 20> io_buffer; // 1 MiB chunk for efficient copy
					for (int tid = 0; tid < used_threads; ++tid)
					{
						std::string part_path = tmp_hash_thread_path(taxid, feature_suffix, tid);
						std::ifstream part(part_path, std::ios::binary);
						part.tie(nullptr);
						if (!part.is_open()) {
							continue;
						}

						while (part)
						{
							part.read(io_buffer.data(), static_cast<std::streamsize>(io_buffer.size()));
							std::streamsize got = part.gcount();
							if (got > 0)
							{
								base.write(io_buffer.data(), got);
							}
						}
						part.close();
						std::error_code ec;
						std::filesystem::remove(part_path, ec);
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
	 * This function builds the IMCF by iterating over groups of taxids, reading feature-hash files, and inserting the hashes
	 * into the IMCF. Each group represents a set of taxids that are processed together, with each hash being inserted
	 * into the corresponding index in the filter.
	 *
	 * @param imcf A reference to the Interleaved Merged Cuckoo Filter object to be built.
	 * @param groups A vector of Group objects, where each Group contains a collection of taxids.
	 *
	 * @return A vector of vectors, where each sub-vector maps the index in the IMCF to the corresponding taxid strings.
	 *
	 * @details
	 * The function first resizes `indexToTaxid` to match the number of groups and resizes each sub-vector to the number of taxids
	 * in the group. It then iterates over each group, reads the feature hashes from a file named in the format "tmp/{taxid}<suffix>",
	 * and inserts each hash into the IMCF using the `insertTag` method.
	 *
	 * The function is parallelized using OpenMP to process each group concurrently, speeding up the construction process.
	 *
	 * If the feature file for a taxid cannot be opened, an error message is printed to `std::cerr`.
	 *
	 * @note
	 * Ensure that the temporary files containing feature hashes ("tmp/{taxid}<suffix>") are pre-generated and accessible.
	 * The function assumes that `imcf.insertTag` is implemented and capable of inserting 64-bit feature hashes into the filter.
	 */
	std::vector<std::vector<std::string>> buildIMCF(
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		const std::vector<chimera::imcf::Group>& groups,
		const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount,
		std::string_view featureSuffix)
	{
		std::vector<std::vector<std::string>> indexToTaxid(groups.size());
		const std::string suffix(featureSuffix);
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
					std::cerr << "Warning: feature count mismatch for taxid " << taxid
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
			std::ifstream ifile(tmp_hash_path(taxid, suffix), std::ios::binary);
			if (!ifile) {
#pragma omp critical(imcf_log)
				std::cerr << "Failed to open feature hash file: " << taxid << suffix << std::endl;
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
						// Extra syncmers beyond planned quota
	#pragma omp critical(imcf_log)
						std::cerr << "Warning: taxid " << taxid
							<< " produced more syncmers than expected; extra entries ignored" << std::endl;
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
				std::cerr << "Warning: failed to read syncmer file completely for taxid " << taxid << std::endl;
			}
			ifile.close();

			for (size_t shardIdx = 0; shardIdx < plans.size(); ++shardIdx) {
				flushSlot(shardIdx);
			}

			uint64_t missing = totalRemaining;
			if (missing != 0) {
#pragma omp critical(imcf_log)
				std::cerr << "Warning: taxid " << taxid
					<< " ended with " << missing << " syncmers short of expected quota" << std::endl;
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
		IMCFConfig& imcfConfig,
		bool needIndexStructures)
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

	std::future<void> indexFuture;
	Clock::time_point indexLaunchTime{};
	if (needIndexStructures && runAsyncIndexBuild) {
		indexLaunchTime = Clock::now();
		indexFuture = std::async(std::launch::async, [&imcf]() {
			imcf.buildActiveGroups();
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

	if (needIndexStructures) {
		const char *indexTitle = runAsyncIndexBuild ?
		    "Building index structures (.idx) [async]" :
		    "Building index structures (.idx)";
		std::cout << "  - " << indexTitle << "... " << std::flush;
		long long indexElapsedMs = 0;
		if (indexFuture.valid()) {
			indexFuture.get();
			auto end = Clock::now();
			indexElapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(end - indexLaunchTime).count();
		} else {
			auto start = Clock::now();
			imcf.buildActiveGroups();
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
	}

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
		const char* kTaxonomyMetaFilename = "taxonomy.meta";
		auto trim = [](std::string& value) {
			while (!value.empty() && std::isspace(static_cast<unsigned char>(value.front()))) {
				value.erase(value.begin());
			}
			while (!value.empty() && std::isspace(static_cast<unsigned char>(value.back()))) {
				value.pop_back();
			}
		};
		if ((config.taxonomy_kind == "auto" || config.taxonomy_kind.empty()) ||
			(config.taxonomy_version == "auto" || config.taxonomy_version.empty())) {
			std::filesystem::path metaPath = std::filesystem::path(config.input_file).parent_path() / kTaxonomyMetaFilename;
			if (std::filesystem::exists(metaPath)) {
				std::ifstream metaStream(metaPath);
				std::string line;
				while (std::getline(metaStream, line)) {
					auto pos = line.find('=');
					if (pos == std::string::npos) {
						continue;
					}
					std::string key = line.substr(0, pos);
					std::string value = line.substr(pos + 1);
					trim(key);
					trim(value);
					if (key == "taxonomy_kind") {
						if (config.taxonomy_kind == "auto" || config.taxonomy_kind.empty()) {
							std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
								return static_cast<char>(std::tolower(ch));
							});
							if (!value.empty()) {
								config.taxonomy_kind = value;
							}
						}
					}
					else if (key == "taxonomy_version") {
						if (config.taxonomy_version == "auto" || config.taxonomy_version.empty()) {
							if (!value.empty()) {
								config.taxonomy_version = value;
							}
						}
					}
				}
			}
		}
		if (config.verbose) {
			std::cout << config << std::endl;
		}
		omp_set_num_threads(config.threads);
		auto build_start = std::chrono::high_resolution_clock::now();
		auto read_start = std::chrono::high_resolution_clock::now();
		std::cout << "Reading input files..." << std::endl;
		FileInfo fileInfo;
		if (config.smer_size == 0) {
			std::cerr << "Syncmer s-mer size must be greater than 0." << std::endl;
			return;
		}
		if (config.smer_size >= config.kmer_size) {
			std::cerr << "Syncmer s-mer size must be smaller than k-mer size." << std::endl;
			return;
		}
		const uint16_t syncmer_span = static_cast<uint16_t>(config.kmer_size - config.smer_size + 1);
		if (config.syncmer_position >= syncmer_span) {
			std::cerr << "Syncmer offset must satisfy 0 <= pos < k - s + 1." << std::endl;
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

		HashFrequencyContext hashFreqContext;
		hashFreqContext.mode = config.hash_freq_mode;
		hashFreqContext.quantile = config.hash_filter_quantile;
		if (hashFreqContext.mode != HashFreqMode::Off) {
			build_hash_frequency_sketch(config, inputFiles, hashFreqContext);
		}

		std::filesystem::path dir = "tmp";
		createOrResetDirectory(dir, config);
		auto calculate_start = std::chrono::high_resolution_clock::now();
		std::cout << "Calculating feature hashes..." << std::endl;
		syncmer_count(config, inputFiles, hashCount, fileInfo, hashFreqContext.enabled() ? &hashFreqContext : nullptr);
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
		if (hashFreqContext.enabled()) {
			const uint64_t total_checked = hashFreqContext.passB_total_hashes.load(std::memory_order_relaxed);
			const uint64_t filtered = hashFreqContext.passB_filtered_hashes.load(std::memory_order_relaxed);
			std::cout << "Hash frequency filter dropped " << filtered << " / " << total_checked << " hashes";
			if (total_checked > 0) {
				std::ostringstream ratio_stream;
				double ratio = static_cast<double>(filtered) * 100.0 / static_cast<double>(total_checked);
				ratio_stream << std::fixed << std::setprecision(2) << ratio;
				std::cout << " (" << ratio_stream.str() << "%)";
			}
			std::cout << " with df >= " << hashFreqContext.stats.df_high_threshold << std::endl;
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

		auto imcf_build_start = std::chrono::high_resolution_clock::now();
		std::cout << "Building IMCF..." << std::endl;
		uint64_t imcf_feature_seed = 0;
		chimera::feature::Method imcf_feature_method{};
		make_feature_params(config, imcf_feature_method, imcf_feature_seed);
		const std::string feature_suffix = (imcf_feature_method == chimera::feature::Method::Strobemer) ? ".strb" : ".sync";
		IMCFConfig imcfConfig;
		imcfConfig.loadFactor = config.load_factor;
		imcfConfig.kmerSize = config.kmer_size;
		imcfConfig.smerSize = config.smer_size;
		imcfConfig.syncmerPosition = config.syncmer_position;
		imcfConfig.seed64 = imcf_feature_seed;
		imcfConfig.fpSalt = IMCFConfig::DefaultFingerprintSalt;
		imcfConfig.hashVersion = IMCFConfig::CurrentHashVersion;
		if (imcf_feature_method == chimera::feature::Method::Strobemer) {
			imcfConfig.featureMethod = 1;
			imcfConfig.strobeOrder = config.strobemer_order;
			imcfConfig.strobeWmin = config.strobemer_w_min;
			imcfConfig.strobeWmax = config.strobemer_w_max;
			imcfConfig.strobeK = config.strobemer_k;
		} else {
			imcfConfig.featureMethod = 0;
			imcfConfig.strobeOrder = 0;
			imcfConfig.strobeWmin = 0;
			imcfConfig.strobeWmax = 0;
			imcfConfig.strobeK = 0;
		}
		if (config.taxonomy_kind == "auto" || config.taxonomy_kind.empty()) {
			imcfConfig.taxonomyKind = "ncbi";
		} else {
			imcfConfig.taxonomyKind = config.taxonomy_kind;
		}
		if (config.taxonomy_version == "auto" || config.taxonomy_version.empty()) {
			imcfConfig.taxonomyVersion = "ncbi-taxdump";
		} else {
			imcfConfig.taxonomyVersion = config.taxonomy_version;
		}
		chimera::imcf::InterleavedMergedCuckooFilter imcf(groups, imcfConfig);
		std::vector<std::vector<std::string>> indexToTaxid = buildIMCF(imcf, groups, hashCount, feature_suffix);
		auto imcf_build_end = std::chrono::high_resolution_clock::now();
		auto imcf_build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(imcf_build_end - imcf_build_start).count();
		if (config.verbose) {
			std::cout << "Build time: ";
			print_build_time(imcf_build_total_time);
			std::cout << std::endl;
			std::cout << imcf << std::endl;
		}

		auto save_start = std::chrono::high_resolution_clock::now();
		std::cout << "Saving IMCF..." << std::endl;
		saveIMCF(imcf, config.output_file, indexToTaxid, imcfConfig, true);
		auto save_end = std::chrono::high_resolution_clock::now();
		auto save_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(save_end - save_start).count();
		if (config.verbose) {
			std::cout << "Save time: ";
			print_build_time(save_total_time);
			std::cout << std::endl;
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
