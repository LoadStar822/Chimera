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
#include <array>
#include <cctype>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <string_view>
#include <chrono>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <queue>
#include <mutex>
#include <span>
#include <cstdint>

#include <utils/Syncmer.hpp>

namespace fs = std::filesystem;

namespace ChimeraBuild {
	static inline uint64_t mix64(uint64_t value) {
		return XXH3_64bits(&value, sizeof(value));
	}

	struct TaxidShardPlan {
		size_t groupIndex;
		size_t slotIndex;
		uint64_t expectedCount;
	};

	class BottomKSampler {
	public:
		BottomKSampler() = default;
		explicit BottomKSampler(size_t capacity)
		{
			configure(capacity);
		}

		void configure(size_t capacity)
		{
			capacity_ = capacity;
			heap_.clear();
			present_.clear();
		}

		void insert(uint64_t value)
		{
			if (capacity_ == 0)
			{
				return;
			}

			auto [_, inserted] = present_.insert(value);
			if (!inserted)
			{
				return;
			}

			if (heap_.size() < capacity_)
			{
				heap_.push_back(value);
				std::push_heap(heap_.begin(), heap_.end());
				return;
			}

			// Heap is full; keep the new value only if it is smaller than the current maximum.
			if (value >= heap_.front())
			{
				present_.erase(value);
				return;
			}

			const uint64_t removed = heap_.front();
			std::pop_heap(heap_.begin(), heap_.end());
			heap_.back() = value;
			std::push_heap(heap_.begin(), heap_.end());
			present_.erase(removed);
		}

		void merge(const BottomKSampler& other)
		{
			if (capacity_ == 0)
			{
				return;
			}
			for (uint64_t value : other.heap_)
			{
				insert(value);
			}
		}

		size_t size() const
		{
			return heap_.size();
		}

		std::vector<uint64_t> release_sorted()
		{
			std::sort(heap_.begin(), heap_.end());
			present_.clear();
			return std::move(heap_);
		}

	private:
		size_t capacity_{ 0 };
		std::vector<uint64_t> heap_;
		robin_hood::unordered_flat_set<uint64_t> present_;
	};

	class SpaceSaving {
	public:
		struct Entry {
			uint64_t key{ 0 };
			uint64_t count{ 0 };
			uint64_t error{ 0 };
		};

		SpaceSaving() = default;

		explicit SpaceSaving(size_t capacity)
		{
			configure(capacity);
		}

		void configure(size_t capacity)
		{
			capacity_ = capacity;
			total_weight_ = 0;
			counters_.clear();
			index_.clear();
			min_heap_ = {};
			counters_.reserve(capacity);
			index_.reserve(capacity);
		}

		void add(uint64_t key, uint64_t weight = 1)
		{
			if (capacity_ == 0 || weight == 0)
			{
				return;
			}

			total_weight_ += weight;

			auto it = index_.find(key);
			if (it != index_.end())
			{
				const size_t idx = it->second;
				auto& slot = counters_[idx];
				slot.count += weight;
				++slot.version;
				min_heap_.push({ slot.count, slot.version, idx });
				maybe_rebuild_heap();
				return;
			}

			if (counters_.size() < capacity_)
			{
				const size_t idx = counters_.size();
				counters_.push_back({ key, weight, 0, 1 });
				index_[key] = idx;
				min_heap_.push({ counters_[idx].count, counters_[idx].version, idx });
				maybe_rebuild_heap();
				return;
			}

			prune_heap();
			if (min_heap_.empty())
			{
				return;
			}

			const auto node = min_heap_.top();
			auto& slot = counters_[node.index];
			index_.erase(slot.key);
			const uint64_t previous = slot.count;
			slot.key = key;
			slot.error = previous;
			slot.count = previous + weight;
			++slot.version;
			index_[key] = node.index;
			min_heap_.push({ slot.count, slot.version, node.index });
			maybe_rebuild_heap();
		}

		void merge(const SpaceSaving& other)
		{
			if (capacity_ == 0)
			{
				return;
			}
			for (const auto& slot : other.counters_)
			{
				if (slot.count == 0)
				{
					continue;
				}
				add(slot.key, slot.count);
			}
		}

		uint64_t total_weight() const noexcept
		{
			return total_weight_;
		}

		std::vector<Entry> entries() const
		{
			std::vector<Entry> result;
			result.reserve(counters_.size());
			for (const auto& slot : counters_)
			{
				if (slot.count == 0)
				{
					continue;
				}
				result.push_back({ slot.key, slot.count, slot.error });
			}
			return result;
		}

	private:
		struct Counter {
			uint64_t key{ 0 };
			uint64_t count{ 0 };
			uint64_t error{ 0 };
			uint64_t version{ 0 };
		};

		struct HeapNode {
			uint64_t count;
			uint64_t version;
			size_t index;
		};

		struct MinCompare {
			bool operator()(const HeapNode& lhs, const HeapNode& rhs) const
			{
				if (lhs.count != rhs.count)
				{
					return lhs.count > rhs.count;
				}
				if (lhs.version != rhs.version)
				{
					return lhs.version > rhs.version;
				}
				return lhs.index > rhs.index;
			}
		};

		void prune_heap()
		{
			while (!min_heap_.empty())
			{
				const auto node = min_heap_.top();
				if (node.index >= counters_.size())
				{
					min_heap_.pop();
					continue;
				}
				const auto& slot = counters_[node.index];
				if (slot.version != node.version || slot.count != node.count)
				{
					min_heap_.pop();
					continue;
				}
				break;
			}
		}

		void maybe_rebuild_heap()
		{
			if (counters_.empty())
			{
				return;
			}
			constexpr size_t HeapGrowthFactor = 8;
			if (min_heap_.size() <= counters_.size() * HeapGrowthFactor)
			{
				return;
			}
			rebuild_heap();
		}

		void rebuild_heap()
		{
			std::priority_queue<HeapNode, std::vector<HeapNode>, MinCompare> fresh;
			for (size_t i = 0; i < counters_.size(); ++i)
			{
				const auto& slot = counters_[i];
				if (slot.count == 0)
				{
					continue;
				}
				fresh.push({ slot.count, slot.version, i });
			}
			min_heap_.swap(fresh);
		}

		size_t capacity_{ 0 };
		uint64_t total_weight_{ 0 };
		std::vector<Counter> counters_;
		robin_hood::unordered_flat_map<uint64_t, size_t> index_;
		std::priority_queue<HeapNode, std::vector<HeapNode>, MinCompare> min_heap_;
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
	 * Count syncmers for each taxid and its associated files.
	 *
	 * @param config The build configuration.
	 * @param inputFiles The map of input files.
	 * @param hashCount The map to store the hash count.
	 * @param sampledHashFiles 临时文件表，用于存放各 taxid 的采样 hash。
	 * @param fileInfo The struct to store file information.
	 * @param scratchDir 输出参数，指向创建的临时目录。
	 */
	void syncmer_count(
		BuildConfig& config,
		robin_hood::unordered_flat_map<std::string, std::vector<std::string>>& inputFiles,
		robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount,
		robin_hood::unordered_flat_map<std::string, std::filesystem::path>& sampledHashFiles,
		FileInfo& fileInfo,
		std::filesystem::path& scratchDir)
{
	if (config.max_hashes_per_taxid == 0)
	{
		throw std::runtime_error("bottom-k 采样要求 --max-hashes > 0");
	}
	if (config.adaptive_cutoff)
	{
		std::cerr << "Warning: adaptive cutoff 已被 bottom-k 采样忽略" << std::endl;
	}

	sampledHashFiles.clear();

	fs::path outputPath{ config.output_file };
	fs::path baseDir = outputPath.has_parent_path() ? outputPath.parent_path() : fs::current_path();
	std::string stem = outputPath.has_filename() ? outputPath.filename().string() : "ChimeraDB";
	if (stem.empty()) {
		stem = "ChimeraDB";
	}
	std::string scratchName = stem + ".syncmer.tmp";
	scratchDir = baseDir / scratchName;
	std::error_code ec;
	fs::remove_all(scratchDir, ec);
	fs::create_directories(scratchDir);

	auto make_temp_path = [&](const std::string& taxid, size_t ordinal) {
		std::string safe = taxid;
		for (char& ch : safe) {
			unsigned char uch = static_cast<unsigned char>(ch);
			if (!std::isalnum(uch) && ch != '-' && ch != '_') {
				ch = '_';
			}
		}
		std::ostringstream oss;
		oss << ordinal << '_' << safe << ".syncbin";
		return scratchDir / oss.str();
	};

	const std::array<size_t, 1> syncmer_positions{ static_cast<size_t>(config.syncmer_position) };
	const std::span<const size_t> syncmer_pos_span(syncmer_positions);
	const uint64_t syncmer_seed = ChimeraBuild::adjust_seed(config.kmer_size);
	const size_t max_hashes = config.max_hashes_per_taxid;
	const bool has_top_n = config.toxic_top_n > 0;
	const bool has_quantile = config.toxic_quantile < 1.0;
	const double toxic_fraction = has_quantile ? (1.0 - config.toxic_quantile) : 0.0;
	constexpr size_t kMinHeavyCapacity = 1024;
	constexpr size_t kMaxHeavyCapacity = 1'000'000;
	size_t heavy_capacity = 0;
	if (has_top_n)
	{
		uint64_t scaled = config.toxic_top_n;
		if (scaled > std::numeric_limits<uint64_t>::max() / 4ULL)
		{
			scaled = std::numeric_limits<uint64_t>::max();
		}
		else
		{
			scaled *= 4ULL;
		}
		const uint64_t capped = std::min<uint64_t>(scaled, static_cast<uint64_t>(kMaxHeavyCapacity));
		heavy_capacity = static_cast<size_t>(std::max<uint64_t>(capped, kMinHeavyCapacity));
	}
	else if (has_quantile && toxic_fraction > 0.0)
	{
		double estimate = std::ceil((1.0 / toxic_fraction) * 4.0);
		if (std::isfinite(estimate) && estimate > 0.0)
		{
			if (estimate < static_cast<double>(kMinHeavyCapacity))
			{
				estimate = static_cast<double>(kMinHeavyCapacity);
			}
			if (estimate > static_cast<double>(kMaxHeavyCapacity))
			{
				estimate = static_cast<double>(kMaxHeavyCapacity);
			}
			heavy_capacity = static_cast<size_t>(estimate);
		}
	}
	const bool enable_toxic_filter = heavy_capacity > 0;

	std::vector<std::string> index_to_taxid;
	index_to_taxid.reserve(inputFiles.size());
	for (auto& [taxid, _] : inputFiles) {
		index_to_taxid.push_back(taxid);
	}

	std::mutex fileInfo_mutex;
	std::mutex heavy_mutex;
	SpaceSaving globalHeavy(heavy_capacity);
	std::mutex hash_mutex;
	std::mutex raw_mutex;
	robin_hood::unordered_flat_map<std::string, fs::path> rawSampleFiles;

	auto write_vector_to_file = [](const fs::path& path, const std::vector<uint64_t>& data) {
		std::ofstream ofs(path, std::ios::binary | std::ios::trunc);
		if (!ofs.is_open()) {
			throw std::runtime_error("无法写入临时 syncmer 文件: " + path.string());
		}
		if (!data.empty()) {
			ofs.write(reinterpret_cast<const char*>(data.data()),
			          static_cast<std::streamsize>(data.size() * sizeof(uint64_t)));
		}
		if (!ofs) {
			throw std::runtime_error("写入临时 syncmer 文件失败: " + path.string());
		}
	};

#pragma omp parallel for schedule(dynamic)
	for (size_t taxIdx = 0; taxIdx < index_to_taxid.size(); ++taxIdx) {
		const auto& taxid = index_to_taxid[taxIdx];
		auto filesIt = inputFiles.find(taxid);
		if (filesIt == inputFiles.end()) {
			continue;
		}
		const auto& files = filesIt->second;

		FileInfo localFileInfo{};
		BottomKSampler sampler(max_hashes);
		robin_hood::unordered_flat_set<uint64_t> perReadSeen;
		SpaceSaving localHeavy(enable_toxic_filter ? heavy_capacity : 0);

		for (const auto& filename : files) {
			try {
				seqan3::sequence_file_input<raptor::dna4_traits,
				                            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
				    fin{filename};
				for (auto& record : fin) {
					auto& seq = record.sequence();
					if (seq.size() < config.min_length) {
						localFileInfo.skippedSeqNum++;
						continue;
					}
					localFileInfo.sequenceNum++;
					localFileInfo.bpLength += seq.size();

					perReadSeen.clear();
					auto hashes = chimera::syncmer::compute_hashes(
					    seq, config.smer_size, config.kmer_size, syncmer_pos_span, syncmer_seed, true);
					for (uint64_t hash : hashes) {
						if (perReadSeen.insert(hash).second) {
							sampler.insert(hash);
						}
						if (enable_toxic_filter) {
							localHeavy.add(hash);
						}
					}
				}
			} catch (const std::exception& ex) {
#pragma omp critical(syncmer_log)
				std::cerr << "读取序列文件失败 " << filename << ": " << ex.what() << std::endl;
				std::lock_guard<std::mutex> lock(fileInfo_mutex);
				fileInfo.skippedNum++;
			}
		}

		auto values = sampler.release_sorted();
		fs::path rawPath = make_temp_path(taxid, taxIdx);
		write_vector_to_file(rawPath, values);

		{
			std::lock_guard<std::mutex> lock(hash_mutex);
			hashCount[taxid] = static_cast<uint64_t>(values.size());
		}
		{
			std::lock_guard<std::mutex> lock(raw_mutex);
			rawSampleFiles.emplace(taxid, std::move(rawPath));
		}

		if (enable_toxic_filter && localHeavy.total_weight() > 0) {
			std::lock_guard<std::mutex> heavy_lock(heavy_mutex);
			globalHeavy.merge(localHeavy);
		}

		{
			std::lock_guard<std::mutex> lock(fileInfo_mutex);
			fileInfo.skippedSeqNum += localFileInfo.skippedSeqNum;
			fileInfo.sequenceNum += localFileInfo.sequenceNum;
			fileInfo.bpLength += localFileInfo.bpLength;
		}
	}

	robin_hood::unordered_flat_set<uint64_t> toxic_hashes;
	uint64_t min_count_threshold = 0;
	size_t candidate_count = 0;
	uint64_t target_weight = 0;
	uint64_t removed_weight = 0;
	bool used_topn = false;
	bool used_quantile = false;
	if (enable_toxic_filter)
	{
		auto entries = globalHeavy.entries();
		if (!entries.empty())
		{
			std::sort(entries.begin(), entries.end(), [](const SpaceSaving::Entry& lhs, const SpaceSaving::Entry& rhs)
			{
				if (lhs.count != rhs.count)
				{
					return lhs.count > rhs.count;
				}
				return lhs.key < rhs.key;
			});

			const uint64_t total_seen = globalHeavy.total_weight();
			if (config.toxic_min_fraction > 0.0 && total_seen > 0)
			{
				const double scaled = static_cast<double>(total_seen) * config.toxic_min_fraction;
				min_count_threshold = static_cast<uint64_t>(std::ceil(scaled));
				if (min_count_threshold == 0)
				{
					min_count_threshold = 1;
				}
			}

			std::vector<SpaceSaving::Entry> candidates;
			candidates.reserve(entries.size());
			for (const auto& entry : entries)
			{
				if (entry.count >= min_count_threshold)
				{
					candidates.push_back(entry);
				}
			}

			if (!candidates.empty())
			{
				candidate_count = candidates.size();
				if (has_top_n)
				{
					const uint64_t top_n = config.toxic_top_n;
					const size_t limit = static_cast<size_t>(std::min<uint64_t>(top_n, static_cast<uint64_t>(candidates.size())));
					if (limit > 0)
					{
						toxic_hashes.reserve(limit);
						for (size_t i = 0; i < limit; ++i)
						{
							toxic_hashes.insert(candidates[i].key);
							removed_weight += candidates[i].count;
						}
						used_topn = true;
					}
				}
				else if (has_quantile)
				{
					const double removal_fraction = std::max(0.0, 1.0 - config.toxic_quantile);
					if (removal_fraction > 0.0 && total_seen > 0)
					{
						target_weight = static_cast<uint64_t>(std::ceil(removal_fraction * static_cast<double>(total_seen)));
						if (target_weight == 0)
						{
							target_weight = 1;
						}
						std::sort(candidates.begin(), candidates.end(), [](const SpaceSaving::Entry& lhs, const SpaceSaving::Entry& rhs)
						{
							const uint64_t lhs_base = lhs.count > lhs.error ? lhs.count - lhs.error : 0;
							const uint64_t rhs_base = rhs.count > rhs.error ? rhs.count - rhs.error : 0;
							if (lhs_base != rhs_base)
							{
								return lhs_base > rhs_base;
							}
							return lhs.key < rhs.key;
						});
						toxic_hashes.reserve(candidates.size());
						for (const auto& entry : candidates)
						{
							if (entry.count == 0)
							{
								continue;
							}
							toxic_hashes.insert(entry.key);
							removed_weight += entry.count;
							if (removed_weight >= target_weight)
							{
								break;
							}
						}
						if (removed_weight > 0)
						{
							used_quantile = true;
						}
					}
				}
			}
		}
	}

	const bool has_toxic_hashes = !toxic_hashes.empty();
	size_t safety_base = 0;
	if (has_toxic_hashes)
	{
		size_t fraction_floor = 0;
		if (config.toxic_safety_fraction > 0.0)
		{
			const double scaled = static_cast<double>(max_hashes) * config.toxic_safety_fraction;
			fraction_floor = static_cast<size_t>(std::ceil(scaled));
		}
		safety_base = std::max(fraction_floor, static_cast<size_t>(config.toxic_safety_min));
		if (safety_base > max_hashes)
		{
			safety_base = max_hashes;
		}
	}

sampledHashFiles.reserve(index_to_taxid.size());
	size_t taxa_filtered_count = 0;
	size_t taxa_safety_hits = 0;
	uint64_t hashes_removed_total = 0;
	auto read_vector_from_file = [](const fs::path& path) -> std::vector<uint64_t> {
		std::ifstream ifs(path, std::ios::binary);
		if (!ifs.is_open()) {
			throw std::runtime_error("无法读取临时 syncmer 文件: " + path.string());
		}
		ifs.seekg(0, std::ios::end);
		std::streamsize bytes = ifs.tellg();
		ifs.seekg(0, std::ios::beg);
		if (bytes <= 0) {
			return {};
		}
		size_t count = static_cast<size_t>(bytes) / sizeof(uint64_t);
		std::vector<uint64_t> data(count);
		ifs.read(reinterpret_cast<char*>(data.data()),
		         static_cast<std::streamsize>(count * sizeof(uint64_t)));
		if (!ifs) {
			throw std::runtime_error("读取临时 syncmer 文件失败: " + path.string());
		}
		return data;
	};

	for (size_t i = 0; i < index_to_taxid.size(); ++i) {
		const std::string& taxid = index_to_taxid[i];
		auto rawIt = rawSampleFiles.find(taxid);
		if (rawIt == rawSampleFiles.end()) {
			continue;
		}
		fs::path tmpFile = rawIt->second;
		auto values = read_vector_from_file(tmpFile);
		const size_t original_size = values.size();

		if (has_toxic_hashes) {
			std::vector<uint64_t> filtered;
			filtered.reserve(values.size());
			for (uint64_t hash : values) {
				if (toxic_hashes.find(hash) == toxic_hashes.end()) {
					filtered.push_back(hash);
				}
			}
			size_t min_keep = config.toxic_safety_min;
			if (config.toxic_safety_fraction > 0.0 && original_size > 0) {
				const double scaled =
				    static_cast<double>(original_size) * config.toxic_safety_fraction;
				const size_t fractional_keep = static_cast<size_t>(std::ceil(scaled));
				if (fractional_keep > min_keep) {
					min_keep = fractional_keep;
				}
			}
			if (min_keep > original_size) {
				min_keep = original_size;
			}
			if (filtered.size() < min_keep) {
				{
					std::lock_guard<std::mutex> lock(hash_mutex);
					hashCount[taxid] = static_cast<uint64_t>(values.size());
				}
				write_vector_to_file(tmpFile, values);
				++taxa_safety_hits;
			} else {
				const size_t removed_here = original_size - filtered.size();
				if (removed_here > 0) {
					hashes_removed_total += static_cast<uint64_t>(removed_here);
					++taxa_filtered_count;
				}
				{
					std::lock_guard<std::mutex> lock(hash_mutex);
					hashCount[taxid] = static_cast<uint64_t>(filtered.size());
				}
				write_vector_to_file(tmpFile, filtered);
			}
		} else {
			{
				std::lock_guard<std::mutex> lock(hash_mutex);
				hashCount[taxid] = static_cast<uint64_t>(values.size());
			}
			write_vector_to_file(tmpFile, values);
		}

		{
			std::lock_guard<std::mutex> lock(raw_mutex);
			sampledHashFiles.emplace(taxid, std::move(tmpFile));
		}
	}

	if (config.verbose && has_toxic_hashes)
	{
		std::cerr << "[chimera] 毒 syncmer 过滤: 黑名单 " << toxic_hashes.size()
			<< " 个 hash（候选 " << candidate_count
			<< "，容量 " << heavy_capacity
			<< "，频率阈值 " << min_count_threshold;
		if (used_topn)
		{
			std::cerr << "，模式 topN";
		}
		else if (used_quantile)
		{
			const double weight_ratio = (globalHeavy.total_weight() > 0)
				? static_cast<double>(removed_weight) / static_cast<double>(globalHeavy.total_weight())
				: 0.0;
			std::ostringstream ratio_stream;
			ratio_stream << std::fixed << std::setprecision(4) << weight_ratio;
			std::cerr << "，模式 quantile，目标权重 "
				<< target_weight << "，实际裁剪权重 " << removed_weight
				<< " (" << ratio_stream.str() << ")";
		}
		std::cerr << "），被过滤 taxid 数 " << taxa_filtered_count
			<< "，合计剔除 " << hashes_removed_total
			<< " 个采样 hash，安全网回退 " << taxa_safety_hits << " 次" << std::endl;
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
	 * @brief 构建 IMCF 并返回槽位对应的 taxid 列表。
	 *
	 * 该函数根据分组计划遍历每个 taxid 的 bottom-k 采样结果，并将哈希值直接插入 IMCF。
	 * @param imcf 已初始化的 IMCF 对象。
	 * @param groups shard 规划结果。
	 * @param hashCount 每个 taxid 对应的采样数量。
	 * @param sampledHashFiles bottom-k 采样出的哈希集合在磁盘上的临时文件。
	 */
	std::vector<std::vector<std::string>> buildIMCF(
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		const std::vector<chimera::imcf::Group>& groups,
		const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount,
		const robin_hood::unordered_flat_map<std::string, std::filesystem::path>& sampledHashFiles)
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
					std::cerr << "Warning: syncmer count mismatch for taxid " << taxid
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

		const auto shardEntriesView = std::span<const ShardEntry>(shardEntries.data(), shardEntries.size());

#pragma omp parallel for schedule(dynamic)
		for (size_t entryIdx = 0; entryIdx < shardEntriesView.size(); ++entryIdx) {
			const auto& entry = shardEntriesView[entryIdx];
			const std::string& taxid = entry.taxid;
			const auto& plans = *entry.plans;
			auto fileIt = sampledHashFiles.find(taxid);
			if (fileIt == sampledHashFiles.end()) {
#pragma omp critical(imcf_log)
				std::cerr << "Warning: 无采样数据可用于 taxid " << taxid << std::endl;
				continue;
			}
			const fs::path& hashFile = fileIt->second;
			std::ifstream hashStream(hashFile, std::ios::binary);
			if (!hashStream.is_open()) {
#pragma omp critical(imcf_log)
				std::cerr << "Warning: 无法打开 taxid " << taxid << " 的 syncmer 临时文件: " << hashFile << std::endl;
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

			bool overflow = false;
			std::vector<uint64_t> ioBuffer(1u << 14);
			while (!overflow && hashStream) {
				hashStream.read(reinterpret_cast<char*>(ioBuffer.data()),
					static_cast<std::streamsize>(ioBuffer.size() * sizeof(uint64_t)));
				std::streamsize bytesRead = hashStream.gcount();
				if (bytesRead <= 0) {
					break;
				}
				size_t count = static_cast<size_t>(bytesRead) / sizeof(uint64_t);
				for (size_t idx = 0; idx < count; ++idx) {
					uint64_t hash = ioBuffer[idx];
					size_t shard = pickShard(hash);
					if (shard >= plans.size()) {
#pragma omp critical(imcf_log)
						std::cerr << "Warning: taxid " << taxid
							<< " produced more syncmers than expected; extra entries ignored" << std::endl;
						overflow = true;
						break;
					}
					slotBuffers[shard].push_back(hash);
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
			}
			hashStream.close();
			std::error_code removeEc;
			fs::remove(hashFile, removeEc);

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
		uint32_t magic = chimera::imcf::LayoutHeaderV2::Magic;
		if (!os.write(reinterpret_cast<const char *>(&magic), sizeof(magic))) {
			os.close();
			std::cerr << "写入 IMCF 布局头失败" << std::endl;
			return false;
		}
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
		robin_hood::unordered_flat_map<std::string, std::filesystem::path> sampledHashFiles;
		robin_hood::unordered_flat_map<std::string, std::vector<std::string>> inputFiles;
		parseInputFile(config.input_file, inputFiles, hashCount, fileInfo);
		auto read_end = std::chrono::high_resolution_clock::now();
		auto read_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(read_end - read_start).count();
		if (config.verbose) {
			std::cout << "Read time: ";
			print_build_time(read_total_time);
			std::cout << std::endl;
		}

		if (config.filter == "imcf")
		{
			auto calculate_start = std::chrono::high_resolution_clock::now();
			std::cout << "Calculating syncmers..." << std::endl;
			std::filesystem::path scratchDir;
			syncmer_count(config, inputFiles, hashCount, sampledHashFiles, fileInfo, scratchDir);
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
			imcfConfig.smerSize = config.smer_size;
			imcfConfig.syncmerPosition = config.syncmer_position;
			imcfConfig.seed64 = ChimeraBuild::adjust_seed(config.kmer_size);
			imcfConfig.fpSalt = IMCFConfig::DefaultFingerprintSalt;
			imcfConfig.hashVersion = IMCFConfig::CurrentHashVersion;
			chimera::imcf::InterleavedMergedCuckooFilter imcf(groups, imcfConfig);
			std::vector<std::vector<std::string>> indexToTaxid = buildIMCF(imcf, groups, hashCount, sampledHashFiles);
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
		saveIMCF(imcf, config.output_file, indexToTaxid, imcfConfig,
		        config.filter == "imcf");
			auto save_end = std::chrono::high_resolution_clock::now();
			auto save_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(save_end - save_start).count();
			if (config.verbose) {
				std::cout << "Save time: ";
				print_build_time(save_total_time);
				std::cout << std::endl;
			}

			if (!scratchDir.empty()) {
				std::error_code cleanupEc;
				fs::remove_all(scratchDir, cleanupEc);
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
		}
	}
}
