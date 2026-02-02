#include "ChimeraBuildCommon.hpp"

#include <dna4_traits.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <thread>

namespace ChimeraBuild {

static inline bool file_is_compressed(const std::filesystem::path& filepath)
	{
		std::string extension = filepath.extension().string();
		return extension == ".gz" || extension == ".bgzf" || extension == ".bz2";
	}

	std::string tmp_hash_path(const std::string &taxid, std::string_view suffix)
	{
		std::string path = "tmp/";
		path.append(taxid);
		path.append(suffix);
		return path;
	}

	std::string tmp_hash_thread_path(const std::string &taxid,
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
		HashFrequencyContext* hashFreqContext,
		robin_hood::unordered_flat_map<std::string, uint64_t>* bpCount)
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
		std::vector<std::atomic<uint64_t>> bpCounters(index);
		for (auto &bp : bpCounters) {
			bp.store(0, std::memory_order_relaxed);
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

			// Open the sequence file
			seqan3::sequence_file_input<raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq >> fin{ filename };
			std::vector<uint64_t> hashes;
			hashes.reserve(4096);

			for (auto &record : fin)
			{
					auto &seq = record.sequence();

					const uint64_t seq_len = seq.size();
					bpCounters[taxid_index].fetch_add(seq_len, std::memory_order_relaxed);
					if (seq_len < min_required)
					{
						localFileInfo.skippedSeqNum++;
						continue;
					}
					localFileInfo.sequenceNum++;
					localFileInfo.bpLength += seq_len;

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

		if (bpCount) {
			for (size_t idx = 0; idx < index_to_taxid.size(); ++idx) {
				const std::string &taxid = index_to_taxid[idx];
				(*bpCount)[taxid] = bpCounters[idx].load(std::memory_order_relaxed);
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

	chimera::presence::CoverageMeta compute_presence_meta(
	    const robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
	    std::string_view featureSuffix,
	    const HashFrequencyContext *hashFreqContext,
	    const robin_hood::unordered_flat_map<std::string, uint64_t> *bpCount,
	    uint16_t effectiveSpan,
	    uint16_t refReadLen,
	    uint32_t uniqueDegThreshold,
	    uint16_t threads) {
		chimera::presence::CoverageMeta meta;
		meta.unique_deg_threshold = std::max<uint32_t>(1, uniqueDegThreshold);
		meta.ref_read_length = refReadLen;
		meta.effective_span = effectiveSpan;
		if (hashCount.empty()) {
			return meta;
		}
		const bool hasSketch = (hashFreqContext != nullptr) &&
		                       hashFreqContext->enabled() &&
		                       static_cast<bool>(hashFreqContext->sketch);
		std::vector<std::string> taxids;
		taxids.reserve(hashCount.size());
		for (const auto &kv : hashCount) {
			taxids.push_back(kv.first);
		}

#pragma omp parallel
		{
			std::vector<chimera::presence::CoverageEntry> local;
			std::vector<uint64_t> buffer;
			buffer.resize(1u << 14);

#pragma omp for schedule(dynamic)
			for (size_t i = 0; i < taxids.size(); ++i) {
				const std::string &taxid = taxids[i];
				std::string path = tmp_hash_path(taxid, featureSuffix);
				std::ifstream is(path, std::ios::binary);
				if (!is.is_open()) {
#pragma omp critical(imcf_log)
					std::cerr << "Coverage meta: failed to open feature file "
					          << path << std::endl;
					continue;
				}
				uint64_t unique = 0;
				uint64_t total = 0;
				while (is) {
					is.read(reinterpret_cast<char *>(buffer.data()),
					        static_cast<std::streamsize>(buffer.size() *
					                                     sizeof(uint64_t)));
					std::streamsize bytes = is.gcount();
					if (bytes <= 0) {
						break;
					}
					size_t got = static_cast<size_t>(bytes) / sizeof(uint64_t);
					total += got;
					for (size_t k = 0; k < got; ++k) {
						uint32_t df_est = 1;
						if (hasSketch) {
							df_est = hashFreqContext->sketch->estimate(buffer[k]);
						}
						if (df_est <= meta.unique_deg_threshold) {
							++unique;
						}
					}
				}
				if (unique > total) {
					unique = total;
				}
				uint64_t genome_bp = 0;
				if (bpCount) {
					auto itbp = bpCount->find(taxid);
					if (itbp != bpCount->end()) {
						genome_bp = itbp->second;
					}
				}
				const uint16_t span_used = (effectiveSpan > 0) ? effectiveSpan : static_cast<uint16_t>(1);
				const double window = std::max<int64_t>(
				    1, static_cast<int64_t>(refReadLen) - static_cast<int64_t>(span_used) + 1);
				double density = 0.0;
				if (genome_bp > 0 && unique > 0) {
					density = static_cast<double>(unique) /
					          static_cast<double>(genome_bp);
				}
				double expected_ref = (density > 0.0) ? density * window : 0.0;
				chimera::presence::CoverageEntry entry{};
				entry.taxid = taxid;
				entry.unique_signatures = unique;
				entry.total_signatures = total;
				entry.genome_length = genome_bp;
				entry.unique_density = density;
				entry.expected_unique_per_ref_read = expected_ref;
				local.push_back(entry);
			}

#pragma omp critical
			{
				meta.entries.insert(meta.entries.end(), local.begin(),
				                    local.end());
			}
		}

		if (hasSketch && hashFreqContext->sketch) {
			meta.freq_model.depth = hashFreqContext->sketch->depth();
			meta.freq_model.width = hashFreqContext->sketch->width();
			meta.freq_model.quantile = hashFreqContext->quantile;
			meta.freq_model.stats = hashFreqContext->stats;
			meta.freq_model.total_hashes = hashFreqContext->passA_total_hashes.load(std::memory_order_relaxed);
			meta.freq_model.filtered_hashes = hashFreqContext->passB_filtered_hashes.load(std::memory_order_relaxed);
			hashFreqContext->sketch->exportCounts(meta.freq_model.counters);
		}

		return meta;
	}

	/**
	* Get the maximum value from the hashCount map.
	*
	* @param hashCount The map containing the values.
	* @return The maximum value.
	*/
	
}
