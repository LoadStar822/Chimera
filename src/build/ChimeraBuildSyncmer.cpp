#include "ChimeraBuildCommon.hpp"

#include <dna4_traits.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <sstream>
#include <thread>

namespace ChimeraBuild {

	std::string tmp_hash_path(const std::string &taxid, std::string_view suffix)
	{
		std::string path = "tmp/";
		path.append(taxid);
		path.append(suffix);
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
	const uint64_t rank_salt = feature_seed ^ 0xD1B54A32D192ED03ULL;
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

	const size_t feature_min_length = chimera::feature::min_required_length(feature_params);
	const size_t min_required = std::max<size_t>(config.min_length, feature_min_length);
	if (config.verbose) {
		std::cout << "Feature hash minimum read length: " << min_required << std::endl;
	}
	const bool freq_filter_enabled = (hashFreqContext != nullptr) && hashFreqContext->enabled();
	const uint32_t freq_threshold = freq_filter_enabled ? hashFreqContext->stats.df_high_threshold
		: std::numeric_limits<uint32_t>::max();
	auto hash_df_est = [&](uint64_t hash) -> uint32_t {
		if (!freq_filter_enabled) {
			return 1u;
		}
		uint32_t df_est = hashFreqContext->sketch->estimate(hash);
		return df_est ? df_est : 1u;
	};
	auto keep_hash = [&](uint64_t hash, uint32_t df_est) -> bool {
		if (!freq_filter_enabled) {
			return true;
		}
		hashFreqContext->passB_total_hashes.fetch_add(1, std::memory_order_relaxed);
		if (df_est >= freq_threshold) {
			hashFreqContext->passB_filtered_hashes.fetch_add(1, std::memory_order_relaxed);
			return false;
		}
		return true;
	};

	struct TaxonTask {
		const std::string* taxid;
		const std::vector<std::string>* files;
	};
	std::vector<TaxonTask> tasks;
	tasks.reserve(inputFiles.size());
	for (const auto& kv : inputFiles) {
		tasks.push_back({ &kv.first, &kv.second });
	}
	const uint64_t taxa_total = static_cast<uint64_t>(tasks.size());
	const double core_alpha = std::max(0.0, config.core_alpha);
	const double core_beta = std::max(0.0, config.core_beta);
	const size_t file_cap =
	    (config.taxid_file_cap > 0) ? static_cast<size_t>(config.taxid_file_cap) : 0u;
	const double sig_oversample = std::max(1.0, config.sig_oversample);
	const size_t sig_s_min = std::max<size_t>(1, config.sig_s_min);
	const size_t sig_s_max = std::max<size_t>(sig_s_min, config.sig_s_max);
	const size_t k_base = std::max<size_t>(1, config.k_base);
	const size_t k_min = std::max<size_t>(1, config.k_min);
	const uint64_t k_max = (config.max_hashes_per_taxid > 0)
	                           ? static_cast<uint64_t>(config.max_hashes_per_taxid)
	                           : 2'000'000ull;
	constexpr double kGenomeRefBp = 5'000'000.0;
	constexpr double kBuildIdfMax = 8.0; // keep in sync with classify default
	constexpr double kCoverageReadLen = 150.0;
	constexpr double kCoverageHitsSingle = 8.0;
	constexpr double kCoverageHitsMulti = 6.0;
	if (config.verbose) {
		std::cout << "Core-IDF coverage target: read_len=" << static_cast<int>(kCoverageReadLen)
		          << ", hits_single=" << kCoverageHitsSingle
		          << ", hits_multi=" << kCoverageHitsMulti
		          << std::endl;
	}

	std::vector<uint64_t> finalCounts(tasks.size(), 0);
	std::vector<uint64_t> finalBp(tasks.size(), 0);

#pragma omp parallel
	{
		FileInfo localFileInfo{};
		std::vector<uint64_t> hashes;
		hashes.reserve(4096);
		std::vector<const std::string*> file_ptrs;
		std::vector<const std::string*> files_used;
		std::vector<uint64_t> file_hashes;
		struct SigCand {
			uint32_t df;
			uint64_t rank;
			uint64_t hash;
		};
		auto worse_df = [](const SigCand& a, const SigCand& b) {
			if (a.df != b.df) {
				return a.df < b.df; // max-heap by df (worst on top)
			}
			if (a.rank != b.rank) {
				return a.rank < b.rank; // max-heap by rank (worst on top)
			}
			return a.hash < b.hash;
		};
		auto worse_rank = [](const SigCand& a, const SigCand& b) {
			if (a.rank != b.rank) {
				return a.rank < b.rank; // max-heap by rank (worst on top)
			}
			if (a.df != b.df) {
				return a.df < b.df;
			}
			return a.hash < b.hash;
		};
		std::vector<SigCand> rare_heap;
		std::vector<SigCand> cov_heap;
		robin_hood::unordered_flat_set<uint64_t> rare_seen;
		robin_hood::unordered_flat_set<uint64_t> cov_seen;
		robin_hood::unordered_flat_map<uint64_t, uint16_t> within_df;
		struct Scored {
			double score;
			uint64_t hash;
		};
		auto heap_comp = [](const Scored& a, const Scored& b) {
			if (a.score != b.score) {
				return a.score > b.score; // min-heap by score
			}
			return a.hash > b.hash;
		};
		std::vector<Scored> best;
		struct RankedHash {
			uint64_t rank;
			uint64_t hash;
		};
		auto rank_comp = [](const RankedHash& a, const RankedHash& b) {
			if (a.rank != b.rank) {
				return a.rank < b.rank; // max-heap by rank (worst on top)
			}
			return a.hash < b.hash;
		};
		std::vector<RankedHash> cov_best;
		robin_hood::unordered_flat_set<uint64_t> selected_set;
		std::vector<uint64_t> selected;

#pragma omp for schedule(dynamic)
		for (size_t idx = 0; idx < tasks.size(); ++idx)
		{
			const auto& task = tasks[idx];
			if (task.taxid == nullptr || task.files == nullptr || task.files->empty()) {
				continue;
			}
			const std::string& taxid = *task.taxid;
			const std::string out_path = tmp_hash_path(taxid, feature_suffix);

			std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
			out.tie(nullptr);
			if (!out.is_open())
			{
#pragma omp critical(imcf_log)
				std::cerr << "Unable to open the feature hash file: " << out_path << std::endl;
				finalCounts[idx] = 0;
				finalBp[idx] = 0;
				continue;
			}

			file_ptrs.clear();
			file_ptrs.reserve(task.files->size());
			for (const auto& filename : *task.files) {
				file_ptrs.push_back(&filename);
			}
			std::sort(file_ptrs.begin(), file_ptrs.end(),
			          [](const std::string* a, const std::string* b) {
				          return *a < *b;
			          });
			files_used.clear();
			if (file_cap == 0 || file_ptrs.size() <= file_cap) {
				files_used = file_ptrs;
			} else if (file_cap == 1) {
				files_used.push_back(file_ptrs.front());
			} else {
				files_used.reserve(file_cap);
				const size_t n = file_ptrs.size();
				for (size_t i = 0; i < file_cap; ++i) {
					size_t pick = (i * (n - 1)) / (file_cap - 1);
					files_used.push_back(file_ptrs[pick]);
				}
				std::sort(files_used.begin(), files_used.end(),
				          [](const std::string* a, const std::string* b) {
					          return *a < *b;
				          });
				files_used.erase(std::unique(files_used.begin(), files_used.end()), files_used.end());
			}
			const size_t f_used = files_used.size();
			if (f_used == 0) {
				out.close();
				finalCounts[idx] = 0;
				finalBp[idx] = 0;
				continue;
			}

			const bool exact_within_df = (f_used <= 4);
			const size_t sig_size_raw = static_cast<size_t>(
			    std::ceil((static_cast<double>(k_base) * sig_oversample) / static_cast<double>(f_used)));
			const size_t sig_size = std::clamp(sig_size_raw, sig_s_min, sig_s_max);

			within_df.clear();
			within_df.reserve(static_cast<size_t>(sig_size * std::min<size_t>(f_used, 64)));

			uint64_t bp_sum = 0;
			for (const std::string* filename_ptr : files_used)
			{
				if (filename_ptr == nullptr) {
					continue;
				}
				const std::string& filename = *filename_ptr;
				try
				{
					seqan3::sequence_file_input<raptor::dna4_traits,
					                            seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{ filename };
					const size_t sig_rare = exact_within_df ? 0 : (sig_size / 4);
					const size_t sig_cov = exact_within_df ? 0 : (sig_size - sig_rare);
					if (exact_within_df) {
						file_hashes.clear();
					} else {
						rare_heap.clear();
						cov_heap.clear();
						rare_seen.clear();
						cov_seen.clear();
						rare_heap.reserve(sig_rare);
						cov_heap.reserve(sig_cov);
						rare_seen.reserve(sig_rare * 2 + 8);
						cov_seen.reserve(sig_cov * 2 + 8);
					}
					for (auto& record : fin)
					{
						auto& seq = record.sequence();
						const uint64_t seq_len = seq.size();
						bp_sum += seq_len;
						if (seq_len < min_required)
						{
							localFileInfo.skippedSeqNum++;
							continue;
						}
						localFileInfo.sequenceNum++;
						localFileInfo.bpLength += seq_len;

						hashes.clear();
						chimera::feature::compute_hashes_append(seq, feature_params, hashes);
						if (hashes.empty()) {
							continue;
						}

						for (uint64_t hash : hashes)
						{
							const uint32_t df_est = hash_df_est(hash);
							if (!keep_hash(hash, df_est)) {
								continue;
							}
							if (exact_within_df) {
								file_hashes.push_back(hash);
								continue;
							}
							const uint64_t rank = mix64(hash ^ rank_salt);
							if (sig_rare > 0 && !rare_seen.contains(hash)) {
								SigCand cand{df_est, rank, hash};
								if (rare_heap.size() < sig_rare) {
									rare_heap.push_back(cand);
									std::push_heap(rare_heap.begin(), rare_heap.end(), worse_df);
									rare_seen.insert(hash);
								} else {
									const SigCand& worst = rare_heap.front();
									if (df_est < worst.df || (df_est == worst.df && rank < worst.rank)) {
										std::pop_heap(rare_heap.begin(), rare_heap.end(), worse_df);
										rare_seen.erase(rare_heap.back().hash);
										rare_heap.back() = cand;
										std::push_heap(rare_heap.begin(), rare_heap.end(), worse_df);
										rare_seen.insert(hash);
									}
								}
							}
							if (sig_cov > 0 && !cov_seen.contains(hash)) {
								SigCand cand{df_est, rank, hash};
								if (cov_heap.size() < sig_cov) {
									cov_heap.push_back(cand);
									std::push_heap(cov_heap.begin(), cov_heap.end(), worse_rank);
									cov_seen.insert(hash);
								} else {
									const SigCand& worst = cov_heap.front();
									if (rank < worst.rank || (rank == worst.rank && df_est < worst.df)) {
										std::pop_heap(cov_heap.begin(), cov_heap.end(), worse_rank);
										cov_seen.erase(cov_heap.back().hash);
										cov_heap.back() = cand;
										std::push_heap(cov_heap.begin(), cov_heap.end(), worse_rank);
										cov_seen.insert(hash);
									}
								}
							}
						}
					}
					if (exact_within_df) {
						if (file_hashes.empty()) {
							continue;
						}
						std::sort(file_hashes.begin(), file_hashes.end());
						file_hashes.erase(std::unique(file_hashes.begin(), file_hashes.end()), file_hashes.end());
						for (uint64_t h : file_hashes) {
							auto it = within_df.find(h);
							if (it == within_df.end()) {
								within_df.emplace(h, static_cast<uint16_t>(1));
							} else if (it->second < std::numeric_limits<uint16_t>::max()) {
								++it->second;
							}
						}
					} else {
						if (rare_heap.empty() && cov_heap.empty()) {
							continue;
						}
						file_hashes.clear();
						file_hashes.reserve(rare_heap.size() + cov_heap.size());
						for (const auto& cand : rare_heap) {
							file_hashes.push_back(cand.hash);
						}
						for (const auto& cand : cov_heap) {
							file_hashes.push_back(cand.hash);
						}
						std::sort(file_hashes.begin(), file_hashes.end());
						file_hashes.erase(std::unique(file_hashes.begin(), file_hashes.end()), file_hashes.end());
						for (uint64_t h : file_hashes) {
							auto it = within_df.find(h);
							if (it == within_df.end()) {
								within_df.emplace(h, static_cast<uint16_t>(1));
							} else if (it->second < std::numeric_limits<uint16_t>::max()) {
								++it->second;
							}
						}
					}
				}
				catch (const std::exception& ex)
				{
#pragma omp critical(imcf_log)
					{
						std::cerr << "Failed to read sequence file: " << filename
						          << " taxid=" << taxid
						          << " (" << ex.what() << ")" << std::endl;
					}
					localFileInfo.invalidNum++;
				}
			}

			uint64_t rep_len = (f_used > 0) ? (bp_sum / static_cast<uint64_t>(f_used)) : bp_sum;
			if (rep_len == 0) {
				rep_len = bp_sum;
			}

			size_t k_target = k_min;
			if (rep_len > 0) {
				const double l_rep = static_cast<double>(rep_len);
				const double scaled = static_cast<double>(k_base) * (l_rep / kGenomeRefBp);
				if (std::isfinite(scaled) && scaled > 0.0) {
					k_target = static_cast<size_t>(std::llround(scaled));
				}
			}
			k_target = std::clamp<size_t>(k_target, k_min, static_cast<size_t>(k_max));
			if (within_df.size() < k_target) {
				k_target = within_df.size();
			}

			const double cov_hits = (f_used <= 1) ? kCoverageHitsSingle : kCoverageHitsMulti;
			size_t k_cov_min = 0;
			if (rep_len > 0 && cov_hits > 0.0) {
				const double expect = (static_cast<double>(rep_len) * cov_hits) / kCoverageReadLen;
				if (std::isfinite(expect) && expect > 0.0) {
					k_cov_min = static_cast<size_t>(std::ceil(expect));
				}
			}
			if (k_cov_min > k_target) {
				k_target = std::min<size_t>(k_cov_min, static_cast<size_t>(k_max));
			}
			if (within_df.size() < k_target) {
				k_target = within_df.size();
			}
			const size_t k_cov = std::min<size_t>(k_target, k_cov_min);
			const size_t k_score = (k_target > k_cov) ? (k_target - k_cov) : 0;

			best.clear();
			best.reserve(k_score);
			selected.clear();
			selected.reserve(k_target);

			if (k_score > 0 && !within_df.empty()) {
				for (const auto& kv : within_df) {
					const uint64_t h = kv.first;
					const uint32_t df_est = hash_df_est(h);
					const double df = static_cast<double>(std::max<uint32_t>(1u, df_est));
					const double n_taxa = static_cast<double>(taxa_total);
					double idf = 0.0;
					if (n_taxa > 0.0) {
						idf = std::log2((n_taxa + 1.0) / (df + 1.0));
					}
					idf = std::clamp(idf, 0.0, kBuildIdfMax);
					const double numerator = static_cast<double>(kv.second) + core_alpha;
					const double denom = static_cast<double>(f_used) + core_alpha;
					double core = 0.0;
					if (denom > 0.0) {
						core = numerator / denom;
					}
					if (core_beta != 1.0) {
						core = std::pow(core, core_beta);
					}
					double score = idf * core;
					if (!std::isfinite(score)) {
						continue;
					}
					if (score < 0.0) {
						score = 0.0;
					}
					if (best.size() < k_score) {
						best.push_back({score, h});
						std::push_heap(best.begin(), best.end(), heap_comp);
						continue;
					}
					const Scored& worst = best.front();
					if (score > worst.score || (score == worst.score && h > worst.hash)) {
						std::pop_heap(best.begin(), best.end(), heap_comp);
						best.back() = {score, h};
						std::push_heap(best.begin(), best.end(), heap_comp);
					}
				}
				selected.reserve(best.size());
				for (const auto& item : best) {
					selected.push_back(item.hash);
				}
			}

			if (k_cov > 0 && !within_df.empty()) {
				selected_set.clear();
				selected_set.reserve(selected.size() * 2 + 8);
				for (uint64_t h : selected) {
					selected_set.insert(h);
				}
				cov_best.clear();
				cov_best.reserve(k_cov);
				for (const auto& kv : within_df) {
					const uint64_t h = kv.first;
					if (selected_set.contains(h)) {
						continue;
					}
					const uint64_t rank = mix64(h ^ rank_salt);
					if (cov_best.size() < k_cov) {
						cov_best.push_back({rank, h});
						std::push_heap(cov_best.begin(), cov_best.end(), rank_comp);
						continue;
					}
					const RankedHash& worst = cov_best.front();
					if (rank < worst.rank || (rank == worst.rank && h > worst.hash)) {
						std::pop_heap(cov_best.begin(), cov_best.end(), rank_comp);
						cov_best.back() = {rank, h};
						std::push_heap(cov_best.begin(), cov_best.end(), rank_comp);
					}
				}
				for (const auto& item : cov_best) {
					selected.push_back(item.hash);
				}
			}

			if (!selected.empty()) {
				std::sort(selected.begin(), selected.end());
				selected.erase(std::unique(selected.begin(), selected.end()), selected.end());
			}

			if (config.verbose) {
				const uint64_t tid_key =
				    static_cast<uint64_t>(std::hash<std::string>{}(taxid));
				const bool log_sample = ((mix64(tid_key) % 1000ull) == 0ull);
				if (log_sample) {
					uint16_t max_within = 0;
					uint64_t ones = 0;
					for (const auto& kv : within_df) {
						if (kv.second > max_within) {
							max_within = kv.second;
						}
						if (kv.second == 1) {
							++ones;
						}
					}
					double pct_one = 0.0;
					if (!within_df.empty()) {
						pct_one = 100.0 * static_cast<double>(ones) /
						          static_cast<double>(within_df.size());
					}
					std::vector<uint32_t> df_vals;
					df_vals.reserve(selected.size());
					for (uint64_t h : selected) {
						df_vals.push_back(hash_df_est(h));
					}
					std::sort(df_vals.begin(), df_vals.end());
					auto pick = [&](double q) -> uint32_t {
						if (df_vals.empty()) {
							return 0;
						}
						q = std::clamp(q, 0.0, 1.0);
						size_t pos = static_cast<size_t>(
						    std::floor(q * static_cast<double>(df_vals.size() - 1)));
						if (pos >= df_vals.size()) {
							pos = df_vals.size() - 1;
						}
						return df_vals[pos];
					};
					std::ostringstream ss;
					ss << "[core-idf] taxid=" << taxid << " files=" << f_used
					   << " sig=" << sig_size << " cand=" << within_df.size()
					   << " max_within=" << max_within << " within_df==1="
					   << std::fixed << std::setprecision(1) << pct_one << "%"
					   << std::defaultfloat << " sel=" << selected.size()
					   << " df_p50=" << pick(0.50) << " df_p90=" << pick(0.90)
					   << " df_p99=" << pick(0.99) << '\n';
#pragma omp critical(imcf_log)
					std::cout << ss.str();
				}
			}

			if (!selected.empty()) {
				out.write(reinterpret_cast<const char*>(selected.data()),
				          static_cast<std::streamsize>(selected.size() * sizeof(uint64_t)));
			}
			if (!out.good())
			{
#pragma omp critical(imcf_log)
				std::cerr << "Failed to write feature hashes for taxid " << taxid << std::endl;
				selected.clear();
			}

			finalCounts[idx] = static_cast<uint64_t>(selected.size());
			finalBp[idx] = rep_len;
		}

#pragma omp critical
		{
			fileInfo.invalidNum += localFileInfo.invalidNum;
			fileInfo.skippedNum += localFileInfo.skippedNum;
			fileInfo.skippedSeqNum += localFileInfo.skippedSeqNum;
			fileInfo.sequenceNum += localFileInfo.sequenceNum;
			fileInfo.bpLength += localFileInfo.bpLength;
		}
	}

	for (size_t idx = 0; idx < tasks.size(); ++idx)
	{
		const auto& task = tasks[idx];
		if (!task.taxid) {
			continue;
		}
		const std::string& taxid = *task.taxid;
		auto it = hashCount.find(taxid);
		if (it != hashCount.end()) {
			it->second = finalCounts[idx];
		} else {
			hashCount.emplace(taxid, finalCounts[idx]);
		}
		if (bpCount) {
			(*bpCount)[taxid] = finalBp[idx];
		}
	}
	}

} // namespace ChimeraBuild
