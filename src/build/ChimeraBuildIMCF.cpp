#include "ChimeraBuildCommon.hpp"

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <limits>
#include <mutex>
#include <span>
#include <sstream>
#include <thread>

namespace ChimeraBuild {

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
		std::string_view featureSuffix,
		const HashFrequencyContext* hashFreqContext,
		const robin_hood::unordered_flat_map<std::string, uint64_t>* bpCount,
		uint16_t effectiveSpan,
		uint16_t refReadLen,
		uint32_t uniqueDegThreshold,
		chimera::presence::CoverageMeta* coverageMeta)
	{
		std::vector<std::vector<std::string>> indexToTaxid(groups.size());
		const std::string suffix(featureSuffix);
		robin_hood::unordered_flat_map<std::string, std::vector<TaxidShardPlan>> shardPlan;
		shardPlan.reserve(hashCount.size());

		const bool collectCoverage = coverageMeta != nullptr;
		const uint32_t uniqueDeg =
		    std::max<uint32_t>(1, uniqueDegThreshold);
		const uint16_t spanUsed =
		    (effectiveSpan > 0) ? effectiveSpan : static_cast<uint16_t>(1);
		const uint16_t refLenUsed =
		    (refReadLen == 0) ? static_cast<uint16_t>(150) : refReadLen;
		const bool hasSketch =
		    collectCoverage && hashFreqContext && hashFreqContext->enabled() &&
		    static_cast<bool>(hashFreqContext->sketch);
		if (collectCoverage) {
			coverageMeta->unique_deg_threshold = uniqueDeg;
			coverageMeta->ref_read_length = refLenUsed;
			coverageMeta->effective_span = spanUsed;
			coverageMeta->entries.clear();
		}

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
		uint64_t startHash{0};
		uint64_t nHashes{0};
	};
	std::vector<ShardEntry> shardEntries;
	shardEntries.reserve(shardPlan.size() * 4);
#ifdef _OPENMP
	const uint16_t threadCount = static_cast<uint16_t>(omp_get_max_threads());
#else
	const uint16_t threadCount = static_cast<uint16_t>(std::thread::hardware_concurrency());
#endif
	const uint64_t targetChunk = 1'000'000ull; // 1e6 hashes per chunk baseline

	for (auto& kv : shardPlan) {
		const std::string& taxid = kv.first;
		const auto* plans = &kv.second;
		auto it = hashCount.find(taxid);
		uint64_t totalHashes = (it == hashCount.end()) ? 0 : it->second;
		if (totalHashes == 0) {
			continue;
		}
		uint64_t chunkSize = std::max<uint64_t>(
		    targetChunk,
		    threadCount > 0 ? (totalHashes / (static_cast<uint64_t>(threadCount) * 4ull)) : targetChunk);
		chunkSize = std::max<uint64_t>(chunkSize, targetChunk);
		for (uint64_t start = 0; start < totalHashes; start += chunkSize) {
			uint64_t len = std::min<uint64_t>(chunkSize, totalHashes - start);
			shardEntries.push_back({taxid, plans, start, len});
		}
	}
	// 优先分配重负载 taxid，避免长尾串行
	std::sort(shardEntries.begin(), shardEntries.end(),
		      [&](const ShardEntry& a, const ShardEntry& b) {
			          if (a.nHashes == b.nHashes) {
				          return a.taxid < b.taxid;
			          }
			          return a.nHashes > b.nHashes;
		          });

	std::vector<std::mutex> groupLocks(groups.size());
	std::unique_ptr<std::atomic<uint64_t>[]> globalUnique;
	std::unique_ptr<std::atomic<uint64_t>[]> globalTotal;
	robin_hood::unordered_flat_map<std::string, size_t> taxidToIndex;
	if (collectCoverage) {
		size_t idx = 0;
		taxidToIndex.reserve(shardPlan.size());
		for (const auto& kv : shardPlan) {
			taxidToIndex[kv.first] = idx++;
		}
		size_t n = taxidToIndex.size();
		globalUnique.reset(new std::atomic<uint64_t>[n]);
		globalTotal.reset(new std::atomic<uint64_t>[n]);
		for (size_t i = 0; i < n; ++i) {
			globalUnique[i].store(0, std::memory_order_relaxed);
			globalTotal[i].store(0, std::memory_order_relaxed);
		}
	}

		std::vector<uint64_t> groupSuccess(groups.size(), 0);
		std::vector<uint64_t> groupFailures(groups.size(), 0);
		std::vector<std::vector<uint64_t>> slotSuccess(groups.size());
		std::vector<std::vector<uint64_t>> slotFailures(groups.size());
		for (size_t g = 0; g < groups.size(); ++g) {
			size_t slotCount = groups[g].taxids.size();
			slotSuccess[g].assign(slotCount, 0);
			slotFailures[g].assign(slotCount, 0);
		}

		constexpr size_t kFlushThreshold = 131072;     // 128K hashes
		constexpr size_t kReadBlockBytes = 4 * 1024 * 1024; // 4MB

		const auto shardEntriesView = std::span<const ShardEntry>(shardEntries.data(), shardEntries.size());

	#pragma omp parallel
		{
#pragma omp for schedule(dynamic, 1)
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

				const uint64_t totalWeight = [&]() {
					uint64_t sum = 0;
					for (const auto& p : plans) {
						sum += p.expectedCount;
					}
					return sum;
				}();
				if (totalWeight == 0) {
					continue;
				}

				std::vector<std::vector<uint64_t>> slotBuffers(plans.size());
				for (size_t i = 0; i < plans.size(); ++i) {
					uint64_t expect = plans[i].expectedCount;
					size_t reserveSize = static_cast<size_t>(std::min<uint64_t>(expect, static_cast<uint64_t>(kFlushThreshold)));
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

				auto pickShard = [&](uint64_t hash) -> size_t {
					uint64_t r = mix64(hash) % totalWeight;
					for (size_t idx = 0; idx < plans.size(); ++idx) {
						uint64_t quota = plans[idx].expectedCount;
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
				uint64_t uniqueCount = 0;
				uint64_t totalCount = 0;

				uint64_t remaining = entry.nHashes;
				if (entry.startHash > 0) {
					ifile.seekg(static_cast<std::streamoff>(entry.startHash * sizeof(uint64_t)),
					            std::ios::beg);
				}

				while (ifile && remaining > 0) {
					uint64_t toRead = std::min<uint64_t>(remaining,
					                                     static_cast<uint64_t>(readBlock.size()));
					ifile.read(reinterpret_cast<char*>(readBlock.data()),
					           static_cast<std::streamsize>(toRead * sizeof(uint64_t)));
					std::streamsize bytesRead = ifile.gcount();
					if (bytesRead <= 0) {
						break;
					}
					size_t hashesRead = static_cast<size_t>(bytesRead) / sizeof(uint64_t);
					if (hashesRead == 0) {
						break;
					}
					bool overflow = false;
					for (size_t h = 0; h < hashesRead; ++h) {
						uint64_t hash = readBlock[h];
						if (collectCoverage) {
							++totalCount;
							uint32_t df_est = 1;
							if (hasSketch) {
								df_est = hashFreqContext->sketch->estimate(hash);
							}
							if (df_est <= uniqueDeg) {
								++uniqueCount;
							}
						}
						size_t shard = pickShard(hash);
						if (shard >= plans.size()) {
							// Extra syncmers beyond planned quota
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
					}
					if (overflow) {
						break;
					}
					if (hashesRead > remaining) {
						remaining = 0;
					} else {
						remaining -= hashesRead;
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

				if (collectCoverage) {
					if (uniqueCount > totalCount) {
						uniqueCount = totalCount;
					}
					chimera::presence::CoverageEntry entry{};
					entry.taxid = taxid;
					entry.unique_signatures = uniqueCount;
					entry.total_signatures = totalCount;
					uint64_t genome_bp = 0;
					if (bpCount) {
						auto itbp = bpCount->find(taxid);
						if (itbp != bpCount->end()) {
							genome_bp = itbp->second;
						}
					}
					auto itIdx = taxidToIndex.find(taxid);
					if (itIdx != taxidToIndex.end()) {
						size_t aggIdx = itIdx->second;
						globalUnique[aggIdx].fetch_add(uniqueCount, std::memory_order_relaxed);
						globalTotal[aggIdx].fetch_add(totalCount, std::memory_order_relaxed);
					}
				}
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

		if (collectCoverage) {
			coverageMeta->entries.clear();
			coverageMeta->entries.reserve(taxidToIndex.size());
			for (const auto& kv : taxidToIndex) {
				const std::string& taxid = kv.first;
				size_t idx = kv.second;
				chimera::presence::CoverageEntry entry{};
				entry.taxid = taxid;
				entry.unique_signatures = globalUnique[idx].load(std::memory_order_relaxed);
				entry.total_signatures = globalTotal[idx].load(std::memory_order_relaxed);
				uint64_t genome_bp = 0;
				if (bpCount) {
					auto itbp = bpCount->find(taxid);
					if (itbp != bpCount->end()) {
						genome_bp = itbp->second;
					}
				}
				entry.genome_length = genome_bp;
				const double window = std::max<int64_t>(
				    1, static_cast<int64_t>(refLenUsed) - static_cast<int64_t>(spanUsed) + 1);
				entry.unique_density = (genome_bp > 0 && entry.unique_signatures > 0)
					? static_cast<double>(entry.unique_signatures) / static_cast<double>(genome_bp)
					: 0.0;
				entry.expected_unique_per_ref_read =
					(entry.unique_density > 0.0) ? entry.unique_density * window : 0.0;
				coverageMeta->entries.push_back(entry);
			}
			if (hasSketch) {
				coverageMeta->freq_model.depth = hashFreqContext->sketch->depth();
				coverageMeta->freq_model.width = hashFreqContext->sketch->width();
				coverageMeta->freq_model.quantile = hashFreqContext->quantile;
				coverageMeta->freq_model.stats = hashFreqContext->stats;
				coverageMeta->freq_model.total_hashes =
				    hashFreqContext->passA_total_hashes.load(std::memory_order_relaxed);
				coverageMeta->freq_model.filtered_hashes =
				    hashFreqContext->passB_filtered_hashes.load(std::memory_order_relaxed);
				hashFreqContext->sketch->exportCounts(
				    coverageMeta->freq_model.counters);
			}
		}

		return indexToTaxid;
	}

	void saveIMCF(chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		const std::string& output_file,
		std::vector<std::vector<std::string>>& indexToTaxid,
		IMCFConfig& imcfConfig,
		bool needIndexStructures,
		const chimera::presence::CoverageMeta* presenceMeta)
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
		if (presenceMeta) {
			archive(*presenceMeta);
		}
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

}
