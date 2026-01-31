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
		robin_hood::unordered_flat_map<std::string, uint64_t> shardOffset;
		shardOffset.reserve(hashCount.size());

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
				uint64_t count = group.assignedHashes[slot];
				if (count == 0) {
					continue;
				}
				auto& plans = shardPlan[group.taxids[slot]];
				uint64_t offset = 0;
				auto itOff = shardOffset.find(group.taxids[slot]);
				if (itOff != shardOffset.end()) {
					offset = itOff->second;
				}
				plans.push_back({ groupIdx, slot, count, offset });
				shardOffset[group.taxids[slot]] = offset + count;
			}
		}

		const double kShardToleranceRatio = 0.001; // 0.1%
		for (const auto& [taxid, plans] : shardPlan) {
			uint64_t shardTotal = 0;
			for (const auto& plan : plans) {
				shardTotal += plan.count;
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
		size_t groupIndex;
		size_t slotIndex;
		uint64_t fileOffset;
		uint64_t nHashes;
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
		const auto& plans = kv.second;
		auto it = hashCount.find(taxid);
		uint64_t totalHashes = (it == hashCount.end()) ? 0 : it->second;
		if (totalHashes == 0) {
			continue;
		}
		uint64_t chunkSize = std::max<uint64_t>(
		    targetChunk,
		    threadCount > 0 ? (totalHashes / (static_cast<uint64_t>(threadCount) * 4ull)) : targetChunk);
		chunkSize = std::max<uint64_t>(chunkSize, targetChunk);

		for (const auto& plan : plans) {
			if (plan.count == 0) {
				continue;
			}
			uint64_t start = plan.fileOffsetStart;
			uint64_t remaining = plan.count;
			while (remaining > 0) {
				uint64_t len = std::min<uint64_t>(chunkSize, remaining);
				shardEntries.push_back({taxid, plan.groupIndex, plan.slotIndex, start, len});
				start += len;
				remaining -= len;
			}
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

	if (shardEntries.empty()) {
		return indexToTaxid;
	}

	constexpr size_t kReadBlockBytes = 4 * 1024 * 1024; // 4MB
	using InsertionPair = std::pair<uint64_t, uint16_t>;

	const size_t kDefaultBatches = 10;
	const uint64_t kMemoryBudgetBytes = 15ull * 1024ull * 1024ull * 1024ull; // ~15GB target for batch buffers
	const uint64_t bytesPerPair = sizeof(InsertionPair);
	uint64_t totalShardHashes = 0;
	for (const auto& entry : shardEntries) {
		totalShardHashes += entry.nHashes;
	}
	uint64_t minBatches = 1;
	if (bytesPerPair > 0 && kMemoryBudgetBytes > 0) {
		minBatches = (totalShardHashes * bytesPerPair + kMemoryBudgetBytes - 1) / kMemoryBudgetBytes;
		if (minBatches == 0) {
			minBatches = 1;
		}
	}
	size_t numBatches = std::max<size_t>(kDefaultBatches, static_cast<size_t>(minBatches));
	numBatches = std::min<size_t>(numBatches, shardEntries.size());
	if (numBatches == 0) {
		numBatches = 1;
	}
	uint64_t targetBatchSize = (totalShardHashes + numBatches - 1) / numBatches;
	std::vector<std::pair<size_t, size_t>> batchRanges;
	batchRanges.reserve(numBatches);
	size_t batchStart = 0;
	uint64_t batchAccum = 0;
	for (size_t idx = 0; idx < shardEntries.size(); ++idx) {
		batchAccum += shardEntries[idx].nHashes;
		bool split = (batchRanges.size() + 1 < numBatches) && (batchAccum >= targetBatchSize);
		if (split) {
			batchRanges.emplace_back(batchStart, idx + 1);
			batchStart = idx + 1;
			batchAccum = 0;
		}
	}
	batchRanges.emplace_back(batchStart, shardEntries.size());

	std::unique_ptr<std::atomic<uint64_t>[]> globalUnique;
	std::unique_ptr<std::atomic<uint64_t>[]> globalTotal;
	robin_hood::unordered_flat_map<std::string, size_t> taxidToIndex;
	if (collectCoverage) {
		size_t idx = 0;
		taxidToIndex.reserve(shardPlan.size());
		for (const auto& kvPlan : shardPlan) {
			taxidToIndex[kvPlan.first] = idx++;
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

	const auto shardEntriesView = std::span<const ShardEntry>(shardEntries.data(), shardEntries.size());

	const int numThreads = static_cast<int>(threadCount == 0 ? 1 : threadCount);
	std::vector<std::vector<std::vector<InsertionPair>>> threadBuffers(static_cast<size_t>(numThreads));

	for (const auto& range : batchRanges) {
		size_t batchStartIdx = range.first;
		size_t batchEndIdx = range.second;

		for (auto& tb : threadBuffers) {
			tb.resize(groups.size());
			for (auto& buf : tb) {
				buf.clear();
			}
		}

	#pragma omp parallel
		{
			int tid = 0;
#ifdef _OPENMP
			tid = omp_get_thread_num();
#endif
			auto& myBuffers = threadBuffers[static_cast<size_t>(tid)];
			if (myBuffers.size() < groups.size()) {
				myBuffers.resize(groups.size());
			}
			std::vector<uint64_t> readBlock(kReadBlockBytes / sizeof(uint64_t));

		#pragma omp for schedule(guided, 4)
			for (size_t entryIdx = batchStartIdx; entryIdx < batchEndIdx; ++entryIdx) {
				const auto& entry = shardEntriesView[entryIdx];
				const std::string& taxid = entry.taxid;
				std::ifstream ifile(tmp_hash_path(taxid, suffix), std::ios::binary);
				if (!ifile) {
#pragma omp critical(imcf_log)
					std::cerr << "Failed to open feature hash file: " << taxid << suffix << std::endl;
					continue;
				}

				ifile.seekg(static_cast<std::streamoff>(entry.fileOffset * sizeof(uint64_t)), std::ios::beg);
				uint64_t remaining = entry.nHashes;
				uint64_t uniqueCount = 0;
				uint64_t totalCount = 0;

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
						myBuffers[entry.groupIndex].emplace_back(hash, static_cast<uint16_t>(entry.slotIndex));
					}
					if (hashesRead > remaining) {
						remaining = 0;
					} else {
						remaining -= hashesRead;
					}
				}

				if (remaining != 0) {
#pragma omp critical(imcf_log)
					std::cerr << "Warning: shard read underflow for taxid " << taxid
						<< ", expected " << entry.nHashes << " hashes, missing " << remaining << std::endl;
				}

				if (!ifile.eof() && ifile.fail()) {
#pragma omp critical(imcf_log)
					std::cerr << "Warning: failed to read syncmer file completely for taxid " << taxid << std::endl;
				}
				ifile.close();

				if (collectCoverage) {
					auto itIdx = taxidToIndex.find(taxid);
					if (itIdx != taxidToIndex.end()) {
						size_t aggIdx = itIdx->second;
						globalUnique[aggIdx].fetch_add(uniqueCount, std::memory_order_relaxed);
						globalTotal[aggIdx].fetch_add(totalCount, std::memory_order_relaxed);
					}
				}
			}
		}

	#pragma omp parallel for schedule(dynamic)
		for (size_t groupIdx = 0; groupIdx < groups.size(); ++groupIdx) {
			size_t totalSize = 0;
			for (const auto& tb : threadBuffers) {
				totalSize += tb[groupIdx].size();
			}
			if (totalSize == 0) {
				continue;
			}
			std::vector<InsertionPair> aggregatedBuffer;
			aggregatedBuffer.reserve(totalSize);
			for (auto& tb : threadBuffers) {
				auto& buf = tb[groupIdx];
				aggregatedBuffer.insert(aggregatedBuffer.end(),
				                      std::make_move_iterator(buf.begin()),
				                      std::make_move_iterator(buf.end()));
				std::vector<InsertionPair>().swap(buf); // release memory for this group asap
			}

			uint64_t success = 0;
			uint64_t failure = 0;
			for (const auto& pair : aggregatedBuffer) {
				if (imcf.insertTag(groupIdx, pair.first, pair.second)) {
					++success;
					++slotSuccess[groupIdx][static_cast<size_t>(pair.second)];
				} else {
					++failure;
					++slotFailures[groupIdx][static_cast<size_t>(pair.second)];
				}
			}
			groupSuccess[groupIdx] += success;
			groupFailures[groupIdx] += failure;
			std::vector<InsertionPair>().swap(aggregatedBuffer);
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
