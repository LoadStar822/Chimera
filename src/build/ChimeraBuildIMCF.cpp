#include "ChimeraBuildCommon.hpp"

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <algorithm>
#include <atomic>
#include <cstdlib>
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

	namespace {
	struct SpoolSegmentRef {
		uint16_t threadId{0};
		uint64_t spoolOffset{0};
		uint64_t count{0};
	};

	struct ShardEntry {
		uint32_t taxidId{0};
		size_t slotIndex{0};
		uint64_t fileOffset{0};
		uint64_t nHashes{0};
		std::vector<SpoolSegmentRef> segments;
	};

	constexpr size_t kMinReadBlockBytes = 1 * 1024 * 1024;
	constexpr size_t kMaxReadBlockBytes = 4 * 1024 * 1024;
	constexpr size_t kTargetTotalReadBlockBytes =
	    256ull * 1024ull * 1024ull;

	static size_t resolve_read_block_bytes(size_t workerCount)
	{
		return std::clamp(
		    kTargetTotalReadBlockBytes / std::max<size_t>(1, workerCount),
		    kMinReadBlockBytes, kMaxReadBlockBytes);
	}

	static bool map_virtual_range_to_chunks(const TaxidFeatureLayout& layout,
	                                       uint64_t virtualOffset,
	                                       uint64_t count,
	                                       std::vector<SpoolSegmentRef>& out,
	                                       size_t& chunkCursor,
	                                       uint64_t& chunkBase)
	{
		out.clear();
		uint64_t remaining = count;
		if (remaining == 0) {
			return true;
		}
		if (virtualOffset < chunkBase) {
			chunkCursor = 0;
			chunkBase = 0;
		}
		while (chunkCursor < layout.chunkRefs.size() &&
		       virtualOffset >= chunkBase + layout.chunkRefs[chunkCursor].count) {
			chunkBase += layout.chunkRefs[chunkCursor].count;
			++chunkCursor;
		}
		size_t cursor = chunkCursor;
		uint64_t base = chunkBase;
		while (remaining != 0 && cursor < layout.chunkRefs.size()) {
			const auto& chunk = layout.chunkRefs[cursor];
			const uint64_t offset = virtualOffset > base ? virtualOffset - base : 0;
			const uint64_t available = chunk.count - offset;
			const uint64_t take = std::min<uint64_t>(remaining, available);
			out.push_back(
			    {chunk.threadId, chunk.spoolOffset + offset, take});
			remaining -= take;
			virtualOffset += take;
			if (offset + take >= chunk.count) {
				base += chunk.count;
				++cursor;
			}
		}
		chunkCursor = cursor;
		chunkBase = base;
		return remaining == 0;
	}

	static std::string format_imcf_failure_rate(double rate)
	{
		std::ostringstream oss;
		oss << std::fixed << std::setprecision(4) << rate * 100.0 << '%';
		return oss.str();
	}

	struct ImcfInsertionGroupSummary {
		size_t index;
		uint64_t attempts;
		uint64_t failures;
		double rate;
		size_t worstSlot;
		uint64_t worstSlotAttempts;
		uint64_t worstSlotFailures;
		double worstSlotRate;
	};

	static void report_imcf_insertion_failures(
	    const chimera::imcf::InterleavedMergedCuckooFilter& imcf,
	    const std::vector<std::vector<std::string>>& indexToTaxid,
	    const std::vector<uint64_t>& groupSuccess,
	    const std::vector<uint64_t>& groupFailures,
	    const std::vector<std::vector<uint64_t>>& slotSuccess,
	    const std::vector<std::vector<uint64_t>>& slotFailures)
	{
		uint64_t totalInserted = 0;
		uint64_t totalFailed = 0;
		for (size_t i = 0; i < groupSuccess.size(); ++i) {
			totalInserted += groupSuccess[i];
			totalFailed += groupFailures[i];
		}
		uint64_t totalAttempts = totalInserted + totalFailed;
		if (totalAttempts == 0) {
			return;
		}

		double failureRate =
		    static_cast<double>(totalFailed) / static_cast<double>(totalAttempts);
		std::cout << "  - IMCF insertion result: " << totalFailed << "/"
		          << totalAttempts << " failed ("
		          << format_imcf_failure_rate(failureRate) << ")" << std::endl;
		constexpr double kFailureWarningThreshold = 0.001; // 0.1%
		constexpr uint64_t kFailureMinPrint = 100;
		size_t groupsWithFailures = 0;
		size_t groupsAboveWarn = 0;
		double maxGroupRate = 0.0;
		size_t maxGroupIndex = std::numeric_limits<size_t>::max();
		std::vector<ImcfInsertionGroupSummary> summaries;
		summaries.reserve(groupSuccess.size());

		for (size_t i = 0; i < groupSuccess.size(); ++i) {
			uint64_t attempts = groupSuccess[i] + groupFailures[i];
			if (attempts == 0) {
				continue;
			}
			double groupFailureRate =
			    static_cast<double>(groupFailures[i]) / static_cast<double>(attempts);
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
			if (groupFailureRate <= kFailureWarningThreshold &&
			    groupFailures[i] < kFailureMinPrint) {
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
				double slotRate = static_cast<double>(slotFailures[i][slot]) /
				                  static_cast<double>(slotAttempts);
				if (slotRate > worstSlotRate) {
					worstSlotRate = slotRate;
					worstSlot = slot;
					worstSlotAttempts = slotAttempts;
					worstSlotFailures = slotFailures[i][slot];
				}
			}

			ImcfInsertionGroupSummary summary{
			    i, attempts, groupFailures[i], groupFailureRate, worstSlot,
			    worstSlotAttempts, worstSlotFailures, worstSlotRate};
			summaries.push_back(summary);
		}

		std::cout << "    groups with failures: " << groupsWithFailures << "/"
		          << groupSuccess.size();
		if (maxGroupIndex != std::numeric_limits<size_t>::max()) {
			std::cout << ", worst group=" << maxGroupIndex << " ("
			          << format_imcf_failure_rate(maxGroupRate) << ")";
		}
		std::cout << std::endl;
		if (groupsAboveWarn > 0) {
			std::cout << "    groups above warning threshold ("
			          << format_imcf_failure_rate(kFailureWarningThreshold)
			          << "): " << groupsAboveWarn << std::endl;
		}

		std::sort(summaries.begin(), summaries.end(),
		          [](const ImcfInsertionGroupSummary& a,
		             const ImcfInsertionGroupSummary& b) {
			          return a.rate > b.rate;
		          });

		const size_t maxDetails = 5;
		for (size_t idx = 0; idx < summaries.size() && idx < maxDetails; ++idx) {
			const auto& summary = summaries[idx];
			std::ostringstream line;
			line << "    group " << summary.index << ": " << summary.failures
			     << "/" << summary.attempts << " failed ("
			     << format_imcf_failure_rate(summary.rate) << ")";
			if (summary.rate > kFailureWarningThreshold) {
				line << " - consider lowering load_factor or increasing MaxCuckooCount/binSize";
			}
			std::cout << line.str() << std::endl;
			if (summary.worstSlot != std::numeric_limits<size_t>::max() &&
			    summary.worstSlot < indexToTaxid[summary.index].size() &&
			    summary.worstSlotFailures > 0 &&
			    summary.worstSlotRate > kFailureWarningThreshold) {
				std::cout << "      worst slot " << summary.worstSlot
				          << " (taxid=" << indexToTaxid[summary.index][summary.worstSlot]
				          << ") : " << summary.worstSlotFailures << "/"
				          << summary.worstSlotAttempts << " failed ("
				          << format_imcf_failure_rate(summary.worstSlotRate) << ")"
				          << std::endl;
			}
			auto bucketStats = imcf.computeBucketStats(summary.index, 256);
			std::ostringstream statsLine;
			statsLine << std::fixed << std::setprecision(4)
			          << "      bucket load mean=" << bucketStats.meanLoad
			          << ", stdev=" << bucketStats.stddevLoad
			          << ", full=" << bucketStats.occupancy[4] << "/"
			          << bucketStats.bucketCount << " ("
			          << format_imcf_failure_rate(bucketStats.percentFull) << ")";
			statsLine << std::setprecision(3)
			          << ", low-bit load ratio range="
			          << bucketStats.lowBitMinRatio << " - "
			          << bucketStats.lowBitMaxRatio;
			if (bucketStats.stashEntries > 0) {
				statsLine << std::setprecision(0)
				          << ", stash_species=" << bucketStats.stashEntries
				          << " (entries=" << bucketStats.stashEntryCount
				          << ", max_per_entry=" << bucketStats.stashMaxPerEntry
				          << ")";
			}
			std::cout << statsLine.str() << std::endl;
		}
		uint64_t imcfFailureTotal = imcf.getInsertFailureTotal();
		uint64_t imcfFailureSaturated =
		    imcf.getInsertFailureSaturatedFingerprint();
		if (imcfFailureTotal > 0) {
			double saturatedRatio = static_cast<double>(imcfFailureSaturated) /
			                        static_cast<double>(imcfFailureTotal);
			std::cout << "  - IMCF failure counters: total=" << imcfFailureTotal
			          << ", saturated_fp=" << imcfFailureSaturated << " ("
			          << format_imcf_failure_rate(saturatedRatio) << ")"
			          << std::endl;
		}
	}
	} // namespace

	std::vector<std::vector<std::string>> buildIMCF(
		chimera::imcf::InterleavedMergedCuckooFilter& imcf,
		const std::vector<chimera::imcf::Group>& groups,
		const robin_hood::unordered_flat_map<std::string, uint64_t>& hashCount,
		const HashFrequencyContext* hashFreqContext,
		const FeatureBuildLayout* featureLayout,
		uint16_t effectiveSpan,
		uint16_t refReadLen,
		uint32_t uniqueDegThreshold,
		chimera::presence::CoverageMeta* coverageMeta,
		size_t groupIndexOffset,
		bool verifyShardTotals,
		const robin_hood::unordered_flat_map<std::string, uint64_t>*
		    shardOffsetBase,
		size_t groupRangeBegin,
		size_t groupRangeEnd)
	{
		if (groupRangeEnd == 0) {
			groupRangeEnd = groups.size();
		}
		if (groupRangeBegin > groupRangeEnd || groupRangeEnd > groups.size()) {
			throw std::runtime_error("IMCF build: invalid group range");
		}
		const size_t groupCount = groupRangeEnd - groupRangeBegin;
		std::vector<std::vector<std::string>> indexToTaxid(groupCount);
		robin_hood::unordered_flat_map<std::string, std::vector<TaxidShardPlan>> shardPlan;
		shardPlan.reserve(hashCount.size());
		robin_hood::unordered_flat_map<std::string, uint64_t> shardOffset;
		shardOffset.reserve(hashCount.size());
		if (shardOffsetBase != nullptr) {
			for (const auto& [taxid, offset] : *shardOffsetBase) {
				if (offset != 0) {
					shardOffset.emplace(taxid, offset);
				}
			}
		}
		robin_hood::unordered_flat_map<std::string, size_t> featureLayoutIndex;
		if (featureLayout) {
			featureLayoutIndex.reserve(featureLayout->taxids.size());
			for (size_t idx = 0; idx < featureLayout->taxids.size(); ++idx) {
				featureLayoutIndex.emplace(featureLayout->taxids[idx], idx);
			}
		}

		for (size_t localGroupIdx = 0; localGroupIdx < groupCount; ++localGroupIdx) {
			const auto& group = groups[groupRangeBegin + localGroupIdx];
			indexToTaxid[localGroupIdx] = group.taxids;
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
				plans.push_back({ localGroupIdx, slot, count, offset });
				shardOffset[group.taxids[slot]] = offset + count;
			}
		}

		if (verifyShardTotals) {
			const double kShardToleranceRatio = 0.001;
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
				int64_t delta = static_cast<int64_t>(shardTotal) -
				                static_cast<int64_t>(expected);
				uint64_t absDelta =
				    delta >= 0 ? static_cast<uint64_t>(delta)
				               : static_cast<uint64_t>(-delta);
				double ratio =
				    expected == 0
				        ? (delta == 0 ? 0.0 : std::numeric_limits<double>::infinity())
				        : static_cast<double>(absDelta) /
				              static_cast<double>(expected);
				if (ratio > kShardToleranceRatio) {
					std::cerr << "Warning: feature count mismatch for taxid " << taxid
					          << ": expected=" << expected << ", planned=" << shardTotal
					          << ", delta=" << delta << " (" << std::fixed
					          << std::setprecision(4) << ratio * 100.0 << "%)"
					          << std::endl;
				}
			}
		}

		std::vector<std::string> taxidNames;
		taxidNames.reserve(shardPlan.size());
		robin_hood::unordered_flat_map<std::string, uint32_t> taxidNameToId;
		taxidNameToId.reserve(shardPlan.size());
		for (const auto& kv : shardPlan) {
			const uint32_t taxidId = static_cast<uint32_t>(taxidNames.size());
			taxidNameToId.emplace(kv.first, taxidId);
			taxidNames.push_back(kv.first);
		}
		std::vector<std::vector<ShardEntry>> groupShardEntries(groupCount);

		for (const auto& kv : shardPlan) {
			const std::string& taxid = kv.first;
			const auto& plans = kv.second;
			auto it = hashCount.find(taxid);
			uint64_t totalHashes = (it == hashCount.end()) ? 0 : it->second;
			if (totalHashes == 0) {
				continue;
			}
			auto itTaxidId = taxidNameToId.find(taxid);
			if (itTaxidId == taxidNameToId.end()) {
				continue;
			}
			auto itLayout = featureLayoutIndex.find(taxid);
			if (itLayout == featureLayoutIndex.end()) {
				throw std::runtime_error("IMCF build: missing feature layout for taxid " + taxid);
			}
			const auto& layout = featureLayout->perTaxid[itLayout->second];
			const uint32_t taxidId = itTaxidId->second;
			size_t chunkCursor = 0;
			uint64_t chunkBase = 0;
			for (const auto& plan : plans) {
				if (plan.count == 0) {
					continue;
				}
				ShardEntry entry;
				entry.taxidId = taxidId;
				entry.slotIndex = plan.slotIndex;
				entry.fileOffset = plan.fileOffsetStart;
				entry.nHashes = plan.count;
					if (!map_virtual_range_to_chunks(layout, plan.fileOffsetStart, plan.count,
					                                entry.segments, chunkCursor, chunkBase)) {
						throw std::runtime_error("IMCF build: chunk layout underflow for taxid " + taxid);
					}
					if (plan.groupIndex < groupShardEntries.size()) {
						groupShardEntries[plan.groupIndex].push_back(std::move(entry));
					}
			}
		}

		bool hasAnyShards = false;
		for (auto &groupEntries : groupShardEntries) {
			if (!groupEntries.empty()) {
				hasAnyShards = true;
			}
			std::sort(groupEntries.begin(), groupEntries.end(),
			          [](const ShardEntry &a, const ShardEntry &b) {
				          if (a.taxidId != b.taxidId) {
					          return a.taxidId < b.taxidId;
				          }
				          return a.fileOffset < b.fileOffset;
			          });
		}
		if (!hasAnyShards) {
			return indexToTaxid;
		}

#ifdef _OPENMP
		const size_t workerCount =
		    std::max<size_t>(1, static_cast<size_t>(omp_get_max_threads()));
#else
		const size_t workerCount = std::max<size_t>(
		    1, static_cast<size_t>(std::thread::hardware_concurrency()));
#endif
		const size_t readBlockBytes = resolve_read_block_bytes(workerCount);

		std::vector<uint64_t> groupSuccess(groupCount, 0);
		std::vector<uint64_t> groupFailures(groupCount, 0);
		std::vector<std::vector<uint64_t>> slotSuccess(groupCount);
		std::vector<std::vector<uint64_t>> slotFailures(groupCount);
		for (size_t g = 0; g < groupCount; ++g) {
			size_t slotCount = groups[groupRangeBegin + g].taxids.size();
			slotSuccess[g].assign(slotCount, 0);
			slotFailures[g].assign(slotCount, 0);
		}

#pragma omp parallel
		{
			std::vector<uint64_t> readBlock(readBlockBytes / sizeof(uint64_t));

#pragma omp for schedule(dynamic, 1)
			for (size_t groupIdx = 0; groupIdx < groupCount; ++groupIdx) {
				const auto &entries = groupShardEntries[groupIdx];
				if (entries.empty()) {
					continue;
				}
					uint64_t localSuccess = 0;
					uint64_t localFailure = 0;
					auto &localSlotSuccess = slotSuccess[groupIdx];
					auto &localSlotFailures = slotFailures[groupIdx];
					uint16_t openedThreadId = std::numeric_limits<uint16_t>::max();
					std::ifstream ifile;

					for (const auto &entry : entries) {
						const uint32_t taxidId = entry.taxidId;
					const std::string &taxid = taxidNames[taxidId];
						uint64_t remainingEntry = entry.nHashes;
						for (const auto& segment : entry.segments) {
							if (openedThreadId != segment.threadId) {
								ifile.close();
								ifile.clear();
								if (static_cast<size_t>(segment.threadId) >=
								    featureLayout->threadSpoolPaths.size()) {
#pragma omp critical(imcf_log)
									std::cerr << "Missing feature spool path for thread "
									          << segment.threadId << std::endl;
									remainingEntry = 0;
									break;
								}
								const std::string& spoolPath =
								    featureLayout->threadSpoolPaths[segment.threadId];
								ifile.open(spoolPath, std::ios::binary);
								openedThreadId = segment.threadId;
								if (!ifile) {
#pragma omp critical(imcf_log)
									std::cerr << "Failed to open feature spool file: "
									          << spoolPath << std::endl;
									remainingEntry = 0;
									break;
								}
							}
							ifile.clear();
							ifile.seekg(
							    static_cast<std::streamoff>(segment.spoolOffset * sizeof(uint64_t)),
							    std::ios::beg);
							if (!ifile.good()) {
#pragma omp critical(imcf_log)
								std::cerr << "Warning: failed to seek feature spool for taxid "
								          << taxid << " thread " << segment.threadId << std::endl;
								remainingEntry = 0;
								break;
							}
							uint64_t remainingSegment = segment.count;
						while (ifile && remainingSegment > 0) {
							uint64_t toRead = std::min<uint64_t>(
							    remainingSegment, static_cast<uint64_t>(readBlock.size()));
							ifile.read(reinterpret_cast<char *>(readBlock.data()),
							           static_cast<std::streamsize>(toRead * sizeof(uint64_t)));
							std::streamsize bytesRead = ifile.gcount();
							if (bytesRead <= 0) {
								break;
							}
							size_t hashesRead =
							    static_cast<size_t>(bytesRead) / sizeof(uint64_t);
							if (hashesRead == 0) {
								break;
							}
							for (size_t h = 0; h < hashesRead; ++h) {
								const uint64_t hash = readBlock[h];
								const size_t slotIndex = entry.slotIndex;
								if (imcf.insertTag(groupIdx, groupIndexOffset + groupIdx,
								                   hash, slotIndex)) {
									++localSuccess;
									if (slotIndex < localSlotSuccess.size()) {
										++localSlotSuccess[slotIndex];
									}
								} else {
									++localFailure;
									if (slotIndex < localSlotFailures.size()) {
										++localSlotFailures[slotIndex];
									}
								}
							}
							const uint64_t consumed = static_cast<uint64_t>(hashesRead);
							remainingSegment = (consumed >= remainingSegment)
							                       ? 0
							                       : (remainingSegment - consumed);
							remainingEntry = (consumed >= remainingEntry)
							                     ? 0
							                     : (remainingEntry - consumed);
						}
							if (remainingSegment != 0) {
#pragma omp critical(imcf_log)
								std::cerr << "Warning: spool read underflow for taxid " << taxid
								          << " thread " << segment.threadId
								          << ", expected " << segment.count << " hashes, missing "
								          << remainingSegment << std::endl;
								break;
							}
							if (!ifile.eof() && ifile.fail()) {
#pragma omp critical(imcf_log)
								std::cerr << "Warning: failed to read feature spool completely for taxid "
								          << taxid << " thread " << segment.threadId << std::endl;
							}
						}
					if (remainingEntry != 0) {
#pragma omp critical(imcf_log)
						std::cerr << "Warning: shard read underflow for taxid " << taxid
						          << ", expected " << entry.nHashes << " hashes, missing "
						          << remainingEntry << std::endl;
					}
				}
				groupSuccess[groupIdx] = localSuccess;
				groupFailures[groupIdx] = localFailure;
			}
		}

		report_imcf_insertion_failures(imcf, indexToTaxid, groupSuccess,
		                               groupFailures, slotSuccess, slotFailures);

		if (coverageMeta) {
			populateCoverageMeta(hashFreqContext, featureLayout, effectiveSpan,
			                     refReadLen, uniqueDegThreshold, *coverageMeta);
		}

		return indexToTaxid;
	}

	void populateCoverageMeta(const HashFrequencyContext* hashFreqContext,
	                          const FeatureBuildLayout* featureLayout,
	                          uint16_t effectiveSpan,
	                          uint16_t refReadLen,
	                          uint32_t uniqueDegThreshold,
	                          chimera::presence::CoverageMeta& coverageMeta)
	{
		const uint32_t uniqueDeg = std::max<uint32_t>(1, uniqueDegThreshold);
		const uint16_t spanUsed =
		    (effectiveSpan > 0) ? effectiveSpan : static_cast<uint16_t>(1);
		const uint16_t refLenUsed =
		    (refReadLen == 0) ? static_cast<uint16_t>(150) : refReadLen;
		const bool hasSketch = hashFreqContext && hashFreqContext->enabled();
		coverageMeta.unique_deg_threshold = uniqueDeg;
		coverageMeta.ref_read_length = refLenUsed;
		coverageMeta.effective_span = spanUsed;
		coverageMeta.entries.clear();
		if (!featureLayout) {
			return;
		}
		coverageMeta.entries.reserve(featureLayout->taxids.size());
		for (size_t layoutIdx = 0; layoutIdx < featureLayout->taxids.size();
		     ++layoutIdx) {
			const auto& taxid = featureLayout->taxids[layoutIdx];
			const auto& layout = featureLayout->perTaxid[layoutIdx];
			chimera::presence::CoverageEntry entry{};
			entry.taxid = taxid;
			entry.unique_signatures = layout.uniqueSignatures;
			entry.total_signatures = layout.totalSignatures;
			entry.genome_length = layout.genomeBp;
			const double window = std::max<int64_t>(
			    1, static_cast<int64_t>(refLenUsed) - static_cast<int64_t>(spanUsed) + 1);
			entry.unique_density = (layout.genomeBp > 0 && entry.unique_signatures > 0)
			                           ? static_cast<double>(entry.unique_signatures) /
			                                 static_cast<double>(layout.genomeBp)
			                           : 0.0;
			entry.expected_unique_per_ref_read =
			    (entry.unique_density > 0.0) ? entry.unique_density * window : 0.0;
			coverageMeta.entries.push_back(entry);
		}
		if (hasSketch) {
			coverageMeta.freq_model.depth = hashFreqContext->sketch->depth();
			coverageMeta.freq_model.width = hashFreqContext->sketch->width();
			coverageMeta.freq_model.quantile = hashFreqContext->quantile;
			coverageMeta.freq_model.stats = hashFreqContext->stats;
			coverageMeta.freq_model.total_hashes =
			    hashFreqContext->passA_total_hashes.load(std::memory_order_relaxed);
			coverageMeta.freq_model.filtered_hashes =
			    hashFreqContext->passB_filtered_hashes.load(std::memory_order_relaxed);
			hashFreqContext->sketch->exportCounts(coverageMeta.freq_model.counters);
		}
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

	try {
		logStep("Writing filter archive (.imcf)", [&]() {
			cereal::BinaryOutputArchive archive(os);
			imcf.save_for_archive(archive);
			archive(indexToTaxid);
			archive(imcfConfig);
			if (presenceMeta) {
				archive(*presenceMeta);
			}
			os.close();
			return true;
		});
	} catch (...) {
		imcf.cleanup_qidx_spool_files();
		os.close();
		throw;
	}
	imcf.cleanup_qidx_spool_files();

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
