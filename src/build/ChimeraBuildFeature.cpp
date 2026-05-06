#include "ChimeraBuildCommon.hpp"

#include <dna4_traits.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <filesystem>
#include <fcntl.h>
#include <iostream>
#include <limits>
#include <mutex>
#include <unistd.h>

namespace ChimeraBuild {

namespace {
std::filesystem::path g_tmp_work_dir = "tmp";

struct SpoolChunkRecord {
  uint32_t taxidIndex{0};
  uint64_t spoolOffset{0};
  uint64_t count{0};
};

struct PendingBufferChunk {
  uint32_t taxidIndex{0};
  uint64_t bufferOffset{0};
  uint64_t count{0};
};

constexpr size_t kMinThreadBufferFlushBytes = 1 * 1024 * 1024;
constexpr size_t kMaxThreadBufferFlushBytes = 4 * 1024 * 1024;
constexpr size_t kTargetTotalThreadBufferBytes = 384ull * 1024ull * 1024ull;

inline size_t resolve_thread_buffer_flush_bytes(size_t workerCount) {
  return std::clamp(
      kTargetTotalThreadBufferBytes / std::max<size_t>(1, workerCount),
      kMinThreadBufferFlushBytes, kMaxThreadBufferFlushBytes);
}

class DirectSpoolWriter {
public:
  explicit DirectSpoolWriter(const std::vector<std::string> &spoolPathsIn)
      : spoolPaths(spoolPathsIn), fds(spoolPathsIn.size(), -1) {}

  DirectSpoolWriter(const DirectSpoolWriter &) = delete;
  DirectSpoolWriter &operator=(const DirectSpoolWriter &) = delete;

  ~DirectSpoolWriter() { close_all_fds(); }

  void finish() {
    close_all_fds();
    if (failedFlag) {
      throw std::runtime_error(failureMessage.empty()
                                   ? "Feature spool direct writer failed"
                                   : failureMessage);
    }
  }

  bool failed() const { return failedFlag.load(std::memory_order_relaxed); }

  bool write_hashes(uint16_t threadId, uint64_t spoolOffset,
                    const std::vector<uint64_t> &hashes) {
    if (hashes.empty()) {
      return true;
    }
    if (failedFlag.load(std::memory_order_relaxed)) {
      return false;
    }
    const size_t tid = static_cast<size_t>(threadId);
    if (tid >= fds.size()) {
      set_failure("Feature spool thread id out of range");
      return false;
    }
    int fd = -1;
    if (!ensure_fd(tid, fd)) {
      return false;
    }
    const char *data = reinterpret_cast<const char *>(hashes.data());
    const size_t bytes = hashes.size() * sizeof(uint64_t);
    size_t written = 0;
    while (written < bytes) {
      const ssize_t rc =
          ::pwrite(fd, data + written, bytes - written,
                   static_cast<off_t>(spoolOffset * sizeof(uint64_t) +
                                      written));
      if (rc < 0) {
        if (errno == EINTR) {
          continue;
        }
        const int err = errno;
        set_failure("Failed to write feature spool file: " +
                    spoolPaths[tid] + " (" + std::strerror(err) +
                    ")");
        return false;
      }
      if (rc == 0) {
        set_failure("Failed to write feature spool file: " +
                    spoolPaths[tid] + " (short write)");
        return false;
      }
      written += static_cast<size_t>(rc);
    }
    return true;
  }

private:
  void set_failure(std::string message) {
    std::lock_guard<std::mutex> lock(failureMutex);
    if (failureMessage.empty()) {
      failureMessage = std::move(message);
    }
    failedFlag.store(true, std::memory_order_relaxed);
  }

  bool ensure_fd(size_t threadId, int &fd) {
    fd = fds[threadId];
    if (fd >= 0) {
      return true;
    }
    fd = ::open(spoolPaths[threadId].c_str(),
                O_CREAT | O_TRUNC | O_WRONLY | O_CLOEXEC, 0644);
    if (fd < 0) {
      const int err = errno;
      set_failure("Unable to open feature spool file: " + spoolPaths[threadId] +
                  " (" + std::strerror(err) + ")");
      return false;
    }
    fds[threadId] = fd;
    return true;
  }

  void close_all_fds() {
    for (int &fd : fds) {
      if (fd >= 0) {
        ::close(fd);
        fd = -1;
      }
    }
  }

  const std::vector<std::string> &spoolPaths;
  std::vector<int> fds;
  std::atomic<bool> failedFlag{false};
  std::mutex failureMutex;
  std::string failureMessage;
};

inline void append_kept_hash(uint32_t taxidIndex, uint64_t hash,
                             bool isUnique,
                             std::vector<uint64_t> &threadBuffer,
                             std::vector<PendingBufferChunk> &pendingChunks,
                             uint64_t *totalSignatures,
                             uint64_t *uniqueSignatures) {
  if (pendingChunks.empty() || pendingChunks.back().taxidIndex != taxidIndex) {
    pendingChunks.push_back(
        {taxidIndex, static_cast<uint64_t>(threadBuffer.size()), 0});
  }
  threadBuffer.push_back(hash);
  ++pendingChunks.back().count;
  ++totalSignatures[taxidIndex];
  if (isUnique) {
    ++uniqueSignatures[taxidIndex];
  }
}

inline void flush_feature_thread_buffer(
    DirectSpoolWriter &spoolWriter, uint16_t threadId,
    std::vector<uint64_t> &threadBuffer,
    std::vector<PendingBufferChunk> &pendingChunks,
    uint64_t &localSpoolOffset, std::vector<SpoolChunkRecord> &localChunks,
    uint64_t *localHashCounts) {
  if (threadBuffer.empty()) {
    return;
  }
  if (pendingChunks.empty()) {
    threadBuffer.clear();
    return;
  }
  if (spoolWriter.failed()) {
    threadBuffer.clear();
    pendingChunks.clear();
    return;
  }

  const size_t buffer_size = threadBuffer.size();
  const uint64_t chunk_offset = localSpoolOffset;
  std::vector<PendingBufferChunk> flush_chunks;
  flush_chunks.swap(pendingChunks);
  if (!spoolWriter.write_hashes(threadId, chunk_offset, threadBuffer)) {
    threadBuffer.clear();
    pendingChunks.clear();
    return;
  }
  for (const auto &chunk : flush_chunks) {
    localChunks.push_back(
        {chunk.taxidIndex, chunk_offset + chunk.bufferOffset, chunk.count});
    localHashCounts[chunk.taxidIndex] += chunk.count;
  }
  localSpoolOffset += static_cast<uint64_t>(buffer_size);
  threadBuffer.clear();
}
} // namespace

void set_tmp_work_dir(const std::filesystem::path &dir) {
  if (dir.empty()) {
    g_tmp_work_dir = "tmp";
    return;
  }
  g_tmp_work_dir = dir;
}

const std::filesystem::path &tmp_work_dir() { return g_tmp_work_dir; }

std::string tmp_hash_spool_path(std::string_view suffix, int thread_id) {
  std::filesystem::path path = tmp_work_dir();
  path /=
      "feature" + std::string(suffix) + ".spool." + std::to_string(thread_id);
  return path;
}

void feature_count(
    BuildConfig &config,
    robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    FileInfo &fileInfo, HashFrequencyContext *hashFreqContext,
    robin_hood::unordered_flat_map<std::string, uint64_t> *bpCount,
    FeatureBuildLayout *featureLayout) {
  uint64_t feature_seed = 0;
  chimera::feature::Params feature_params =
      make_feature_params(config, feature_seed);
  const std::string feature_suffix = ".strb";
  if (config.verbose) {
    std::cout << "Feature method: strobemer (k="
              << static_cast<int>(feature_params.strobe.k)
              << ", order=" << static_cast<int>(feature_params.strobe.order)
              << ", w=[" << feature_params.strobe.w_min << ','
              << feature_params.strobe.w_max << "], seed="
              << static_cast<unsigned long long>(feature_params.strobe.seed)
              << ")" << std::endl;
  }

  struct TaxidFileTask {
    uint32_t taxidIndex;
    const std::string *filePath;
  };

  size_t totalFiles = 0;
  for (const auto &[taxid, files] : inputFiles) {
    (void)taxid;
    totalFiles += files.size();
  }

  std::vector<TaxidFileTask> taxid_file_tasks;
  taxid_file_tasks.reserve(totalFiles);
  std::vector<std::string> index_to_taxid;
  index_to_taxid.reserve(inputFiles.size());
  for (auto &[taxid, files] : inputFiles) {
    const uint32_t taxid_index = static_cast<uint32_t>(index_to_taxid.size());
    index_to_taxid.push_back(taxid);
    for (auto &file : files) {
      taxid_file_tasks.push_back({taxid_index, &file});
    }
  }

  const size_t worker_hint =
      std::max<size_t>(1, static_cast<size_t>(config.threads));
  const size_t hash_buffer_flush_bytes =
      resolve_thread_buffer_flush_bytes(worker_hint);
  const size_t taxidCount = index_to_taxid.size();

  std::vector<uint64_t> threadHashCounts(taxidCount * worker_hint, 0u);
  std::vector<uint64_t> threadBpCounts(worker_hint * taxidCount, 0u);
  std::vector<uint64_t> threadTotalSignatures(worker_hint * taxidCount, 0u);
  std::vector<uint64_t> threadUniqueSignatures(worker_hint * taxidCount, 0u);
  std::vector<FileInfo> threadFileInfos(worker_hint);
  std::vector<uint64_t> threadPassBTotal(worker_hint, 0u);
  std::vector<uint64_t> threadPassBFiltered(worker_hint, 0u);
  std::vector<std::vector<SpoolChunkRecord>> threadChunkManifests(worker_hint);
  std::vector<std::string> threadSpoolPaths(worker_hint);
  for (size_t tid = 0; tid < worker_hint; ++tid) {
    threadSpoolPaths[tid] =
        tmp_hash_spool_path(feature_suffix, static_cast<int>(tid));
  }
  DirectSpoolWriter spoolWriter(threadSpoolPaths);

  const size_t feature_min_length =
      chimera::feature::min_required_length(feature_params);
  const size_t min_required =
      std::max<size_t>(config.min_length, feature_min_length);
  if (config.verbose) {
    std::cout << "Feature hash minimum read length: " << min_required
              << std::endl;
  }
  const bool has_freq_sketch =
      (hashFreqContext != nullptr) && hashFreqContext->enabled();

  int used_threads = 1;
#pragma omp parallel
  {
    std::vector<uint64_t> thread_buffer;
    thread_buffer.reserve(8192);
    std::vector<PendingBufferChunk> pending_buffer_chunks;
    pending_buffer_chunks.reserve(256);
    const int tid = omp_get_thread_num();
    uint64_t *local_bp_counts =
        threadBpCounts.data() + static_cast<size_t>(tid) * taxidCount;
    uint64_t *local_hash_counts =
        threadHashCounts.data() + static_cast<size_t>(tid) * taxidCount;
    uint64_t *local_total_signatures =
        threadTotalSignatures.data() + static_cast<size_t>(tid) * taxidCount;
    uint64_t *local_unique_signatures =
        threadUniqueSignatures.data() + static_cast<size_t>(tid) * taxidCount;
    FileInfo localThreadInfo{};
    uint64_t localPassBTotal = 0;
    uint64_t localPassBFiltered = 0;
    uint64_t localSpoolOffset = 0;
    std::vector<SpoolChunkRecord> localChunks;
#ifdef _OPENMP
#pragma omp single
    { used_threads = omp_get_num_threads(); }
#endif

#pragma omp for schedule(dynamic)
    for (size_t idx = 0; idx < taxid_file_tasks.size(); ++idx) {
      if (spoolWriter.failed()) {
        continue;
      }
      const TaxidFileTask &task = taxid_file_tasks[idx];
      const size_t taxid_index = static_cast<size_t>(task.taxidIndex);
      const std::string &filename = *(task.filePath);

      seqan3::sequence_file_input<raptor::dna4_traits,
                                  seqan3::fields<seqan3::field::id,
                                                  seqan3::field::seq>>
          fin{filename};
      std::vector<uint64_t> hashes;
      hashes.reserve(4096);

      for (auto &record : fin) {
        if (spoolWriter.failed()) {
          break;
        }
        auto &seq = record.sequence();
        const uint64_t seq_len = seq.size();
        local_bp_counts[taxid_index] += seq_len;
        if (seq_len < min_required) {
          localThreadInfo.skippedSeqNum++;
          continue;
        }
        localThreadInfo.sequenceNum++;
        localThreadInfo.bpLength += seq_len;

        hashes.clear();
        chimera::feature::compute_hashes_append(seq, feature_params, hashes);
        if (hashes.empty()) {
          continue;
        }

        size_t appended = 0;
        for (uint64_t hash : hashes) {
          if (has_freq_sketch) {
            ++localPassBTotal;
            if (hashFreqContext->should_filter_hash(hash)) {
              ++localPassBFiltered;
              continue;
            }
          }
          append_kept_hash(task.taxidIndex, hash,
                           !has_freq_sketch ||
                               hashFreqContext->is_unique_signature(hash),
                           thread_buffer, pending_buffer_chunks,
                           local_total_signatures, local_unique_signatures);
          ++appended;
        }
        if (appended == 0) {
          continue;
        }
        if (thread_buffer.size() * sizeof(uint64_t) >=
            hash_buffer_flush_bytes) {
          flush_feature_thread_buffer(
              spoolWriter, static_cast<uint16_t>(tid), thread_buffer,
              pending_buffer_chunks, localSpoolOffset, localChunks,
              local_hash_counts);
        }
      }
    }

    if (!thread_buffer.empty()) {
      flush_feature_thread_buffer(spoolWriter, static_cast<uint16_t>(tid),
                                  thread_buffer, pending_buffer_chunks,
                                  localSpoolOffset, localChunks,
                                  local_hash_counts);
    }
    threadFileInfos[static_cast<size_t>(tid)] = localThreadInfo;
    threadPassBTotal[static_cast<size_t>(tid)] = localPassBTotal;
    threadPassBFiltered[static_cast<size_t>(tid)] = localPassBFiltered;
    threadChunkManifests[static_cast<size_t>(tid)] = std::move(localChunks);
  }
  spoolWriter.finish();

  uint64_t passBTotal = 0;
  uint64_t passBFiltered = 0;
  for (size_t tid = 0; tid < static_cast<size_t>(used_threads); ++tid) {
    fileInfo.skippedNum += threadFileInfos[tid].skippedNum;
    fileInfo.skippedSeqNum += threadFileInfos[tid].skippedSeqNum;
    fileInfo.sequenceNum += threadFileInfos[tid].sequenceNum;
    fileInfo.bpLength += threadFileInfos[tid].bpLength;
    passBTotal += threadPassBTotal[tid];
    passBFiltered += threadPassBFiltered[tid];
  }
  if (hashFreqContext) {
    hashFreqContext->passB_total_hashes.store(passBTotal,
                                              std::memory_order_relaxed);
    hashFreqContext->passB_filtered_hashes.store(passBFiltered,
                                                 std::memory_order_relaxed);
  }

  if (featureLayout) {
    featureLayout->taxids = index_to_taxid;
    featureLayout->perTaxid.clear();
    featureLayout->perTaxid.resize(taxidCount);
    featureLayout->threadSpoolPaths.assign(
        threadSpoolPaths.begin(),
        threadSpoolPaths.begin() + static_cast<size_t>(used_threads));
    std::vector<size_t> chunkCounts(taxidCount, 0u);
    for (size_t tid = 0; tid < static_cast<size_t>(used_threads); ++tid) {
      for (const auto &chunk : threadChunkManifests[tid]) {
        ++chunkCounts[chunk.taxidIndex];
      }
    }
    for (size_t idx = 0; idx < taxidCount; ++idx) {
      featureLayout->perTaxid[idx].chunkRefs.reserve(chunkCounts[idx]);
    }
  }

  for (size_t idx = 0; idx < index_to_taxid.size(); ++idx) {
    const std::string &taxid = index_to_taxid[idx];
    uint64_t total_hashes = 0;
    uint64_t genome_bp = 0;
    uint64_t total_signatures = 0;
    uint64_t unique_signatures = 0;
    for (size_t tid = 0; tid < static_cast<size_t>(used_threads); ++tid) {
      total_hashes += threadHashCounts[idx + tid * taxidCount];
      const size_t thread_offset = tid * taxidCount + idx;
      genome_bp += threadBpCounts[thread_offset];
      total_signatures += threadTotalSignatures[thread_offset];
      unique_signatures += threadUniqueSignatures[thread_offset];
    }
    hashCount[taxid] = total_hashes;
    if (bpCount) {
      (*bpCount)[taxid] = genome_bp;
    }
    if (featureLayout) {
      auto &layout = featureLayout->perTaxid[idx];
      layout.totalSignatures = total_signatures;
      layout.uniqueSignatures = unique_signatures;
      layout.genomeBp = genome_bp;
    }
  }

  if (featureLayout) {
    for (size_t tid = 0; tid < static_cast<size_t>(used_threads); ++tid) {
      for (const auto &chunk : threadChunkManifests[tid]) {
        featureLayout->perTaxid[chunk.taxidIndex].chunkRefs.push_back(
            {static_cast<uint16_t>(tid), chunk.spoolOffset, chunk.count});
      }
    }
  }
}

} // namespace ChimeraBuild
