#include "ChimeraClassifyCommon.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <stdexcept>

#include <dna4_traits.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/utility/views/chunk.hpp>

namespace ChimeraClassify {

namespace {

constexpr char kSpoolMagic[] = {'C', 'H', 'S', 'P', '2', '\0', '\0', '\0'};

template <typename T>
void write_pod(std::ostream &os, const T &value) {
  os.write(reinterpret_cast<const char *>(&value), sizeof(T));
  if (!os) {
    throw std::runtime_error("Failed to write classify spool");
  }
}

template <typename T>
bool read_pod(std::istream &is, T &value) {
  is.read(reinterpret_cast<char *>(&value), sizeof(T));
  return static_cast<bool>(is);
}

void write_string(std::ostream &os, const std::string &value) {
  const uint32_t len = static_cast<uint32_t>(value.size());
  write_pod(os, len);
  if (len > 0) {
    os.write(value.data(), static_cast<std::streamsize>(len));
    if (!os) {
      throw std::runtime_error("Failed to write classify spool string");
    }
  }
}

void read_string(std::istream &is, std::string &value,
                 const std::string &path) {
  uint32_t len = 0;
  if (!read_pod(is, len)) {
    throw std::runtime_error("Truncated classify spool string: " + path);
  }
  value.resize(len);
  if (len > 0) {
    is.read(value.data(), static_cast<std::streamsize>(len));
    if (!is) {
      throw std::runtime_error("Truncated classify spool string: " + path);
    }
  }
}

} // namespace

void write_spool_header(std::ostream &os) {
  os.write(kSpoolMagic, sizeof(kSpoolMagic));
  if (!os) {
    throw std::runtime_error("Failed to write classify spool header");
  }
}

void read_spool_header(std::istream &is, const std::string &path) {
  char magic[sizeof(kSpoolMagic)]{};
  is.read(magic, sizeof(magic));
  if (!is || !std::equal(std::begin(magic), std::end(magic),
                         std::begin(kSpoolMagic))) {
    throw std::runtime_error("Invalid classify spool header: " + path);
  }
}

void write_spool_record(std::ostream &os, const SpoolReadRecord &record) {
  write_string(os, record.id);
  write_pod(os, record.evaluated);
  write_pod(os, record.best_taxid_hint);
  write_string(os, record.reject_reason);
  const uint32_t cand_len = static_cast<uint32_t>(record.candidates.size());
  write_pod(os, cand_len);
  if (cand_len > 0) {
    os.write(reinterpret_cast<const char *>(record.candidates.data()),
             static_cast<std::streamsize>(cand_len * sizeof(SpoolCandidate)));
    if (!os) {
      throw std::runtime_error("Failed to write classify spool candidates");
    }
  }
  const uint32_t abundance_len =
      static_cast<uint32_t>(record.abundance_candidates.size());
  write_pod(os, abundance_len);
  if (abundance_len > 0) {
    os.write(reinterpret_cast<const char *>(record.abundance_candidates.data()),
             static_cast<std::streamsize>(abundance_len *
                                          sizeof(SpoolCandidate)));
    if (!os) {
      throw std::runtime_error(
          "Failed to write classify spool abundance candidates");
    }
  }
}

void write_spool_record(std::ostream &os,
                        const CompactClassifyResult &record) {
  write_string(os, record.id);
  write_pod(os, record.evaluated);
  write_pod(os, record.best_taxid_hint);
  write_string(os, record.reject_reason);
  const uint32_t cand_len = static_cast<uint32_t>(record.candidates.size());
  write_pod(os, cand_len);
  if (cand_len > 0) {
    os.write(reinterpret_cast<const char *>(record.candidates.data()),
             static_cast<std::streamsize>(cand_len * sizeof(SpoolCandidate)));
    if (!os) {
      throw std::runtime_error("Failed to write classify spool candidates");
    }
  }
  const uint32_t abundance_len =
      static_cast<uint32_t>(record.abundance_candidates.size());
  write_pod(os, abundance_len);
  if (abundance_len > 0) {
    os.write(reinterpret_cast<const char *>(record.abundance_candidates.data()),
             static_cast<std::streamsize>(abundance_len *
                                          sizeof(SpoolCandidate)));
    if (!os) {
      throw std::runtime_error(
          "Failed to write classify spool abundance candidates");
    }
  }
}

bool read_spool_record(std::istream &is, SpoolReadRecord &record) {
  uint32_t id_len = 0;
  is.read(reinterpret_cast<char *>(&id_len), sizeof(id_len));
  if (!is) {
    return false;
  }
  record.id.resize(id_len);
  if (id_len > 0) {
    is.read(record.id.data(), static_cast<std::streamsize>(id_len));
    if (!is) {
      throw std::runtime_error("Truncated classify spool record id");
    }
  }
  if (!read_pod(is, record.evaluated) ||
      !read_pod(is, record.best_taxid_hint)) {
    throw std::runtime_error("Truncated classify spool record header");
  }
  uint32_t reject_len = 0;
  if (!read_pod(is, reject_len)) {
    throw std::runtime_error("Truncated classify spool reject length");
  }
  record.reject_reason.resize(reject_len);
  if (reject_len > 0) {
    is.read(record.reject_reason.data(),
            static_cast<std::streamsize>(reject_len));
    if (!is) {
      throw std::runtime_error("Truncated classify spool reject reason");
    }
  }
  uint32_t cand_len = 0;
  if (!read_pod(is, cand_len)) {
    throw std::runtime_error("Truncated classify spool candidate length");
  }
  record.candidates.resize(cand_len);
  if (cand_len > 0) {
    is.read(reinterpret_cast<char *>(record.candidates.data()),
            static_cast<std::streamsize>(cand_len * sizeof(SpoolCandidate)));
    if (!is) {
      throw std::runtime_error("Truncated classify spool candidates");
    }
  }
  uint32_t abundance_len = 0;
  if (!read_pod(is, abundance_len)) {
    throw std::runtime_error("Truncated classify spool abundance length");
  }
  record.abundance_candidates.resize(abundance_len);
  if (abundance_len > 0) {
    is.read(reinterpret_cast<char *>(record.abundance_candidates.data()),
            static_cast<std::streamsize>(abundance_len *
                                         sizeof(SpoolCandidate)));
    if (!is) {
      throw std::runtime_error("Truncated classify spool abundance candidates");
    }
  }
  return true;
}

void parseReads(std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
                ClassifyConfig config, FileInfo &fileInfo,
                size_t max_reads,
                std::vector<QueueThrottle> *queueThrottles) {
  if (readQueues.empty()) {
    throw std::runtime_error("readQueues is empty");
  }

  const size_t shardCount = readQueues.size();
  std::hash<std::string> hasher;

  auto init_batch = [&](batchReads &batch, bool paired) {
    batch.ids.reserve(config.batchSize);
    batch.seqs.reserve(config.batchSize);
    if (paired) {
      batch.seqs2.reserve(config.batchSize);
    }
  };

  auto flush_batch = [&](size_t shard, batchReads &batch, bool paired) {
    if (batch.ids.empty()) {
      return;
    }
    QueueThrottle *throttle =
        (queueThrottles != nullptr && shard < queueThrottles->size())
            ? &(*queueThrottles)[shard]
            : nullptr;
    acquire_queue_slot(throttle, estimate_batch_bytes(batch));
    readQueues[shard].enqueue(std::move(batch));
    batchReads fresh;
    init_batch(fresh, paired);
    batch = std::move(fresh);
  };

  size_t totalSequences = 0;
  size_t totalFiles = 0;
  bool reached_limit = false;

  if (!config.singleFiles.empty()) {
    std::vector<batchReads> pending(shardCount);
    for (auto &b : pending) {
      init_batch(b, false);
    }

    for (const auto &file : config.singleFiles) {
      ++totalFiles;

      seqan3::sequence_file_input<
          raptor::dna4_traits,
          seqan3::fields<seqan3::field::id, seqan3::field::seq>>
          fin{file};

      for (auto &&r : fin) {
        std::string id = std::move(r.id());
        size_t shard = hasher(id) % shardCount;
        auto &batch = pending[shard];
        batch.ids.emplace_back(std::move(id));
        batch.seqs.emplace_back(std::move(r.sequence()));
        ++totalSequences;
        if (max_reads > 0 && totalSequences >= max_reads) {
          reached_limit = true;
          break;
        }
        if (batch.ids.size() >= config.batchSize) {
          flush_batch(shard, batch, false);
        }
      }
      if (reached_limit) {
        break;
      }
    }

    for (size_t shard = 0; shard < shardCount; ++shard) {
      flush_batch(shard, pending[shard], false);
    }
  } else if (!config.pairedFiles.empty()) {
    std::vector<batchReads> pending(shardCount);
    for (auto &b : pending) {
      init_batch(b, true);
    }

    for (size_t i = 0; i < config.pairedFiles.size(); i += 2) {
      totalFiles += 2;

      seqan3::sequence_file_input<
          raptor::dna4_traits,
          seqan3::fields<seqan3::field::id, seqan3::field::seq>>
          fin1{config.pairedFiles[i]};
      seqan3::sequence_file_input<
          raptor::dna4_traits,
          seqan3::fields<seqan3::field::id, seqan3::field::seq>>
          fin2{config.pairedFiles[i + 1]};

      auto it1 = fin1.begin();
      auto end1 = fin1.end();
      auto it2 = fin2.begin();
      auto end2 = fin2.end();

      for (; it1 != end1 && it2 != end2; ++it1, ++it2) {
        auto rec1 = *it1;
        auto rec2 = *it2;
        std::string id = std::move(rec1.id());
        size_t shard = hasher(id) % shardCount;
        auto &batch = pending[shard];
        batch.ids.emplace_back(std::move(id));
        batch.seqs.emplace_back(std::move(rec1.sequence()));
        batch.seqs2.emplace_back(std::move(rec2.sequence()));
        ++totalSequences;
        if (max_reads > 0 && totalSequences >= max_reads) {
          reached_limit = true;
          break;
        }
        if (batch.ids.size() >= config.batchSize) {
          flush_batch(shard, batch, true);
        }
      }
      if (reached_limit) {
        break;
      }
    }

    for (size_t shard = 0; shard < shardCount; ++shard) {
      flush_batch(shard, pending[shard], true);
    }
  } else {
    throw std::runtime_error("No input files specified");
  }

  fileInfo.sequenceNum += totalSequences;
  fileInfo.fileNum += totalFiles;

}

void saveResult(const std::vector<classifyResult> &classifyResults,
                const ClassifyConfig &config) {
  auto resolve_tsv_path = [](std::string path) {
    if (std::filesystem::path(path).extension() != ".tsv") {
      path += ".tsv";
    }
    return path;
  };

  const std::string outputFile = resolve_tsv_path(config.outputFile);

  std::ofstream os(outputFile, std::ios::out);
  if (!os.is_open()) {
    throw std::runtime_error("Failed to open file: " + outputFile);
  }
  std::vector<char> outputBuffer(1 << 20, '\0');
  os.rdbuf()->pubsetbuf(outputBuffer.data(),
                        static_cast<std::streamsize>(outputBuffer.size()));
  std::ostringstream postTopkOss;

  for (const auto &result : classifyResults) {
    writeResultRecord(os, result, postTopkOss);
  }
  os.close();
}

void writeResultRecord(std::ostream &os, const classifyResult &result,
                       std::ostringstream &postTopkOss) {
  constexpr size_t kDumpPostTopK = 16;
  os << result.id << '\t';
  bool handled = false;
  if (!result.taxidCount.empty() &&
      result.taxidCount.front().first == "unclassified") {
    os << "unclassified";
    if (!result.reject_reason.empty()) {
      os << "\tREJECT=" << result.reject_reason;
    }
    if (!result.best_taxid_hint.empty()) {
      os << "\tHINT=" << result.best_taxid_hint;
    }
    handled = true;
  }
  if (!handled) {
    for (const auto &[taxid, count] : result.taxidCount) {
      if (taxid == "unclassified") {
        os << taxid;
        continue;
      }
      os << taxid << ':' << count << '\t';
    }
  }
  if (!result.posteriors.empty()) {
    if (handled) {
      os << '\t';
    }
    const size_t k = std::min(kDumpPostTopK, result.posteriors.size());
    postTopkOss.str("");
    postTopkOss.clear();
    postTopkOss.setf(std::ios::fixed);
    postTopkOss << std::setprecision(6);
    postTopkOss << "POST_TOPK=";
    for (size_t i = 0; i < k; ++i) {
      if (i) {
        postTopkOss << ',';
      }
      postTopkOss << result.posteriors[i].first << ':'
                  << result.posteriors[i].second;
    }
    os << postTopkOss.str() << '\t';
  }
  os << '\n';
}

} // namespace ChimeraClassify
