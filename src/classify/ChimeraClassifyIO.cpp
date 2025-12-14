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
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/chunk.hpp>

namespace ChimeraClassify {

void ReadStats::update(size_t len) {
  if (len == 0)
    return;
  ++count;
  total_len += len;
  if (len < min_len)
    min_len = len;
  if (len > max_len)
    max_len = len;
}

ReadStats sample_read_stats(const ClassifyConfig &config, size_t max_reads) {
  ReadStats stats;
  if (!config.singleFiles.empty()) {
    for (const auto &file : config.singleFiles) {
      try {
        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin{file};
        for (auto &record : fin) {
          stats.update(record.sequence().size());
          if (stats.count >= max_reads)
            break;
        }
      } catch (const std::exception &e) {
        std::cerr << "Warning: 采样读取文件 " << file
                  << " 失败: " << e.what() << std::endl;
      }
      if (stats.count >= max_reads)
        break;
    }
  } else if (!config.pairedFiles.empty()) {
    for (size_t i = 0; i + 1 < config.pairedFiles.size(); i += 2) {
      try {
        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin1{config.pairedFiles[i]};
        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin2{config.pairedFiles[i + 1]};
        auto it1 = fin1.begin();
        auto it2 = fin2.begin();
        for (; it1 != fin1.end() && it2 != fin2.end(); ++it1, ++it2) {
          stats.update((*it1).sequence().size());
          stats.update((*it2).sequence().size());
          if (stats.count >= max_reads)
            break;
        }
      } catch (const std::exception &e) {
        std::cerr << "Warning: 采样读取文件 " << config.pairedFiles[i]
                  << " 或 " << config.pairedFiles[i + 1]
                  << " 失败: " << e.what() << std::endl;
      }
      if (stats.count >= max_reads)
        break;
    }
  }
  if (stats.count == 0) {
    stats.min_len = 0;
    stats.max_len = 0;
  }
  return stats;
}

void parseReads(std::vector<moodycamel::ConcurrentQueue<batchReads>> &readQueues,
                ClassifyConfig config, FileInfo &fileInfo) {
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
    readQueues[shard].enqueue(std::move(batch));
    batchReads fresh;
    init_batch(fresh, paired);
    batch = std::move(fresh);
  };

  size_t totalSequences = 0;
  size_t totalFiles = 0;

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
        if (batch.ids.size() >= config.batchSize) {
          flush_batch(shard, batch, false);
        }
      }
    }

    for (size_t shard = 0; shard < shardCount; ++shard) {
      flush_batch(shard, pending[shard], false);
    }
  } else if (!config.pairedFiles.empty()) {
    if (config.pairedFiles.size() % 2 != 0) {
      throw std::runtime_error("Paired input requires an even number of files");
    }

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
        if (batch.ids.size() >= config.batchSize) {
          flush_batch(shard, batch, true);
        }
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

  if (config.verbose) {
    std::cout << "Number of files: " << fileInfo.fileNum << std::endl;
    std::cout << "Number of sequences: " << fileInfo.sequenceNum << std::endl
              << std::endl;
  }
}

void print_classify_time(long long milliseconds) {
  long long total_seconds = milliseconds / 1000;
  long long seconds = total_seconds % 60;
  long long total_minutes = total_seconds / 60;
  long long minutes = total_minutes % 60;
  long long hours = total_minutes / 60;

  if (hours > 0) {
    std::cout << hours << "h " << minutes << "min " << seconds << "s "
              << milliseconds % 1000 << "ms" << std::endl;
  } else if (minutes > 0) {
    std::cout << minutes << "min " << seconds << "s "
              << milliseconds % 1000 << "ms" << std::endl;
  } else {
    std::cout << seconds << "s " << milliseconds % 1000 << "ms" << std::endl;
  }
}

void saveResult(std::vector<classifyResult> classifyResults,
                ClassifyConfig config) {
  std::string outputFile = config.outputFile;
  if (std::filesystem::path(outputFile).extension() != ".tsv") {
    outputFile += ".tsv";
  }

  std::ofstream os(outputFile, std::ios::out);
  if (!os.is_open()) {
    throw std::runtime_error("Failed to open file: " + config.outputFile);
  }

  for (const auto &result : classifyResults) {
    os << result.id << '\t';
    bool handled = false;
    if (!result.taxidCount.empty() &&
        result.taxidCount.front().first == "unclassified") {
      os << "unclassified";
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
    os << '\n';
  }
  os.close();
}

} // namespace ChimeraClassify
