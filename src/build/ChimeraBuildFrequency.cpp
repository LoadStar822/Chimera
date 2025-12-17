#include "ChimeraBuildCommon.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>

#include <dna4_traits.hpp>
#include <robin_hood.h>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace ChimeraBuild {

static chimera::presence::HashFreqStats
compute_hash_freq_stats(const CountMinSketch &cms, double quantile) {
  chimera::presence::HashFreqStats stats{};
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
  uint64_t target = static_cast<uint64_t>(
      std::ceil(tail_fraction * static_cast<double>(stats.nonzero_counters)));
  if (target == 0) {
    target = 1;
  }
  uint64_t accumulated = 0;
  for (int bucket = static_cast<int>(histogram.size()) - 1; bucket >= 0;
       --bucket) {
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

void build_hash_frequency_sketch(
    const BuildConfig &config,
    const robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    HashFrequencyContext &context) {
  constexpr uint32_t kSketchDepth = 4;
  constexpr uint32_t kSketchWidth = 1u << 23;
  constexpr size_t kSketchFileCap = 16;
  constexpr size_t kDocSigSize = 4096;
  context.passA_total_hashes.store(0, std::memory_order_relaxed);
  context.passB_total_hashes.store(0, std::memory_order_relaxed);
  context.passB_filtered_hashes.store(0, std::memory_order_relaxed);
  context.passB_sampled_out_hashes.store(0, std::memory_order_relaxed);
  context.stats = {};
  context.quantile = 0.999;
  std::cout << "Building hash frequency sketch (depth=" << kSketchDepth
            << ", width=" << kSketchWidth << ", taxon-docsig)";
  std::cout << "..." << std::endl;
  context.sketch = std::make_unique<CountMinSketch>(kSketchDepth, kSketchWidth);
  if (inputFiles.empty()) {
    context.stats = {};
    return;
  }
  struct TaxonFiles {
    const std::string *taxid;
    const std::vector<std::string> *files;
  };
  std::vector<TaxonFiles> taxa;
  taxa.reserve(inputFiles.size());
  for (const auto &entry : inputFiles) {
    taxa.push_back({&entry.first, &entry.second});
  }
  chimera::feature::Method feature_method{};
  uint64_t feature_seed = 0;
  auto feature_params =
      make_feature_params(config, feature_method, feature_seed);
  const size_t feature_min_length =
      chimera::feature::min_required_length(feature_params);
  const size_t min_required =
      std::max<size_t>(config.min_length, feature_min_length);
  const uint64_t rank_salt = feature_seed ^ 0xD1B54A32D192ED03ULL;
  auto sketch_start = std::chrono::high_resolution_clock::now();

  const uint16_t sketch_threads =
      std::max<uint16_t>(1, std::min<uint16_t>(config.threads, 32));
#pragma omp parallel num_threads(sketch_threads)
  {
    std::vector<const std::string *> file_ptrs;
    std::vector<const std::string *> files_used;
    std::vector<uint64_t> taxon_hashes;
    taxon_hashes.reserve(kSketchFileCap * kDocSigSize);
    struct SigItem {
      uint64_t rank;
      uint64_t hash;
    };
    auto worse_rank = [](const SigItem &a, const SigItem &b) {
      if (a.rank != b.rank) {
        return a.rank < b.rank; // max-heap by rank (worst on top)
      }
      return a.hash < b.hash;
    };
    std::vector<SigItem> heap;
    heap.reserve(kDocSigSize);
    robin_hood::unordered_flat_set<uint64_t> heap_seen;
    heap_seen.reserve(kDocSigSize * 2 + 8);
    std::vector<uint64_t> file_sig;
    file_sig.reserve(kDocSigSize);
    std::vector<uint64_t> hashes;
    hashes.reserve(4096);
    uint64_t local_hash_total = 0;

#pragma omp for schedule(dynamic)
    for (size_t idx = 0; idx < taxa.size(); ++idx) {
      const auto &taxon = taxa[idx];
      if (taxon.files == nullptr || taxon.files->empty()) {
        continue;
      }
      file_ptrs.clear();
      file_ptrs.reserve(taxon.files->size());
      for (const auto &filename : *taxon.files) {
        file_ptrs.push_back(&filename);
      }
      std::sort(file_ptrs.begin(), file_ptrs.end(),
                [](const std::string *a, const std::string *b) {
                  return *a < *b;
                });
      files_used.clear();
      if (file_ptrs.size() <= kSketchFileCap) {
        files_used = file_ptrs;
      } else if (kSketchFileCap == 1) {
        files_used.push_back(file_ptrs.front());
      } else {
        files_used.reserve(kSketchFileCap);
        const size_t n = file_ptrs.size();
        for (size_t i = 0; i < kSketchFileCap; ++i) {
          size_t pick = (i * (n - 1)) / (kSketchFileCap - 1);
          files_used.push_back(file_ptrs[pick]);
        }
        std::sort(files_used.begin(), files_used.end(),
                  [](const std::string *a, const std::string *b) {
                    return *a < *b;
                  });
        files_used.erase(std::unique(files_used.begin(), files_used.end()),
                         files_used.end());
      }

      taxon_hashes.clear();
      for (const std::string *filename_ptr : files_used) {
        if (!filename_ptr) {
          continue;
        }
        const std::string &filename = *filename_ptr;
        try {
          seqan3::sequence_file_input<
              raptor::dna4_traits,
              seqan3::fields<seqan3::field::id, seqan3::field::seq>>
              fin{filename};
          heap.clear();
          heap_seen.clear();
          for (auto &record : fin) {
            auto &seq = record.sequence();
            if (seq.size() < min_required) {
              continue;
            }
            hashes.clear();
            chimera::feature::compute_hashes_append(seq, feature_params, hashes);
            if (hashes.empty()) {
              continue;
            }
            for (uint64_t hash : hashes) {
              if (heap_seen.contains(hash)) {
                continue;
              }
              uint64_t rank = mix64(hash ^ rank_salt);
              SigItem cand{rank, hash};
              if (heap.size() < kDocSigSize) {
                heap.push_back(cand);
                std::push_heap(heap.begin(), heap.end(), worse_rank);
                heap_seen.insert(hash);
                continue;
              }
              const SigItem &worst = heap.front();
              if (rank < worst.rank || (rank == worst.rank && hash < worst.hash)) {
                std::pop_heap(heap.begin(), heap.end(), worse_rank);
                heap_seen.erase(heap.back().hash);
                heap.back() = cand;
                std::push_heap(heap.begin(), heap.end(), worse_rank);
                heap_seen.insert(hash);
              }
            }
          }
          if (!heap.empty()) {
            file_sig.clear();
            file_sig.reserve(heap.size());
            for (const auto &item : heap) {
              file_sig.push_back(item.hash);
            }
            std::sort(file_sig.begin(), file_sig.end());
            file_sig.erase(std::unique(file_sig.begin(), file_sig.end()),
                           file_sig.end());
            taxon_hashes.insert(taxon_hashes.end(), file_sig.begin(),
                                file_sig.end());
          }
        } catch (const std::exception &ex) {
#pragma omp critical(sketch_log)
          {
            std::cerr << "Failed to read sequence file for sketch: "
                      << filename;
            if (taxon.taxid) {
              std::cerr << " taxid=" << *taxon.taxid;
            }
            std::cerr << " (" << ex.what() << ")" << std::endl;
          }
        }
      }

      if (!taxon_hashes.empty()) {
        std::sort(taxon_hashes.begin(), taxon_hashes.end());
        taxon_hashes.erase(std::unique(taxon_hashes.begin(), taxon_hashes.end()),
                           taxon_hashes.end());
        for (uint64_t hash : taxon_hashes) {
          context.sketch->add(hash);
        }
        local_hash_total += static_cast<uint64_t>(taxon_hashes.size());
      }
    }

    if (local_hash_total > 0) {
      context.passA_total_hashes.fetch_add(local_hash_total,
                                           std::memory_order_relaxed);
    }
  }
  auto sketch_end = std::chrono::high_resolution_clock::now();
  auto sketch_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                         sketch_end - sketch_start)
                         .count();
  context.stats = compute_hash_freq_stats(*context.sketch, context.quantile);
  std::cout << "Hash frequency sketch time: " << sketch_time / 1000 << "s "
            << sketch_time % 1000 << "ms" << std::endl;
  const auto hashed =
      context.passA_total_hashes.load(std::memory_order_relaxed);
  std::cout << "Sketch summary:" << std::endl;
  std::cout << "  Streamed hashes: " << hashed << std::endl;
  std::cout << "  Active counters: " << context.stats.nonzero_counters
            << std::endl;
  std::cout << "  Max df estimate: " << context.stats.df_max_observed
            << std::endl;
  if (context.stats.df_high_threshold ==
      std::numeric_limits<uint32_t>::max()) {
    std::cout << "  BasicFilter threshold: disabled (insufficient data)"
              << std::endl;
  } else {
    std::cout << "  BasicFilter threshold: df >= "
              << context.stats.df_high_threshold
              << " (quantile=" << context.quantile << ")" << std::endl;
  }
}

} // namespace ChimeraBuild
