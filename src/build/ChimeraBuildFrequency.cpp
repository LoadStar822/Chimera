#include "ChimeraBuildCommon.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>

#include <dna4_traits.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace ChimeraBuild {

static chimera::presence::HashFreqStats
compute_hash_freq_stats(const std::vector<uint32_t> &counters, double quantile) {
  chimera::presence::HashFreqStats stats{};
  std::array<uint64_t, 32> histogram{};
  for (uint32_t value : counters) {
    if (value == 0) {
      continue;
    }
    ++stats.nonzero_counters;
    if (value > stats.df_max_observed) {
      stats.df_max_observed = value;
    }
    int bucket = 0;
    if (value > 1) {
      bucket = static_cast<int>(std::log2(static_cast<double>(value)));
      if (bucket > 31) {
        bucket = 31;
      }
    }
    histogram[static_cast<size_t>(bucket)]++;
  }
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

static void build_threshold_masks(const std::vector<uint32_t> &counters,
                                  uint32_t depth, uint32_t width,
                                  uint32_t high_threshold,
                                  uint32_t unique_threshold,
                                  CmsThresholdBitmask &high_mask,
                                  CmsThresholdBitmask &unique_mask) {
  const bool high_enabled =
      high_threshold != std::numeric_limits<uint32_t>::max();
  if (!high_enabled) {
    high_mask.clear();
  } else {
    high_mask.reset(depth, width);
  }
  unique_mask.reset(depth, width);
  for (uint32_t row = 0; row < depth; ++row) {
    const size_t row_base =
        static_cast<size_t>(row) * static_cast<size_t>(width);
    for (uint32_t column = 0; column < width; ++column) {
      const uint32_t value = counters[row_base + column];
      if (high_enabled && value >= high_threshold) {
        high_mask.set(row, column);
      }
      if (value > unique_threshold) {
        unique_mask.set(row, column);
      }
    }
  }
}

void build_hash_frequency_sketch(
    const BuildConfig &config,
    const robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    HashFrequencyContext &context) {
  constexpr uint32_t kSketchDepth = 4;
  constexpr uint32_t kSketchWidth = 1u << 22;
  context.passA_total_hashes.store(0, std::memory_order_relaxed);
  context.passB_total_hashes.store(0, std::memory_order_relaxed);
  context.passB_filtered_hashes.store(0, std::memory_order_relaxed);
  context.stats = {};
  context.quantile = 0.999;
  context.unique_deg_threshold =
      std::max<uint32_t>(1u, config.presence_unique_deg);
  context.high_df_mask.clear();
  context.gt_unique_mask.clear();
  std::cout << "Building hash frequency sketch (depth=" << kSketchDepth
            << ", width=" << kSketchWidth << ")..." << std::endl;
  context.sketch = std::make_unique<CountMinSketch>(kSketchDepth, kSketchWidth);
  std::vector<std::string> all_files;
  for (const auto &entry : inputFiles) {
    for (const auto &filename : entry.second) {
      all_files.push_back(filename);
    }
  }
  if (all_files.empty()) {
    context.stats = {};
    return;
  }
  uint64_t feature_seed = 0;
  auto feature_params = make_feature_params(config, feature_seed);
  const size_t feature_min_length =
      chimera::feature::min_required_length(feature_params);
  const size_t min_required =
      std::max<size_t>(config.min_length, feature_min_length);
  auto sketch_start = std::chrono::high_resolution_clock::now();

#pragma omp parallel
  {
    std::vector<uint64_t> hashes;
    hashes.reserve(4096);
    chimera::feature::FeatureHashScratch featureScratch;
    uint64_t local_hash_total = 0;

#pragma omp for schedule(dynamic)
    for (size_t idx = 0; idx < all_files.size(); ++idx) {
      const std::string &filename = all_files[idx];
      try {
        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            fin{filename};
        for (auto &record : fin) {
          auto &seq = record.sequence();
          if (seq.size() < min_required) {
            continue;
          }
          hashes.clear();
          chimera::feature::compute_hashes_append(seq, feature_params, hashes,
                                                  featureScratch);
          local_hash_total += hashes.size();
          context.sketch->add_many(hashes);
        }
      } catch (const std::exception &ex) {
#pragma omp critical(sketch_log)
        {
          std::cerr << "Failed to read sequence file for sketch: " << filename
                    << " (" << ex.what() << ")" << std::endl;
        }
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
  context.sketch->freeze();
  std::vector<uint32_t> sketch_counters;
  context.sketch->exportCounts(sketch_counters);
  context.stats = compute_hash_freq_stats(sketch_counters, context.quantile);
  build_threshold_masks(sketch_counters, kSketchDepth, kSketchWidth,
                        context.stats.df_high_threshold,
                        context.unique_deg_threshold, context.high_df_mask,
                        context.gt_unique_mask);
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
