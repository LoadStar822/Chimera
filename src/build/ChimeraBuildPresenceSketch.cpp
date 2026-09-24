#include "ChimeraBuildPresenceSketch.hpp"

#include "BuildTaxonomy.hpp"

#include <utils/LocalResolutionManifest.hpp>
#include <utils/PresenceSketch.hpp>
#include <dna4_traits.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

namespace ChimeraBuild {

namespace {

namespace psk = chimera::presence_sketch;

struct GenomeTask {
  std::string path;
  uint32_t taxid{0};
  uint32_t species{0};
  uintmax_t bytes{0};
};

struct GenomeSketch {
  size_t task{0};
  uint64_t bases{0};
  std::vector<uint64_t> keys; // sorted, unique
};

std::string reference_label(const std::string &path) {
  std::string name = std::filesystem::path(path).filename().string();
  if (name.rfind("GCF_", 0) == 0 || name.rfind("GCA_", 0) == 0) {
    const size_t second = name.find('_', 4);
    if (second != std::string::npos) {
      return name.substr(0, second);
    }
  }
  const size_t dot = name.find('.');
  return dot == std::string::npos ? name : name.substr(0, dot);
}

std::vector<GenomeTask> parse_input(const std::filesystem::path &input) {
  std::ifstream in(input);
  if (!in) {
    throw std::runtime_error("failed to open build input: " + input.string());
  }
  std::vector<GenomeTask> tasks;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    std::istringstream iss(line);
    GenomeTask task;
    std::string taxid_text;
    if (!(iss >> task.path >> taxid_text)) {
      continue;
    }
    unsigned long long taxid = 0;
    try {
      taxid = std::stoull(taxid_text);
    } catch (const std::exception &) {
      continue;
    }
    if (taxid > std::numeric_limits<uint32_t>::max()) {
      throw std::runtime_error("taxid " + taxid_text +
                               " is larger than 4294967295; taxids must fit "
                               "in 32 bits, renumber custom taxa below this "
                               "limit");
    }
    task.taxid = static_cast<uint32_t>(taxid);
    std::error_code ec;
    task.bytes = std::filesystem::file_size(task.path, ec);
    if (ec) {
      task.bytes = 0;
    }
    tasks.push_back(std::move(task));
  }
  return tasks;
}

GenomeSketch sketch_genome(const GenomeTask &task, size_t task_index,
                           const psk::Params &params) {
  GenomeSketch sketch;
  sketch.task = task_index;
  seqan3::sequence_file_input<raptor::dna4_traits,
                              seqan3::fields<seqan3::field::seq>>
      input{task.path};
  for (auto &record : input) {
    const auto &seq = record.sequence();
    sketch.bases += seq.size();
    psk::sample_hashes(seq, params, sketch.keys);
  }
  psk::sort_unique(sketch.keys);
  return sketch;
}

// number of keys in `a` that are not in sorted `b`
size_t count_new_keys(const std::vector<uint64_t> &a,
                      const std::vector<uint64_t> &b) {
  size_t i = 0;
  size_t j = 0;
  size_t novel = 0;
  while (i < a.size()) {
    if (j >= b.size() || a[i] < b[j]) {
      ++novel;
      ++i;
    } else if (a[i] == b[j]) {
      ++i;
      ++j;
    } else {
      ++j;
    }
  }
  return novel;
}

std::vector<uint64_t> merge_keys(const std::vector<uint64_t> &a,
                                 const std::vector<uint64_t> &b) {
  std::vector<uint64_t> out;
  out.reserve(a.size() + b.size());
  std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));
  return out;
}

// greedy selection of up to max_refs genomes by covered union (0 = all)
std::vector<size_t> select_references(const std::vector<GenomeSketch> &sketches,
                                      uint32_t max_refs) {
  std::vector<size_t> selected;
  if (max_refs == 0 || sketches.size() <= max_refs) {
    selected.resize(sketches.size());
    for (size_t i = 0; i < sketches.size(); ++i) {
      selected[i] = i;
    }
    return selected;
  }
  std::vector<bool> used(sketches.size(), false);
  size_t first = 0;
  for (size_t i = 1; i < sketches.size(); ++i) {
    if (sketches[i].keys.size() > sketches[first].keys.size()) {
      first = i;
    }
  }
  selected.push_back(first);
  used[first] = true;
  std::vector<uint64_t> covered = sketches[first].keys;
  while (selected.size() < max_refs) {
    size_t best = sketches.size();
    size_t best_gain = 0;
    for (size_t i = 0; i < sketches.size(); ++i) {
      if (used[i]) {
        continue;
      }
      const size_t gain = count_new_keys(sketches[i].keys, covered);
      if (gain > best_gain) {
        best_gain = gain;
        best = i;
      }
    }
    if (best == sketches.size() || best_gain == 0) {
      break;
    }
    used[best] = true;
    selected.push_back(best);
    covered = merge_keys(covered, sketches[best].keys);
  }
  return selected;
}

psk::SpeciesMarkers build_species_markers(
    uint32_t species, const std::vector<GenomeTask> &tasks,
    const std::vector<GenomeSketch> &sketches, const std::vector<size_t> &selected) {
  psk::SpeciesMarkers markers;
  markers.species = species;
  size_t total = 0;
  for (size_t ref : selected) {
    total += sketches[ref].keys.size();
  }
  // sorted union of all selected sketches
  markers.keys.reserve(total);
  for (size_t ref : selected) {
    const auto &keys = sketches[ref].keys;
    markers.keys.insert(markers.keys.end(), keys.begin(), keys.end());
  }
  psk::sort_unique(markers.keys);
  // per-reference index lists (both sides sorted -> linear merge)
  markers.index.reserve(total);
  markers.ref_begin.assign(1, 0);
  markers.ref_begin.reserve(selected.size() + 1);
  for (size_t ref : selected) {
    const GenomeSketch &sketch = sketches[ref];
    psk::RefInfo info;
    info.bases = sketch.bases;
    info.sketch_size = static_cast<uint32_t>(sketch.keys.size());
    info.name = reference_label(tasks[sketch.task].path);
    markers.refs.push_back(std::move(info));
    size_t u = 0;
    for (uint64_t key : sketch.keys) {
      while (markers.keys[u] < key) {
        ++u;
      }
      markers.index.push_back(static_cast<uint32_t>(u));
    }
    markers.ref_begin.push_back(markers.index.size());
  }
  return markers;
}

} // namespace

std::string default_presence_taxonomy_dir(const std::filesystem::path &db_path) {
  std::vector<std::filesystem::path> candidates;
  if (std::filesystem::is_directory(db_path)) {
    candidates.push_back(db_path / "taxonomy" / "taxdump");
    candidates.push_back(db_path / "taxonomy");
  } else {
    std::filesystem::path core = db_path;
    if (core.extension() != ".imcf") {
      core.replace_extension(".imcf");
    }
    const auto parent = core.parent_path().empty() ? std::filesystem::path(".")
                                                   : core.parent_path();
    candidates.push_back(parent / (core.stem().string() + ".profiledb") /
                         "taxonomy" / "taxdump");
    candidates.push_back(parent / "taxonomy" / "taxdump");
  }
  for (const auto &candidate : candidates) {
    if (std::filesystem::exists(candidate / "nodes.dmp")) {
      return candidate.string();
    }
  }
  return {};
}

PresenceSketchBuildStats
build_presence_sketch(const PresenceSketchBuildOptions &options) {
  const auto started = std::chrono::steady_clock::now();
  PresenceSketchBuildStats stats;
  if (options.scaled == 0) {
    throw std::runtime_error("presence sketch: --scaled must be >= 1");
  }
  BuildConfig taxonomy_config;
  taxonomy_config.taxonomy_dir = options.taxonomy_dir;
  const BuildTaxonomy taxonomy =
      BuildTaxonomy::load_required(taxonomy_config, "presence sketch build");

  std::vector<GenomeTask> tasks = parse_input(options.input_file);
  if (tasks.empty()) {
    throw std::runtime_error("presence sketch build: no genomes in input " +
                             options.input_file.string());
  }
  stats.genomes = tasks.size();
  std::unordered_map<uint32_t, std::vector<size_t>> by_species;
  by_species.reserve(tasks.size());
  for (size_t i = 0; i < tasks.size(); ++i) {
    tasks[i].species = taxonomy.to_species(tasks[i].taxid);
    by_species[tasks[i].species].push_back(i);
  }
  stats.species = by_species.size();

  psk::Params params;
  params.scaled = options.scaled;
  psk::SketchWriter writer(options.output_path, params, options.max_refs,
                           tasks.size());

  std::mutex error_mutex;
  std::exception_ptr first_error;
  std::atomic<bool> stop{false};
  auto record_error = [&](std::exception_ptr error) {
    std::lock_guard<std::mutex> lock(error_mutex);
    if (!first_error) {
      first_error = error;
    }
    stop.store(true);
  };
  auto run_workers = [&](size_t count, const std::function<void()> &body) {
    std::vector<std::thread> threads;
    threads.reserve(count);
    for (size_t i = 0; i < count; ++i) {
      threads.emplace_back(body);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    if (first_error) {
      std::rethrow_exception(first_error);
    }
  };

  // phase 1: sketch every genome, largest files first
  std::vector<size_t> genome_order(tasks.size());
  for (size_t i = 0; i < tasks.size(); ++i) {
    genome_order[i] = i;
  }
  std::sort(genome_order.begin(), genome_order.end(), [&](size_t a, size_t b) {
    return tasks[a].bytes > tasks[b].bytes;
  });
  std::vector<GenomeSketch> genome_sketches(tasks.size());
  std::vector<char> readable(tasks.size(), 0);
  std::atomic<size_t> next_genome{0};
  std::atomic<size_t> done_genomes{0};
  std::atomic<uint64_t> without_markers{0};
  std::atomic<uint64_t> total_bases{0};
  const size_t genome_workers =
      std::max<size_t>(1, std::min<size_t>(options.threads, tasks.size()));
  run_workers(genome_workers, [&]() {
    while (!stop.load(std::memory_order_relaxed)) {
      const size_t position = next_genome.fetch_add(1);
      if (position >= genome_order.size()) {
        break;
      }
      const size_t idx = genome_order[position];
      try {
        GenomeSketch sketch = sketch_genome(tasks[idx], idx, params);
        if (sketch.bases == 0 || sketch.keys.empty()) {
          without_markers.fetch_add(1);
        } else {
          total_bases.fetch_add(sketch.bases);
          genome_sketches[idx] = std::move(sketch);
          readable[idx] = 1;
        }
      } catch (const std::exception &ex) {
        record_error(std::make_exception_ptr(std::runtime_error(
            "presence sketch build: cannot read reference " + tasks[idx].path +
            ": " + ex.what())));
      } catch (...) {
        record_error(std::current_exception());
      }
      const size_t finished = done_genomes.fetch_add(1) + 1;
      if (options.verbose && (finished % 10000 == 0)) {
        std::cout << "  presence sketch: " << finished << "/" << tasks.size()
                  << " genomes" << std::endl;
      }
    }
  });

  // phase 2: per species, choose representatives and write the marker block
  std::vector<uint32_t> species_order;
  species_order.reserve(by_species.size());
  for (const auto &[species, members] : by_species) {
    species_order.push_back(species);
  }
  std::sort(species_order.begin(), species_order.end(),
            [&](uint32_t a, uint32_t b) {
              return by_species[a].size() > by_species[b].size();
            });
  std::atomic<size_t> next_species{0};
  std::atomic<uint64_t> capped{0};
  const size_t species_workers =
      std::max<size_t>(1, std::min<size_t>(options.threads, species_order.size()));
  run_workers(species_workers, [&]() {
    while (!stop.load(std::memory_order_relaxed)) {
      const size_t index = next_species.fetch_add(1);
      if (index >= species_order.size()) {
        break;
      }
      const uint32_t species = species_order[index];
      try {
        std::vector<GenomeSketch> sketches;
        for (size_t idx : by_species[species]) {
          if (readable[idx]) {
            sketches.push_back(std::move(genome_sketches[idx]));
          }
        }
        if (sketches.empty()) {
          continue;
        }
        if (options.max_refs != 0 && sketches.size() > options.max_refs) {
          capped.fetch_add(1);
        }
        const std::vector<size_t> selected =
            select_references(sketches, options.max_refs);
        writer.add_species(
            build_species_markers(species, tasks, sketches, selected));
      } catch (...) {
        record_error(std::current_exception());
      }
    }
  });
  writer.finish();
  stats.genomes_without_markers = without_markers.load();
  stats.bases = total_bases.load();
  stats.species_capped = capped.load();
  stats.keys = writer.keys_written();
  stats.references = writer.refs_written();
  stats.species = writer.species_written();
  if (!options.database.empty()) {
    stats.registered = chimera::local_resolution::stamp_presence_sketch(
        chimera::local_resolution::core_archive_path_for(options.database),
        options.output_path);
  }
  stats.seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  if (options.verbose) {
    std::cout << "presence sketch written: " << options.output_path.string()
              << "\n  genomes      " << stats.genomes
              << (stats.genomes_without_markers > 0
                      ? " (" + std::to_string(stats.genomes_without_markers) +
                            " without markers)"
                      : "")
              << "\n  species      " << stats.species
              << (stats.species_capped > 0
                      ? " (" + std::to_string(stats.species_capped) +
                            " capped at " + std::to_string(options.max_refs) +
                            " references)"
                      : "")
              << "\n  references   " << stats.references
              << "\n  markers      " << stats.keys
              << "\n  scaled       " << options.scaled
              << (options.database.empty()
                      ? ""
                      : stats.registered
                            ? "\n  manifest     registered"
                            : "\n  manifest     not registered (database has "
                              "no manifest)")
              << "\n  time         " << static_cast<uint64_t>(stats.seconds)
              << "s" << std::endl;
  }
  return stats;
}

} // namespace ChimeraBuild
