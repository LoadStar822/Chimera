#include "ChimeraBuildNativeBounded.hpp"

#include <dna4_traits.hpp>
#include <utils/LocalResolutionMetadata.hpp>
#include <utils/NativeBoundedIndex.hpp>

#include <algorithm>
#include <atomic>
#include <array>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <unordered_set>
#include <unordered_map>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <stdexcept>
#include <string>
#include <sstream>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace ChimeraBuild {
namespace {

struct InputTask {
  uint32_t species{};
  std::string filename;
};

std::string make_target_name(size_t source_id, uint32_t contig,
                             uint32_t species) {
  return "src=" + std::to_string(source_id) + "|contig=" +
         std::to_string(contig) + "|sp=" + std::to_string(species);
}

uint32_t parse_taxid_u32(const std::string &taxid) {
  try {
    return static_cast<uint32_t>(std::stoul(taxid));
  } catch (...) {
    return 0;
  }
}

struct Taxdump {
  std::vector<uint32_t> parent;
  std::vector<uint8_t> is_species;
  std::vector<uint8_t> is_genus;

  uint32_t to_species(uint32_t tid) const {
    if (tid == 0 || tid >= is_species.size()) {
      return tid;
    }
    if (is_species[tid]) {
      return tid;
    }
    uint32_t cur = tid;
    for (int steps = 0; steps < 128; ++steps) {
      if (cur == 0 || cur >= parent.size()) {
        break;
      }
      const uint32_t p = parent[cur];
      if (p == 0 || p == cur) {
        break;
      }
      cur = p;
      if (cur < is_species.size() && is_species[cur]) {
        return cur;
      }
    }
    return tid;
  }

  uint32_t to_genus(uint32_t tid) const {
    if (tid == 0 || tid >= parent.size()) {
      return 0;
    }
    uint32_t cur = tid;
    for (int steps = 0; steps < 128; ++steps) {
      if (cur == 0 || cur >= parent.size()) {
        break;
      }
      if (cur < is_genus.size() && is_genus[cur]) {
        return cur;
      }
      const uint32_t p = parent[cur];
      if (p == 0 || p == cur) {
        break;
      }
      cur = p;
    }
    return 0;
  }
};

std::optional<Taxdump> load_taxdump(const BuildConfig &config) {
  std::vector<std::filesystem::path> roots;
  if (!config.taxonomy_dir.empty()) {
    roots.emplace_back(config.taxonomy_dir);
  }
  if (const char *env_dir = std::getenv("CHIMERA_NCBI_TAXDUMP_DIR")) {
    if (*env_dir) {
      roots.emplace_back(env_dir);
    }
  }

  std::filesystem::path nodes;
  for (const auto &root : roots) {
    const auto candidate = root / "nodes.dmp";
    if (std::filesystem::exists(candidate)) {
      nodes = candidate;
      break;
    }
  }
  if (nodes.empty()) {
    return std::nullopt;
  }
  std::ifstream in(nodes);
  if (!in) {
    throw std::runtime_error("failed to open taxonomy nodes.dmp: " +
                             nodes.string());
  }

  Taxdump tax;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, '\t')) {
      while (!tok.empty() &&
             std::isspace(static_cast<unsigned char>(tok.front()))) {
        tok.erase(tok.begin());
      }
      while (!tok.empty() &&
             std::isspace(static_cast<unsigned char>(tok.back()))) {
        tok.pop_back();
      }
      if (tok.empty() || tok == "|") {
        continue;
      }
      fields.push_back(tok);
      if (fields.size() >= 3) {
        break;
      }
    }
    if (fields.size() < 3) {
      continue;
    }
    const uint32_t tid = parse_taxid_u32(fields[0]);
    const uint32_t parent = parse_taxid_u32(fields[1]);
    if (tid == 0) {
      continue;
    }
    if (tid >= tax.parent.size()) {
      tax.parent.resize(static_cast<size_t>(tid) + 1, 0);
      tax.is_species.resize(static_cast<size_t>(tid) + 1, 0);
      tax.is_genus.resize(static_cast<size_t>(tid) + 1, 0);
    }
    tax.parent[tid] = parent;
    tax.is_species[tid] = (fields[2] == "species") ? 1 : 0;
    tax.is_genus[tid] = (fields[2] == "genus") ? 1 : 0;
  }
  if (tax.parent.empty()) {
    return std::nullopt;
  }
  return tax;
}

Taxdump load_required_taxdump_for_local_resolution(const BuildConfig &config) {
  auto taxdump = load_taxdump(config);
  if (!taxdump.has_value()) {
    throw std::runtime_error(
        "local read resolution build requires taxonomy nodes.dmp; provide "
        "--taxonomy-dir or CHIMERA_NCBI_TAXDUMP_DIR");
  }
  return std::move(*taxdump);
}

std::vector<InputTask> make_tasks(
    const robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles) {
  std::vector<InputTask> tasks;
  for (const auto &[taxid, files] : inputFiles) {
    const uint32_t species = parse_taxid_u32(taxid);
    if (species == 0) {
      continue;
    }
    for (const auto &file : files) {
      tasks.push_back({species, file});
    }
  }
  std::sort(tasks.begin(), tasks.end(), [](const auto &lhs, const auto &rhs) {
    if (lhs.species != rhs.species) {
      return lhs.species < rhs.species;
    }
    return lhs.filename < rhs.filename;
  });
  return tasks;
}

std::string shard_filename(uint32_t genus) {
  return "g" + std::to_string(genus) + ".nbcidx";
}

struct TargetPlan {
  size_t source_index{};
  uint32_t species{};
  uint32_t genus{};
  uint32_t contig_index{};
  uint32_t len{};
  uint32_t anchor_count{};
  uint64_t anchor_offset{};
  uint64_t anchor_bytes{};
  std::string name;
};

bool target_plan_rank_less(const TargetPlan &lhs, const TargetPlan &rhs) {
  if (lhs.anchor_count != rhs.anchor_count) {
    return lhs.anchor_count > rhs.anchor_count;
  }
  if (lhs.len != rhs.len) {
    return lhs.len > rhs.len;
  }
  return lhs.name < rhs.name;
}

std::string local_resolution_source_key(const std::string &targetName) {
  const size_t pos = targetName.find("|contig=");
  if (pos == std::string::npos) {
    return targetName;
  }
  return targetName.substr(0, pos);
}

struct EncodedTarget {
  TargetPlan plan;
  std::vector<char> anchors;
};

size_t worst_target_index(const std::vector<EncodedTarget> &targets) {
  size_t worst = 0;
  for (size_t i = 1; i < targets.size(); ++i) {
    if (target_plan_rank_less(targets[worst].plan, targets[i].plan)) {
      worst = i;
    }
  }
  return worst;
}

bool target_would_be_retained(const std::vector<EncodedTarget> &targets,
                              const TargetPlan &plan, uint32_t cap,
                              size_t &replaceIndex) {
  if (cap == 0 || targets.size() < cap) {
    replaceIndex = targets.size();
    return true;
  }
  replaceIndex = worst_target_index(targets);
  return target_plan_rank_less(plan, targets[replaceIndex].plan);
}

void retain_target(std::vector<EncodedTarget> &targets,
                   EncodedTarget target, uint32_t cap) {
  size_t replaceIndex = 0;
  if (!target_would_be_retained(targets, target.plan, cap, replaceIndex)) {
    return;
  }
  if (replaceIndex == targets.size()) {
    targets.push_back(std::move(target));
  } else {
    targets[replaceIndex] = std::move(target);
  }
}

std::vector<EncodedTarget>
finalize_species_targets(std::vector<EncodedTarget> targets,
                         uint32_t sourceCap, uint32_t targetCap) {
  std::sort(targets.begin(), targets.end(),
            [](const EncodedTarget &lhs, const EncodedTarget &rhs) {
              return target_plan_rank_less(lhs.plan, rhs.plan);
            });
  if (sourceCap > 0) {
    std::unordered_set<std::string> selectedSources;
    selectedSources.reserve(sourceCap);
    targets.erase(
        std::remove_if(
            targets.begin(), targets.end(), [&](const EncodedTarget &target) {
              const std::string source =
                  local_resolution_source_key(target.plan.name);
              if (selectedSources.contains(source)) {
                return false;
              }
              if (selectedSources.size() >= sourceCap) {
                return true;
              }
              selectedSources.insert(source);
              return false;
            }),
        targets.end());
  }
  if (targetCap > 0) {
    std::unordered_map<std::string, uint32_t> sourceCounts;
    targets.erase(
        std::remove_if(
            targets.begin(), targets.end(), [&](const EncodedTarget &target) {
              const std::string source =
                  local_resolution_source_key(target.plan.name);
              return sourceCounts[source]++ >= targetCap;
            }),
        targets.end());
  }
  std::sort(targets.begin(), targets.end(),
            [](const EncodedTarget &lhs, const EncodedTarget &rhs) {
              if (lhs.plan.source_index != rhs.plan.source_index) {
                return lhs.plan.source_index < rhs.plan.source_index;
              }
              return lhs.plan.contig_index < rhs.plan.contig_index;
            });
  return targets;
}

struct ShardOutput {
  std::filesystem::path path;
  std::atomic<uint64_t> next_offset{0};

  explicit ShardOutput(std::filesystem::path outputPath)
      : path(std::move(outputPath)) {}
};

void write_all_at(int fd, const char *data, size_t size, uint64_t offset,
                  const std::filesystem::path &path) {
  size_t written = 0;
  while (written < size) {
    const ssize_t rc = ::pwrite(
        fd, data + written, size - written,
        static_cast<off_t>(offset + static_cast<uint64_t>(written)));
    if (rc < 0) {
      if (errno == EINTR) {
        continue;
      }
      throw std::runtime_error("failed to write native bounded shard " +
                               path.string() + " errno=" +
                               std::to_string(errno) + " " +
                               std::strerror(errno));
    }
    if (rc == 0) {
      throw std::runtime_error("short native bounded shard write: " +
                               path.string());
    }
    written += static_cast<size_t>(rc);
  }
}

class ShardFdCache {
public:
  ShardFdCache(
      const std::unordered_map<uint32_t, std::unique_ptr<ShardOutput>> &outputs,
      size_t max_open)
      : outputs_(outputs), max_open_(std::max<size_t>(1, max_open)) {}

  int fd_for(uint32_t genus) {
    auto found = fds_.find(genus);
    if (found != fds_.end()) {
      return found->second;
    }
    if (fds_.size() >= max_open_) {
      close_oldest();
    }
    auto output = outputs_.find(genus);
    if (output == outputs_.end()) {
      throw std::runtime_error(
          "native bounded direct build missing shard output");
    }
    const int fd =
        ::open(output->second->path.c_str(), O_WRONLY | O_CREAT, 0644);
    if (fd < 0) {
      throw std::runtime_error("failed to open native bounded shard for "
                               "direct write: " +
                               output->second->path.string() + " errno=" +
                               std::to_string(errno) + " " +
                               std::strerror(errno));
    }
    fds_.emplace(genus, fd);
    order_.push_back(genus);
    return fd;
  }

  void close_all() noexcept {
    for (auto &[_, fd] : fds_) {
      if (fd >= 0) {
        ::close(fd);
      }
    }
    fds_.clear();
    order_.clear();
  }

  ~ShardFdCache() { close_all(); }

private:
  void close_oldest() {
    while (!order_.empty()) {
      const uint32_t genus = order_.front();
      order_.pop_front();
      auto found = fds_.find(genus);
      if (found == fds_.end()) {
        continue;
      }
      ::close(found->second);
      fds_.erase(found);
      return;
    }
  }

  const std::unordered_map<uint32_t, std::unique_ptr<ShardOutput>> &outputs_;
  size_t max_open_;
  std::unordered_map<uint32_t, int> fds_;
  std::deque<uint32_t> order_;
};

uint64_t write_representative_pool_for_genus(
    std::ostream *out, uint32_t genus,
    const std::vector<chimera::native_bounded::TargetMeta> &targets,
    uint32_t cap, uint64_t anchorDataOffset, uint32_t k,
    std::vector<chimera::local_resolution::TargetRep> &metadataRows) {
  if (targets.empty()) {
    return 0;
  }
  std::unordered_map<uint32_t, std::vector<size_t>> bySpecies;
  bySpecies.reserve(targets.size());
  for (size_t i = 0; i < targets.size(); ++i) {
    bySpecies[targets[i].species].push_back(i);
  }
  std::vector<uint32_t> speciesIds;
  speciesIds.reserve(bySpecies.size());
  for (const auto &[species, _] : bySpecies) {
    speciesIds.push_back(species);
  }
  std::sort(speciesIds.begin(), speciesIds.end());
  for (uint32_t species : speciesIds) {
    auto &group = bySpecies[species];
    std::sort(group.begin(), group.end(), [&](size_t lhs, size_t rhs) {
      const auto &lt = targets[lhs];
      const auto &rt = targets[rhs];
      if (lt.anchor_count != rt.anchor_count) {
        return lt.anchor_count > rt.anchor_count;
      }
      if (lt.len != rt.len) {
        return lt.len > rt.len;
      }
      return lt.name < rt.name;
    });
  }

  const uint64_t limit =
      cap == 0 ? std::numeric_limits<uint64_t>::max() : cap;
  std::unordered_map<uint32_t, size_t> cursor;
  cursor.reserve(speciesIds.size());
  uint64_t written = 0;
  while (written < limit) {
    bool advanced = false;
    for (uint32_t species : speciesIds) {
      auto &idx = cursor[species];
      const auto &group = bySpecies[species];
      if (idx >= group.size()) {
        continue;
      }
      const size_t targetIndex = group[idx++];
      const auto &target = targets[targetIndex];
      const uint64_t anchorByteOffset = anchorDataOffset + target.anchor_offset;
      if (out != nullptr) {
        *out << genus << '\t' << target.species << '\t' << targetIndex << '\t'
             << target.name << '\t' << target.len << '\t'
             << target.anchor_count << '\t' << target.anchor_offset << '\t'
             << anchorByteOffset << '\t' << target.anchor_bytes << '\t'
             << (written + 1) << '\t'
             << "species_round_robin_anchor_desc" << '\n';
        if (!out->good()) {
          throw std::runtime_error("failed to write native bounded rep pool");
        }
      }
      metadataRows.push_back(chimera::local_resolution::TargetRep{
          genus,
          target.species,
          target.name,
          target.len,
          target.anchor_count,
          anchorByteOffset,
          target.anchor_bytes,
          static_cast<uint32_t>(written + 1),
      });
      ++written;
      advanced = true;
      if (written >= limit) {
        break;
      }
    }
    if (!advanced) {
      break;
    }
  }
  return written;
}

struct SpeciesBuildState {
  uint32_t species{};
  uint32_t genus{};
  size_t remaining_tasks{};
  std::mutex mutex;
  std::vector<EncodedTarget> candidates;
  std::vector<TargetPlan> selected;
};

struct ResolvedInputTask {
  size_t species_state{};
  uint32_t species{};
  uint32_t genus{};
};

NativeBoundedBuildStats build_native_bounded_index_fused(
    const BuildConfig &config, const std::vector<InputTask> &tasks,
    const NativeBoundedOutputPaths &paths) {
  const auto anchorBuildStarted = std::chrono::steady_clock::now();
  std::filesystem::create_directories(paths.metadata_index.parent_path().empty()
                                          ? std::filesystem::path(".")
                                          : paths.metadata_index.parent_path());
  std::filesystem::create_directories(paths.rep_metadata.parent_path().empty()
                                          ? std::filesystem::path(".")
                                          : paths.rep_metadata.parent_path());
  std::filesystem::create_directories(paths.shard_manifest.parent_path().empty()
                                          ? std::filesystem::path(".")
                                          : paths.shard_manifest.parent_path());
  std::filesystem::remove_all(paths.shard_dir);
  std::filesystem::create_directories(paths.shard_dir);

  const auto taxdump = load_required_taxdump_for_local_resolution(config);
  NativeBoundedBuildStats stats;
  const size_t workerCount =
      std::max<size_t>(1, std::min<size_t>(config.threads, tasks.size()));

  std::vector<std::unique_ptr<SpeciesBuildState>> speciesStates;
  std::unordered_map<uint32_t, size_t> speciesStateIndex;
  std::vector<ResolvedInputTask> resolvedTasks(tasks.size());
  speciesStateIndex.reserve(tasks.size());
  for (size_t taskIndex = 0; taskIndex < tasks.size(); ++taskIndex) {
    const uint32_t species = taxdump.to_species(tasks[taskIndex].species);
    uint32_t genus = taxdump.to_genus(species);
    if (genus == 0) {
      genus = species;
    }
    auto [found, inserted] =
        speciesStateIndex.emplace(species, speciesStates.size());
    if (inserted) {
      auto state = std::make_unique<SpeciesBuildState>();
      state->species = species;
      state->genus = genus;
      speciesStates.push_back(std::move(state));
    } else if (speciesStates[found->second]->genus != genus) {
      throw std::runtime_error(
          "native bounded species maps to multiple genera");
    }
    auto &state = *speciesStates[found->second];
    ++state.remaining_tasks;
    resolvedTasks[taskIndex] =
        ResolvedInputTask{found->second, species, genus};
  }

  std::unordered_map<uint32_t, std::unique_ptr<ShardOutput>> shardOutputs;
  shardOutputs.reserve(speciesStates.size());
  for (const auto &state : speciesStates) {
    if (shardOutputs.contains(state->genus)) {
      continue;
    }
    const std::filesystem::path path =
        paths.shard_dir / shard_filename(state->genus);
    std::filesystem::remove(path);
    shardOutputs.emplace(state->genus,
                         std::make_unique<ShardOutput>(path));
  }

  std::atomic<size_t> nextTask{0};
  std::atomic<bool> stopWorkers{false};
  std::exception_ptr workerError;
  std::mutex errorMutex;
  const size_t shardFdCacheMaxOpen =
      std::min<size_t>(64, std::max<size_t>(8, shardOutputs.size()));
  auto buildWorker = [&]() {
    ShardFdCache fdCache(shardOutputs, shardFdCacheMaxOpen);
    while (true) {
      if (stopWorkers.load(std::memory_order_relaxed)) {
        break;
      }
      const size_t taskIndex = nextTask.fetch_add(1);
      if (taskIndex >= tasks.size()) {
        break;
      }
      try {
        const auto &task = tasks[taskIndex];
        const auto &resolved = resolvedTasks[taskIndex];
        std::vector<EncodedTarget> taskCandidates;
        taskCandidates.reserve(config.native_bounded_targets_per_species == 0
                                   ? 32
                                   : config.native_bounded_targets_per_species);
        seqan3::sequence_file_input<
            raptor::dna4_traits,
            seqan3::fields<seqan3::field::id, seqan3::field::seq>>
            input{task.filename};
        uint32_t contig = 0;
        for (auto &record : input) {
          auto &seq = record.sequence();
          (void)record.id();
          TargetPlan plan;
          plan.source_index = taskIndex;
          plan.species = resolved.species;
          plan.genus = resolved.genus;
          plan.contig_index = contig;
          plan.len = static_cast<uint32_t>(seq.size());
          plan.name =
              make_target_name(taskIndex, contig, resolved.species);
          auto anchors = chimera::native_bounded::extract_minimizers(
              seq, config.native_bounded_k, config.native_bounded_w);
          plan.anchor_count = static_cast<uint32_t>(anchors.size());

          size_t replaceIndex = 0;
          if (target_would_be_retained(
                  taskCandidates, plan,
                  config.native_bounded_targets_per_species, replaceIndex)) {
            auto encoded = chimera::native_bounded::encode_anchor_block(
                anchors, config.native_bounded_k);
            plan.anchor_bytes = encoded.size();
            retain_target(
                taskCandidates,
                EncodedTarget{std::move(plan), std::move(encoded)},
                config.native_bounded_targets_per_species);
          }
          ++contig;
        }

        auto &speciesState = *speciesStates[resolved.species_state];
        std::vector<EncodedTarget> speciesTargets;
        bool finalizeSpecies = false;
        {
          std::lock_guard<std::mutex> lock(speciesState.mutex);
          for (auto &candidate : taskCandidates) {
            retain_target(speciesState.candidates, std::move(candidate),
                          config.native_bounded_targets_per_species);
          }
          if (speciesState.remaining_tasks == 0) {
            throw std::runtime_error(
                "native bounded species task accounting underflow");
          }
          --speciesState.remaining_tasks;
          if (speciesState.remaining_tasks == 0) {
            speciesTargets = std::move(speciesState.candidates);
            finalizeSpecies = true;
          }
        }

        if (finalizeSpecies) {
          speciesTargets = finalize_species_targets(
              std::move(speciesTargets),
              config.native_bounded_sources_per_species,
              config.native_bounded_targets_per_source);
          uint64_t speciesBytes = 0;
          for (const auto &target : speciesTargets) {
            if (target.anchors.size() >
                std::numeric_limits<uint64_t>::max() - speciesBytes) {
              throw std::runtime_error(
                  "native bounded species anchor size overflow");
            }
            speciesBytes += target.anchors.size();
          }
          auto output = shardOutputs.find(resolved.genus);
          if (output == shardOutputs.end()) {
            throw std::runtime_error(
                "native bounded direct build missing shard output");
          }
          const uint64_t speciesOffset =
              output->second->next_offset.fetch_add(
                  speciesBytes, std::memory_order_relaxed);
          if (speciesBytes >
              std::numeric_limits<uint64_t>::max() - speciesOffset) {
            throw std::runtime_error(
                "native bounded shard anchor size overflow");
          }
          uint64_t relativeOffset = 0;
          int fd = -1;
          if (speciesBytes > 0) {
            fd = fdCache.fd_for(resolved.genus);
          }
          speciesState.selected.reserve(speciesTargets.size());
          for (auto &target : speciesTargets) {
            target.plan.anchor_offset = speciesOffset + relativeOffset;
            if (!target.anchors.empty()) {
              write_all_at(fd, target.anchors.data(), target.anchors.size(),
                           target.plan.anchor_offset,
                           output->second->path);
            }
            relativeOffset += target.anchors.size();
            speciesState.selected.push_back(std::move(target.plan));
          }
          if (relativeOffset != speciesBytes) {
            throw std::runtime_error(
                "native bounded species anchor size mismatch");
          }
        }
      } catch (...) {
        std::lock_guard<std::mutex> lock(errorMutex);
        if (workerError == nullptr) {
          workerError = std::current_exception();
          stopWorkers.store(true, std::memory_order_relaxed);
        }
        break;
      }
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(workerCount);
  for (size_t worker = 0; worker < workerCount; ++worker) {
    workers.emplace_back(buildWorker);
  }
  for (auto &worker : workers) {
    worker.join();
  }
  if (workerError != nullptr) {
    std::rethrow_exception(workerError);
  }
  const auto anchorBuildFinished = std::chrono::steady_clock::now();
  stats.anchor_build_seconds =
      std::chrono::duration<double>(anchorBuildFinished - anchorBuildStarted)
          .count();

  std::unordered_map<uint32_t, std::vector<TargetPlan *>> plansByGenus;
  std::vector<chimera::native_bounded::TargetMeta> rootTargets;
  std::vector<TargetPlan *> selectedPlans;
  for (auto &state : speciesStates) {
    if (state->remaining_tasks != 0) {
      throw std::runtime_error(
          "native bounded species task accounting mismatch");
    }
    for (auto &plan : state->selected) {
      selectedPlans.push_back(&plan);
    }
  }
  std::sort(selectedPlans.begin(), selectedPlans.end(),
            [](const TargetPlan *lhs, const TargetPlan *rhs) {
              if (lhs->source_index != rhs->source_index) {
                return lhs->source_index < rhs->source_index;
              }
              return lhs->contig_index < rhs->contig_index;
            });
  rootTargets.reserve(selectedPlans.size());
  for (auto *plan : selectedPlans) {
    rootTargets.push_back({plan->name, plan->species, plan->len, 0,
                           plan->anchor_count, plan->anchor_bytes});
    plansByGenus[plan->genus].push_back(plan);
  }
  if (rootTargets.empty()) {
    throw std::runtime_error("native bounded direct build selected no targets");
  }

  const auto layoutStarted = std::chrono::steady_clock::now();
  chimera::native_bounded::write_index_metadata_only(
      paths.metadata_index, config.native_bounded_k, config.native_bounded_w,
      rootTargets);
  std::ofstream manifestOut(paths.shard_manifest, std::ios::trunc);
  if (!manifestOut) {
    throw std::runtime_error("failed to open native bounded shard manifest");
  }
  manifestOut << "genus\tpath\ttargets\tanchors\n";

  std::vector<uint32_t> shardGenera;
  shardGenera.reserve(plansByGenus.size());
  for (const auto &[genus, _] : plansByGenus) {
    shardGenera.push_back(genus);
  }
  std::sort(shardGenera.begin(), shardGenera.end());

  std::vector<chimera::local_resolution::TargetRep> repMetadataRows;
  for (uint32_t genus : shardGenera) {
    auto &plans = plansByGenus[genus];
    std::vector<chimera::native_bounded::TargetMeta> mergedTargets;
    mergedTargets.reserve(plans.size());
    uint64_t mergedAnchors = 0;
    for (auto *plan : plans) {
      mergedTargets.push_back({plan->name, plan->species, plan->len,
                               plan->anchor_offset, plan->anchor_count,
                               plan->anchor_bytes});
      mergedAnchors += plan->anchor_count;
    }
    auto output = shardOutputs.find(genus);
    if (output == shardOutputs.end()) {
      throw std::runtime_error(
          "native bounded direct build missing shard output");
    }
    if (!std::filesystem::exists(output->second->path)) {
      std::ofstream emptyShard(output->second->path,
                               std::ios::binary | std::ios::trunc);
      if (!emptyShard) {
        throw std::runtime_error(
            "failed to create empty native bounded shard");
      }
    }
    const uint64_t expectedShardBytes =
        output->second->next_offset.load(std::memory_order_relaxed);
    if (std::filesystem::file_size(output->second->path) !=
        expectedShardBytes) {
      throw std::runtime_error(
          "native bounded shard size does not match allocated anchor data");
    }
    manifestOut << genus << '\t' << output->second->path.filename().string()
                << '\t' << mergedTargets.size() << '\t' << mergedAnchors
                << '\n';
    stats.representativeRecords += write_representative_pool_for_genus(
        nullptr, genus, mergedTargets, config.native_bounded_rep_pool_cap, 0,
        config.native_bounded_k, repMetadataRows);
  }
  chimera::local_resolution::write_rep_metadata(paths.rep_metadata,
                                                std::move(repMetadataRows));
  const auto layoutFinished = std::chrono::steady_clock::now();
  stats.layout_seconds =
      std::chrono::duration<double>(layoutFinished - layoutStarted).count();

  for (const auto *plan : selectedPlans) {
    ++stats.targets;
    ++stats.sequences;
    stats.bp += plan->len;
    stats.anchors += plan->anchor_count;
  }
  return stats;
}

} // namespace

NativeBoundedBuildStats build_native_bounded_index(
    const BuildConfig &config,
    const robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    const NativeBoundedOutputPaths &paths) {
  const auto tasks = make_tasks(inputFiles);
  if (tasks.empty()) {
    throw std::runtime_error("native bounded index build has no valid inputs");
  }
  return build_native_bounded_index_fused(config, tasks, paths);
}

} // namespace ChimeraBuild
