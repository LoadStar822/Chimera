#include <build/filter/interleaved-merged-cuckoo-filter.h>

#include <algorithm>
#include <bitset>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <filesystem>

#include <robin_hood.h>

namespace {
template <typename T, typename U>
void expect_equal(const T &lhs, const U &rhs, const std::string &message) {
  if (!(lhs == rhs)) {
    throw std::runtime_error(message + "：实际值=" + std::to_string(lhs) +
                             "，期望值=" + std::to_string(rhs));
  }
}

void expect_true(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

constexpr uint64_t kGoldenIncrement = 0x9E3779B97F4A7C15ULL;

static std::mt19937_64 rng(0xC0FFEE1234ULL);

uint64_t rand64() { return rng(); }

struct TestCase {
  std::string name;
  std::function<void()> body;
};

struct StageResult {
  std::string filter;
  std::string stage;
  double ms;
  double throughputMops;
};

class TestRunner {
public:
  void add(std::string name, std::function<void()> body) {
    tests_.push_back(TestCase{std::move(name), std::move(body)});
  }

  int run() const {
    auto format_duration = [](double ms) {
      std::ostringstream oss;
      oss << std::fixed << std::setprecision(3) << ms;
      return oss.str();
    };

    struct Result {
      std::string name;
      double elapsed_ms;
      bool success;
      std::string message;
    };

    std::vector<Result> results;
    results.reserve(tests_.size());

    int failures = 0;
    double total_ms = 0.0;
    for (const auto &test : tests_) {
      auto start = std::chrono::steady_clock::now();
      bool success = true;
      std::string message;
      try {
        test.body();
      } catch (const std::exception &ex) {
        success = false;
        message = ex.what();
        ++failures;
      }
      auto end = std::chrono::steady_clock::now();
      double elapsed_ms =
          std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
              end - start)
              .count();
      total_ms += elapsed_ms;
      results.push_back(Result{test.name, elapsed_ms, success, message});
    }

    size_t nameWidth = std::string("用例").size();
    for (const auto &result : results) {
      nameWidth = std::max(nameWidth, result.name.size());
    }
    const size_t statusWidth = std::string("结果").size();
    const size_t timeWidth = std::string("耗时 (ms)").size();
    const size_t totalWidth =
        nameWidth + statusWidth + timeWidth + 8; // padding and separators
    std::string divider(totalWidth, '-');

    std::cout << "\n测试结果\n" << divider << '\n';
    std::cout << std::left << std::setw(nameWidth) << "用例"
              << " | " << std::setw(statusWidth) << "结果"
              << " | " << std::right << std::setw(timeWidth) << "耗时 (ms)"
              << '\n';
    std::cout << divider << '\n';

    for (const auto &result : results) {
      std::string timeString = format_duration(result.elapsed_ms);
      std::cout << std::left << std::setw(nameWidth) << result.name << " | "
                << std::setw(statusWidth) << (result.success ? "通过" : "失败")
                << " | " << std::right << std::setw(timeWidth) << timeString
                << '\n';
      std::cout << std::left; // reset alignment for next iteration
    }
    std::cout << divider << '\n';

    if (failures > 0) {
      std::cout << "\n失败用例详情：" << '\n';
      for (const auto &result : results) {
        if (!result.success) {
          std::cout << " - " << result.name << "：" << result.message << '\n';
        }
      }
    }

    std::cout << "\n共 " << results.size() << " 个用例，其中 "
              << (results.size() - failures) << " 个通过，用时 "
              << format_duration(total_ms) << " ms" << std::endl;

    return failures == 0 ? 0 : 1;
  }

private:
  std::vector<TestCase> tests_;
};

chimera::imcf::InterleavedMergedCuckooFilter
make_filter_with_hashes(ChimeraBuild::IMCFConfig &config,
                        const std::vector<uint64_t> &groupHashes) {
  std::vector<chimera::imcf::Group> groups;
  groups.reserve(groupHashes.size());
  for (size_t i = 0; i < groupHashes.size(); ++i) {
    chimera::imcf::Group group;
    group.taxids = {"taxid_" + std::to_string(i)};
    group.totalHash = groupHashes[i];
    groups.push_back(std::move(group));
  }

  if (config.loadFactor <= 0.0) {
    config.loadFactor = 0.95;
  }

  return chimera::imcf::InterleavedMergedCuckooFilter(groups, config);
}

chimera::imcf::InterleavedMergedCuckooFilter
make_filter(ChimeraBuild::IMCFConfig &config, std::size_t groupHashes) {
  return make_filter_with_hashes(config, {groupHashes});
}

std::unordered_set<uint16_t>
collect_fingerprints(chimera::imcf::InterleavedMergedCuckooFilter &filter,
                     const std::vector<size_t> &values) {
  std::unordered_set<uint16_t> fps;
  for (auto value : values) {
    fps.insert(filter.reduceTo12bit(value));
  }
  return fps;
}

size_t next_absent_with_unique_fingerprint(
    chimera::imcf::InterleavedMergedCuckooFilter &filter,
    const std::unordered_set<uint16_t> &used, size_t seed) {
  size_t candidate = seed;
  for (int attempt = 0; attempt < 4096; ++attempt) {
    if (!used.contains(filter.reduceTo12bit(candidate))) {
      return candidate;
    }
    candidate += kGoldenIncrement;
  }
  return candidate;
}

bool bitset_vectors_equal(const std::vector<std::bitset<16>> &lhs,
                          const std::vector<std::bitset<16>> &rhs) {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (size_t i = 0; i < lhs.size(); ++i) {
    if (lhs[i] != rhs[i]) {
      return false;
    }
  }
  return true;
}

std::vector<chimera::imcf::HashChunk> build_hash_chunks(
    const robin_hood::unordered_flat_map<std::string, uint64_t> &counts) {
  if (counts.empty()) {
    return {};
  }

  std::vector<uint64_t> values;
  values.reserve(counts.size());
  for (const auto &kv : counts) {
    values.push_back(kv.second);
  }
  std::sort(values.begin(), values.end());
  uint64_t median = values[values.size() / 2];
  uint64_t threshold = median * 64;

  std::vector<chimera::imcf::HashChunk> chunks;
  chunks.reserve(counts.size());
  for (const auto &kv : counts) {
    if (threshold != 0 && kv.second > threshold) {
      int numChunks = static_cast<int>(std::ceil(
          static_cast<double>(kv.second) / static_cast<double>(threshold)));
      uint64_t chunkSize = kv.second / numChunks;
      for (int i = 0; i < numChunks; ++i) {
        uint64_t current = (i == numChunks - 1)
                               ? (kv.second - chunkSize * (numChunks - 1))
                               : chunkSize;
        chunks.push_back({kv.first, current});
      }
    } else {
      chunks.push_back({kv.first, kv.second});
    }
  }
  return chunks;
}

uint64_t total_hash(const std::vector<chimera::imcf::Group> &groups) {
  uint64_t sum = 0;
  for (const auto &group : groups) {
    sum += group.totalHash;
  }
  return sum;
}

std::vector<size_t> generate_values(size_t binNum, size_t perBin) {
  std::vector<size_t> values;
  values.reserve(binNum * perBin);
  uint64_t seed = 42;
  for (size_t bin = 0; bin < binNum; ++bin) {
    for (size_t i = 0; i < perBin; ++i) {
      seed += kGoldenIncrement;
      values.push_back((static_cast<uint64_t>(bin) << 48) ^ seed);
    }
  }
  return values;
}
} // namespace

static void test_partition_default_small() {
  robin_hood::unordered_flat_map<std::string, uint64_t> counts{
      {"tax1", 120}, {"tax2", 80}, {"tax3", 20}, {"tax4", 10}};

  auto groups = chimera::imcf::partitionHashCount(counts, 2);
  expect_equal(groups.size(), static_cast<std::size_t>(2),
               "默认策略（小样本）应得到 2 个分组");

  uint64_t expectedTotal = 0;
  for (const auto &kv : counts) {
    expectedTotal += kv.second;
  }
  expect_equal(total_hash(groups), expectedTotal,
               "默认策略（小样本）总哈希量应守恒");

  std::unordered_set<std::string> covered;
  for (const auto &group : groups) {
    expect_true(!group.taxids.empty(), "默认策略（小样本）不允许空分组");
    expect_equal(group.taxids.size(), group.assignedHashes.size(),
                 "默认策略（小样本）每个 taxid 应有匹配的配额");
    uint64_t shardSum = 0;
    for (uint64_t assigned : group.assignedHashes) {
      shardSum += assigned;
    }
    expect_equal(shardSum, group.totalHash,
                 "默认策略（小样本）组内配额总和应匹配 totalHash");
    covered.insert(group.taxids.begin(), group.taxids.end());
  }
  expect_equal(covered.size(), counts.size(),
               "默认策略（小样本）应覆盖全部 taxid");
}

static void test_partition_default_heavy_tail() {
  robin_hood::unordered_flat_map<std::string, uint64_t> counts;
  counts["tax_big"] = 1'000'000ULL;
  for (int i = 0; i < 50; ++i) {
    counts["tax_mid_" + std::to_string(i)] = 1'000ULL;
  }
  for (int i = 0; i < 200; ++i) {
    counts["tax_small_" + std::to_string(i)] = 100ULL;
  }
  for (int i = 0; i < 30; ++i) {
    counts["tax_tail_" + std::to_string(i)] =
        20ULL + static_cast<uint64_t>(i % 5);
  }

  auto groups = chimera::imcf::partitionHashCount(counts, 4);
  expect_true(!groups.empty(), "重尾数据分组结果不应为空");

  uint64_t expectedTotal = 0;
  for (const auto &kv : counts) {
    expectedTotal += kv.second;
  }
  expect_equal(total_hash(groups), expectedTotal, "重尾数据总哈希量应守恒");

  std::unordered_set<std::string> covered;
  uint64_t maxTotal = 0;
  uint64_t minTotal = std::numeric_limits<uint64_t>::max();
  for (const auto &group : groups) {
    expect_true(!group.taxids.empty(), "重尾数据每个分组都应包含 taxid");
    expect_equal(group.taxids.size(), group.assignedHashes.size(),
                 "重尾数据每个分组都应保留配额信息");
    uint64_t shardSum = 0;
    for (uint64_t assigned : group.assignedHashes) {
      shardSum += assigned;
    }
    expect_equal(shardSum, group.totalHash,
                 "重尾数据分组的配额总和应守恒");
    covered.insert(group.taxids.begin(), group.taxids.end());
    maxTotal = std::max(maxTotal, group.totalHash);
    minTotal = std::min(minTotal, group.totalHash);
  }
  expect_equal(covered.size(), counts.size(), "重尾数据应覆盖全部 taxid");

  auto chunks = build_hash_chunks(counts);
  robin_hood::unordered_flat_map<std::string, size_t> chunksPerTaxid;
  for (const auto &chunk : chunks) {
    ++chunksPerTaxid[chunk.taxid];
  }

  constexpr size_t kGroupCap = 4;
  size_t naiveLowerBound = (chunks.size() + kGroupCap - 1) / kGroupCap;
  size_t taxidLowerBound = 0;
  for (const auto &kv : chunksPerTaxid) {
    taxidLowerBound = std::max(taxidLowerBound, kv.second);
  }

  size_t expectedGroupCount = std::max(naiveLowerBound, taxidLowerBound);
  expect_equal(groups.size(), expectedGroupCount,
               "重尾数据分组数应与切块策略一致");

  if (minTotal == std::numeric_limits<uint64_t>::max()) {
    minTotal = 0;
  }
  double mean =
      static_cast<double>(expectedTotal) / static_cast<double>(groups.size());
  if (mean > 0.0) {
    double spread = static_cast<double>(maxTotal - minTotal) / mean;
    double tolerance = 0.6;
    expect_true(spread <= tolerance,
                "重尾数据的分组差异应控制在 60% 以内（spread=" +
                    std::to_string(spread) + "）");
  }
}

static void test_partition_edge_cases() {
  robin_hood::unordered_flat_map<std::string, uint64_t> single{
      {"tax_only", 42}};
  auto singleGroups = chimera::imcf::partitionHashCount(single, 8);
  expect_equal(singleGroups.size(), static_cast<std::size_t>(1),
               "单个 taxid 应只生成一个分组");
  expect_equal(total_hash(singleGroups), static_cast<uint64_t>(42),
               "单个 taxid 的哈希总量应一致");
  expect_equal(singleGroups[0].taxids.size(), static_cast<std::size_t>(1),
               "单个 taxid 应完整保留");

  robin_hood::unordered_flat_map<std::string, uint64_t> equal;
  for (int i = 0; i < 10; ++i) {
    equal["tax_eq_" + std::to_string(i)] = 500;
  }
  auto equalGroups = chimera::imcf::partitionHashCount(equal, 4);
  expect_equal(total_hash(equalGroups), static_cast<uint64_t>(5000),
               "哈希值相等场景总量应守恒");
  std::unordered_set<std::string> covered;
  for (const auto &group : equalGroups) {
    covered.insert(group.taxids.begin(), group.taxids.end());
  }
  expect_equal(covered.size(), equal.size(), "哈希值相等场景应覆盖全部 taxid");

  auto wideCapGroups = chimera::imcf::partitionHashCount(
      equal, static_cast<int>(equal.size()) + 5);
  expect_equal(wideCapGroups.size(), static_cast<std::size_t>(1),
               "上限足够大时应只生成一个分组");
  expect_equal(total_hash(wideCapGroups), static_cast<uint64_t>(5000),
               "上限足够大时总哈希量应守恒");
}

static void test_insert_and_bulk_contain_basic() {
  ChimeraBuild::IMCFConfig config{};
  auto filter = make_filter(config, 64);
  std::vector<std::bitset<16>> result(config.binNum);

  bool inserted = filter.insertTag(0, 424242ULL, 3);
  expect_true(inserted, "insertTag 应成功插入新条目");

  filter.bulkContain(424242ULL, result);
  expect_true(result[0].test(3), "bulkContain 应命中刚插入的标签");

  std::vector<size_t> insertedValues{424242ULL};
  auto fps = collect_fingerprints(filter, insertedValues);
  size_t absent = next_absent_with_unique_fingerprint(filter, fps, 111111ULL);

  for (auto &bitset : result) {
    bitset.reset();
  }
  filter.bulkContain(absent, result);
  for (const auto &bitset : result) {
    expect_true(!bitset.any(), "bulkContain 不应为缺失标签置位");
  }
}

static void test_bulk_count_repeat() {
  ChimeraBuild::IMCFConfig config{};
  auto filter = make_filter(config, 96);

  filter.insertTag(0, 100ULL, 1);
  filter.insertTag(0, 200ULL, 2);

  std::vector<std::vector<std::size_t>> counters(config.binNum);
  filter.bulkCount(std::vector<std::size_t>{100ULL, 200ULL, 100ULL}, counters);

  expect_true(counters.size() == config.binNum,
              "bulkCount 返回的桶维度应与配置一致");
  expect_true(counters[0].size() > 2, "bulkCount 结果需保留足够的标签槽位");
  expect_equal(counters[0][1], static_cast<std::size_t>(2),
               "bulkCount 应将标签 1 计数为 2");
  expect_equal(counters[0][2], static_cast<std::size_t>(1),
               "bulkCount 应将标签 2 计数为 1");
}

static void test_router_route_and_subset() {
  ChimeraBuild::IMCFConfig config{};
  auto filter = make_filter_with_hashes(config, {128, 128});

  expect_true(config.binNum >= 2, "测试路由需至少两个分组");

  const size_t value = 987654321ULL;
  expect_true(filter.insertTag(0, value, 1), "应能向分组 0 插入标签");
  expect_true(filter.insertTag(1, value, 2), "应能向分组 1 插入标签");

  filter.buildRouterIndex();
  expect_true(filter.hasRouterIndex(), "构建路由索引后应标记为可用");

  std::vector<uint32_t> routed;
  filter.route(value, routed);
  expect_true(routed.size() == 2, "路由应返回两个命中分组");
  expect_true(routed[0] == 0 && routed[1] == 1, "路由结果需升序列出分组");

#ifdef IMCF_MIRROR64
  filter.releaseBitStorage();
  routed.clear();
  filter.route(value, routed);
  expect_true(routed.size() == 2 && routed[0] == 0 && routed[1] == 1,
              "释放 bit_vector 后路由仍应依赖镜像正常工作");
#endif

  auto tmpPath = std::filesystem::temp_directory_path() / "chimera_imcf_router_test.rtr";
  expect_true(filter.saveRouterIndex(tmpPath.string()), "路由索引应成功写入磁盘");

  filter.clearRouterIndex();
  expect_true(!filter.hasRouterIndex(), "clearRouterIndex 应清空索引状态");
  expect_true(filter.loadRouterIndex(tmpPath.string()), "路由索引应可重新载入");

  std::error_code ec;
  std::filesystem::remove(tmpPath, ec);

  routed.clear();
  filter.route(value, routed);
  expect_true(routed.size() == 2 && routed[0] == 0 && routed[1] == 1,
              "重新加载后路由结果应保持一致");

  std::vector<size_t> minimizers{value};
  std::vector<uint32_t> subset{0, 1, 5};
  std::sort(subset.begin(), subset.end());

  std::vector<std::vector<uint32_t>> counters;
  std::vector<std::pair<uint32_t, uint16_t>> touched;
  filter.bulkCount_sparse_subset(minimizers, subset, counters, &touched);

  bool found0 = false;
  bool found1 = false;
  for (auto [bin, sp] : touched) {
    if (bin == 0 && sp == 1) {
      found0 = counters[bin][sp] == 1;
    }
    if (bin == 1 && sp == 2) {
      found1 = counters[bin][sp] == 1;
    }
  }

  expect_true(found0, "路由计数应命中分组 0 的标签 1");
  expect_true(found1, "路由计数应命中分组 1 的标签 2");
}

static void test_serialize_roundtrip() {
  ChimeraBuild::IMCFConfig config{};
  auto filter = make_filter_with_hashes(config, {128, 96});

  std::vector<size_t> values;
  values.reserve(40);
  for (size_t i = 0; i < 40; ++i) {
    size_t value = 1'000'000ULL * (i + 3) + 7;
    size_t bin = i % config.binNum;
    size_t index = i % 16;
    expect_true(filter.insertTag(bin, value, index),
                "序列化往返前的插入应成功");
    values.push_back(value);
  }

  auto fps = collect_fingerprints(filter, values);
  std::vector<size_t> queries = values;
  size_t seed = 777ULL;
  for (int i = 0; i < 5; ++i) {
    seed = next_absent_with_unique_fingerprint(filter, fps,
                                               seed + kGoldenIncrement);
    queries.push_back(seed);
  }

  std::vector<std::vector<std::bitset<16>>> baseline;
  for (auto query : queries) {
    std::vector<std::bitset<16>> result(config.binNum);
    filter.bulkContain(query, result);
    baseline.push_back(result);
  }

  std::stringstream buffer(std::ios::in | std::ios::out | std::ios::binary);
  {
    cereal::BinaryOutputArchive archive(buffer);
    archive(filter);
  }

  chimera::imcf::InterleavedMergedCuckooFilter restored;
  {
    cereal::BinaryInputArchive archive(buffer);
    archive(restored);
  }

  for (size_t i = 0; i < queries.size(); ++i) {
    std::vector<std::bitset<16>> result(config.binNum);
    restored.bulkContain(queries[i], result);
    expect_true(bitset_vectors_equal(result, baseline[i]),
                "序列化往返后的 bulkContain 结果应与基线一致");
  }

  size_t newValue =
      next_absent_with_unique_fingerprint(filter, fps, seed + kGoldenIncrement);
  expect_true(filter.insertTag(0, newValue, 5), "往返后原始过滤器仍应接受新值");
  expect_true(restored.insertTag(0, newValue, 5),
              "反序列化后的过滤器应接受新值插入");

  std::vector<std::bitset<16>> resultOriginal(config.binNum);
  std::vector<std::bitset<16>> resultRestored(config.binNum);
  filter.bulkContain(newValue, resultOriginal);
  restored.bulkContain(newValue, resultRestored);
  expect_true(bitset_vectors_equal(resultOriginal, resultRestored),
              "新值插入后原始与恢复过滤器的结果应保持一致");
}

static void test_insert_alt_bucket_and_find() {
  ChimeraBuild::IMCFConfig config{};
  config.loadFactor = 0.95;
  auto filter = make_filter(config, 8);
  std::vector<size_t> collidingValues;
  size_t baseHash = filter.hashIndex(1ULL);
  for (size_t candidate = 2;
       collidingValues.size() < 6 && candidate < 2'000'000; ++candidate) {
    if (filter.hashIndex(candidate) != baseHash) {
      continue;
    }
    uint16_t fp = filter.reduceTo12bit(candidate);
    bool duplicate = false;
    for (auto existing : collidingValues) {
      if (filter.reduceTo12bit(existing) == fp) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) {
      collidingValues.push_back(candidate);
    }
  }

  expect_equal(collidingValues.size(), static_cast<std::size_t>(6),
               "应找到足够数量的碰撞值");

  for (size_t i = 0; i < collidingValues.size(); ++i) {
    expect_true(filter.insertTag(0, collidingValues[i], i % 16),
                "碰撞值应可成功插入");
  }

  for (size_t i = 0; i < collidingValues.size(); ++i) {
    std::vector<std::bitset<16>> result(config.binNum);
    filter.bulkContain(collidingValues[i], result);
    expect_true(result[0].test(i % 16), "bulkContain 应能定位这些碰撞值");
  }
}

static void test_insert_failure_throw() {
  ChimeraBuild::IMCFConfig config{};
  config.loadFactor = 0.95;
  auto filter = make_filter(config, 4);

  std::unordered_set<uint16_t> fps;
  size_t seed = 10ULL;
  size_t capacity = config.binNum * config.binSize * 4;
  for (size_t i = 0; i < capacity && seed < 10'000ULL; ++i) {
    size_t bin = config.binNum > 0 ? (i % config.binNum) : 0;
    size_t index = i % 16;
    expect_true(filter.insertTag(bin, seed, index),
                "预填充阶段的插入操作应成功");
    fps.insert(filter.reduceTo12bit(seed));
    seed += 1;
  }

  size_t failingValue = next_absent_with_unique_fingerprint(filter, fps, seed);
  bool ok = true;
  bool caught = false;
  try {
    ok = filter.insertTag(0, failingValue, 5);
  } catch (const std::runtime_error &ex) {
    caught = std::string(ex.what()).find("Filter is full") != std::string::npos;
  }
  expect_true(caught || !ok, "超出容量的插入应抛出异常或返回失败");
}

static void test_bulk_contain_negative() {
  ChimeraBuild::IMCFConfig config{};
  auto filter = make_filter_with_hashes(config, {64, 64});

  std::vector<size_t> inserted;
  inserted.reserve(64);
  for (size_t i = 0; i < 64; ++i) {
    size_t value = 500'000ULL + i * 97ULL;
    size_t bin = i % config.binNum;
    size_t index = i % 16;
    expect_true(filter.insertTag(bin, value, index),
                "负例基准数据的插入应成功");
    inserted.push_back(value);
  }

  auto fps = collect_fingerprints(filter, inserted);
  size_t seed = 123456789ULL;
  for (int attempt = 0; attempt < 32; ++attempt) {
    size_t absent = next_absent_with_unique_fingerprint(
        filter, fps, seed + attempt * kGoldenIncrement);
    std::vector<std::bitset<16>> result(config.binNum);
    filter.bulkContain(absent, result);
    for (const auto &bitset : result) {
      expect_true(!bitset.any(), "查询缺失值时 bulkContain 结果应保持为空");
    }
  }
}

static void benchmark_imcf_large_scale() {
  const size_t binNum = 128;
  const size_t valuesPerBin = 20000;
  const size_t containSamples = 100000;
  const size_t countSamples = 50000;
  const size_t fprSamples = 100000;

  std::vector<size_t> values = generate_values(binNum, valuesPerBin);
  std::vector<size_t> queryValues(
      values.begin(), values.begin() + std::min(containSamples, values.size()));
  std::vector<size_t> countValues(
      values.begin(), values.begin() + std::min(countSamples, values.size()));

  std::vector<uint64_t> groupHashes(binNum, valuesPerBin);
  ChimeraBuild::IMCFConfig config{};
  config.loadFactor = 0.3;
  auto filter = make_filter_with_hashes(config, groupHashes);

  std::vector<std::pair<size_t, size_t>> entries;
  entries.reserve(values.size());
  for (size_t bin = 0; bin < binNum; ++bin) {
    for (size_t i = 0; i < valuesPerBin; ++i) {
      size_t value = values[bin * valuesPerBin + i];
      entries.emplace_back(value, bin);
    }
  }

  std::vector<size_t> perBinInsert(binNum, 0);
  std::vector<StageResult> stages;

  // Insert benchmark
  {
    auto start = std::chrono::steady_clock::now();
    for (const auto &entry : entries) {
      size_t bin = entry.second;
      size_t idx = perBinInsert[bin]++ % 16;
      filter.insertTag(bin, entry.first, idx);
    }
    auto end = std::chrono::steady_clock::now();
    double ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            end - start)
            .count();
    double throughput =
        entries.empty() ? 0.0 : (entries.size() / (ms / 1000.0) / 1e6);
    stages.push_back({"IMCF", "插入", ms, throughput});
  }

  for (size_t bin = 0; bin < binNum; ++bin) {
    expect_equal(perBinInsert[bin], valuesPerBin,
                 "每个桶都应插入预期数量的条目");
  }

  // Contain benchmark
  {
    std::vector<std::bitset<16>> result(config.binNum);
    auto start = std::chrono::steady_clock::now();
    for (size_t i = 0; i < queryValues.size(); ++i) {
      for (auto &bits : result) {
        bits.reset();
      }
      filter.bulkContain(queryValues[i], result);
    }
    auto end = std::chrono::steady_clock::now();
    double ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            end - start)
            .count();
    double throughput =
        queryValues.empty() ? 0.0 : (queryValues.size() / (ms / 1000.0) / 1e6);
    stages.push_back({"IMCF", "查找", ms, throughput});
  }

  // Count benchmark
  {
    std::vector<std::vector<size_t>> counters(config.binNum);
    auto start = std::chrono::steady_clock::now();
    filter.bulkCount(countValues, counters);
    auto end = std::chrono::steady_clock::now();
    double ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            end - start)
            .count();
    double throughput =
        countValues.empty() ? 0.0 : (countValues.size() / (ms / 1000.0) / 1e6);
    stages.push_back({"IMCF", "计数", ms, throughput});
  }

  // Pretty print benchmark results
  std::cout << "\n[IMCF 基准] 大规模测试" << std::endl;
  std::cout << "桶数量=" << binNum << "，每桶元素数=" << valuesPerBin
            << "，contain 样本数=" << queryValues.size()
            << "，count 样本数=" << countValues.size()
            << "，FPR 采样数=" << fprSamples << std::endl;

  std::cout << std::string(78, '-') << '\n';
  std::cout << std::left << std::setw(12) << "过滤器" << std::setw(12) << "阶段"
            << std::right << std::setw(16) << "耗时(ms)" << ' ' << std::setw(18)
            << "吞吐量(M/s)" << '\n';
  std::cout << std::string(78, '-') << '\n';
  for (const auto &stage : stages) {
    std::cout << std::left << std::setw(12) << stage.filter << std::setw(12)
              << stage.stage << std::right << std::setw(16) << std::fixed
              << std::setprecision(3) << stage.ms << ' ' << std::setw(18)
              << std::fixed << std::setprecision(3) << stage.throughputMops
              << '\n';
  }
  std::cout << std::string(78, '-') << '\n';

  // False positive sampling
  std::unordered_set<uint16_t> usedFingerprints;
  usedFingerprints.reserve(values.size() * 2);
  for (auto value : values) {
    usedFingerprints.insert(filter.reduceTo12bit(value));
  }

  std::vector<std::bitset<16>> fprResult(config.binNum);
  size_t binsWithHit = 0;
  size_t bitHitCount = 0;
  uint64_t seed = 123456789ULL;
  auto fprStart = std::chrono::steady_clock::now();
  for (size_t sample = 0; sample < fprSamples; ++sample) {
    seed = next_absent_with_unique_fingerprint(filter, usedFingerprints,
                                               seed + kGoldenIncrement);
    usedFingerprints.insert(filter.reduceTo12bit(seed));
    for (auto &bits : fprResult) {
      bits.reset();
    }
    filter.bulkContain(seed, fprResult);
    for (const auto &bits : fprResult) {
      if (bits.any()) {
        ++binsWithHit;
        bitHitCount += bits.count();
      }
    }
  }
  auto fprEnd = std::chrono::steady_clock::now();
  double fprMs =
      std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
          fprEnd - fprStart)
          .count();
  double avgBinsPerQuery =
      fprSamples == 0
          ? 0.0
          : static_cast<double>(binsWithHit) / static_cast<double>(fprSamples);
  double perBinRate = fprSamples == 0
                          ? 0.0
                          : static_cast<double>(binsWithHit) /
                                (static_cast<double>(fprSamples) *
                                 static_cast<double>(config.binNum));
  double avgBitsPerQuery =
      fprSamples == 0
          ? 0.0
          : static_cast<double>(bitHitCount) / static_cast<double>(fprSamples);

  std::cout << std::fixed << std::setprecision(3) << "假阳性采样耗时: " << fprMs
            << " ms" << '\n'
            << "  每次查询触发的桶平均数: " << avgBinsPerQuery << " ("
            << std::setprecision(4) << perBinRate * 100.0 << "%/桶)\n"
            << std::setprecision(3)
            << "  每次查询命中的位平均数: " << avgBitsPerQuery << std::endl;
}

int main() {
  TestRunner runner;
  runner.add("分组-默认策略-小样本", test_partition_default_small);
  runner.add("分组-默认策略-重尾", test_partition_default_heavy_tail);
  runner.add("分组-边界情况", test_partition_edge_cases);
  runner.add("插入与 bulkContain-基础", test_insert_and_bulk_contain_basic);
  runner.add("bulkCount-重复计数", test_bulk_count_repeat);
  runner.add("路由索引与子集计数", test_router_route_and_subset);
  runner.add("序列化往返", test_serialize_roundtrip);
  runner.add("插入备用桶并命中", test_insert_alt_bucket_and_find);
  runner.add("超容量插入抛异常", test_insert_failure_throw);
  runner.add("bulkContain-缺失值", test_bulk_contain_negative);
  runner.add("IMCF 基准-大规模", benchmark_imcf_large_scale);
  return runner.run();
}
