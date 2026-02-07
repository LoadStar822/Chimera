#pragma once

inline InterleavedMergedCuckooFilter::QidxLookupState
InterleavedMergedCuckooFilter::make_qidx_lookup_state(uint64_t value) const {
  QidxLookupState state;
  const uint16_t fingerprint = reduceTo12bit(value);
  state.hash1 = hashIndex(value);
  state.hash2 = altHash(state.hash1, fingerprint);
  state.group = static_cast<uint32_t>(fingerprint & (qidx->groupSize - 1));
  state.hi = static_cast<uint16_t>(fingerprint >> qidx->g);
  state.hiMask =
      (qidx->fpHiBits >= 32) ? 0xFFFFFFFFu : ((1u << qidx->fpHiBits) - 1u);
  return state;
}

inline void InterleavedMergedCuckooFilter::qidx_group_range(
    size_t bucket, uint32_t group, uint64_t &start, uint64_t &end) const {
  const uint64_t base = qidx->bucketBase[bucket];
  start = base + qidx_prefix(bucket, group);
  end = base + qidx_prefix(bucket, group + 1);
}

inline void InterleavedMergedCuckooFilter::route_qidx(
    uint64_t value, std::vector<uint32_t> &bins) const {
  bins.clear();
  if (!has_qidx()) {
    return;
  }
  const QidxLookupState state = make_qidx_lookup_state(value);

  auto append_bucket = [&](size_t bucket) {
    uint64_t s = 0;
    uint64_t e = 0;
    qidx_group_range(bucket, state.group, s, e);
    uint32_t last_bin = std::numeric_limits<uint32_t>::max();
    for (uint64_t i = s; i < e; ++i) {
      uint64_t code = qidx->entries[i];
      uint16_t hi2 = static_cast<uint16_t>((code >> 4) & state.hiMask);
      if (hi2 != state.hi) {
        continue;
      }
      uint32_t bin = static_cast<uint32_t>(code >> (qidx->fpHiBits + 4));
      if (bin != last_bin) {
        bins.push_back(bin);
        last_bin = bin;
      }
    }
  };

  append_bucket(state.hash1);
  if (state.hash2 != state.hash1) {
    append_bucket(state.hash2);
  }

  if (bins.empty()) {
    return;
  }
  std::sort(bins.begin(), bins.end());
  bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
}
inline void InterleavedMergedCuckooFilter::build_query_index(
    bool include_stash, bool verify, bool low_peak_mode,
    bool drop_classic_before_materialize) {
    if (binNum == 0 || hashSize == 0) {
      qidx = std::make_unique<QueryIndex>();
      qidx->refresh();
      qidx->bucketBase.assign(hashSize + 1, 0);
      qidx->prefix_bits = 1;
      qidx->entry_bits = 1;
      qidx->prefix = sdsl::int_vector<0>(0, 0, qidx->prefix_bits);
      qidx->entries = sdsl::int_vector<0>(0, 0, qidx->entry_bits);
      return;
    }
    if (data.size() == 0) {
      throw std::runtime_error(
          "QIMCF build requires classic in-memory data");
    }

    std::vector<StashFlat> stashFlat;
    std::vector<size_t> stashBucketStart;
    std::vector<size_t> stashBucketEnd;
    if (include_stash) {
      size_t totalStash = 0;
      for (const auto &entries : stash) {
        totalStash += entries.size();
      }
      stashFlat.reserve(totalStash);
      for (size_t bin = 0; bin < binNum; ++bin) {
        if (bin >= stash.size()) {
          break;
        }
        for (const auto &entry : stash[bin]) {
          if (entry.bucket >= hashSize) {
            continue;
          }
          StashFlat flat;
          flat.bucket = static_cast<uint32_t>(entry.bucket);
          flat.fp = entry.fingerprint;
          flat.mask = entry.speciesMask;
          flat.bin = static_cast<uint32_t>(bin);
          stashFlat.push_back(flat);
        }
      }
      std::sort(stashFlat.begin(), stashFlat.end(),
                [](const StashFlat &a, const StashFlat &b) {
                  if (a.bucket != b.bucket) {
                    return a.bucket < b.bucket;
                  }
                  return a.bin < b.bin;
                });
      stashBucketStart.assign(hashSize + 1, 0);
      stashBucketEnd.assign(hashSize + 1, 0);
      size_t idx = 0;
      for (size_t b = 0; b < hashSize; ++b) {
        stashBucketStart[b] = idx;
        while (idx < stashFlat.size() && stashFlat[idx].bucket == b) {
          ++idx;
        }
        stashBucketEnd[b] = idx;
      }
      stashBucketStart[hashSize] = stashFlat.size();
      stashBucketEnd[hashSize] = stashFlat.size();
    } else {
      stashBucketStart.assign(hashSize + 1, 0);
      stashBucketEnd.assign(hashSize + 1, 0);
    }

    qidx = std::make_unique<QueryIndex>();
    constexpr uint8_t kQidxGroupBits = 8;
    qidx->g = kQidxGroupBits;
    qidx->refresh();

    const uint32_t groupMask = qidx->groupSize - 1;
    const uint64_t prefixSize = static_cast<uint64_t>(hashSize) * qidx->stride;
    std::vector<uint32_t> bucketTotal(hashSize, 0);
    uint32_t maxBucketTotal = 0;

#pragma omp parallel
    {
      std::vector<uint32_t> counts(qidx->groupSize, 0);
      std::vector<uint64_t> bucketWords(binNum, 0);
      uint32_t localMaxBucketTotal = 0;

#pragma omp for schedule(dynamic, 4)
      for (int64_t bucket_i = 0; bucket_i < static_cast<int64_t>(hashSize);
           ++bucket_i) {
        size_t bucket = static_cast<size_t>(bucket_i);
        std::fill(counts.begin(), counts.end(), 0);
        loadBucketWords(bucket, bucketWords);
        for (size_t bin = 0; bin < binNum; ++bin) {
          uint64_t q = bucketWords[bin];
          if (q == 0) {
            continue;
          }
          for (int lane = 0; lane < static_cast<int>(tagNum); ++lane) {
            uint16_t tag = static_cast<uint16_t>((q >> (lane * 16)) & 0xFFFFu);
            if (tag == 0u) {
              continue;
            }
            uint16_t fp = static_cast<uint16_t>(tag & 0x0FFFu);
            uint32_t lo = static_cast<uint32_t>(fp & groupMask);
            ++counts[lo];
          }
        }
        if (include_stash) {
          size_t s = stashBucketStart[bucket];
          size_t e = stashBucketEnd[bucket];
          for (size_t idx = s; idx < e; ++idx) {
            uint16_t fp = stashFlat[idx].fp;
            uint32_t lo = static_cast<uint32_t>(fp & groupMask);
            counts[lo] += popcount16(stashFlat[idx].mask);
          }
        }

        uint32_t running = 0;
        for (uint32_t lo = 0; lo < qidx->groupSize; ++lo) {
          running += counts[lo];
        }
        bucketTotal[bucket] = running;
        if (running > localMaxBucketTotal) {
          localMaxBucketTotal = running;
        }
      }

#pragma omp critical
      {
        if (localMaxBucketTotal > maxBucketTotal) {
          maxBucketTotal = localMaxBucketTotal;
        }
      }
    }

    qidx->bucketBase.resize(hashSize + 1);
    uint64_t acc = 0;
    for (size_t b = 0; b < hashSize; ++b) {
      qidx->bucketBase[b] = acc;
      acc += bucketTotal[b];
    }
    qidx->bucketBase[hashSize] = acc;

    qidx->prefix_bits = bits_required_u64(maxBucketTotal);
    const bool use_prefix_spool = low_peak_mode;
    if (!use_prefix_spool) {
      qidx->prefix = sdsl::int_vector<0>(prefixSize, 0, qidx->prefix_bits);
    } else {
      qidx->prefix =
          sdsl::int_vector<0>(0, 0, std::max<uint8_t>(1, qidx->prefix_bits));
    }

    uint8_t bin_bits = bits_required_u64(binNum > 0 ? (binNum - 1) : 0);
    qidx->entry_bits = static_cast<uint8_t>(bin_bits + qidx->fpHiBits + 4);

    std::vector<uint32_t> counts(qidx->groupSize, 0);
    std::vector<uint32_t> pos(qidx->groupSize, 0);
    std::vector<uint32_t> localPrefix(qidx->stride, 0);
    std::vector<uint32_t> bucketCodes;
    bucketCodes.reserve(1u << 16);
    std::vector<uint64_t> payloadBucketWords(binNum, 0);
    std::filesystem::path entriesSpoolPath;

    auto build_bucket_payload = [&](size_t bucket) {
      std::fill(counts.begin(), counts.end(), 0u);
      loadBucketWords(bucket, payloadBucketWords);

      size_t stashStart = include_stash ? stashBucketStart[bucket] : 0;
      size_t stashEnd = include_stash ? stashBucketEnd[bucket] : 0;

      // Pass A: count per fp_lo group and materialize prefix row.
      for (size_t bin = 0; bin < binNum; ++bin) {
        uint64_t q = payloadBucketWords[bin];
        if (q != 0) {
          for (int lane = 0; lane < static_cast<int>(tagNum); ++lane) {
            uint16_t tag = static_cast<uint16_t>((q >> (lane * 16)) & 0xFFFFu);
            if (tag == 0u) {
              continue;
            }
            uint16_t fp = static_cast<uint16_t>(tag & 0x0FFFu);
            uint32_t lo = static_cast<uint32_t>(fp & groupMask);
            ++counts[lo];
          }
        }
      }
      if (include_stash) {
        for (size_t idx = stashStart; idx < stashEnd; ++idx) {
          uint16_t fp = stashFlat[idx].fp;
          uint32_t lo = static_cast<uint32_t>(fp & groupMask);
          counts[lo] += popcount16(stashFlat[idx].mask);
        }
      }

      uint32_t running = 0;
      const uint64_t prefixOffset = static_cast<uint64_t>(bucket) * qidx->stride;
      localPrefix[0] = 0u;
      if (!use_prefix_spool) {
        qidx->prefix[prefixOffset] = 0u;
      }
      for (uint32_t lo = 0; lo < qidx->groupSize; ++lo) {
        running += counts[lo];
        localPrefix[lo + 1] = running;
        if (!use_prefix_spool) {
          qidx->prefix[prefixOffset + lo + 1] = running;
        }
      }
      if (running != bucketTotal[bucket]) {
        throw std::runtime_error(
            "QIMCF build: bucket total mismatch during prefix materialization");
      }

      bucketCodes.assign(running, 0u);
      for (uint32_t lo = 0; lo < qidx->groupSize; ++lo) {
        pos[lo] = localPrefix[lo];
      }

      // Pass B: emit per-group payload using the same scan order as classic
      // table.
      size_t stashIdx = stashStart;
      for (size_t bin = 0; bin < binNum; ++bin) {
        uint64_t q = payloadBucketWords[bin];
        if (q != 0) {
          for (int lane = 0; lane < static_cast<int>(tagNum); ++lane) {
            uint16_t tag = static_cast<uint16_t>((q >> (lane * 16)) & 0xFFFFu);
            if (tag == 0u) {
              continue;
            }
            uint16_t fp = static_cast<uint16_t>(tag & 0x0FFFu);
            uint32_t lo = static_cast<uint32_t>(fp & groupMask);
            uint16_t hi = static_cast<uint16_t>(fp >> qidx->g);
            uint16_t sp = static_cast<uint16_t>(tag >> 12);
            uint32_t code =
                (static_cast<uint32_t>(bin) << (qidx->fpHiBits + 4)) |
                (static_cast<uint32_t>(hi) << 4) | sp;
            uint32_t local = pos[lo]++;
            bucketCodes[local] = code;
          }
        }
        if (include_stash) {
          while (stashIdx < stashEnd && stashFlat[stashIdx].bin == bin) {
            uint16_t fp = stashFlat[stashIdx].fp;
            uint32_t lo = static_cast<uint32_t>(fp & groupMask);
            uint16_t hi = static_cast<uint16_t>(fp >> qidx->g);
            uint16_t mask = stashFlat[stashIdx].mask;
            while (mask) {
              uint16_t sp = static_cast<uint16_t>(std::countr_zero(mask));
              mask &= static_cast<uint16_t>(mask - 1);
              uint32_t code =
                  (static_cast<uint32_t>(bin) << (qidx->fpHiBits + 4)) |
                  (static_cast<uint32_t>(hi) << 4) | sp;
              uint32_t local = pos[lo]++;
              bucketCodes[local] = code;
            }
            ++stashIdx;
          }
        }
      }

      for (uint32_t lo = 0; lo < qidx->groupSize; ++lo) {
        if (pos[lo] != localPrefix[lo + 1]) {
          throw std::runtime_error(
              "QIMCF build: bucket total mismatch during payload build");
        }
      }
    };

    if (!low_peak_mode) {
      qidx->entries = sdsl::int_vector<0>(acc, 0, qidx->entry_bits);
      for (size_t bucket = 0; bucket < hashSize; ++bucket) {
        build_bucket_payload(bucket);
        const uint64_t base = qidx->bucketBase[bucket];
        for (uint64_t i = 0; i < bucketCodes.size(); ++i) {
          qidx->entries[base + i] = bucketCodes[i];
        }
      }
    } else {
      const auto nonce = static_cast<uint64_t>(
          std::chrono::steady_clock::now().time_since_epoch().count());
      entriesSpoolPath =
          std::filesystem::path("/tmp") /
          ("chimera_qidx_entries_" + std::to_string(nonce) + ".bin");
      std::filesystem::path prefixSpoolPath =
          std::filesystem::path("/tmp") /
          ("chimera_qidx_prefix_" + std::to_string(nonce) + ".bin");

      auto cleanup_spool = [&]() {
        if (!entriesSpoolPath.empty()) {
          std::error_code ec;
          std::filesystem::remove(entriesSpoolPath, ec);
        }
        if (!prefixSpoolPath.empty()) {
          std::error_code ec;
          std::filesystem::remove(prefixSpoolPath, ec);
        }
      };

      try {
        {
          std::ofstream entries_spool(entriesSpoolPath, std::ios::binary);
          if (!entries_spool.is_open()) {
            throw std::runtime_error("QIMCF low-peak: failed to open spool file");
          }
          std::ofstream prefix_spool(prefixSpoolPath, std::ios::binary);
          if (!prefix_spool.is_open()) {
            throw std::runtime_error(
                "QIMCF low-peak: failed to open prefix spool file");
          }

          for (size_t bucket = 0; bucket < hashSize; ++bucket) {
            build_bucket_payload(bucket);
            prefix_spool.write(reinterpret_cast<const char *>(localPrefix.data()),
                              static_cast<std::streamsize>(qidx->stride *
                                                           sizeof(uint32_t)));
            if (!prefix_spool.good()) {
              throw std::runtime_error(
                  "QIMCF low-peak: failed to write prefix spool");
            }
            if (!bucketCodes.empty()) {
              entries_spool.write(
                  reinterpret_cast<const char *>(bucketCodes.data()),
                  static_cast<std::streamsize>(bucketCodes.size() *
                                               sizeof(uint32_t)));
              if (!entries_spool.good()) {
                throw std::runtime_error(
                    "QIMCF low-peak: failed to write entries spool");
              }
            }
          }
          entries_spool.flush();
          if (!entries_spool.good()) {
            throw std::runtime_error(
                "QIMCF low-peak: failed to flush entries spool");
          }
          entries_spool.close();
          if (!entries_spool.good()) {
            throw std::runtime_error(
                "QIMCF low-peak: failed to close entries spool");
          }
          prefix_spool.flush();
          if (!prefix_spool.good()) {
            throw std::runtime_error(
                "QIMCF low-peak: failed to flush prefix spool");
          }
          prefix_spool.close();
          if (!prefix_spool.good()) {
            throw std::runtime_error(
                "QIMCF low-peak: failed to close prefix spool");
          }
        }
        const uint64_t expectedBytes = acc * sizeof(uint32_t);
        const uint64_t entries_spool_bytes =
            std::filesystem::file_size(entriesSpoolPath);
        if (entries_spool_bytes != expectedBytes) {
          throw std::runtime_error(
              "QIMCF low-peak: entries spool size mismatch (expected " +
              std::to_string(expectedBytes) + ", got " +
              std::to_string(entries_spool_bytes) + ")");
        }
        const uint64_t expectedPrefixBytes = prefixSize * sizeof(uint32_t);
        const uint64_t prefix_spool_bytes =
            std::filesystem::file_size(prefixSpoolPath);
        if (prefix_spool_bytes != expectedPrefixBytes) {
          throw std::runtime_error(
              "QIMCF low-peak: prefix spool size mismatch (expected " +
              std::to_string(expectedPrefixBytes) + ", got " +
              std::to_string(prefix_spool_bytes) + ")");
        }

        if (drop_classic_before_materialize && !verify) {
          release_classic_storage();
        }
        std::vector<StashFlat>().swap(stashFlat);
        std::vector<size_t>().swap(stashBucketStart);
        std::vector<size_t>().swap(stashBucketEnd);
        std::vector<uint32_t>().swap(bucketTotal);
        std::vector<uint32_t>().swap(counts);
        std::vector<uint32_t>().swap(pos);
        std::vector<uint32_t>().swap(localPrefix);
        std::vector<uint32_t>().swap(bucketCodes);
        std::vector<uint64_t>().swap(payloadBucketWords);

        constexpr size_t kSpoolChunkU32 = 1u << 20;
        auto load_u32_spool =
            [&](const std::filesystem::path &path, uint64_t totalCount,
                const char *readError, const char *truncateError,
                auto &&sink) {
              if (totalCount == 0) {
                return;
              }
              std::ifstream spool(path, std::ios::binary);
              if (!spool.is_open()) {
                throw std::runtime_error(readError);
              }
              std::vector<uint32_t> buf(kSpoolChunkU32, 0u);
              uint64_t offset = 0;
              while (offset < totalCount) {
                uint64_t need =
                    std::min<uint64_t>(kSpoolChunkU32, totalCount - offset);
                auto needBytes =
                    static_cast<std::streamsize>(need * sizeof(uint32_t));
                spool.read(reinterpret_cast<char *>(buf.data()), needBytes);
                if (spool.gcount() != needBytes) {
                  throw std::runtime_error(truncateError);
                }
                for (uint64_t i = 0; i < need; ++i) {
                  sink(offset + i, buf[i]);
                }
                offset += need;
              }
            };

        qidx->entries = sdsl::int_vector<0>(acc, 0, qidx->entry_bits);
        load_u32_spool(
            entriesSpoolPath, acc, "QIMCF low-peak: failed to read spool file",
            "QIMCF low-peak: truncated entries spool",
            [&](uint64_t idx, uint32_t value) { qidx->entries[idx] = value; });

        qidx->prefix = sdsl::int_vector<0>(prefixSize, 0, qidx->prefix_bits);
        load_u32_spool(
            prefixSpoolPath, prefixSize,
            "QIMCF low-peak: failed to read prefix spool",
            "QIMCF low-peak: truncated prefix spool",
            [&](uint64_t idx, uint32_t value) { qidx->prefix[idx] = value; });

        cleanup_spool();
      } catch (...) {
        cleanup_spool();
        throw;
      }
    }

    if (!verify) {
      return;
    }

    std::mt19937_64 rng{0xC0FFEEULL};
    std::vector<uint64_t> out_scan;
    std::vector<uint64_t> out_qidx;
    out_scan.reserve(256);
    out_qidx.reserve(256);

    auto collect_full = [&](uint64_t value,
                            std::vector<uint64_t> &out,
                            auto &&fn) {
      out.clear();
      fn(value, [&](uint32_t bin, uint16_t sp) {
        out.push_back((static_cast<uint64_t>(bin) << 16) |
                      static_cast<uint64_t>(sp));
      });
      std::sort(out.begin(), out.end());
    };

    for (size_t i = 0; i < 1000; ++i) {
      uint64_t value = rng();
      collect_full(value, out_scan,
                   [&](uint64_t v, auto &&emit) {
                     bulkContain_events_scan(v, std::forward<decltype(emit)>(emit));
                   });
      collect_full(value, out_qidx,
                   [&](uint64_t v, auto &&emit) {
                     bulkContain_events_qidx(v, std::forward<decltype(emit)>(emit));
                   });
      if (out_scan != out_qidx) {
        throw std::runtime_error("QIMCF verify failed: bulkContain_events mismatch");
      }
    }

    size_t subsetSize = std::min<size_t>(binNum, 256);
    std::vector<uint32_t> subset;
    if (subsetSize > 0) {
      std::vector<uint32_t> bins;
      bins.resize(binNum);
      std::iota(bins.begin(), bins.end(), 0u);
      std::shuffle(bins.begin(), bins.end(), rng);
      subset.assign(bins.begin(), bins.begin() + subsetSize);
      std::sort(subset.begin(), subset.end());
    }

    if (!subset.empty()) {
      for (size_t i = 0; i < 200; ++i) {
        uint64_t value = rng();
        collect_full(value, out_scan,
                     [&](uint64_t v, auto &&emit) {
                       bulkContain_events_subset_scan(
                           v, subset, std::forward<decltype(emit)>(emit));
                     });
        collect_full(value, out_qidx,
                     [&](uint64_t v, auto &&emit) {
                       bulkContain_events_subset_qidx(
                           v, subset, std::forward<decltype(emit)>(emit));
                     });
        if (out_scan != out_qidx) {
          throw std::runtime_error("QIMCF verify failed: bulkContain_events_subset mismatch");
        }
      }
    }

    for (size_t i = 0; i < 200; ++i) {
      uint64_t value = rng();
      std::vector<uint32_t> route_scan_bins;
      std::vector<uint32_t> route_qidx_bins;
      route_scan(value, route_scan_bins);
      route_qidx(value, route_qidx_bins);
      if (route_scan_bins != route_qidx_bins) {
        throw std::runtime_error("QIMCF verify failed: route mismatch");
      }
    }
  }

template <class EmitFn>
inline void InterleavedMergedCuckooFilter::bulkContain_events_qidx(
    uint64_t value, EmitFn &&emit) {
  if (!has_qidx()) {
    return;
  }
  const QidxLookupState state = make_qidx_lookup_state(value);

  auto process_bucket = [&](size_t bucket) {
    uint64_t s = 0;
    uint64_t e = 0;
    qidx_group_range(bucket, state.group, s, e);
    for (uint64_t i = s; i < e; ++i) {
      uint64_t code = qidx->entries[i];
      uint16_t hi2 = static_cast<uint16_t>((code >> 4) & state.hiMask);
      if (hi2 != state.hi) {
        continue;
      }
      uint32_t bin = static_cast<uint32_t>(code >> (qidx->fpHiBits + 4));
      uint16_t sp = static_cast<uint16_t>(code & 0xF);
      emit(bin, sp);
    }
  };

  process_bucket(state.hash1);
  if (state.hash2 != state.hash1) {
    process_bucket(state.hash2);
  }
}

template <class EmitFn>
inline void InterleavedMergedCuckooFilter::bulkContain_events_subset_qidx(
    uint64_t value, const std::vector<uint32_t> &binSubset, EmitFn &&emit) {
  if (!has_qidx() || binSubset.empty()) {
    return;
  }
  const QidxLookupState state = make_qidx_lookup_state(value);

  auto process_bucket = [&](size_t bucket) {
    uint64_t s = 0;
    uint64_t e = 0;
    qidx_group_range(bucket, state.group, s, e);
    size_t j = 0;
    for (uint64_t i = s; i < e && j < binSubset.size(); ++i) {
      uint64_t code = qidx->entries[i];
      uint32_t bin = static_cast<uint32_t>(code >> (qidx->fpHiBits + 4));
      while (j < binSubset.size() && binSubset[j] < bin) {
        ++j;
      }
      if (j >= binSubset.size()) {
        break;
      }
      if (binSubset[j] == bin) {
        uint16_t hi2 = static_cast<uint16_t>((code >> 4) & state.hiMask);
        if (hi2 != state.hi) {
          continue;
        }
        uint16_t sp = static_cast<uint16_t>(code & 0xF);
        emit(bin, sp);
      }
    }
  };

  process_bucket(state.hash1);
  if (state.hash2 != state.hash1) {
    process_bucket(state.hash2);
  }
}

template <class EmitFn>
inline void InterleavedMergedCuckooFilter::bulkContain_events_subset_qidx_mask(
    uint64_t value, const std::vector<uint8_t> &binMask, EmitFn &&emit) {
  if (!has_qidx() || binMask.empty()) {
    return;
  }
  const QidxLookupState state = make_qidx_lookup_state(value);
  const uint8_t *mask = binMask.data();
  const uint32_t maskSize = static_cast<uint32_t>(binMask.size());

  auto process_bucket = [&](size_t bucket) {
    uint64_t s = 0;
    uint64_t e = 0;
    qidx_group_range(bucket, state.group, s, e);
    for (uint64_t i = s; i < e; ++i) {
      uint64_t code = qidx->entries[i];
      uint32_t bin = static_cast<uint32_t>(code >> (qidx->fpHiBits + 4));
      if (bin >= maskSize || mask[bin] == 0u) {
        continue;
      }
      uint16_t hi2 = static_cast<uint16_t>((code >> 4) & state.hiMask);
      if (hi2 != state.hi) {
        continue;
      }
      uint16_t sp = static_cast<uint16_t>(code & 0xF);
      emit(bin, sp);
    }
  };

  process_bucket(state.hash1);
  if (state.hash2 != state.hash1) {
    process_bucket(state.hash2);
  }
}

template <class EmitFn>
inline void InterleavedMergedCuckooFilter::bulkContain_events_subset_qidx_marked(
    uint64_t value, const std::vector<uint32_t> &binMarks, uint32_t activeMark,
    EmitFn &&emit) {
  if (!has_qidx() || binMarks.empty() || activeMark == 0u) {
    return;
  }
  const QidxLookupState state = make_qidx_lookup_state(value);
  const uint32_t *marks = binMarks.data();
  const uint32_t marksSize = static_cast<uint32_t>(binMarks.size());

  auto process_bucket = [&](size_t bucket) {
    uint64_t s = 0;
    uint64_t e = 0;
    qidx_group_range(bucket, state.group, s, e);
    for (uint64_t i = s; i < e; ++i) {
      uint64_t code = qidx->entries[i];
      uint32_t bin = static_cast<uint32_t>(code >> (qidx->fpHiBits + 4));
      if (bin >= marksSize || marks[bin] != activeMark) {
        continue;
      }
      uint16_t hi2 = static_cast<uint16_t>((code >> 4) & state.hiMask);
      if (hi2 != state.hi) {
        continue;
      }
      uint16_t sp = static_cast<uint16_t>(code & 0xF);
      emit(bin, sp);
    }
  };

  process_bucket(state.hash1);
  if (state.hash2 != state.hash1) {
    process_bucket(state.hash2);
  }
}
