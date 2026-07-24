#include "ChimeraClassifyCommon.hpp"

#include <utils/LocalResolutionManifest.hpp>

#include <chrono>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cereal/archives/binary.hpp>

namespace ChimeraClassify {

namespace {

template <typename T>
void load_required_archive_section(cereal::BinaryInputArchive &archive,
                                   T &value, const std::string &label) {
  try {
    archive(value);
  } catch (const cereal::Exception &exc) {
    throw std::runtime_error("Failed to load IMCF " + label + ": " +
                             exc.what());
  }
}

void load_optional_coverage_meta(cereal::BinaryInputArchive &archive,
                                 chimera::presence::CoverageMeta *coverageMeta) {
  if (coverageMeta) {
    try {
      archive(*coverageMeta);
    } catch (const cereal::Exception &) {
      coverageMeta->entries.clear();
      coverageMeta->unique_deg_threshold = 1;
    }
    return;
  }
  try {
    chimera::presence::CoverageMeta tmp;
    archive(tmp);
  } catch (const cereal::Exception &) {
  }
}

class MappedArchive {
public:
  explicit MappedArchive(const std::filesystem::path &path) {
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
      throw std::runtime_error("Cannot open IMCF archive: " + path.string());
    }

    struct stat statBuf {};
    if (::fstat(fd, &statBuf) != 0 || statBuf.st_size <= 0) {
      ::close(fd);
      throw std::runtime_error("Cannot stat IMCF archive: " + path.string());
    }
    bytes_ = static_cast<size_t>(statBuf.st_size);
    base_ = ::mmap(nullptr, bytes_, PROT_READ, MAP_SHARED, fd, 0);
    ::close(fd);
    if (base_ == MAP_FAILED) {
      base_ = nullptr;
      bytes_ = 0;
      throw std::runtime_error("Cannot map IMCF archive: " + path.string());
    }
  }

  ~MappedArchive() {
    if (base_ != nullptr) {
      ::munmap(base_, bytes_);
    }
  }

  MappedArchive(const MappedArchive &) = delete;
  MappedArchive &operator=(const MappedArchive &) = delete;

  const uint8_t *data() const {
    return static_cast<const uint8_t *>(base_);
  }
  size_t size() const { return bytes_; }

  void *release() {
    void *base = base_;
    base_ = nullptr;
    bytes_ = 0;
    return base;
  }

private:
  void *base_{nullptr};
  size_t bytes_{0};
};

class ArchiveCursor {
public:
  ArchiveCursor(const uint8_t *data, size_t bytes)
      : data_(data), bytes_(bytes) {}

  template <typename T> T read(const char *label) {
    const uint8_t *source = take(sizeof(T), label);
    T value{};
    std::memcpy(&value, source, sizeof(T));
    return value;
  }

  const uint8_t *take(uint64_t count, const char *label) {
    if (offset_ > bytes_ ||
        count > std::numeric_limits<size_t>::max() ||
        static_cast<size_t>(count) > bytes_ - offset_) {
      throw std::runtime_error(std::string("Truncated IMCF ") + label);
    }
    const uint8_t *result = data_ + offset_;
    offset_ += static_cast<size_t>(count);
    return result;
  }

  size_t offset() const { return offset_; }

private:
  const uint8_t *data_;
  size_t bytes_;
  size_t offset_{0};
};

struct MappedQidxLayout {
  uint8_t storageMode{};
  size_t binNum{};
  size_t binSize{};
  size_t tagNum{};
  int maxCuckooCount{};
  size_t hashSize{};
  const uint8_t *bucketBase{};
  uint64_t bucketCount{};
  const uint8_t *prefixWords{};
  uint64_t prefixBitSize{};
  uint8_t prefixBits{};
  const uint8_t *entryWords{};
  uint64_t entryBitSize{};
  uint8_t entryBits{};
  uint8_t groupBits{};
  size_t metadataOffset{};
};

uint64_t packed_word_bytes(uint64_t bitSize, const char *label) {
  if (bitSize > std::numeric_limits<uint64_t>::max() - 63u) {
    throw std::runtime_error(std::string("Invalid IMCF ") + label +
                             " bit size");
  }
  const uint64_t wordCount = (bitSize + 63u) >> 6;
  if (wordCount > std::numeric_limits<uint64_t>::max() / sizeof(uint64_t)) {
    throw std::runtime_error(std::string("Invalid IMCF ") + label +
                             " bit size");
  }
  return wordCount * sizeof(uint64_t);
}

uint64_t read_unaligned_u64(const uint8_t *data) {
  uint64_t value = 0;
  std::memcpy(&value, data, sizeof(value));
  return value;
}

MappedQidxLayout parse_mapped_qidx(const uint8_t *data, size_t bytes) {
  static_assert(sizeof(size_t) == sizeof(uint64_t),
                "Mapped QIDX archives require 64-bit size_t");

  // Locate the three packed qidx arrays without materializing them. The
  // trailing taxid/configuration sections remain decoded by Cereal.
  ArchiveCursor cursor(data, bytes);
  MappedQidxLayout layout;
  layout.storageMode = cursor.read<uint8_t>("storage mode");
  layout.binNum = cursor.read<size_t>("bin count");
  layout.binSize = cursor.read<size_t>("bin size");
  layout.tagNum = cursor.read<size_t>("tag count");
  layout.maxCuckooCount = cursor.read<int>("cuckoo limit");
  layout.hashSize = cursor.read<size_t>("hash size");

  layout.bucketCount = cursor.read<cereal::size_type>("bucket-base count");
  if (layout.bucketCount >
      std::numeric_limits<uint64_t>::max() / sizeof(uint64_t)) {
    throw std::runtime_error("Invalid IMCF bucket-base count");
  }
  layout.bucketBase =
      cursor.take(layout.bucketCount * sizeof(uint64_t), "bucket-base data");

  layout.prefixBits = cursor.read<uint8_t>("prefix width");
  cursor.read<float>("prefix growth factor");
  layout.prefixBitSize = cursor.read<cereal::size_type>("prefix bit size");
  layout.prefixWords =
      cursor.take(packed_word_bytes(layout.prefixBitSize, "prefix"),
                  "prefix data");

  layout.entryBits = cursor.read<uint8_t>("entry width");
  cursor.read<float>("entry growth factor");
  layout.entryBitSize = cursor.read<cereal::size_type>("entry bit size");
  layout.entryWords =
      cursor.take(packed_word_bytes(layout.entryBitSize, "entry"),
                  "entry data");

  layout.groupBits = cursor.read<uint8_t>("group width");
  const uint8_t serializedPrefixBits =
      cursor.read<uint8_t>("serialized prefix width");
  const uint8_t serializedEntryBits =
      cursor.read<uint8_t>("serialized entry width");
  layout.metadataOffset = cursor.offset();

  if (layout.storageMode != 1 || layout.hashSize == 0 ||
      layout.hashSize == std::numeric_limits<size_t>::max() ||
      layout.binSize != layout.hashSize ||
      layout.bucketCount != layout.hashSize + 1 ||
      layout.prefixBits == 0 || layout.prefixBits > 64 ||
      layout.entryBits == 0 || layout.entryBits > 64 ||
      layout.groupBits > 12 || serializedPrefixBits != layout.prefixBits ||
      serializedEntryBits != layout.entryBits ||
      layout.prefixBitSize % layout.prefixBits != 0 ||
      layout.entryBitSize % layout.entryBits != 0) {
    throw std::runtime_error(
        "Incompatible or malformed qidx-only IMCF archive");
  }

  const uint64_t prefixCount = layout.prefixBitSize / layout.prefixBits;
  const uint64_t stride = (uint64_t{1} << layout.groupBits) + 1u;
  if (layout.hashSize > std::numeric_limits<uint64_t>::max() / stride ||
      prefixCount != static_cast<uint64_t>(layout.hashSize) * stride) {
    throw std::runtime_error("IMCF prefix index dimensions are inconsistent");
  }
  const uint64_t entryCount = layout.entryBitSize / layout.entryBits;
  const uint64_t finalBucketBase =
      read_unaligned_u64(layout.bucketBase +
                         layout.hashSize * sizeof(uint64_t));
  if (finalBucketBase != entryCount) {
    throw std::runtime_error("IMCF entry index dimensions are inconsistent");
  }
  return layout;
}

} // namespace

void rebuild_bin_slot_rep_lookup(
    TaxDict &tax, const std::vector<uint32_t> *tid2speciesRep) {
  tax.binSlotRepTid.assign(tax.idx2id.size() * kTaxSlotCount, kInvalidTidId);
  for (size_t bin = 0; bin < tax.idx2id.size(); ++bin) {
    const auto &speciesVec = tax.idx2id[bin];
    const size_t slotLimit =
        std::min<size_t>(speciesVec.size(), kTaxSlotCount);
    const size_t base = bin * kTaxSlotCount;
    for (size_t slot = 0; slot < slotLimit; ++slot) {
      uint32_t tid = speciesVec[slot];
      if (tid >= tax.id2str.size()) {
        continue;
      }
      if (tid2speciesRep && tid < tid2speciesRep->size()) {
        tid = (*tid2speciesRep)[tid];
      }
      tax.binSlotRepTid[base + slot] = tid;
    }
  }
}

TaxDict build_tax_dict(const std::vector<std::vector<std::string>> &idx2tax) {
  robin_hood::unordered_flat_map<std::string, uint32_t> dict;
  dict.reserve(1ull << 20);

  TaxDict td;
  td.idx2id.resize(idx2tax.size());
  for (size_t b = 0; b < idx2tax.size(); ++b) {
    td.idx2id[b].resize(idx2tax[b].size());
    for (size_t s = 0; s < idx2tax[b].size(); ++s) {
      const std::string &t = idx2tax[b][s];
      auto it = dict.find(t);
      uint32_t id;
      if (it == dict.end()) {
        id = static_cast<uint32_t>(td.id2str.size());
        dict.emplace(t, id);
        td.id2str.push_back(t);
        td.tid2bin.emplace_back();
      } else {
        id = it->second;
      }
      td.idx2id[b][s] = id;
      td.tid2bin[id].push_back(static_cast<uint32_t>(b));
    }
  }
  for (auto &bins : td.tid2bin) {
    std::sort(bins.begin(), bins.end());
    bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
  }
  td.str2id = std::move(dict);
  rebuild_bin_slot_rep_lookup(td);
  return td;
}

void loadFilter(
    const std::string &input_file,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid,
    chimera::presence::CoverageMeta *coverageMeta) {
  namespace fs = std::filesystem;

  fs::path archivePath =
      chimera::local_resolution::core_archive_path_for(input_file);
  if (!fs::exists(archivePath)) {
    throw std::runtime_error("Cannot find IMCF archive: " + input_file);
  }


  MappedArchive mappedArchive(archivePath);
  const MappedQidxLayout layout =
      parse_mapped_qidx(mappedArchive.data(), mappedArchive.size());

  std::ifstream is(archivePath, std::ios::binary);
  if (!is.is_open()) {
    throw std::runtime_error("Cannot open IMCF archive: " + archivePath.string());
  }
  is.seekg(static_cast<std::streamoff>(layout.metadataOffset), std::ios::beg);
  if (!is.good()) {
    throw std::runtime_error("Cannot seek IMCF archive metadata: " +
                             archivePath.string());
  }

  cereal::BinaryInputArchive archive(is);
  load_required_archive_section(archive, indexToTaxid, "taxid index");
  load_required_archive_section(archive, imcfConfig, "configuration");
  load_optional_coverage_meta(archive, coverageMeta);
  is.close();

  const size_t mappedBytes = mappedArchive.size();
  void *mappedBase = mappedArchive.release();
  try {
    imcf.initialize_mapped_qidx(
        mappedBase, mappedBytes, layout.binNum, layout.binSize, layout.tagNum,
        layout.maxCuckooCount, layout.hashSize, layout.bucketBase,
        layout.prefixWords, layout.prefixBits, layout.entryWords,
        layout.entryBits, layout.groupBits);
  } catch (...) {
    ::munmap(mappedBase, mappedBytes);
    throw;
  }

  if (imcfConfig.hashVersion != ChimeraBuild::IMCFConfig::CurrentHashVersion) {
    throw std::runtime_error(
        "Incompatible IMCF database version: hash_version=" +
        std::to_string(imcfConfig.hashVersion) +
        ". Rebuild the database with the current Chimera.");
  }
  if (imcf.get_storage_mode() != 1) {
    throw std::runtime_error(
        "Incompatible IMCF storage mode: only qidx-only (storageMode=1) is supported. "
        "Rebuild the database with the current Chimera.");
  }
  if (imcfConfig.seed64 == 0) {
    throw std::runtime_error(
        "The IMCF database is missing feature seed metadata. Please rerun Chimera build.");
  }
  if (imcfConfig.fpSalt != ChimeraBuild::IMCFConfig::DefaultFingerprintSalt) {
    throw std::runtime_error(
        "IMCF fingerprint salt does not match the current implementation, detected fp_salt=" +
        std::to_string(imcfConfig.fpSalt) +
        ". Rebuild the database or upgrade the program.");
  }
}

} // namespace ChimeraClassify
