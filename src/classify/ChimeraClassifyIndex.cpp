#include "ChimeraClassifyCommon.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <cereal/archives/binary.hpp>

namespace ChimeraClassify {

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
  return td;
}

void loadFilter(
    const std::string &input_file,
    chimera::imcf::InterleavedMergedCuckooFilter &imcf,
    ChimeraBuild::IMCFConfig &imcfConfig,
    std::vector<std::vector<std::string>> &indexToTaxid,
    chimera::presence::CoverageMeta *coverageMeta) {
  namespace fs = std::filesystem;

  fs::path archivePath{input_file};
  if (!fs::exists(archivePath)) {
    archivePath = fs::path{input_file + ".imcf"};
  }
  if (!fs::exists(archivePath)) {
    throw std::runtime_error("Cannot find IMCF archive: " + input_file);
  }


  std::ifstream is(archivePath, std::ios::binary);
  if (!is.is_open()) {
    throw std::runtime_error("Cannot open IMCF archive: " + archivePath.string());
  }

  cereal::BinaryInputArchive archive(is);
  try {
    archive(imcf);
  } catch (const cereal::Exception &exc) {
    throw std::runtime_error(std::string("Failed to load IMCF archive: ") + exc.what());
  }
  try {
    archive(indexToTaxid);
  } catch (const cereal::Exception &exc) {
    throw std::runtime_error(std::string("Failed to load IMCF taxid index: ") + exc.what());
  }
  try {
    archive(imcfConfig);
  } catch (const cereal::Exception &exc) {
    throw std::runtime_error(std::string("Failed to load IMCF configuration: ") + exc.what());
  }
  if (coverageMeta) {
    try {
      archive(*coverageMeta);
    } catch (const cereal::Exception &) {
      coverageMeta->entries.clear();
      coverageMeta->unique_deg_threshold = 1;
    }
  } else {
    try {
      chimera::presence::CoverageMeta tmp;
      archive(tmp);
    } catch (const cereal::Exception &) {
    }
  }
  is.close();

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
