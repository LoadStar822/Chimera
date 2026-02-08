#include "ChimeraBuildCommon.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace ChimeraBuild {

void parseInputFile(
    const std::string &filePath,
    robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    FileInfo &fileInfo) {
  std::ifstream inputFile(filePath);
  if (!inputFile.is_open()) {
    std::cerr << "Failed to open input file: " << filePath << std::endl;
    return;
  }
  std::string line;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    std::string filePathEntry;
    std::string taxidStr;
    if (!(iss >> filePathEntry >> taxidStr)) {
      std::cerr << "Failed to parse line: " << line << std::endl;
      fileInfo.invalidNum++;
      continue;
    }
    hashCount[taxidStr] = 0;
    inputFiles[taxidStr].push_back(filePathEntry);
    fileInfo.fileNum++;
  }
  inputFile.close();
}

void createOrResetDirectory(const std::string &dir,
                            const BuildConfig &config) {
  try {
    std::filesystem::path directoryPath = dir;

    if (std::filesystem::exists(directoryPath)) {
      if (std::filesystem::is_directory(directoryPath)) {
        std::filesystem::remove_all(directoryPath);
        if (config.verbose) {
          std::cout << "Directory '" << dir << "' existed and was removed."
                    << std::endl;
        }
      } else {
        if (config.verbose) {
          std::cerr << "'" << dir
                    << "' exists but is not a directory, can't be replaced."
                    << std::endl;
        }
        return;
      }
    }

    if (std::filesystem::create_directories(directoryPath)) {
      if (config.verbose) {
        std::cout << "Directory '" << dir << "' created successfully."
                  << std::endl;
      }
    } else {
      if (config.verbose) {
        std::cerr << "Failed to create directory '" << dir << "'." << std::endl;
      }
    }
  } catch (const std::filesystem::filesystem_error &e) {
    if (config.verbose) {
      std::cerr << "Error: " << e.what() << std::endl;
    }
  }
}

} // namespace ChimeraBuild
