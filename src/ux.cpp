// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/ux.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace SparseLib::UX {

namespace fs = std::filesystem;

fs::path projectRoot(const std::string& execPath) {
  return fs::canonical(fs::path(execPath).parent_path() / "..");
}

std::vector<std::string> listFiles(const std::string& execPath,
                                   const std::string& extension) {
  fs::path dataDir = projectRoot(execPath) / "data";
  std::vector<std::string> files;
  for (const auto& entry : fs::directory_iterator(dataDir)) {
    if (entry.path().extension() == extension)
      files.push_back(entry.path().string());
  }
  std::sort(files.begin(), files.end());
  return files;
}

std::string selectFileFromDirectory(const std::string& execPath,
                                    const std::string& extension) {
  auto files = listFiles(execPath, extension);

  if (files.empty())
    throw std::runtime_error(
        "No " + extension + " files found in data/"
    );

  std::cout << "\nAvailable files:\n";
  for (size_t i = 0; i < files.size(); ++i)
    std::cout << "  [" << i + 1 << "] " << files[i] << "\n";

  int choice = 0;
  std::cout << "\nSelect a file (1-" << files.size() << "): ";
  std::cin >> choice;

  if (choice < 1 || choice > static_cast<int>(files.size()))
    throw std::out_of_range("Invalid selection: " + std::to_string(choice));

  return files[choice - 1];
}

} // namespace SparseLib::UX
