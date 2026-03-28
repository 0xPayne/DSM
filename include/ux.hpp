// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#pragma once
#include <filesystem>
#include <string>
#include <vector>

namespace SparseLib::UX {

    // Resolve the project root directory from the executable path.
    std::filesystem::path projectRoot(const std::string& execPath);

    // List all files in data/ matching the given extension, sorted by name.
    std::vector<std::string> listFiles(const std::string& execPath,
                                       const std::string& extension = ".mtx");

    // Presents a numbered menu of matching files and returns the selected path.
    std::string selectFileFromDirectory(const std::string& execPath,
                                        const std::string& extension = ".mtx");

} // namespace SparseLib::UX
