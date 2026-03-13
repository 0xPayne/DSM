#include "../include/ux.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace SparseLib::UX {

namespace fs = std::filesystem;

std::string selectFileFromDirectory(const std::string& execPath,
                                    const std::string& extension) {

    const fs::path dataDir = fs::canonical(
        fs::path(execPath).parent_path() / "../data"
    );

    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(dataDir)) {
        if (entry.path().extension() == extension)
            files.push_back(entry.path().string());
    }

    if (files.empty())
        throw std::runtime_error(
            "No " + extension + " files found in " + dataDir.string()
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
} 
