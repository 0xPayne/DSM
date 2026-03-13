#pragma once
#include <string>

namespace SparseLib::UX {

    // Scans the /data directory (relative) for files matching .mtx
    // presents a numbered list to the user, and returns the path of the selected file.
    std::string selectFileFromDirectory(const std::string& execPath,
                                        const std::string& extension = ".mtx");
}
