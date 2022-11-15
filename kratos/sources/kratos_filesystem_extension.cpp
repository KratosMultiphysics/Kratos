//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <algorithm>
#include <thread>
#include <chrono>
#include <set>

// External includes

// Project includes
#include "includes/kratos_filesystem_extension.h"

namespace Kratos {

namespace FilesystemExtensions {

std::string CurrentWorkingDirectory()
{
    return std::filesystem::current_path().string();
}

std::string JoinPaths(const std::vector<std::string>& rPaths)
{
    auto paths(rPaths); // create local copy

    // first remove empty paths
    paths.erase(std::remove_if(paths.begin(), paths.end(),
                         [](const std::string& s)
                         { return s.empty(); }), paths.end());

    const std::size_t num_paths = paths.size();

    if (num_paths == 0) { return ""; }

    std::string full_path = paths[0];
    if (num_paths > 1) {
        for(std::size_t i=1; i<num_paths; ++i) {
            full_path += "/" + paths[i]; // using portable separator "/"
        }
    }

    return full_path;
}

std::vector<std::filesystem::path> ListDirectory(const std::filesystem::path& rPath)
{
    std::vector<std::filesystem::path> result;
    for (const auto& current_directory : std::filesystem::directory_iterator(rPath)) {
        result.push_back(current_directory.path());
    }
    return result;
}

void MPISafeCreateDirectories(const std::filesystem::path& rPath)
{
    if (!std::filesystem::exists(rPath)) {
        std::filesystem::create_directories(rPath);
    }
    if (!std::filesystem::exists(rPath)) { // wait for the path to appear in the filesystem
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
}


std::filesystem::path ResolveSymlinks(const std::filesystem::path& rPath)
{
    auto status = std::filesystem::symlink_status(rPath);
    KRATOS_ERROR_IF(status.type() == std::filesystem::file_type::not_found) << "File not found: " << rPath;

    std::filesystem::path path = rPath;
    std::set<std::filesystem::path> symlinks;

    while (status.type() == std::filesystem::file_type::symlink) {
        const auto insert_result = symlinks.insert(path);
        KRATOS_ERROR_IF_NOT(insert_result.second) << rPath << " leads to cyclic symlinks";
        path = std::filesystem::read_symlink(path);
        status = std::filesystem::symlink_status(path);
    }

    return path;
}

} // namespace FilesystemExtensions
} // namespace Kratos
