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

// External includes
#include "ghc/filesystem.hpp" // TODO after moving to C++17 this can be removed since the functions can be used directly

// Project includes
#include "includes/kratos_filesystem.h"

namespace Kratos {
namespace filesystem {

bool exists(const std::string& rPath)
{
    return ghc::filesystem::exists(rPath);
}


bool is_regular_file(const std::string& rPath)
{
    return ghc::filesystem::is_regular_file(rPath);
}


bool is_directory(const std::string& rPath)
{
    return ghc::filesystem::is_directory(rPath);
}


bool create_directory(const std::string& rPath)
{
    return ghc::filesystem::create_directory(rPath);
}


bool create_directories(const std::string& rPath)
{
    return ghc::filesystem::create_directories(rPath);
}


bool remove(const std::string& rPath)
{
    return ghc::filesystem::remove(rPath);
}


std::uintmax_t remove_all(const std::string& rPath)
{
    return ghc::filesystem::remove_all(rPath);
}


void rename(const std::string& rPathFrom, const std::string& rPathTo)
{
    return ghc::filesystem::rename(rPathFrom, rPathTo);
}
    
std::string parent_path(const std::string& rPath)
{
    return ghc::filesystem::path(rPath).parent_path();
}

std::string filename(const std::string& rPath)
{
    return ghc::filesystem::path(rPath).filename();
}

} // namespace filesystem


namespace FilesystemExtensions {

std::string CurrentWorkingDirectory()
{
    return ghc::filesystem::current_path().string();
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
    
std::vector<std::string> ListDirectory(const std::string& rPath)
{
    std::vector<std::string> result;
    for (const auto& current_directory : ghc::filesystem::directory_iterator(rPath)) {
        result.push_back(current_directory.path());
    }
    return result;
}

} // namespace FilesystemExtensions
} // namespace Kratos
