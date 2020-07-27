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

} // namespace filesystem


namespace FilesystemExtensions {

std::string CurrentWorkingDirectory()
{
    return ghc::filesystem::current_path().string();
}

std::string JoinPaths(const std::vector<std::string>& rPaths)
{
    const std::size_t num_paths = rPaths.size();

    if (num_paths == 0) { return ""; }

    std::string full_path = rPaths[0];
    if (num_paths > 1) {
        for(std::size_t i=1; i<num_paths; ++i) {
            full_path += "/" + rPaths[i]; // using portable separator "/"
        }
    }

    return full_path;
}

} // namespace FilesystemExtensions
} // namespace Kratos