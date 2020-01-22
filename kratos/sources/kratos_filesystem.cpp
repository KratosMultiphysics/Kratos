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
// #include "ghc/filesystem.hpp" // TODO after moving to C++17 this can be exchanged with std::filesystem

// Project includes
#include "includes/kratos_filesystem.h"

namespace Kratos {
namespace FileSystem {

bool Exists(const std::string& rPath)
{
    return false;
    // return ghc::filesystem::exists(rPath);
}


bool IsRegularFile(const std::string& rPath)
{
    return false;
    // return ghc::filesystem::is_regular_file(rPath);
}


bool IsDirectory(const std::string& rPath)
{
    return false;
    // return ghc::filesystem::is_directory(rPath);
}


bool CreateDirectory(const std::string& rPath)
{
    return false;
    // return ghc::filesystem::create_directory(rPath);
}


bool CreateDirectories(const std::string& rPath)
{
    return false;
    // return ghc::filesystem::create_directories(rPath);
}


std::string CurrentPath()
{
    return "";
    // return ghc::filesystem::current_path().string();
}


bool Remove(const std::string& rPath)
{
    return false;
    // return ghc::filesystem::remove(rPath);
}


std::uintmax_t RemoveAll(const std::string& rPath)
{
    return 0;
    // return ghc::filesystem::remove_all(rPath);
}


void Rename(const std::string& rPathFrom, const std::string& rPathTo)
{
    // return ghc::filesystem::rename(rPathFrom, rPathTo);
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

} // namespace Kratos
} // namespace FileSystem