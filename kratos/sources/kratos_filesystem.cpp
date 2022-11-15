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
#include "includes/kratos_filesystem.h"

namespace Kratos {
namespace filesystem {

bool exists(const std::string& rPath)
{
    return std::filesystem::exists(rPath);
}


bool is_regular_file(const std::string& rPath)
{
    return std::filesystem::is_regular_file(rPath);
}


bool is_directory(const std::string& rPath)
{
    return std::filesystem::is_directory(rPath);
}


bool create_directory(const std::string& rPath)
{
    return std::filesystem::create_directory(rPath);
}


bool create_directories(const std::string& rPath)
{
    return std::filesystem::create_directories(rPath);
}


bool remove(const std::string& rPath)
{
    return std::filesystem::remove(rPath);
}


std::uintmax_t remove_all(const std::string& rPath)
{
    return std::filesystem::remove_all(rPath);
}


void rename(const std::string& rPathFrom, const std::string& rPathTo)
{
    return std::filesystem::rename(rPathFrom, rPathTo);
}

std::string parent_path(const std::string& rPath)
{
    return std::filesystem::path(rPath).parent_path().string();
}

std::string filename(const std::string& rPath)
{
    return std::filesystem::path(rPath).filename().string();
}

} // namespace filesystem
} // namespace Kratos
