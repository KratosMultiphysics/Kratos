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

#if !defined(KRATOS_FILESYSTEM)
#define KRATOS_FILESYSTEM

// System includes
#include <string>
#include <vector>

// External includes

// Project includes

namespace Kratos {
namespace FileSystem {

bool Exists(const std::string& rPath);

bool IsRegularFile(const std::string& rPath);

bool IsDirectory(const std::string& rPath);

bool CreateDirectory(const std::string& rPath);

bool CreateDirectories(const std::string& rPath);

std::string CurrentPath();

bool Remove(const std::string& rPath);

bool RemoveAll(const std::string& rPath);

void Rename(const std::string& rPathFrom, const std::string& rPathTo);

std::string JoinPaths(const std::vector<std::string>& rPaths);

} // namespace Kratos
} // namespace FileSystem

#endif // KRATOS_FILESYSTEM defined