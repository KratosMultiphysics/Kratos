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
#include "includes/define.h"


namespace Kratos {
// Functions that provide some basic functionalities that are provided by std::filesystem (part of C++17)
// please check the documentation of std::filesystem for the function documentation
namespace FileSystem {

bool KRATOS_API(KRATOS_CORE) Exists(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) IsRegularFile(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) IsDirectory(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) CreateDirectory2(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) CreateDirectories(const std::string& rPath);

std::string KRATOS_API(KRATOS_CORE) CurrentPath();

bool KRATOS_API(KRATOS_CORE) Remove(const std::string& rPath);

std::uintmax_t KRATOS_API(KRATOS_CORE) RemoveAll(const std::string& rPath);

void KRATOS_API(KRATOS_CORE) Rename(const std::string& rPathFrom, const std::string& rPathTo);

std::string KRATOS_API(KRATOS_CORE) JoinPaths(const std::vector<std::string>& rPaths);

} // namespace Kratos
} // namespace FileSystem

#endif // KRATOS_FILESYSTEM defined