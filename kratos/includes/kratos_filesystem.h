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
namespace filesystem {

bool KRATOS_API(KRATOS_CORE) exists(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) is_regular_file(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) is_directory(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) create_directory(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) create_directories(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) remove(const std::string& rPath);

std::uintmax_t KRATOS_API(KRATOS_CORE) remove_all(const std::string& rPath);

void KRATOS_API(KRATOS_CORE) rename(const std::string& rPathFrom, const std::string& rPathTo);

} // namespace filesystem


namespace FilesystemExtensions {

std::string KRATOS_API(KRATOS_CORE) CurrentWorkingDirectory();

std::string KRATOS_API(KRATOS_CORE) JoinPaths(const std::vector<std::string>& rPaths);

} // namespace FilesystemExtensions
} // namespace Kratos

#endif // KRATOS_FILESYSTEM defined