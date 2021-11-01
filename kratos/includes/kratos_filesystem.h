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
// wrapper functions for std::filesystem (part of C++17)
// the function signatures are identical, hence after moving to C++17 Kratos::filesystem can be replaced with std::filesystem
// please check the documentation of std::filesystem for the function documentation

// Note: the filesystem functinos have a filesystem::path as input, but currently std::string is used as filesystem::path is not available
// this should not be a problem for upgrading to std::filesystem, since filesystem::path has a constructor accepting a string
namespace filesystem {

bool KRATOS_API(KRATOS_CORE) exists(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) is_regular_file(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) is_directory(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) create_directory(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) create_directories(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) remove(const std::string& rPath);

std::uintmax_t KRATOS_API(KRATOS_CORE) remove_all(const std::string& rPath);

void KRATOS_API(KRATOS_CORE) rename(const std::string& rPathFrom, const std::string& rPathTo);

std::string KRATOS_API(KRATOS_CORE) parent_path(const std::string& rPath);

std::string KRATOS_API(KRATOS_CORE) filename(const std::string& rPath);

} // namespace filesystem


namespace FilesystemExtensions {
// helper functions related to filesystem

/**
 * @brief Returns current working directory
 *
 * @return std::string
 */
std::string KRATOS_API(KRATOS_CORE) CurrentWorkingDirectory();

/**
 * @brief Join paths
 *
 * @param rPaths                        List of strings to be joined to get final path
 * @return std::string                  Final joined path
 */
std::string KRATOS_API(KRATOS_CORE) JoinPaths(const std::vector<std::string>& rPaths);

/**
 * @brief Returns list of files and directories in rPath
 *
 * @param rPath                         Path
 * @return std::vector<std::string>     List of files and folders in rPath
 */
std::vector<std::string> KRATOS_API(KRATOS_CORE) ListDirectory(const std::string& rPath);

/**
 * @brief Create directories in MPI, when sometimes filesystems are slow. Intended to be called by all ranks (that make use of this directory). It returns only after the folder exists
 *
 * @param rPath                         Path
 */
void KRATOS_API(KRATOS_CORE) MPISafeCreateDirectories(const std::string& rPath);

} // namespace FilesystemExtensions
} // namespace Kratos

#endif // KRATOS_FILESYSTEM defined
