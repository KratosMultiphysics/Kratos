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
// We haven't checked which filesystem to include yet
#ifndef INCLUDE_STD_FILESYSTEM_EXPERIMENTAL
// Check for feature test macro for <filesystem>
#   if defined(__cpp_lib_filesystem)
#     define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 0
// Check for feature test macro for <experimental/filesystem>
#   elif defined(__cpp_lib_experimental_filesystem)
#     define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 1
// We can't check if headers exist...
// Let's assume experimental to be safe
#   elif !defined(__has_include)
#     define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 1
// Check if the header "<filesystem>" exists
#   elif __has_include(<filesystem>)
// If we're compiling on Visual Studio and are not compiling with C++17, we need to use experimental
#     ifdef _MSC_VER
// Check and include header that defines "_HAS_CXX17"
#       if __has_include(<yvals_core.h>)
#         include <yvals_core.h>
// Check for enabled C++17 support
#         if defined(_HAS_CXX17) && _HAS_CXX17
// We're using C++17, so let's use the normal version
#           define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 0
#         endif
#       endif
// If the marco isn't defined yet, that means any of the other VS specific checks failed, so we need to use experimental
#       ifndef INCLUDE_STD_FILESYSTEM_EXPERIMENTAL
#         define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 1
#       endif
// Not on Visual Studio. Let's use the normal version
#     else // #ifdef _MSC_VER
#       define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 0
#     endif
// Check if the header "<filesystem>" exists
#   elif __has_include(<experimental/filesystem>)
#     define INCLUDE_STD_FILESYSTEM_EXPERIMENTAL 1
// Fail if neither header is available with a nice error message
#   else
#     error Could not find system header "<filesystem>" or "<experimental/filesystem>"
#   endif
// We priously determined that we need the exprimental version
#   if INCLUDE_STD_FILESYSTEM_EXPERIMENTAL
#     include <experimental/filesystem>
// We need the alias from std::experimental::filesystem to std::filesystem
namespace std {
  namespace filesystem = experimental::filesystem;
}
// We have a decent compiler and can use the normal version
#   else
#     include <filesystem>
#   endif
#endif // #ifndef INCLUDE_STD_FILESYSTEM_EXPERIMENTAL

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
