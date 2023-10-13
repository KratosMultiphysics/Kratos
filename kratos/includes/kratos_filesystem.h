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

#pragma once

// System includes
#include <filesystem>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos {

// deprecated namespaces wrongly issue a warning with GCC < 10, hence disabling until removed
#if defined(__GNUG__) && __GNUC__ < 10 && !defined(__clang__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#endif

namespace KRATOS_DEPRECATED_MESSAGE("Please use std::filesystem directly") filesystem {

#if defined(__GNUG__) && __GNUC__ < 10 && !defined(__clang__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic pop
#endif

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


class KRATOS_API(KRATOS_CORE) FilesystemExtensions
{

public:
    /// Default constructor.
    FilesystemExtensions() = delete;

    /// Copy constructor.
    FilesystemExtensions(FilesystemExtensions const& rOther) = delete;

    /// Assignment operator.
    FilesystemExtensions& operator=(FilesystemExtensions const& rOther) = delete;

    // helper functions related to filesystem

    /**
     * @brief Returns current working directory
     *
     * @return std::string
     */
    KRATOS_DEPRECATED_MESSAGE("Please use std::filesystem directly")
    static std::string CurrentWorkingDirectory();

    /**
     * @brief Join paths
     *
     * @param rPaths                        List of strings to be joined to get final path
     * @return std::string                  Final joined path
     */
    KRATOS_DEPRECATED_MESSAGE("Please use the /-operator directly")
    static std::string JoinPaths(const std::vector<std::string>& rPaths);

    /**
     * @brief Returns list of files and directories in rPath
     *
     * @param rPath                               Path
     * @return std::vector<std::filesystem::path> List of files and folders in rPath
     */
    [[nodiscard]] static std::vector<std::filesystem::path> ListDirectory(const std::filesystem::path& rPath);

    /**
     * @brief Create directories in MPI, when sometimes filesystems are slow. Intended to be called by all ranks (that make use of this directory). It returns only after the folder exists
     *
     * @param rPath                         Path
     */
    static void MPISafeCreateDirectories(const std::filesystem::path& rPath);

    /** @brief Resolve symlinks recursively.
     *
     *  @param rPath: path to a symbolic link.
     *  @return The result of the recursive dereferencing.
     *  @throws If the input path does not exist or the symlink is cyclic.
     *  @note The existence of the final result is not checked and is up to the user.
     *  @note The input is returned if it is not a symlink.
     */
    [[nodiscard]] static std::filesystem::path ResolveSymlinks(const std::filesystem::path& rPath);

}; // class FilesystemExtensions

} // namespace Kratos
