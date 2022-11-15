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
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos {

namespace FilesystemExtensions {
// helper functions related to filesystem

/**
 * @brief Returns current working directory
 *
 * @return std::string
 */
KRATOS_DEPRECATED_MESSAGE("Please use std::filesystem directly")
std::string KRATOS_API(KRATOS_CORE) CurrentWorkingDirectory();

/**
 * @brief Join paths
 *
 * @param rPaths                        List of strings to be joined to get final path
 * @return std::string                  Final joined path
 */
KRATOS_DEPRECATED_MESSAGE("Please use the /-operator directly")
std::string KRATOS_API(KRATOS_CORE) JoinPaths(const std::vector<std::string>& rPaths);

/**
 * @brief Returns list of files and directories in rPath
 *
 * @param rPath                               Path
 * @return std::vector<std::filesystem::path> List of files and folders in rPath
 */
std::vector<std::filesystem::path> KRATOS_API(KRATOS_CORE) ListDirectory(const std::filesystem::path& rPath);

/**
 * @brief Create directories in MPI, when sometimes filesystems are slow. Intended to be called by all ranks (that make use of this directory). It returns only after the folder exists
 *
 * @param rPath                         Path
 */
void KRATOS_API(KRATOS_CORE) MPISafeCreateDirectories(const std::filesystem::path& rPath);

/** @brief Resolve symlinks recursively.
 *
 *  @param rPath: path to a symbolic link.
 *  @return The result of the recursive dereferencing.
 *  @throws If the input path does not exist or the symlink is cyclic.
 *  @note The existence of the final result is not checked and is up to the user.
 *  @note The input is returned if it is not a symlink.
 */
std::filesystem::path ResolveSymlinks(const std::filesystem::path& rPath);

} // namespace FilesystemExtensions
} // namespace Kratos
