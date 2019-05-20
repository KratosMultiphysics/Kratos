//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// Project includes
#include "utilities/os_utilities.h" // Has to be included before using KRATOS_COMPILED_IN_WINDOWS

// System includes
#include <sys/stat.h>
#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef KRATOS_COMPILED_IN_WINDOWS
#include <direct.h>
#include "Shlwapi.h"
#ifndef S_ISDIR
#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
#endif
#ifndef S_ISREG
#define S_ISREG(mode)  (((mode) & S_IFMT) == S_IFREG)
#endif
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

// External includes

namespace Kratos
{
namespace OSUtilities
{
std::string GetCurrentWorkingDir()
{
    KRATOS_TRY

    char buff[FILENAME_MAX];
    GetCurrentDir( buff, FILENAME_MAX );
    std::string current_working_dir(buff);
    return current_working_dir;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void Remove(const std::string& rFileName)
{
    const char* name = rFileName.c_str();
    const bool is_dir = DirExist(rFileName);

#ifdef KRATOS_COMPILED_IN_WINDOWS
    if (is_dir) // It is a directory
        RemoveDirectoryA(name);
    else  // It is a file
        remove(name);
#else
    if (is_dir) // It is a directory
        rmdir(name);
    else  // It is a file
        remove(name);
#endif
}

/***********************************************************************************/
/***********************************************************************************/

void RemoveOnCurrentWorkingDir(const std::string& rFileName)
{
    Remove(AppendFolderFile(GetCurrentWorkingDir(), rFileName));
}

/***********************************************************************************/
/***********************************************************************************/

std::string AppendFolderFile(
    const std::string& rFolderName,
    const std::string& rFileName
    )
{
#ifdef KRATOS_COMPILED_IN_WINDOWS
    return rFolderName + "\\" + rFileName;
#else
    return rFolderName + "/" + rFileName;
#endif
}

/***********************************************************************************/
/***********************************************************************************/

bool FileExist(const std::string& rFileName)
{
    struct stat buffer;
    stat(rFileName.c_str(), &buffer);
    return S_ISREG(buffer.st_mode);
}

/***********************************************************************************/
/***********************************************************************************/

bool DirExist(const std::string& rFolderName)
{
    struct stat buffer;
    stat(rFolderName.c_str(), &buffer);
    return S_ISDIR(buffer.st_mode);
}

/***********************************************************************************/
/***********************************************************************************/

int CreateDir(const std::string& rFolderName)
{
#ifdef KRATOS_COMPILED_IN_WINDOWS
    const int status = mkdir(rFolderName.c_str());
#else
    const int status = mkdir(rFolderName.c_str(),0777);
#endif
    return status;
}

} // namespace OSUtilities
} // namespace Kratos
