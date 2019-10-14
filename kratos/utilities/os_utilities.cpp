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
#include "input_output/logger.h"
#include "utilities/os_utilities.h" // Has to be included before using KRATOS_COMPILED_IN_WINDOWS

// System includes
#include <sys/stat.h>
#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef KRATOS_COMPILED_IN_WINDOWS
    #include <direct.h>
    #include "Shlwapi.h"
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

void RemoveFile(const std::string& rFileName)
{
    const bool is_file = IsFile(rFileName);
    if (is_file) { // It is a file
        const char* name = rFileName.c_str();
        remove(name);
    } else {
        KRATOS_WARNING("OSUtilities") << "\"" << rFileName << "\" could not be removed because it is not a file" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RemoveFileInCurrentWorkingDir(const std::string& rFileName)
{
    RemoveFile(JoinPath(GetCurrentWorkingDir(), rFileName));
}

/***********************************************************************************/
/***********************************************************************************/

void RemoveDir(const std::string& rFolderName)
{
    const bool is_dir = IsDir(rFolderName);

    if (is_dir) { // It is a directory
        const char* name = rFolderName.c_str();
#ifdef KRATOS_COMPILED_IN_WINDOWS
        RemoveDirectoryA(name);
#else
        rmdir(name);
#endif
    } else {
        KRATOS_WARNING("OSUtilities") << "\"" << rFolderName << "\" could not be removed because it is not a folder" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RemoveDirInCurrentWorkingDir(const std::string& rFolderName)
{
    RemoveDir(JoinPath(GetCurrentWorkingDir(), rFolderName));
}

/***********************************************************************************/
/***********************************************************************************/

void Remove(const std::string& rFileName)
{
    const bool is_dir = IsDir(rFileName);

    if (is_dir) // It is a directory
        RemoveDir(rFileName);
    else  // It is a file
        RemoveFile(rFileName);
}

/***********************************************************************************/
/***********************************************************************************/

void RemoveInCurrentWorkingDir(const std::string& rFileName)
{
    Remove(JoinPath(GetCurrentWorkingDir(), rFileName));
}

/***********************************************************************************/
/***********************************************************************************/

std::string JoinPath(
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

bool IsFile(const std::string& rFileName)
{
    struct stat buffer;
    stat(rFileName.c_str(), &buffer);
    return (((buffer.st_mode) & S_IFMT) == S_IFREG);
}

/***********************************************************************************/
/***********************************************************************************/

bool IsDir(const std::string& rFolderName)
{
    struct stat buffer;
    stat(rFolderName.c_str(), &buffer);
    return (((buffer.st_mode) & S_IFMT) == S_IFDIR);
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
