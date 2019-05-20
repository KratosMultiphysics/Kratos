//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// Project includes
#include "utilities/os_utilities.h"

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

void RemoveOnCurrentWorkingDir(const std::string& rFileName)
{
#ifdef KRATOS_COMPILED_IN_WINDOWS
	if (IsDirExist(rFileName)) // It is a directory
		RemoveDirectoryA((GetCurrentWorkingDir() + "\\" + rFileName).c_str());
	else  // It is a file
		remove((GetCurrentWorkingDir() + "\\" + rFileName).c_str());
#else
    remove((GetCurrentWorkingDir() + "/" + rFileName).c_str());
#endif
}

/***********************************************************************************/
/***********************************************************************************/

bool IsDirExist(const std::string& rFolderName)
{
    struct stat buffer;
    return (stat (rFolderName.c_str(), &buffer) == 0);
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
