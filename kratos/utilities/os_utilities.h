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

#if !defined(KRATOS_OS_UTILITIES)
#define KRATOS_OS_UTILITIES

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/**
 * @namespace OSUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries to access to OS information
 * @author Vicente Mataix Ferrandiz
 */
namespace OSUtilities
{
    /**
     * @brief This method returns the current working directory
     * @return The current working directory
     */
    std::string KRATOS_API(KRATOS_CORE) GetCurrentWorkingDir();

    /**
     * @brief This method removes a file on the current working directory
     * @param rFileName The name of the file to be removed
     */
    void KRATOS_API(KRATOS_CORE) RemoveOnCurrentWorkingDir(const std::string& rFileName);

    /**
     * @brief This method checks if the directory exists
     * @param rFolderName The name of the directory to be checked
     * @return True if exists, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) IsDirExist(const std::string& rFolderName);

    /**
     * @brief This method creates a new directory
     * @param rFolderName The name of the directory to be created
     * @return The status
     */
    int KRATOS_API(KRATOS_CORE) CreateDir(const std::string& rFolderName);

}; // namespace OSUtilities
}  // namespace Kratos
#endif /* KRATOS_OS_UTILITIES defined */
