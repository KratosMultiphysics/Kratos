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
     * @brief This method removes a file
     * @param rFileName The name of the file to be removed
     */
    void KRATOS_API(KRATOS_CORE) RemoveFile(const std::string& rFileName);
    
    /**
     * @brief This method removes a file in the current working directory
     * @param rFileName The name of the file to be removed
     */
    void KRATOS_API(KRATOS_CORE) RemoveFileInCurrentWorkingDir(const std::string& rFileName);

    /**
     * @brief This method removes a directory
     * @param rFileName The name of the directory to be removed
     */
    void KRATOS_API(KRATOS_CORE) RemoveDir(const std::string& rFolderName);

    /**
     * @brief This method removes a directory in the current working directory
     * @param rFileName The name of the directory to be removed
     */
    void KRATOS_API(KRATOS_CORE) RemoveDirInCurrentWorkingDir(const std::string& rFolderName);

    /**
     * @brief This method removes a file
     * @param rFileName The name of the file to be removed
     */
    void KRATOS_API(KRATOS_CORE) Remove(const std::string& rFileName);

    /**
     * @brief This method removes a file in the current working directory
     * @param rFileName The name of the file to be removed
     */
    void KRATOS_API(KRATOS_CORE) RemoveInCurrentWorkingDir(const std::string& rFileName);

    /**
     * @brief This method appends the name of the folder and the file
     * @param rFolderName The name of the directory to be appended
     * @param rFileName The name of the file to be appended
     */
    std::string KRATOS_API(KRATOS_CORE) JoinPath(
        const std::string& rFolderName,
        const std::string& rFileName
        );

    /**
     * @brief This method checks if the file exists
     * @param rFileName The name of the file to be checked
     * @return True if exists, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) IsFile(const std::string& rFileName);

    /**
     * @brief This method checks if the directory exists
     * @param rFolderName The name of the directory to be checked
     * @return True if exists, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) IsDir(const std::string& rFolderName);

    /**
     * @brief This method creates a new directory
     * @param rFolderName The name of the directory to be created
     * @return The status
     */
    int KRATOS_API(KRATOS_CORE) CreateDir(const std::string& rFolderName);

}; // namespace OSUtilities
}  // namespace Kratos
#endif /* KRATOS_OS_UTILITIES defined */
