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
    KRATOS_DEPRECATED_MESSAGE("OSUtilities is deprecated, please use Kratos::filesystem") std::string KRATOS_API(KRATOS_CORE) GetCurrentWorkingDir();

}; // namespace OSUtilities
}  // namespace Kratos
#endif /* KRATOS_OS_UTILITIES defined */
