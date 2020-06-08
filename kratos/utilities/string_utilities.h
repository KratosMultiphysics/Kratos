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

#if !defined(KRATOS_STRING_UTILITIES)
#define KRATOS_STRING_UTILITIES

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
 * @namespace StringUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for computing string operations
 * @author Vicente Mataix Ferrandiz
 */
namespace StringUtilities
{
    /**
     * @brief This method converts CammelCase to snake_case
     * @param rString The string to be transformed into snake_case
     * @return The string in snake_case
     */
    std::string KRATOS_API(KRATOS_CORE) ConvertCammelCaseToSnakeCase(const std::string& rString);

    /**
     * @brief Erase first occurrence of given  substring from main string.
     * @param rMainString The string to be transformed
     * @param rToErase The string to remove
     * @return The string without the part to remove
     */
    std::string KRATOS_API(KRATOS_CORE) ErasePartialString(
        const std::string& rMainString,
        const std::string& rToErase
        );

    /**
     * @brief Checks the existence of a substring from main string.
     * @param rMainString The string to be transformed
     * @param rToCheck The string to search
     * @return True if the substring is found and false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) ContainsPartialString(
        const std::string& rMainString,
        const std::string& rToCheck
        );

    /**
     * @brief This method removes whitespaces
     * @param rString The string to be transformed
     * @return The string without white spaces
     */
    std::string KRATOS_API(KRATOS_CORE) RemoveWhiteSpaces(const std::string& rString);

}; // namespace StringUtilities
}  // namespace Kratos
#endif /* KRATOS_STRING_UTILITIES defined */
