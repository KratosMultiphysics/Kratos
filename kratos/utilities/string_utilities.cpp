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

// System includes
#include <algorithm>
#include <iostream>
#include <cctype>

// External includes

// Project includes
#include "utilities/string_utilities.h"

namespace Kratos
{
namespace StringUtilities
{
std::string ConvertCammelCaseToSnakeCase(const std::string& rString)
{
    std::string str(1, tolower(rString[0]));

    // First place underscores between contiguous lower and upper case letters.
    // For example, `_LowerCamelCase` becomes `_Lower_Camel_Case`.
    for (auto it = rString.begin() + 1; it != rString.end(); ++it) {
        if (isupper(*it) && *(it-1) != '_' && islower(*(it-1))) {
            str += "_";
        }
        str += *it;
    }

    // Then convert it to lower case.
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);

    return str;
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace StringUtilities
} // namespace Kratos
