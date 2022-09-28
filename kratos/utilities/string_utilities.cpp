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
#include <sstream>

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

std::string ErasePartialString(
    const std::string& rMainString,
    const std::string& rToErase
    )
{
    // Value to return
    std::string sub_string = rMainString;

    // Search for the substring in string
    std::size_t pos = sub_string.find(rToErase);

    if (pos != std::string::npos) {
        // If found then erase it from string
        sub_string.erase(pos, rToErase.length());
    }

    return sub_string;
}

/***********************************************************************************/
/***********************************************************************************/

bool ContainsPartialString(
    const std::string& rMainString,
    const std::string& rToCheck
    )
{
    // Value to return
    std::string sub_string = rMainString;

    // Search for the substring in string
    std::size_t pos = sub_string.find(rToCheck);

    // Return true if found
    if (pos != std::string::npos) {
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

std::string RemoveWhiteSpaces(const std::string& rString)
{
    // Value to return
    std::string output;

    for(char c : rString) {
        if(!std::isspace(c)) {
            output += c ;
        }
    }

    return output;
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> SplitStringByDelimiter(
    const std::string& rString,
    const char Delimiter)
{
    std::istringstream ss(rString);
    std::string token;

    std::vector<std::string> splitted_string;
    while(std::getline(ss, token, Delimiter)) {
        splitted_string.push_back(token);
    }

    return splitted_string;
}

/***********************************************************************************/
/***********************************************************************************/

std::string ReplaceAllSubstrings(
    const std::string& rInputString,
    const std::string& rStringToBeReplaced,
    const std::string& rStringToReplace
    )
{
    std::string output_string(rInputString);
    std::size_t start_pos = 0;
    while((start_pos = output_string.find(rStringToBeReplaced, start_pos)) != std::string::npos) {
        output_string.replace(start_pos, rStringToBeReplaced.length(), rStringToReplace);
        start_pos += rStringToReplace.length(); // Handles case where 'to' is a substring of 'from'
    }
    return output_string;
}

} // namespace StringUtilities
} // namespace Kratos
