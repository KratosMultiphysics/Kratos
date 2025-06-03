//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
#include "input_output/logger.h"
#include "utilities/string_utilities.h"

namespace Kratos::StringUtilities {

std::string ConvertCamelCaseToSnakeCase(const std::string& rString)
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

std::string ConvertSnakeCaseToCamelCase(const std::string& rString)
{
    std::string output;

    if (!rString.empty()) {
        output.reserve(rString.size());
        bool upper_switch = rString[0] == '_' ? false : true;

        for (auto character : rString) {
            KRATOS_ERROR_IF(!(std::isalnum(character) || character == '_') || std::isupper(character))
                << "Invalid character '" << character
                <<"' in snake case string '" << rString << '\'';

            if (character == '_') {
                KRATOS_ERROR_IF(upper_switch)
                    << "Repeated underscores in snake case string '" << rString << '\'';
                upper_switch = true;
            } else { // character != '_'
                // At this point, the character must be in [a-z0-9]
                if (upper_switch) {
                    output.push_back(std::toupper(character));
                    upper_switch = false;
                } else { // !upper_switch
                    output.push_back(character);
                } // else (upper_switch)
            } // else (character == '_')
        } // for character in rString
    } // if rString

    return output;
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
    const std::string& sub_string = rMainString;

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
    const char Delimiter
    )
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

/***********************************************************************************/
/***********************************************************************************/

std::string Trim(
    const std::string& rInputString,
    const bool RemoveNullChar)
{
    return TrimLeft(TrimRight(rInputString, RemoveNullChar), RemoveNullChar);
}

/***********************************************************************************/
/***********************************************************************************/

std::function<bool(std::string::value_type)> TrimChar(const bool RemoveNullChar)
{
    if (RemoveNullChar) {
        return [](auto character) {
            return std::isspace(character) || character == '\0';
        };
    }

    return [](auto character) {
        return std::isspace(character);
    };
}

/***********************************************************************************/
/***********************************************************************************/

std::string TrimLeft(
    const std::string& rInputString,
    const bool RemoveNullChar)
{
    std::string output_string(rInputString);

    const auto trim_char = TrimChar(RemoveNullChar);

    output_string.erase(output_string.begin(), std::find_if(output_string.begin(), output_string.end(),
            [trim_char](std::string::value_type ch) {return !trim_char(ch);}
        )
    );

    return output_string;
}

/***********************************************************************************/
/***********************************************************************************/

std::string TrimRight(
    const std::string& rInputString,
    const bool RemoveNullChar)
{
    std::string output_string(rInputString);
    const auto trim_char = TrimChar(RemoveNullChar);

    output_string.erase(std::find_if(output_string.rbegin(), output_string.rend(),
            [trim_char](std::string::value_type ch) {return !trim_char(ch);}
        ).base(), output_string.end()
    );

    return output_string;
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> SplitStringIntoAVector(const std::string& rText)
{
    std::istringstream iss(rText);
    std::string word;
    std::vector<std::string> words;

    // The stream extraction operator (>>) automatically skips whitespace
    // (including spaces and tabs) by default.
    while (iss >> word) {
        words.push_back(word);
    }
    return words;
};

/***********************************************************************************/
/***********************************************************************************/

template <typename TType>
std::vector<TType> StringToVector(const std::string& rInputString)
{
    std::vector<TType> result;
    std::string str = rInputString; // Make a mutable copy

    // 1. Remove '[' and ']' from the outer string if they exist
    // Trim the input string first to handle cases like "  [1,2,3]  "
    str = Trim(str);
    if (!str.empty() && str.front() == '[') {
        str.erase(0, 1);
    }
    if (!str.empty() && str.back() == ']') {
        str.pop_back();
    }

    // If the string is empty after removing brackets and trimming (e.g., "[]" or "[ ]" or ""), return an empty vector
    if (Trim(str).empty()) {
        return result;
    }

    // 2. Use stringstream to split by comma and parse elements
    std::stringstream ss(str);
    std::string segment;

    // Split the string by comma
    while (std::getline(ss, segment, ',')) {
        // Trim whitespace from the individual segment
        std::string trimmedSegment = Trim(segment);
        if (!trimmedSegment.empty()) {
            TType value;
            std::stringstream converter(trimmedSegment);

            // Try to convert the segment to type TType
            // The "converter >> value" attempts the conversion.
            // The "!(converter >> std::ws).eof()" checks if there's anything left in the stream
            // after the value and any trailing whitespace. If there is, it's an error (e.g. "1.2abc" or "12x").
            if (converter >> value && (converter >> std::ws).eof()) {
                result.push_back(value);
            } else {
                // Handle cases where conversion is not possible (e.g., non-numeric characters)
                KRATOS_WARNING("StringUtilities") << "Invalid argument for conversion: '" << trimmedSegment << "'" << std::endl;
            }
        }
    }
    return result;
}

// Explicit instantiation of the template function for double and integer types
template KRATOS_API(KRATOS_CORE) std::vector<double> StringToVector<double>(const std::string& rInputString);
template KRATOS_API(KRATOS_CORE) std::vector<int> StringToVector<int>(const std::string& rInputString);
template KRATOS_API(KRATOS_CORE) std::vector<unsigned int> StringToVector<unsigned int>(const std::string& rInputString);
template KRATOS_API(KRATOS_CORE) std::vector<std::size_t> StringToVector<std::size_t>(const std::string& rInputString);

/***********************************************************************************/
/***********************************************************************************/

std::size_t CountValuesUntilPrefix(
    const std::string& rInputString,
    const std::string& rStopPrefix
    )
{
    std::istringstream iss(rInputString); // Use a string stream to easily extract tokens
    std::string token;
    std::size_t count = 0;

    while (iss >> token) { // Read tokens one by one
        // Check if the token contains the stop prefix
        // If the token contains the stop prefix, we stop counting
        // We use find to check if the token contains the stop prefix
        // Note: find returns std::string::npos if the substring is not found
        // If the token contains the stop prefix, we stop counting
        // We do not count the token that contains the stop prefix
        if (token.find(rStopPrefix) != std::string::npos) {
            break; // Stop condition met, do not count this token or subsequent ones
        }
        count++; // Increment count for tokens that do not meet the stop condition
    }
    return count;
}

} // namespace Kratos::StringUtilities
