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


} // namespace Kratos::StringUtilities



namespace Kratos {


std::pair<std::string,std::regex> RegexUtility::Integer()
{
    std::pair<std::string,std::regex> output;

    output.first = R"(0|(?:-?[1-9]+[0-9]*))";
    output.second = std::regex(output.first);

    return output;
}


std::pair<std::string,std::regex> RegexUtility::UnsignedInteger()
{
    std::pair<std::string,std::regex> output;

    output.first = R"(0|(?:[1-9]+[0-9]*))";
    output.second = std::regex(output.first);

    return output;
}


std::pair<std::string,std::regex> RegexUtility::FloatingPoint()
{
    std::pair<std::string,std::regex> output;

    // Clutter due to the many uncapturing groups
    output.first = R"(-?(?:(?:(?:[1-9][0-9]*)(?:\.[0-9]*)?)|(?:0(?:\.[0-9]*)?))(?:[eE][\+-]?[0-9]+)?)";
    output.second = std::regex(output.first);

    return output;
}


PlaceholderPattern::PlaceholderPattern(const std::string& rPattern,
                                       const PlaceholderMap& rPlaceholderMap)
    : mPattern(rPattern),
      mPlaceholderGroupMap(),
      mRegexString(FormatRegexLiteral(rPattern)),
      mRegex()
{
    KRATOS_TRY

    using PositionPair = std::pair<std::size_t,PlaceholderMap::const_iterator>;
    const auto number_of_placeholders = rPlaceholderMap.size();

    // Array for tracking position-placeholder pairs for assigning
    // regex groups to placeholders later.
    // {position_in_pattern, it_placeholder_regex_pair}
    std::vector<PositionPair> position_map;
    position_map.reserve(number_of_placeholders);

    // Replace placeholders with their corresponding regex strings
    // and note their positions in the pattern
    PlaceholderMap::const_iterator it_pair = rPlaceholderMap.begin();
    const auto it_end = rPlaceholderMap.end();

    for ( ; it_pair!=it_end; ++it_pair){
        // Make sure that input pattern has no capturing groups of its own.
        KRATOS_ERROR_IF(std::regex(it_pair->second).mark_count())
            << "pattern " << it_pair->second
            << " of placeholder '" << it_pair->first << "'"
            << " has internal capturing group(s) (this is forbidden in PlaceholderPattern)";

        // Wrap the regex in a group
        std::string regex = "(" + it_pair->second + ")";

        const auto placeholder_size = it_pair->first.size();
        const auto regex_size = regex.size();
        const int size_difference = regex_size - placeholder_size;

        while (true) {
            // Find the next instance of the current placeholder
            const auto position_in_pattern = mRegexString.find(it_pair->first);

            // Replace it with its regex (if found)
            if (position_in_pattern != std::string::npos) {
                mRegexString.replace(position_in_pattern, placeholder_size, regex);
            } else {
                break;
            }

            // Update positions
            for (auto& r_pair : position_map) {
                if (position_in_pattern < r_pair.first) r_pair.first += size_difference;
            }

            position_map.emplace_back(position_in_pattern, it_pair);
        } // while placeholder in pattern
    } // for placeholder in rPlaceHolderMap

    // Replace positions with indices in ascending order (based on position)
    // in lieu of std::transform_if
    std::sort(
        position_map.begin(),
        position_map.end(),
        [](const PositionPair& rLeft, const PositionPair& rRight) {return rLeft.first < rRight.first;});

    std::size_t index = 0;
    for (auto& r_pair : position_map) r_pair.first = index++;

    // Populate the placeholder - group index map
    for (auto it_pair=rPlaceholderMap.begin() ; it_pair!=it_end; ++it_pair) {
        // Move the placeholder string and construct an associated empty index array
        auto emplace_result = mPlaceholderGroupMap.emplace(it_pair->first, PlaceholderGroupMap::mapped_type(
            PlaceholderGroupMap::mapped_type::value_type()
        ));

        // Fill the index array with the group indices
        for (const auto& r_pair : position_map) {
            if (r_pair.second == it_pair) emplace_result.first->second.value().push_back(r_pair.first);
        }
    }

    // Disable placeholders that aren't in the pattern
    for (auto& r_pair : mPlaceholderGroupMap) {
        if (r_pair.second.value().empty()) {
            r_pair.second = PlaceholderGroupMap::mapped_type();
        } // if placeholder was not found in the pattern
    } // for key, value in mPlaceholderGroupMap

    // Construct the regex
    mRegexString = "^" + mRegexString + "$";
    mRegex = std::regex(mRegexString);

    KRATOS_CATCH("");
} // PlaceholderPattern::PlaceholderPattern


bool PlaceholderPattern::IsAMatch(const std::string& rString) const
{
    KRATOS_TRY

    return std::regex_match(rString, mRegex);

    KRATOS_CATCH("");
} // PlaceholderPattern::IsAMatch


PlaceholderPattern::MatchType PlaceholderPattern::Match(const std::string& rString) const
{
    KRATOS_TRY

    std::smatch results;
    MatchType output;

    // Perform regex search and extract matches
    if (std::regex_match(rString, results, mRegex)) {
        for (auto& r_pair : mPlaceholderGroupMap) {
            if (r_pair.second.has_value()) {
                // Construct empty group matches
                auto emplace_result = output.emplace(r_pair.first, MatchType::value_type::second_type());

                // Collect matches for the current placeholder
                for (auto i_group : r_pair.second.value()) {
                    // First match (index 0) is irrelevant because it's the entire pattern,
                    // the rest is offset by 1
                    const auto i_group_match = i_group + 1;

                    if (!results.str(i_group_match).empty()) {
                        emplace_result.first->second.push_back(results.str(i_group_match));
                    }
                } // for i_group
            } // for placeholder, group_indices
        } // if group_indices is not None
    } /*if regex_match*/ else {
        KRATOS_ERROR << "'" << rString << "' is not a match for '" << this->GetRegexString() << "'";
    }

    return output;

    KRATOS_CATCH("");
}


std::string PlaceholderPattern::Apply(const PlaceholderMap& rPlaceholderValueMap) const
{
    KRATOS_TRY

    auto output = mPattern;
    const auto it_group_map_end = mPlaceholderGroupMap.end();

    for (const auto& r_pair : rPlaceholderValueMap) {
        auto it_pair = mPlaceholderGroupMap.find(r_pair.first);
        if (it_pair != it_group_map_end) {
            if (it_pair->second.has_value()) {
                while (true) {
                    auto position = output.find(r_pair.first);
                    if (position != output.npos) {
                        output.replace(position, r_pair.first.size(), r_pair.second);
                    } else {
                        break;
                    }
                } // while placeholder in output
            } // if placeholder in pattern
        } else {
            KRATOS_ERROR << r_pair.first << " is not a registered placeholder in " << mPattern;
        } // unrecognized placeholder
    } // for placeholder, value in map

    return output;

    KRATOS_CATCH("");
}


bool PlaceholderPattern::IsConst() const
{
    return std::none_of(mPlaceholderGroupMap.begin(),
                        mPlaceholderGroupMap.end(),
                        [](const auto& rPair) {
                            return rPair.second.has_value();
                        });
}


const std::regex& PlaceholderPattern::GetRegex() const
{
    return mRegex;
}


const std::string& PlaceholderPattern::GetRegexString() const
{
    return mRegexString;
}


const std::string& PlaceholderPattern::GetPatternString() const
{
    return mPattern;
}


ModelPartPattern::PlaceholderMap ModelPartPattern::GetPlaceholderMap()
{
    return PlaceholderMap {
        {"<model_part_name>", ".+"},
        {"<step>", RegexUtility::UnsignedInteger().first},
        {"<time>", RegexUtility::FloatingPoint().first},
        {"<rank>", RegexUtility::Integer().first}
    };
}


std::string PlaceholderPattern::FormatRegexLiteral(const std::string& rLiteral)
{
    KRATOS_TRY

    auto output = rLiteral;

    for (char char_to_escape : R"(!$()*+-?[\]^)") {
        std::size_t position = 0;

        std::string escaped;
        escaped.reserve(2);
        escaped.push_back('\\');
        escaped.push_back(char_to_escape);

        while (true) {
            if (output.size() <= position) break;

            position = output.find(char_to_escape, position);
            if (position == output.npos) {
                break;
            } else {
                // Escape the sensitive character
                if (!position || output[position-1] != '\\') {
                    output.replace(position, 1, escaped);
                    position += 2;
                }
            } // if char_to_escape in output
        } // while True
    } // for char_to_escape

    return output;

    KRATOS_CATCH("");
}


ModelPartPattern::ModelPartPattern(const std::string& rPattern)
    : PlaceholderPattern(rPattern, ModelPartPattern::GetPlaceholderMap())
{
}


std::string ModelPartPattern::Apply(const ModelPart& rModelPart) const
{
    KRATOS_TRY
    ModelPartPattern::PlaceholderMap map;
    this->PopulatePlaceholderMap(map, rModelPart);
    return this->Apply(map);
    KRATOS_CATCH("");
}


void ModelPartPattern::PopulatePlaceholderMap(PlaceholderMap& rMap, const ModelPart& rModelPart) const
{
    // TODO: implement formatting, see the documentation in the header. @matekelemen
    const auto& r_pattern = this->GetPatternString();

    if (r_pattern.find("<model_part_name>") != r_pattern.npos) {
        rMap.emplace("<model_part_name>", rModelPart.Name());
    }

    if (r_pattern.find("<step>") != r_pattern.npos) {
        rMap.emplace("<step>", std::to_string(rModelPart.GetProcessInfo().GetValue(STEP)));
    }

    if (r_pattern.find("<time>") != r_pattern.npos) {
        // Hardcoded formatting - to be changed later
        std::stringstream stream;
        stream << std::scientific << std::setprecision(4) << rModelPart.GetProcessInfo().GetValue(TIME);
        rMap.emplace("<time>", stream.str());
    }

    if (r_pattern.find("<rank>") != r_pattern.npos) {
        rMap.emplace("<rank>", std::to_string(rModelPart.GetCommunicator().MyPID()));
    }
}


ModelPartPattern::ModelPartPattern(const std::string& rPattern, const PlaceholderMap& rPlaceholderMap)
    : PlaceholderPattern(rPattern, rPlaceholderMap)
{
}


} // namespace Kratos
