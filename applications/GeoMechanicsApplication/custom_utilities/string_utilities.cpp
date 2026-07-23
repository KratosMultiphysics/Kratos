// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

#include "string_utilities.h"
#include <algorithm>
#include <cctype>
#include <sstream>

using namespace std::string_literals;

namespace Kratos
{

std::string GeoStringUtilities::ToLower(const std::string& rString)
{
    auto result = rString;
    std::ranges::transform(result, result.begin(), [](auto c) { return std::tolower(c); });
    return result;
}

std::string GeoStringUtilities::Join(const std::vector<std::string>& rStrings, const std::string& rSeparator)
{
    if (rStrings.empty()) return ""s;

    auto oss = std::ostringstream{};
    auto it  = rStrings.begin();
    oss << *it;
    ++it;
    for (; it != rStrings.end(); ++it) {
        oss << rSeparator << *it;
    }

    return oss.str();
}

} // namespace Kratos