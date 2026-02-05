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

namespace Kratos
{

std::string GeoStringUtilities::ToLower(const std::string& rString)
{
    auto result = rString;
    std::ranges::transform(result, result.begin(), [](auto c) { return std::tolower(c); });
    return result;
}

} // namespace Kratos