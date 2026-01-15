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

namespace Kratos
{

std::string GeoStringUtilities::ToLower(const std::string& rString)
{
    auto result = rString;
    std::ranges::transform(result, result.begin(), [](unsigned char c) { return std::tolower(c); });
    return result;
}

} // namespace Kratos