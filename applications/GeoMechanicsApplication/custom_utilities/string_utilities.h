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

#pragma once

#include "includes/kratos_export_api.h"
#include <string>
#include <vector>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoStringUtilities
{
public:
    static std::string ToLower(const std::string& rString);
    static std::string Join(const std::vector<std::string>& rStrings, const std::string& rSeparator);
};

} // namespace Kratos