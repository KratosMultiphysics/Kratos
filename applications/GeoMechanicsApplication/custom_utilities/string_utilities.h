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

#include "includes/kratos_parameters.h"
#include <string>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoStringUtilities
{
public:
    static std::string ToLower(const std::string& rString);
};

} // namespace Kratos