// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include <vector>

#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoMechanicsMathUtilities
{
public:
    [[nodiscard]] static std::vector<double> CalculateDeterminants(const std::vector<Matrix>& rMatrices);

}; // class GeoMechanicsMathUtilities

} // namespace Kratos
