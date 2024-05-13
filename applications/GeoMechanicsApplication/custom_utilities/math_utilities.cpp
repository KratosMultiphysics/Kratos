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

#include "math_utilities.h"

std::vector<double> Kratos::GeoMechanicsMathUtilities::CalculateDeterminants(const std::vector<Matrix>& rMatrices)
{
    std::vector<double> result(rMatrices.size());
    std::transform(rMatrices.cbegin(), rMatrices.cend(), result.begin(),
                   [](const auto& rMatrix) { return MathUtils<>::Det(rMatrix); });

    return result;
}
