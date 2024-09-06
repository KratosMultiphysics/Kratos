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

#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

#include <algorithm>
#include <vector>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoMechanicsMathUtilities
{
public:
    [[nodiscard]] static std::vector<double> CalculateDeterminants(const std::vector<Matrix>& rMatrices);

    template <typename VectorType>
    [[nodiscard]] static VectorType Normalized(const VectorType& rVector)
    {
        KRATOS_ERROR_IF(std::none_of(rVector.begin(), rVector.end(), [](auto component) {
            return component > 0.0;
        })) << "A zero vector cannot be normalized\n";

        return rVector / norm_2(rVector);
    }

}; // class GeoMechanicsMathUtilities

} // namespace Kratos
