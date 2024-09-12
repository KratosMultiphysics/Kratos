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

#include "element_utilities.hpp"

namespace Kratos
{

std::vector<Vector> GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(
    const Geo::IntegrationPointVectorType& rIntegrationPoints, const Geometry<Node>& rGeometry)
{
    auto evaluate_shape_function_values = [&rGeometry](const auto& rIntegrationPoint) {
        auto result = Vector{};
        rGeometry.ShapeFunctionsValues(result, rIntegrationPoint);
        return result;
    };

    auto result = std::vector<Vector>{};
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), std::back_inserter(result),
                   evaluate_shape_function_values);

    return result;
}

} // namespace Kratos