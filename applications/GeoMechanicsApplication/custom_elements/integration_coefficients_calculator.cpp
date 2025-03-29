// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "integration_coefficients_calculator.h"

namespace Kratos
{

std::vector<double> IntegrationCoefficientsCalculator::CalculateIntegrationCoefficients(
    const Geometry<Node>::IntegrationPointsArrayType& rIntegrationPoints,
    const Vector&                                     rDetJs,
    const Geometry<Node>&                             rGeometry) const
{
    auto result = std::vector<double>{};
    result.reserve(rIntegrationPoints.size());
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), rDetJs.begin(),
                   std::back_inserter(result),
                   [this, &rGeometry](const auto& rIntegrationPoint, const auto& rDetJ) {
        return CalculateIntegrationCoefficient(rIntegrationPoint, rDetJ, rGeometry);
    });
    return result;
}

} // namespace Kratos