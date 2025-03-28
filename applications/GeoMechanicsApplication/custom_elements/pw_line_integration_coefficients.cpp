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

#include "pw_line_integration_coefficients.h"

namespace Kratos
{

Vector PwLineIntegrationCoefficients::CalculateIntegrationCoefficients(const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
                                                                       const Vector& rDetJs,
                                                                       std::size_t,
                                                                       double CrossArea) const
{
    auto result = Vector{rIntegrationPoints.size()};
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), rDetJs.begin(),
                   result.begin(), [&CrossArea](const auto& rIntegrationPoint, const auto& rDetJ) {
        return rIntegrationPoint.Weight() * rDetJ * CrossArea;
    });
    return result;
}

void PwLineIntegrationCoefficients::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void PwLineIntegrationCoefficients::load(Serializer&)
{
    // No data members to be loaded (yet)
}

} // namespace Kratos
