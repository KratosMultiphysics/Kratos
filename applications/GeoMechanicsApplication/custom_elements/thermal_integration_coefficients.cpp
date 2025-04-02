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

#include "thermal_integration_coefficients.h"

namespace Kratos
{

Vector ThermalIntegrationCoefficients::CalculateIntegrationCoefficients(const Geometry<Node>::IntegrationPointsArrayType& rIntegrationPoints,
                                                                        const Vector& rDetJs,
                                                                        double        CrossArea,
                                                                        std::size_t LocalDimension) const
{
    auto result = Vector{rIntegrationPoints.size()};
    if (LocalDimension == 1) {
        std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), rDetJs.begin(),
                       result.begin(), [&CrossArea](const auto& rIntegrationPoint, const auto& rDetJ) {
            return rIntegrationPoint.Weight() * rDetJ * CrossArea;
        });
    } else {
        std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), rDetJs.begin(),
                       result.begin(), [](const auto& rIntegrationPoint, const auto& rDetJ) {
            return rIntegrationPoint.Weight() * rDetJ;
        });
    }
    return result;
}

std::unique_ptr<IntegrationCoefficientsCalculator> ThermalIntegrationCoefficients::Clone() const
{
    return std::make_unique<ThermalIntegrationCoefficients>();
}

void ThermalIntegrationCoefficients::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void ThermalIntegrationCoefficients::load(Serializer&)
{
    // No data members to be loaded (yet)
}

} // namespace Kratos
