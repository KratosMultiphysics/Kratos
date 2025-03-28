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
    for (unsigned int integration_point_index = 0;
         integration_point_index < rIntegrationPoints.size(); ++integration_point_index) {
        result[integration_point_index] =
            rIntegrationPoints[integration_point_index].Weight() * rDetJs[integration_point_index];
        if (LocalDimension == 1) {
            result[integration_point_index] *= CrossArea;
        }
    }

    return result;
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
