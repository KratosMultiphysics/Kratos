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

#include "three_dimensional_integration_coefficients.h"

namespace Kratos
{

std::unique_ptr<IntegrationCoefficientsCalculator> ThreeDimensionalIntegrationCoefficients::Clone() const
{
    return std::unique_ptr<ThreeDimensionalIntegrationCoefficients>();
}

double ThreeDimensionalIntegrationCoefficients::CalculateIntegrationCoefficient(
    const Geometry<Node>::IntegrationPointType& rIntegrationPoint, double DetJ, const Geometry<Node>&) const
{
    return rIntegrationPoint.Weight() * DetJ;
}

void ThreeDimensionalIntegrationCoefficients::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void ThreeDimensionalIntegrationCoefficients::load(Serializer&)
{
    // No data members to be loaded (yet)
}

} // namespace Kratos
