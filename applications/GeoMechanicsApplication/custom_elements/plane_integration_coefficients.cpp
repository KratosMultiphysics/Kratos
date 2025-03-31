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

#include "plane_integration_coefficients.h"

namespace Kratos
{

std::unique_ptr<IntegrationCoefficientsCalculator> PlaneIntegrationCoefficients::Clone() const
{
    return std::unique_ptr<PlaneIntegrationCoefficients>();
}

double PlaneIntegrationCoefficients::CalculateIntegrationCoefficient(
    const Geometry<Node>::IntegrationPointType& rIntegrationPoint, double DetJ, const Geometry<Node>&) const
{
    return rIntegrationPoint.Weight() * DetJ;
}

void PlaneIntegrationCoefficients::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void PlaneIntegrationCoefficients::load(Serializer&)
{
    // No data members to be loaded (yet)
}
} // namespace Kratos
