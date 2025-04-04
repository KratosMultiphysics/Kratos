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

#include "interface_integration_coefficients.h"

namespace Kratos
{

std::unique_ptr<IntegrationCoefficientsCalculator> InterfaceIntegrationCoefficients::Clone() const
{
    return std::make_unique<InterfaceIntegrationCoefficients>();
}

double InterfaceIntegrationCoefficients::CalculateIntegrationCoefficient(
    const Geometry<Node>::IntegrationPointType& rIntegrationPoint, double DetJ, const Geometry<Node>&) const
{
    return rIntegrationPoint.Weight() * DetJ;
}

void InterfaceIntegrationCoefficients::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void InterfaceIntegrationCoefficients::load(Serializer&)
{
    // No data members to be loaded (yet)
}

} // namespace Kratos
