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

#include <includes/element.h>
#include <structural_mechanics_application_variables.h>

namespace Kratos
{
double IntegrationCoefficientModifierForThermalElement::operator()(double IntegrationCoefficient,
                                                                   const Geo::IntegrationPointType& rIntegrationPoint,
                                                                   const Element& rElement) const
{
    return rElement.GetGeometry().LocalSpaceDimension() == 1
               ? IntegrationCoefficient * rElement.GetProperties()[CROSS_AREA]
               : IntegrationCoefficient;
}

std::unique_ptr<IntegrationCoefficientModifier> IntegrationCoefficientModifierForThermalElement::Clone() const
{
    return std::make_unique<IntegrationCoefficientModifierForThermalElement>();
}
} // namespace Kratos
