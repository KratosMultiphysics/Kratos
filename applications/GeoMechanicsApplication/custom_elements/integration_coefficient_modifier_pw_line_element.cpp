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

#include "integration_coefficient_modifier_pw_line_element.h"

#include <includes/element.h>
#include <structural_mechanics_application_variables.h>

namespace Kratos
{

double IntegrationCoefficientModifierForPwLineElement::operator()(double IntegrationCoefficient,
                                                                  const Geo::IntegrationPointType& rIntegrationPoint,
                                                                  const Element& rElement) const
{
    return IntegrationCoefficient * rElement.GetProperties()[CROSS_AREA];
}

std::unique_ptr<IntegrationCoefficientModifier> IntegrationCoefficientModifierForPwLineElement::Clone() const
{
    return std::make_unique<IntegrationCoefficientModifierForPwLineElement>();
};

} // namespace Kratos
