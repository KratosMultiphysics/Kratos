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

#include "integration_coefficient_modifier_for_interface_element.h"
#include "geo_mechanics_application_variables.h"
#include "includes/element.h"

namespace Kratos
{

double IntegrationCoefficientModifierForInterfaceElement::operator()(double IntegrationCoefficient,
                                                                     const Geo::IntegrationPointType&,
                                                                     const Element& rElement) const
{
    auto& r_properties = rElement.GetProperties();
    if (r_properties.Has(GEO_CONTACT_SURFACE_PERIMETER)) {
        return IntegrationCoefficient * r_properties[GEO_CONTACT_SURFACE_PERIMETER];
    }

    return IntegrationCoefficient;
}

std::unique_ptr<IntegrationCoefficientModifier> IntegrationCoefficientModifierForInterfaceElement::Clone() const
{
    return std::make_unique<IntegrationCoefficientModifierForInterfaceElement>();
}

} // namespace Kratos
