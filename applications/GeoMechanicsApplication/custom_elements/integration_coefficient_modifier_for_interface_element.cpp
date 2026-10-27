#include "integration_coefficient_modifier_for_interface_element.h"
#include "geo_mechanics_application_variables.h"
#include "includes/element.h"

namespace Kratos
{

double IntegrationCoefficientModifierForInterfaceElement::operator()(
    double                           IntegrationCoefficient,
    const Geo::IntegrationPointType& rIntegrationPoint,
    const Element&                   rElement) const
{
    auto& properties = rElement.GetProperties();
    
    if (!properties.Has(GEO_CONTACT_SURFACE_PERIMETER)) {
        return IntegrationCoefficient;
    }
    
    double perimeter = properties[GEO_CONTACT_SURFACE_PERIMETER];
    return IntegrationCoefficient * perimeter;
}

std::unique_ptr<IntegrationCoefficientModifier> IntegrationCoefficientModifierForInterfaceElement::Clone() const
{
    return std::make_unique<IntegrationCoefficientModifierForInterfaceElement>();
}

} // namespace Kratos
