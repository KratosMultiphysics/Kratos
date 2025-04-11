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
//                   Anne van de Graaf
//

#include "axisymmetric_integration_coefficients.h"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

double IntegrationCoefficientModifierForAxisymmetricElement::operator()(double IntegrationCoefficient,
                                                                        const Geo::IntegrationPointType& rIntegrationPoint,
                                                                        const Element& rElement) const
{
    auto shape_function_values = Vector{};
    shape_function_values      = rElement.GetGeometry().ShapeFunctionsValues(
        shape_function_values, rIntegrationPoint.Coordinates());
    return IntegrationCoefficient * GeoElementUtilities::CalculateAxisymmetricCircumference(
                                        shape_function_values, rElement.GetGeometry());
}

std::unique_ptr<IntegrationCoefficientModifier> IntegrationCoefficientModifierForAxisymmetricElement::Clone() const
{
    return std::make_unique<IntegrationCoefficientModifierForAxisymmetricElement>();
};

} // namespace Kratos
