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

#include "axisymmetric_integration_coefficients.h"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

std::unique_ptr<IntegrationCoefficientsCalculator> AxisymmetricIntegrationCoefficients::Clone() const
{
    return std::make_unique<AxisymmetricIntegrationCoefficients>();
}

double AxisymmetricIntegrationCoefficients::CalculateIntegrationCoefficient(
    const Geometry<Node>::IntegrationPointType& rIntegrationPoint, double DetJ, const Geometry<Node>& rGeometry) const
{
    Vector shape_function_values;
    shape_function_values =
        rGeometry.ShapeFunctionsValues(shape_function_values, rIntegrationPoint.Coordinates());

    const auto radius_weight =
        GeoElementUtilities::CalculateAxisymmetricCircumference(shape_function_values, rGeometry);

    return rIntegrationPoint.Weight() * DetJ * radius_weight;
}

void AxisymmetricIntegrationCoefficients::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void AxisymmetricIntegrationCoefficients::load(Serializer&)
{
    // No data members to be loaded (yet)
}

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

} // namespace Kratos
