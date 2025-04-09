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

#include "pw_line_integration_coefficients.h"

#include <includes/element.h>
#include <structural_mechanics_application_variables.h>

namespace Kratos
{

Vector PwLineIntegrationCoefficients::CalculateIntegrationCoefficients(const Geometry<Node>::IntegrationPointsArrayType& rIntegrationPoints,
                                                                       const Vector& rDetJs,
                                                                       double        CrossArea,
                                                                       std::size_t) const
{
    auto result = Vector{rIntegrationPoints.size()};
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), rDetJs.begin(),
                   result.begin(), [&CrossArea](const auto& rIntegrationPoint, const auto& rDetJ) {
        return rIntegrationPoint.Weight() * rDetJ * CrossArea;
    });
    return result;
}

std::unique_ptr<IntegrationCoefficientsCalculator> PwLineIntegrationCoefficients::Clone() const
{
    return std::make_unique<PwLineIntegrationCoefficients>();
}

void PwLineIntegrationCoefficients::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void PwLineIntegrationCoefficients::load(Serializer&)
{
    // No data members to be loaded (yet)
}

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
