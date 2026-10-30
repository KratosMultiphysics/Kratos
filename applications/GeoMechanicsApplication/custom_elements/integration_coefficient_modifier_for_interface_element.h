#pragma once

#include "integration_coefficients_calculator.hpp"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) IntegrationCoefficientModifierForInterfaceElement : public IntegrationCoefficientModifier
{
public:
    double operator()(double                           IntegrationCoefficient,
                      const Geo::IntegrationPointType& rIntegrationPoint,
                      const Element&                   rElement) const override;

    [[nodiscard]] std::unique_ptr<IntegrationCoefficientModifier> Clone() const override;
};

} // namespace Kratos
