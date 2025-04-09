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

#pragma once

#include "integration_coefficients_calculator.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) AxisymmetricIntegrationCoefficients : public IntegrationCoefficientsCalculator
{
public:
    std::unique_ptr<IntegrationCoefficientsCalculator> Clone() const override;

private:
    [[nodiscard]] double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                         double DetJ,
                                                         const Geometry<Node>& rGeometry) const override;

    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};

class KRATOS_API(GEO_MECHANICS_APPLICATION) IntegrationCoefficientModifierForAxisymmetricElement
    : public IntegrationCoefficientModifier
{
public:
    double operator()(double                           IntegrationCoefficient,
                      const Geo::IntegrationPointType& rIntegrationPoint,
                      const Element&                   rElement) const override;
    [[nodiscard]] std::unique_ptr<IntegrationCoefficientModifier> Clone() const override;
};

} // namespace Kratos
