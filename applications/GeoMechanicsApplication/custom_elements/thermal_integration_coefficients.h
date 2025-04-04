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

class KRATOS_API(GEO_MECHANICS_APPLICATION) ThermalIntegrationCoefficients : public IntegrationCoefficientsCalculator
{
public:
    [[nodiscard]] Vector CalculateIntegrationCoefficients(const Geometry<Node>::IntegrationPointsArrayType& rIntegrationPoints,
                                                          const Vector& rDetJs,
                                                          double        CrossArea,
                                                          std::size_t LocalDimension = 0) const override;

    [[nodiscard]] std::vector<double> CalculateIntegrationCoefficients(const Geometry<Node>::IntegrationPointsArrayType&,
                                                                       const Vector&,
                                                                       const Geometry<Node>&) const override
    {
        KRATOS_ERROR << "PwLineIntegrationCoefficients::CalculateIntegrationCoefficients is called."
                     << std::endl;
    };

    std::unique_ptr<IntegrationCoefficientsCalculator> Clone() const override;

private:
    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};

class KRATOS_API(GEO_MECHANICS_APPLICATION) IntegrationCoefficientModifierForThermalElement : public IntegrationCoefficientModifier
{
public:
    double operator()(double                           IntegrationCoefficient,
                      const Geo::IntegrationPointType& rIntegrationPoint,
                      const Element&                   rElement) const override;
};

} // namespace Kratos
