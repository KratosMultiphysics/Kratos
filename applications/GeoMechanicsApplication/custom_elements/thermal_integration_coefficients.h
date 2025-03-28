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
                                                          std::size_t   LocalDimension,
                                                          double        CrossArea) const override;

private:
    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};

} // namespace Kratos
