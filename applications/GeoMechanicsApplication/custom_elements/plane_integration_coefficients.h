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

class KRATOS_API(GEO_MECHANICS_APPLICATION) PlaneIntegrationCoefficients : public IntegrationCoefficientsCalculator
{
public:
    [[nodiscard]] double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                         double DetJ,
                                                         const Geometry<Node>& rGeometry) const override;

private:
    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};

} // namespace Kratos
