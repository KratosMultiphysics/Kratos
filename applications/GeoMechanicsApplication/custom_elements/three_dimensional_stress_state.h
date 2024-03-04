// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Marjan Fathian
//                   Richard Faasse
//

#pragma once

#include "stress_state_policy.h"

namespace Kratos
{

class ThreeDimensionalStressState : public StressStatePolicy
{
public:
    [[nodiscard]] Matrix CalculateBMatrix(const Matrix& rGradNpT,
                                          const Vector&,
                                          const Geometry<Node>& rGeometry) const override;
    [[nodiscard]] double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                         double DetJ,
                                                         const Geometry<Node>& rGeometry) const override;
    [[nodiscard]] Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const override;
    [[nodiscard]] std::unique_ptr<StressStatePolicy> Clone() const override;
};

} // namespace Kratos
