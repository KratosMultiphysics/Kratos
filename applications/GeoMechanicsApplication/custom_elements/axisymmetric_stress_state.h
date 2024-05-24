// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Marjan Fathian
//

#pragma once

#include "stress_state_policy.h"

namespace Kratos
{

class AxisymmetricStressState : public StressStatePolicy
{
public:
    [[nodiscard]] double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                         double DetJ,
                                                         const Geometry<Node>& rGeometry) const override;
    [[nodiscard]] Matrix CalculateBMatrix(const Matrix&         rDN_DX,
                                          const Vector&         rN,
                                          const Geometry<Node>& rGeometry) const override;
    [[nodiscard]] Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const override;
    [[nodiscard]] std::unique_ptr<StressStatePolicy> Clone() const override;
    [[nodiscard]] const Vector&                      GetVoigtVector() override;
    [[nodiscard]] const SizeType                     GetVoigtSize() override;
};

} // namespace Kratos
