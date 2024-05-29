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
    [[nodiscard]] const Vector&                      GetVoigtVector() const override;
    [[nodiscard]] SizeType                           GetVoigtSize() const override;
    [[nodiscard]] SizeType                           GetStressTensorSize() const override;
};

} // namespace Kratos
