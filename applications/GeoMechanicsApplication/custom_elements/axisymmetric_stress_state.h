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
    double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                           double                detJ,
                                           const Geometry<Node>& rGeometry) const override;
    Matrix CalculateBMatrix(const Matrix& GradNpT, const Vector& Np, const Geometry<Node>& rGeometry) const override;
    unique_ptr<StressStatePolicy> Clone() const override;
    Vector CalculateGreenLagrangeStrain(const Matrix& rTotalDeformationGradient) const override;
};

} // namespace Kratos
