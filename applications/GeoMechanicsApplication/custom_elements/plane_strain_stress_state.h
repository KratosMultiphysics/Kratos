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
//

#pragma once

#include "stress_state_policy.h"

namespace Kratos
{

class PlaneStrainStressState : public StressStatePolicy
{
public:
    Matrix CalculateBMatrix(const Matrix& GradNpT, const Vector& Np, const Geometry<Node>& rGeometry) const override;
    double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                           double                detJ,
                                           const Geometry<Node>& rGeometry) const override;
    unique_ptr<StressStatePolicy> Clone() const override;
};

} // namespace Kratos
