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

class AxisymmetricStressStatePolicy : public StressStatePolicy
{
public:
    Matrix CalculateBMatrix(const Matrix& GradNpT, const Vector& Np, const Geometry<Node>& rGeometry) const override;
};

} // namespace Kratos
