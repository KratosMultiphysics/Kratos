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

#include "axisymmetric_stress_state_policy.h"

namespace Kratos
{

Matrix AxisymmetricStressStatePolicy::CalculateBMatrix(const Matrix& GradNpT, const Vector& Np) const
{
    return ZeroMatrix(0, 0);
}

} // namespace Kratos
