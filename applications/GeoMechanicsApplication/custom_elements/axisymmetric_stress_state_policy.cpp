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

static const SizeType NumberOfUDofPerNode = 2; // Axi-symmetric stress state assumes a 2D problem definition

Matrix AxisymmetricStressStatePolicy::CalculateBMatrix(const Matrix& GradNpT, const Vector& Np, SizeType NumberOfNodes) const
{
    return ZeroMatrix(VOIGT_SIZE_2D_AXISYMMETRIC, NumberOfUDofPerNode * NumberOfNodes);
}

} // namespace Kratos
