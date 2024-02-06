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

static const SizeType voigt_size =
    VOIGT_SIZE_2D_PLANE_STRAIN; // For a axi-symmetric stress state, the voigt size is the same as for plane strain
static const SizeType dimension = 2; // Axi-symmetric stress state assumes a 2D problem definition

Matrix AxisymmetricStressStatePolicy::CalculateBMatrix(const Matrix& GradNpT, const Vector& Np, SizeType NumberOfNodes) const
{
    return ZeroMatrix(voigt_size, dimension * NumberOfNodes);
}

} // namespace Kratos
