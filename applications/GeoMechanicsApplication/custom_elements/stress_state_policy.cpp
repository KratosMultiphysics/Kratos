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

#include "stress_state_policy.h"

namespace Kratos
{

Vector StressStatePolicy::DefineVoigtVector(std::size_t Dimension)
{
    Vector VoigtVector = ZeroVector(GetVoigtSize(Dimension));
    std::fill_n(VoigtVector.begin(), GetStressTensorSize(Dimension), 1.0);
    return VoigtVector;
}

const Vector StressStatePolicy::VoigtVector2D = StressStatePolicy::DefineVoigtVector(2);
const Vector StressStatePolicy::VoigtVector3D = StressStatePolicy::DefineVoigtVector(3);

} // namespace Kratos