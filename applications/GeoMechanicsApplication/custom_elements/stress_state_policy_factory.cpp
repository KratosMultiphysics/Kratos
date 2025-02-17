// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//
#include "stress_state_policy_factory.h"
#include "plane_strain_stress_state.h"

namespace Kratos
{

std::unique_ptr<StressStatePolicy> StressStatePolicyFactory::Create() const
{
    return std::make_unique<PlaneStrainStressState>();
}

} // namespace Kratos
