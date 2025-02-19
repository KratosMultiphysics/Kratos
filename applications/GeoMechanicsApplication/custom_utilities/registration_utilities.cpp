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
#include "registration_utilities.h"
#include "custom_elements/axisymmetric_stress_state.h"
#include "custom_elements/interface_stress_state.h"
#include "custom_elements/plane_strain_stress_state.h"
#include "includes/serializer.h"

namespace Kratos
{

void RegistrationUtilities::RegisterStressStatePolicies()
{
    Serializer::Register("AxisymmetricStressState", AxisymmetricStressState{});
    Serializer::Register("InterfaceStressState", InterfaceStressState{});
    Serializer::Register("PlaneStrainStressState", PlaneStrainStressState{});
}

} // namespace Kratos