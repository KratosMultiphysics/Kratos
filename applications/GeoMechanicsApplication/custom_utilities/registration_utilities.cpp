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
#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_constitutive/tension_cutoff.h"
#include "custom_elements/axisymmetric_stress_state.h"
#include "custom_elements/interface_stress_state.h"
#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/three_dimensional_stress_state.h"
#include "includes/serializer.h"

#include <string>

using namespace std::string_literals;

namespace Kratos
{

void RegistrationUtilities::RegisterStressStatePolicies()
{
    Serializer::Register("AxisymmetricStressState"s, AxisymmetricStressState{});
    Serializer::Register("InterfaceStressState"s, InterfaceStressState{});
    Serializer::Register("PlaneStrainStressState"s, PlaneStrainStressState{});
    Serializer::Register("ThreeDimensionalStressState"s, ThreeDimensionalStressState{});
}

void RegistrationUtilities::DeregisterStressStatePolicies()
{
    Serializer::Deregister("AxisymmetricStressState"s);
    Serializer::Deregister("InterfaceStressState"s);
    Serializer::Deregister("PlaneStrainStressState"s);
    Serializer::Deregister("ThreeDimensionalStressState"s);
}

void RegistrationUtilities::RegisterYieldSurfaces()
{
    Serializer::Register("CoulombYieldSurface"s, CoulombYieldSurface{});
    Serializer::Register("TensionCutoff"s, TensionCutoff{});
}

void RegistrationUtilities::DeregisterYieldSurfaces()
{
    Serializer::Deregister("CoulombYieldSurface"s);
    Serializer::Deregister("TensionCutoff"s);
}

} // namespace Kratos