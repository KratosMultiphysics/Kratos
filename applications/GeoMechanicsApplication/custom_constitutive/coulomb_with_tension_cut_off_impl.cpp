// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohamed Nabi
//                   Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "custom_constitutive/coulomb_with_tension_cut_off_impl.h"

namespace Kratos
{
CoulombWithTensionCutOffImpl::CoulombWithTensionCutOffImpl(double FrictionAngleInRad,
                                                           double Cohesion,
                                                           double DilatationAngleInRad,
                                                           double TensileStrength)
    : mCoulombYieldSurface{FrictionAngleInRad, Cohesion, DilatationAngleInRad}, mTensionCutOff{TensileStrength}
{
}

} // namespace Kratos
