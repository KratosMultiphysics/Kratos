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

#pragma once

#include "coulomb_yield_surface.h"
#include "tension_cutoff.h"

namespace Kratos
{

class CoulombWithTensionCutOffImpl
{
public:
    CoulombWithTensionCutOffImpl() = default;
    CoulombWithTensionCutOffImpl(double FrictionAngleInRad, double Cohesion, double DilatationAngleInRad, double TensileStrength);

private:
    CoulombYieldSurface mCoulombYieldSurface;
    TensionCutoff       mTensionCutOff;
};

} // namespace Kratos
