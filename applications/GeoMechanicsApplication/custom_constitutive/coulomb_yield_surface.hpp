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
//  Main authors:    Mohame Nabi,
//                   Wijtze Pieter Kikstra
//

#pragma once

#include "custom_constitutive/yield_surface.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) CoulombYieldSurface : public YieldSurface
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CoulombYieldSurface);

    CoulombYieldSurface() = default;

    CoulombYieldSurface(double frictionAngle, double cohesion, double dilatationAngle);

    double YieldFunctionValue(const Vector& rPrincipalStress) const override;
    Vector DerivateOfFlowFunction(const Vector& rPrincipalStress) const override;

private:
    // Member Variables
    double mFrictionAngle   = 0.0;
    double mCohesion        = 0.0;
    double mDilatationAngle = 0.0;

}; // Class CoulombYieldFunction

} // namespace Kratos