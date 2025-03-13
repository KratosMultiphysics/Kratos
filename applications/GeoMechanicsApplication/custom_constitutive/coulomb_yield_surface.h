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
//  Main authors:    Mohamed Nabi,
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

    CoulombYieldSurface(double FrictionAngleInRad, double Cohesion, double DilatationAngleInRad);

    [[nodiscard]] double YieldFunctionValue(const Vector& rPrincipalStress) const override;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector& rPrincipalStress) const override;

private:
    double mFrictionAngle   = 0.0;
    double mCohesion        = 0.0;
    double mDilatationAngle = 0.0;

}; // Class CoulombYieldSurface

} // namespace Kratos