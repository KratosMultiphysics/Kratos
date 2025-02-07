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
// //

#pragma once

#include "custom_constitutive/evaluate_yield_function.h"


namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) CoulombYieldFunction : public EvaluateYieldFunction
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CoulombYieldFunction);

    CoulombYieldFunction() = default;

    CoulombYieldFunction(double frictionAngle, double cohesion, double dilatationAngle);

    double CalculateYieldFunction(const Vector& rPrincipalStress) const override;
    Vector CalculateFlowFunctionDerivate(const Vector& rPrincipalStress) const override;

private:
    // Member Variables
    double mFrictionAngle = 0.0;
    double mDilatationAngle = 0.0;
    double mCohesion = 0.0;

}; // Class CoulombYieldFunction

} // namespace Kratos