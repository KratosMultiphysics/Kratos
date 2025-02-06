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

// System includes
#include <cmath>

// Project includes
#include "custom_constitutive/evaluate_yield_function.h"
#include "includes/serializer.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) CoulombYieldFunction : public EvaluateYieldFunction
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CoulombYieldFunction);

    CoulombYieldFunction() = default;

    CoulombYieldFunction(double frictionAngle, double cohesion);

    double operator()(const Vector& rPrincipalStress) const override;

private:
    // Member Variables
    double mFrictionAngle;
    double mCohesion;

}; // Class BilinearCohesive3DLaw

} // namespace Kratos