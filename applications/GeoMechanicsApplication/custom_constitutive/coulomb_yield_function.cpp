// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Wijtze Pieter Kikstra
//

#include "custom_constitutive/coulomb_yield_function.hpp"

namespace Kratos
{
    CoulombYieldFunction::CoulombYieldFunction(double frictionAngle, double cohesion)
        : mFrictionAngle{frictionAngle}, mCohesion{cohesion}
    {
    }

    double CoulombYieldFunction::operator()(const Vector& rPrincipalStress) const
    {
        return 0.5 * (rPrincipalStress(0) - rPrincipalStress(2)) +
               0.5 * (rPrincipalStress(0) + rPrincipalStress(2)) * std::sin(mFrictionAngle) -
               mCohesion * std::cos(mFrictionAngle);
    }

} // Namespace Kratos
