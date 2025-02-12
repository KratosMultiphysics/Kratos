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

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{
CoulombYieldFunction::CoulombYieldFunction(double frictionAngle, double cohesion, double dilatationAngle)
    : mFrictionAngle{frictionAngle}, mCohesion{cohesion}, mDilatationAngle{dilatationAngle}
{
}

double CoulombYieldFunction::CalculateYieldFunction(const Vector& rPrincipalStress) const
{
    return 0.5 * (rPrincipalStress(0) - rPrincipalStress(2)) +
           0.5 * (rPrincipalStress(0) + rPrincipalStress(2)) * std::sin(mFrictionAngle) -
           mCohesion * std::cos(mFrictionAngle);
}

Vector CoulombYieldFunction::CalculateFlowFunctionDerivate(const Vector& rPrincipalStress) const
{
    Vector result(3);
    result <<= 0.5 * (1.0 + std::sin(mDilatationAngle)), 0.0, -0.5 * (1.0 - std::sin(mDilatationAngle));
    return result;
}

} // Namespace Kratos
