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

#include "custom_constitutive/coulomb_yield_surface.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <cmath>

namespace Kratos
{

CoulombYieldSurface::CoulombYieldSurface(double FrictionAngleInRad, double Cohesion, double DilatationAngleInRad)
    : mFrictionAngle{FrictionAngleInRad}, mCohesion{Cohesion}, mDilatationAngle{DilatationAngleInRad}
{
}

double CoulombYieldSurface::YieldFunctionValue(const Vector& rPrincipalStress) const
{
    return 0.5 * (rPrincipalStress(0) - rPrincipalStress(2)) +
           0.5 * (rPrincipalStress(0) + rPrincipalStress(2)) * std::sin(mFrictionAngle) -
           mCohesion * std::cos(mFrictionAngle);
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector& rPrincipalStress) const
{
    Vector result(3);
    result <<= 0.5 * (1.0 + std::sin(mDilatationAngle)), 0.0, -0.5 * (1.0 - std::sin(mDilatationAngle));
    return result;
}

void CoulombYieldSurface::save(Serializer& rSerializer) const {}

void CoulombYieldSurface::load(Serializer& rSerializer) {}

} // Namespace Kratos
