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
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <cmath>

namespace Kratos
{

CoulombYieldSurface::CoulombYieldSurface(double FrictionAngleInRad, double Cohesion, double DilatationAngleInRad)
    : mFrictionAngle{FrictionAngleInRad}, mCohesion{Cohesion}, mDilatationAngle{DilatationAngleInRad}
{
}

double CoulombYieldSurface::YieldFunctionValue(const Vector& rSigmaTau) const
{
    return rSigmaTau[1] + rSigmaTau[0] * std::sin(mFrictionAngle) - mCohesion * std::cos(mFrictionAngle);
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector& rPrincipalStress) const
{
    Vector result(3);
    result <<= 0.5 * (1.0 + std::sin(mDilatationAngle)), 0.0, -0.5 * (1.0 - std::sin(mDilatationAngle));
    return result;
}

void CoulombYieldSurface::save(Serializer& rSerializer) const
{
    rSerializer.save("FrictionAngle", mFrictionAngle);
    rSerializer.save("Cohesion", mCohesion);
    rSerializer.save("DilatationAngle", mDilatationAngle);
}

void CoulombYieldSurface::load(Serializer& rSerializer)
{
    rSerializer.load("FrictionAngle", mFrictionAngle);
    rSerializer.load("Cohesion", mCohesion);
    rSerializer.load("DilatationAngle", mDilatationAngle);
}

} // Namespace Kratos
