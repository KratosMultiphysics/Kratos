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

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector& rSigmaTau) const
{
    return DerivativeOfFlowFunction(rSigmaTau, std::size_t{1});
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector&, std::size_t AveragingType) const
{
    Vector result(2);
    switch (AveragingType) {
    case 0:
        result <<= -(1.0 - 3.0 * std::sin(mDilatationAngle)) / 4.0, (3.0 - std::sin(mDilatationAngle)) / 4.0;
        break;
    case 1:
        result <<= std::sin(mDilatationAngle), 1.0;
        break;
    case 2:
        result <<= (1.0 + 3.0 * std::sin(mDilatationAngle)) / 4.0, (3.0 + std::sin(mDilatationAngle)) / 4.0;
        break;
    default:
        KRATOS_ERROR << "Unsupported Averaging Type: " << AveragingType << "\n";
    }
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
