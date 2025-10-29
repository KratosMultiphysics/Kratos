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

#include "utilities/math_utils.h"
#include <boost/numeric/ublas/assignment.hpp>
#include <cmath>

namespace Kratos
{

CoulombYieldSurface::CoulombYieldSurface(double FrictionAngleInRad, double Cohesion, double DilatationAngleInRad)
    : mFrictionAngle{FrictionAngleInRad}, mCohesion{Cohesion}, mDilatationAngle{DilatationAngleInRad}
{
}

double CoulombYieldSurface::GetFrictionAngleInRadians() const { return mFrictionAngle; }

double CoulombYieldSurface::GetCohesion() const { return mCohesion; }

double CoulombYieldSurface::GetDilatationAngleInRadians() const { return mDilatationAngle; }

double CoulombYieldSurface::YieldFunctionValue(const Vector& rSigmaTau) const
{
    return rSigmaTau[1] + rSigmaTau[0] * std::sin(mFrictionAngle) - mCohesion * std::cos(mFrictionAngle);
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector& rSigmaTau) const
{
    return DerivativeOfFlowFunction(rSigmaTau, CoulombAveragingType::NO_AVERAGING);
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector&, CoulombAveragingType AveragingType) const
{
    Vector result(2);
    switch (AveragingType) {
        using enum CoulombAveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        result <<= -(1.0 - 3.0 * std::sin(mDilatationAngle)) / 4.0, (3.0 - std::sin(mDilatationAngle)) / 4.0;
        break;
    case NO_AVERAGING:
        result <<= std::sin(mDilatationAngle), 1.0;
        break;
    case HIGHEST_PRINCIPAL_STRESSES:
        result <<= (1.0 + 3.0 * std::sin(mDilatationAngle)) / 4.0, (3.0 + std::sin(mDilatationAngle)) / 4.0;
        break;
    default:
        KRATOS_ERROR << "Unsupported Averaging Type: " << static_cast<std::size_t>(AveragingType) << "\n";
    }
    return result;
}

void CoulombYieldSurface::UpdateSurfaceProperties(double InitialFrictionAngle,
                                                  double FrictionAngleStrengthFactor,
                                                  double InitialCohesion,
                                                  double CohesionStrengthFactor,
                                                  double InitialDilatancyAngle,
                                                  double DilatancyAngleStrengthFactor,
                                                  double kappa)
{
    mFrictionAngle =
        this->CalculateUpdatedFrictionAngle(InitialFrictionAngle, FrictionAngleStrengthFactor, kappa);
    mCohesion = this->CalculateUpdatedCohesion(InitialCohesion, CohesionStrengthFactor, kappa);
    mDilatationAngle =
        this->CalculateUpdatedDilatancyAngle(InitialDilatancyAngle, DilatancyAngleStrengthFactor, kappa);
}

double CoulombYieldSurface::CalculateUpdatedFrictionAngle(double InitialFrictionAngle,
                                                          double StrengthFactor,
                                                          double kappa) const
{
    return MathUtils<>::DegreesToRadians(InitialFrictionAngle + StrengthFactor * kappa);
}

double CoulombYieldSurface::CalculateUpdatedCohesion(double InitialCohesion, double StrengthFactor, double kappa) const
{
    return InitialCohesion + StrengthFactor * kappa;
}

double CoulombYieldSurface::CalculateUpdatedDilatancyAngle(double InitialDilatancyAngle,
                                                           double StrengthFactor,
                                                           double kappa) const
{
    return MathUtils<>::DegreesToRadians(InitialDilatancyAngle + StrengthFactor * kappa);
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
