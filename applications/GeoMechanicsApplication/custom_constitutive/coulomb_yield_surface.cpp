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
#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <cmath>

namespace
{

using namespace Kratos;

CoulombYieldSurface::KappaDependentFunction MakeConstantFunction(double Value)
{
    return [Value](double /* unused kappa */) { return Value; };
}

CoulombYieldSurface::KappaDependentFunction MakeKappaDependentFunctionForFrictionAngle(const Properties& rMaterialProperties)
{
    const auto hardening_type = rMaterialProperties[GEO_COULOMB_HARDENING_TYPE];

    if (hardening_type == "None") {
        return MakeConstantFunction(ConstitutiveLawUtilities::GetFrictionAngleInRadians(rMaterialProperties));
    }

    KRATOS_ERROR << "Cannot create a kappa-dependent function for the friction angle: unknown "
                    "hardening type '"
                 << hardening_type << "'\n";
}

} // namespace

namespace Kratos
{

CoulombYieldSurface::CoulombYieldSurface()
{
    mMaterialProperties[GEO_COULOMB_HARDENING_TYPE] = "None";
    mMaterialProperties[GEO_FRICTION_ANGLE]         = 0.0;
    mMaterialProperties[GEO_COHESION]               = 0.0;
    mMaterialProperties[GEO_DILATANCY_ANGLE]        = 0.0;

    InitializeKappaDependentFunctions();
}

CoulombYieldSurface::CoulombYieldSurface(Properties MaterialProperties)
    : mMaterialProperties{std::move(MaterialProperties)}
{
    // For backward compatibility, if no hardening type is given, we assume no hardening at all
    if (!mMaterialProperties.Has(GEO_COULOMB_HARDENING_TYPE)) {
        mMaterialProperties[GEO_COULOMB_HARDENING_TYPE] = "None";
    }

    InitializeKappaDependentFunctions();
}

double CoulombYieldSurface::GetFrictionAngleInRadians() const
{
    return mFrictionAngleCalculator(mKappa);
}

double CoulombYieldSurface::GetCohesion() const { return mCohesionCalculator(mKappa); }

double CoulombYieldSurface::GetDilatancyAngleInRadians() const
{
    return mDilatancyAngleCalculator(mKappa);
}

double CoulombYieldSurface::YieldFunctionValue(const Vector& rSigmaTau) const
{
    return rSigmaTau[1] + rSigmaTau[0] * std::sin(GetFrictionAngleInRadians()) -
           GetCohesion() * std::cos(GetFrictionAngleInRadians());
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector& rSigmaTau) const
{
    return DerivativeOfFlowFunction(rSigmaTau, CoulombAveragingType::NO_AVERAGING);
}

Vector CoulombYieldSurface::DerivativeOfFlowFunction(const Vector&, CoulombAveragingType AveragingType) const
{
    const auto sin_psi = std::sin(GetDilatancyAngleInRadians());
    Vector     result(2);
    switch (AveragingType) {
        using enum CoulombAveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        result <<= -(1.0 - 3.0 * sin_psi) / 4.0, (3.0 - sin_psi) / 4.0;
        break;
    case NO_AVERAGING:
        result <<= sin_psi, 1.0;
        break;
    case HIGHEST_PRINCIPAL_STRESSES:
        result <<= (1.0 + 3.0 * sin_psi) / 4.0, (3.0 + sin_psi) / 4.0;
        break;
    default:
        KRATOS_ERROR << "Unsupported Averaging Type: " << static_cast<std::size_t>(AveragingType) << "\n";
    }
    return result;
}

void CoulombYieldSurface::InitializeKappaDependentFunctions()
{
    // At present, we only support properties that are independent of kappa
    mFrictionAngleCalculator = MakeKappaDependentFunctionForFrictionAngle(mMaterialProperties);
    mCohesionCalculator = MakeConstantFunction(ConstitutiveLawUtilities::GetCohesion(mMaterialProperties));
    mDilatancyAngleCalculator =
        MakeConstantFunction(MathUtils<>::DegreesToRadians(mMaterialProperties[GEO_DILATANCY_ANGLE]));
}

void CoulombYieldSurface::save(Serializer& rSerializer) const
{
    rSerializer.save("Kappa", mKappa);
    rSerializer.save("MaterialProperties", mMaterialProperties);

    // Members `mFrictionAngleCalculator`, `mCohesionCalculator`, and `mDilatancyAngleCalculator`
    // will be reconstructed using data from `mMaterialProperties`. No additional serialization is
    // needed.
}

void CoulombYieldSurface::load(Serializer& rSerializer)
{
    rSerializer.load("Kappa", mKappa);
    rSerializer.load("MaterialProperties", mMaterialProperties);

    InitializeKappaDependentFunctions();
}

} // Namespace Kratos
