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

#include "custom_constitutive/tension_cutoff.h"
#include "custom_constitutive/principal_stresses.hpp"
#include "custom_constitutive/sigma_tau.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{
TensionCutoff::TensionCutoff(const Properties& rMaterialProperties)
    : mMaterialProperties{rMaterialProperties}
{
}

double TensionCutoff::GetTensileStrength() const
{
    return mMaterialProperties[GEO_TENSILE_STRENGTH];
}

// At some point in time we would like to get rid of this API. For now, just forward the request.
double TensionCutoff::YieldFunctionValue(const Vector& rSigmaTau) const
{
    return YieldFunctionValue(Geo::SigmaTau{rSigmaTau});
}

double TensionCutoff::YieldFunctionValue(const Geo::SigmaTau& rSigmaTau) const
{
    // Note that we attempt to calculate the principal stress vector (which consists of three values) from the given
    // sigma and tau (two values). As a result, the second principal stress cannot be uniquely determined.
    const auto principal_stress_vector = Geo::PrincipalStresses{};
    return YieldFunctionValue(StressStrainUtilities::TransformSigmaTauToPrincipalStresses(
        rSigmaTau, principal_stress_vector));
}

double TensionCutoff::YieldFunctionValue(const Geo::PrincipalStresses& rPrincipalStresses) const
{
    return rPrincipalStresses.Values()[0] - mMaterialProperties[GEO_TENSILE_STRENGTH];
}

// At some point in time we would like to get rid of this API. For now, just forward the request.
Vector TensionCutoff::DerivativeOfFlowFunction(const Vector&) const
{
    const auto unused_sigma_tau = Geo::SigmaTau{};
    return DerivativeOfFlowFunction(unused_sigma_tau, YieldSurfaceAveragingType::NO_AVERAGING);
}

Vector TensionCutoff::DerivativeOfFlowFunction(const Geo::SigmaTau&, YieldSurfaceAveragingType AveragingType) const
{
    Vector result(2);
    switch (AveragingType) {
        using enum YieldSurfaceAveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        result <<= 0.5, 0.5;
        break;
    case NO_AVERAGING:
        result <<= 1.0, 1.0;
        break;
    case HIGHEST_PRINCIPAL_STRESSES:
        result <<= 1.0, 1.0;
        break;
    default:
        KRATOS_ERROR << "Unsupported Averaging Type: " << static_cast<std::size_t>(AveragingType) << "\n";
    }
    return result;
}

Vector TensionCutoff::DerivativeOfFlowFunction(const Geo::PrincipalStresses&,
                                               YieldSurfaceAveragingType AveragingType) const
{
    Vector result(3);
    switch (AveragingType) {
        using enum YieldSurfaceAveragingType;
    case LOWEST_PRINCIPAL_STRESSES:
        result <<= 0.5, 0.5, 0.0;
        break;
    case NO_AVERAGING:
        result <<= 1.0, 0.0, 0.0;
        break;
    case HIGHEST_PRINCIPAL_STRESSES:
        result <<= 1.0, 0.0, 0.0;
        break;
    default:
        KRATOS_ERROR << "Unsupported Averaging Type: " << static_cast<std::size_t>(AveragingType) << "\n";
    }
    return result;
}

double TensionCutoff::CalculatePlasticMultiplier(const Geo::SigmaTau& rSigmaTau,
                                                 const Vector& rDerivativeOfFlowFunction) const
{
    const auto numerator = mMaterialProperties[INTERFACE_NORMAL_STIFFNESS] * rDerivativeOfFlowFunction[0] +
                           mMaterialProperties[INTERFACE_SHEAR_STIFFNESS] * rDerivativeOfFlowFunction[1];
    return (mMaterialProperties[GEO_TENSILE_STRENGTH] - rSigmaTau.Sigma() - rSigmaTau.Tau()) / numerator;
}

double TensionCutoff::CalculatePlasticMultiplier(const Geo::PrincipalStresses& rPrincipalStresses,
                                                 const Vector& rDerivativeOfFlowFunction) const
{
    return (mMaterialProperties[GEO_TENSILE_STRENGTH] - rPrincipalStresses.Values()[0]) /
           rDerivativeOfFlowFunction[0];
}

Matrix TensionCutoff::CalculateElasticMatrix(const Geo::PrincipalStresses&) const
{
    const auto youngs_modulus = mMaterialProperties[YOUNG_MODULUS];
    const auto poissons_ratio = mMaterialProperties[POISSON_RATIO];
    const auto c0     = youngs_modulus / ((1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio));
    Matrix     result = ZeroMatrix(3, 3);
    result(0, 0) = result(1, 1) = result(2, 2) = c0 * (1.0 - poissons_ratio);
    result(0, 1) = result(0, 2) = result(1, 0) = result(1, 2) = result(2, 0) = result(2, 1) = c0 * poissons_ratio;
    return result;
}

void TensionCutoff::save(Serializer& rSerializer) const
{
    rSerializer.save("MaterialProperties", mMaterialProperties);
}

void TensionCutoff::load(Serializer& rSerializer)
{
    rSerializer.load("MaterialProperties", mMaterialProperties);
}

} // Namespace Kratos
