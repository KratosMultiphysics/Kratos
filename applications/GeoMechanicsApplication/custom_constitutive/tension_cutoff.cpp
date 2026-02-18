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
    return rSigmaTau.Sigma() + rSigmaTau.Tau() - mMaterialProperties[GEO_TENSILE_STRENGTH];
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

double TensionCutoff::CalculatePlasticMultiplier(const Geo::SigmaTau& rTrialSigmaTau,
                                                 const Vector&        rDerivativeOfFlowFunction,
                                                 const Matrix&        rElasticMatrix) const
{
    const auto numerator = rElasticMatrix(0, 0) * rDerivativeOfFlowFunction[0] +
                           rElasticMatrix(1, 1) * rDerivativeOfFlowFunction[1];
    return -YieldFunctionValue(rTrialSigmaTau) / numerator;
}

double TensionCutoff::CalculatePlasticMultiplier(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                 const Vector& rDerivativeOfFlowFunction,
                                                 const Matrix& rElasticMatrix) const
{
    const auto numerator = inner_prod(row(rElasticMatrix, 0), rDerivativeOfFlowFunction);
    return -YieldFunctionValue(rTrialPrincipalStresses) / numerator;
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
