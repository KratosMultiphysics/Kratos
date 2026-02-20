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
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{
TensionCutoff::TensionCutoff(double TensileStrength) : mTensileStrength{TensileStrength} {}

double TensionCutoff::GetTensileStrength() const { return mTensileStrength; }

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
    return rPrincipalStresses.Values()[0] - mTensileStrength;
}

// At some point in time we would like to get rid of this API. For now, just forward the request.
Vector TensionCutoff::DerivativeOfFlowFunction(const Vector&) const
{
    const auto unused_sigma_tau = Geo::SigmaTau{};
    return DerivativeOfFlowFunction(unused_sigma_tau);
}

Vector TensionCutoff::DerivativeOfFlowFunction(const Geo::SigmaTau&) const
{
    return Vector{ScalarVector{2, 1.0}};
}

void TensionCutoff::save(Serializer& rSerializer) const
{
    rSerializer.save("TensileStrength", mTensileStrength);
}

void TensionCutoff::load(Serializer& rSerializer)
{
    rSerializer.load("TensileStrength", mTensileStrength);
}

} // Namespace Kratos
