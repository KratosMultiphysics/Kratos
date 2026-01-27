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
#include "custom_constitutive/geo_principal_stresses.hpp"
#include "custom_constitutive/geo_sigma_tau.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{
TensionCutoff::TensionCutoff(double TensileStrength) : mTensileStrength{TensileStrength} {}

double TensionCutoff::GetTensileStrength() const { return mTensileStrength; }

double TensionCutoff::YieldFunctionValue(const Vector& rSigmaTau) const
{
    return YieldFunctionValue(Geo::SigmaTau{rSigmaTau});
}

double TensionCutoff::YieldFunctionValue(const Geo::SigmaTau& rSigmaTau) const
{
    // Note that we attempt to calculate the principal stress vector (which consists of three values) from the given
    // sigma and tau (two values). As a result, the second principal stress cannot be uniquely determined.
    auto principal_stress_vector = Vector{ZeroVector{3}};
    return YieldFunctionValue(Geo::PrincipalStresses{StressStrainUtilities::TransformSigmaTauToPrincipalStresses(
        rSigmaTau.CopyTo<Vector>(), principal_stress_vector)});
}

double TensionCutoff::YieldFunctionValue(const Geo::PrincipalStresses& rPrincipalStresses) const
{
    return rPrincipalStresses.values[0] - mTensileStrength;
}

Vector TensionCutoff::DerivativeOfFlowFunction(const Vector&) const
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
