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
#include "custom_utilities/stress_strain_utilities.h"
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{
TensionCutoff::TensionCutoff(double TensileStrength) : mTensileStrength{TensileStrength} {}

double TensionCutoff::GetTensileStrength() const { return mTensileStrength; }

double TensionCutoff::YieldFunctionValue(const Vector& rSigmaTau) const
{
    auto principal_stress_vector = Vector{ZeroVector{3}};
    principal_stress_vector =
        StressStrainUtilities::TransformSigmaTauToPrincipalStresses(rSigmaTau, principal_stress_vector);
    return principal_stress_vector[0] - mTensileStrength;
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
