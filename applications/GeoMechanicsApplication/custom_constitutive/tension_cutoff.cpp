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

#include "custom_constitutive/tension_cutoff.hpp"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{
TensionCutoff::TensionCutoff(double tensileStrength)
    : mTensileStrength{tensileStrength}
{
}

double TensionCutoff::YieldFunctionValue(const Vector& rPrincipalStress) const
{
    return mTensileStrength - rPrincipalStress(2);
}

Vector TensionCutoff::DerivateOfFlowFunction(const Vector& rPrincipalStress) const
{
    Vector result(3);
    result <<= 0.0, 0.0, -1.0;
    return result;
}

} // Namespace Kratos
