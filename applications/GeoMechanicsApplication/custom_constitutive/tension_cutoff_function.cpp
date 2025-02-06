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

#include "custom_constitutive/tension_cutoff_function.hpp"

namespace Kratos
{
    TensionCutoffFunction::TensionCutoffFunction(double tensileStrength)
        : mTensileStrength{tensileStrength}
    {
    }

    double TensionCutoffFunction::operator()(const Vector& rPrincipalStress) const
    {
        return mTensileStrength - rPrincipalStress(2);
    }

} // Namespace Kratos
