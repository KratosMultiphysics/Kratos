// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohame Nabi,
//                   Wijtze Pieter Kikstra
// //

#pragma once

#include "custom_constitutive/yield_surface.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) TensionCutoff : public YieldSurface
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TensionCutoff);

    TensionCutoff() = default;

    explicit TensionCutoff(double tensileStrength);

    double YieldFunctionValue(const Vector& rPrincipalStress) const override;
    Vector DerivateOfFlowFunction(const Vector& rPrincipalStress) const override;

private:
    // Member Variables
    double mTensileStrength = 0.0;

}; // Class TensionCutoffFunction

} // namespace Kratos