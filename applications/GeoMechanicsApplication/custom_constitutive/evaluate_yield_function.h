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

#pragma once

#include "includes/serializer.h"

namespace Kratos
{

class EvaluateYieldFunction
{
public:
    virtual ~EvaluateYieldFunction() = default;

    virtual double CalculateYieldFunction(const Vector& rPrincipalStress) const        = 0;
    virtual Vector CalculateFlowFunctionDerivate(const Vector& rPrincipalStress) const = 0;


};
} // namespace Kratos.
