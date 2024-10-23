// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#pragma once

#include "includes/ublas_interface.h"

namespace Kratos
{

class Calculator
{
public:
    virtual ~Calculator() = default;

    virtual Matrix LHSContribution() = 0;
    virtual Vector RHSContribution() = 0;
    virtual std::pair<Matrix, Vector> CalculateLeftAndRightHandSide() = 0;
};

} // namespace Kratos
