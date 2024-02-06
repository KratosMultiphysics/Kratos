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
//                   Marjan Fathian
//
#pragma once

#include "includes/ublas_interface.h"

namespace Kratos
{

class StressStatePolicy
{
public:
    virtual Matrix CalculateBMatrix(const Matrix& GradNpT, const Vector& Np) const = 0;

    virtual ~StressStatePolicy() = default;
};

}