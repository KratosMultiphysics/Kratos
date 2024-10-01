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

class ConstitutiveDimension
{
public:
    virtual Matrix CreateConstitutiveMatrix(double c1, double c2, double c3) = 0;
    virtual std::unique_ptr<ConstitutiveDimension> Clone() = 0;
};

} // namespace Kratos
