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

#include "containers/flags.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class ConstitutiveType
{
public:
    virtual ~ConstitutiveType() = default;

    virtual Matrix CreateConstitutiveMatrix(double c1, double c2, double c3) = 0;
    virtual std::unique_ptr<ConstitutiveType> Clone()                        = 0;
    virtual std::size_t                       GetStrainSize()                = 0;
    virtual std::size_t                       GetDimension()                 = 0;
    virtual Flags                             GetSpatialType()               = 0;
};

} // namespace Kratos
