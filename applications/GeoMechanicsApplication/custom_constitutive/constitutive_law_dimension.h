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

#include <cstddef>
#include <memory>

namespace Kratos
{

class ConstitutiveLawDimension
{
public:
    virtual ~ConstitutiveLawDimension() = default;

    virtual Matrix FillConstitutiveMatrix(double c1, double c2, double c3) const = 0;
    virtual std::unique_ptr<ConstitutiveLawDimension> Clone() const              = 0;
    virtual std::size_t                               GetStrainSize() const      = 0;
    virtual std::size_t                               GetDimension() const       = 0;
    virtual Flags                                     GetSpatialType() const     = 0;
};

} // namespace Kratos
