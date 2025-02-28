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

#include "constitutive_law_dimension.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) PlaneStrain : public ConstitutiveLawDimension
{
public:
    [[nodiscard]] Matrix CalculateElasticMatrix(double YoungsModulus, double PoissonsRatio) const override;
    [[nodiscard]] std::unique_ptr<ConstitutiveLawDimension> Clone() const override;
    [[nodiscard]] std::size_t                               GetStrainSize() const override;
    [[nodiscard]] std::size_t                               GetDimension() const override;
    [[nodiscard]] Flags                                     GetSpatialType() const override;
};

} // namespace Kratos
