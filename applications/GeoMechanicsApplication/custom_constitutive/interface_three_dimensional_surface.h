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
//                   Gennady Markelov
//

#pragma once

#include "constitutive_law_dimension.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) InterfaceThreeDimensionalSurface : public ConstitutiveLawDimension
{
public:
    [[nodiscard]] Matrix CalculateElasticMatrix(const Properties& rProperties) const override;
    [[nodiscard]] std::unique_ptr<ConstitutiveLawDimension> Clone() const override;
    [[nodiscard]] std::size_t                               GetStrainSize() const override;
    [[nodiscard]] std::size_t                               GetDimension() const override;
    [[nodiscard]] std::size_t GetNumberOfNormalComponents() const override;
    [[nodiscard]] Flags       GetSpatialType() const override;

private:
    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};

} // namespace Kratos
