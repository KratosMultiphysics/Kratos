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

#include "interface_constitutive_law_dimension.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) InterfaceThreeDimensionalSurface : public InterfaceConstitutiveLawDimension
{
public:
    [[nodiscard]] Matrix MakeInterfaceConstitutiveMatrix(double      NormalStiffness,
                                                         double      ShearStiffness,
                                                         std::size_t TractionSize) const override;
    [[nodiscard]] std::unique_ptr<InterfaceConstitutiveLawDimension> Clone() const override;
    [[nodiscard]] std::size_t                                        GetStrainSize() const override;
    [[nodiscard]] std::size_t                                        GetDimension() const override;

private:
    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};

} // namespace Kratos
