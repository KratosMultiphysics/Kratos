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

#include "constitutive_type.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) PlaneStrainType : public ConstitutiveType
{
public:
    Matrix CreateConstitutiveMatrix(double c1, double c2, double c3) override;
    std::unique_ptr<ConstitutiveType> Clone() override;
    std::size_t                       GetStrainSize() override;
    std::size_t                       GetDimension() override;
    Flags                             GetConstitutiveLawType() override;
};

} // namespace Kratos
