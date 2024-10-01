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

#include "constitutive_dimension.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) PlaneStrainDimension : public ConstitutiveDimension
{
public:
    Matrix CreateConstitutiveMatrix(double c1, double c2, double c3) override;
};

} // namespace Kratos
