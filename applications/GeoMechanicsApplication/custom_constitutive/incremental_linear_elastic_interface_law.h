// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#pragma once

#include "includes/constitutive_law.h"
#include "includes/kratos_export_api.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoIncrementalLinearElasticInterfaceLaw : public ConstitutiveLaw
{
public:
    SizeType WorkingSpaceDimension() override;
    SizeType GetStrainSize() const override;
    StressMeasure GetStressMeasure() override;
    bool          IsIncremental() override;
};

} // namespace Kratos
