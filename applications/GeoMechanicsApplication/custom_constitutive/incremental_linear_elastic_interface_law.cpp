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

#include "incremental_linear_elastic_interface_law.h"
#include "geo_mechanics_application_constants.h"

namespace Kratos
{

ConstitutiveLaw::SizeType GeoIncrementalLinearElasticInterfaceLaw::WorkingSpaceDimension()
{
    return 2;
}

ConstitutiveLaw::SizeType GeoIncrementalLinearElasticInterfaceLaw::GetStrainSize() const
{
    return VOIGT_SIZE_2D_INTERFACE;
}

ConstitutiveLaw::StressMeasure GeoIncrementalLinearElasticInterfaceLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

bool GeoIncrementalLinearElasticInterfaceLaw::IsIncremental()
{
    return true;
}

} // namespace Kratos
