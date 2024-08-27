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

namespace Kratos
{

ConstitutiveLaw::SizeType GeoIncrementalLinearElasticInterfaceLaw::WorkingSpaceDimension()
{
    return 2;
}

} // namespace Kratos
