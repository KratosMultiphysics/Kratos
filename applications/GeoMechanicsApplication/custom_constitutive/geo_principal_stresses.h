// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include "includes/ublas_interface.h"

namespace Kratos::Geo
{

struct PrincipalStresses {
    BoundedVector<double, 3> values;
};

} // namespace Kratos::Geo