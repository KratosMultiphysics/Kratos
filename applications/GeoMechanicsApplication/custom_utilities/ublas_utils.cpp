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

#include <algorithm>

#include "ublas_utils.h"

namespace Kratos
{

Vector UBlasUtils::MakeVector(const std::initializer_list<double>& values)
{
    auto result = Vector{values.size()};
    std::copy(values.begin(), values.end(), result.begin());
    return result;
}

}
