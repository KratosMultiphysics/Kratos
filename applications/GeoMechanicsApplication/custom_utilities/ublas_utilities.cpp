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

#include "ublas_utilities.h"
#include <algorithm>

namespace Kratos
{

Vector UblasUtilities::CreateVector(const std::initializer_list<double>& rInitializerList)
{
    Vector result = ZeroVector(rInitializerList.size());
    std::ranges::copy(rInitializerList, result.begin());
    return result;
}
} // namespace Kratos