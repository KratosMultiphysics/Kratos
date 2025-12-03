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
    Vector result(rInitializerList.size(), 0.0);
    std::ranges::copy(rInitializerList, result.begin());
    return result;
}

Matrix UblasUtilities::CreateMatrix(const std::initializer_list<std::initializer_list<double>>& rInitializerList)
{
    if (rInitializerList.size() == 0) return {};

    Matrix result(rInitializerList.size(), rInitializerList.begin()->size(), 0.0);
    for (std::size_t i = 0; const auto& r_row : rInitializerList) {
        std::ranges::copy(r_row, row(result, i).begin());
        ++i;
    }

    return result;
}

} // namespace Kratos