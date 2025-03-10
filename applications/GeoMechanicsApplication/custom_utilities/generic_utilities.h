// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

#pragma once

#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

#include <algorithm>
#include <cstdlib>
#include <vector>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GenericUtilities
{
public:
    template <typename Type>
    [[nodiscard]] static std::vector<Type> ApplyPermutation(const std::vector<Type>&        vec,
                                                            const std::vector<std::size_t>& indices)
    {
        std::vector<Type> result(vec.size());
        for (size_t i = 0; i < vec.size(); i++) {
            result[i] = vec[indices[i]];
        }
        return result;
    }

}; // class GenericUtilities

} // namespace Kratos
