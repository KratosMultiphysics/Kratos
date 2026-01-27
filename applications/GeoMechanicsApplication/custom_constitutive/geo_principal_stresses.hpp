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

#include <algorithm>

namespace Kratos::Geo
{

struct PrincipalStresses {
    BoundedVector<double, 3> values;

    template <typename VectorType>
    explicit PrincipalStresses(const VectorType& rStressVector)
    {
        std::ranges::copy(rStressVector, values.begin());
    }

    template <typename ContainerType>
    ContainerType CopyTo() const
    {
        auto result = ContainerType{values.size()};
        std::ranges::copy(values, result.begin());
        return result;
    }
};

} // namespace Kratos::Geo