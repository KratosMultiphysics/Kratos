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

struct SigmaTau {
    constexpr static auto vector_size = std::size_t{2};

    BoundedVector<double, vector_size> values;
    double&                            sigma;
    double&                            tau;

    template <typename VectorType>
    explicit SigmaTau(const VectorType& rStressVector)
        : values{ZeroVector{vector_size}}, sigma{values[0]}, tau{values[1]}
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