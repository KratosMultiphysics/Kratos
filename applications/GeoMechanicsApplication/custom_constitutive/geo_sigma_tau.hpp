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

    BoundedVector<double, vector_size> values = ZeroVector{vector_size};
    double&                            sigma;
    double&                            tau;

    SigmaTau();

    template <typename VectorType>
    explicit SigmaTau(const VectorType& rStressVector) : sigma{values[0]}, tau{values[1]}
    {
        std::ranges::copy(rStressVector, values.begin());
    }

    ~SigmaTau() = default;
    SigmaTau(const SigmaTau& rOther);
    SigmaTau& operator=(const SigmaTau& rOther);

    // For some unclear reason, I cannot explicitly delete the move operations...
    // SigmaTau(SigmaTau&&) noexcept = delete;
    // SigmaTau& operator=(SigmaTau&&) noexcept = delete;

    template <typename ContainerType>
    ContainerType CopyTo() const
    {
        auto result = ContainerType{values.size()};
        std::ranges::copy(values, result.begin());
        return result;
    }
};

} // namespace Kratos::Geo