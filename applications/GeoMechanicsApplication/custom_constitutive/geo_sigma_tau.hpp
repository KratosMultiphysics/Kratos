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

#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

#include <algorithm>

namespace Kratos::Geo
{

struct KRATOS_API(GEO_MECHANICS_APPLICATION) SigmaTau {
    constexpr static auto vector_size = std::size_t{2};

    BoundedVector<double, vector_size> values;
    double&                            sigma;
    double&                            tau;

    SigmaTau();

    template <typename VectorType>
    explicit SigmaTau(const VectorType& rStressVector)
        : SigmaTau{std::begin(rStressVector), std::end(rStressVector)}
    {
    }

    explicit SigmaTau(const std::initializer_list<double>& rValues);

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

private:
    template <typename InputIt>
    SigmaTau(InputIt First, InputIt Last)
        : values{ZeroVector{vector_size}}, sigma{values[0]}, tau{values[1]}
    {
        std::copy(First, Last, values.begin());
    }
};

} // namespace Kratos::Geo