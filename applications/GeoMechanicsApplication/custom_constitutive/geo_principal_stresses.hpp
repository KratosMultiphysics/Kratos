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

struct KRATOS_API(GEO_MECHANICS_APPLICATION) PrincipalStresses {
    constexpr static auto vector_size = std::size_t{3};

    BoundedVector<double, vector_size> values = ZeroVector{vector_size};

    PrincipalStresses() = default;

    template <typename VectorType>
    explicit PrincipalStresses(const VectorType& rStressVector)
        : PrincipalStresses{std::begin(rStressVector), std::end(rStressVector)}
    {
    }

    explicit PrincipalStresses(const std::initializer_list<double>& rValues);

    template <typename ContainerType>
    ContainerType CopyTo() const
    {
        auto result = ContainerType{values.size()};
        std::ranges::copy(values, result.begin());
        return result;
    }

private:
    template <typename InputIt>
    PrincipalStresses(InputIt First, InputIt Last)
    {
        std::copy(First, Last, values.begin());
    }
};

} // namespace Kratos::Geo