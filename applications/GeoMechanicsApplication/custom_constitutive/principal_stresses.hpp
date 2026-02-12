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

#pragma once

#include "includes/exception.h"
#include "includes/kratos_export_api.h"

#include <algorithm>
#include <array>
#include <initializer_list>
#include <iterator>

namespace Kratos::Geo
{
class KRATOS_API(GEO_MECHANICS_APPLICATION) PrincipalStresses
{
public:
    static constexpr std::size_t msArraySize = 3;
    using InternalArrayType                  = std::array<double, msArraySize>;

    PrincipalStresses() = default;

    template <typename VectorType>
    explicit PrincipalStresses(const VectorType& rStressVector)
        : PrincipalStresses{std::begin(rStressVector), std::end(rStressVector)}
    {
    }

    template <std::forward_iterator InputIt>
    PrincipalStresses(InputIt First, InputIt Last)
    {
        KRATOS_DEBUG_ERROR_IF(std::distance(First, Last) != msArraySize)
            << "Cannot construct a PrincipalStresses instance: expected " << msArraySize
            << " values, but got " << std::distance(First, Last) << " value(s)\n";

        std::copy(First, Last, mValues.begin());
    }

    explicit PrincipalStresses(const std::initializer_list<double>& rValues);

    template <typename VectorType>
    VectorType CopyTo()
    {
        VectorType result(msArraySize);
        std::ranges::copy(mValues, result.begin());
        return result;
    }

    [[nodiscard]] const InternalArrayType& Values() const;
    [[nodiscard]] InternalArrayType&       Values();

private:
    InternalArrayType mValues = {0.0, 0.0, 0.0};
};
} // namespace Kratos::Geo