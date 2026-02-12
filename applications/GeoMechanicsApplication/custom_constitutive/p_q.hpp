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

#include "includes/exception.h"
#include "includes/kratos_export_api.h"

#include <algorithm>
#include <array>
#include <initializer_list>
#include <iterator>

namespace Kratos::Geo
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) PQ
{
public:
    static constexpr std::size_t msArraySize = 2;
    using InternalArrayType                  = std::array<double, msArraySize>;

    PQ() = default;

    template <typename VectorType>
    explicit PQ(const VectorType& rValues) : PQ{std::begin(rValues), std::end(rValues)}
    {
    }

    explicit PQ(const std::initializer_list<double>& rValues);

    [[nodiscard]] const InternalArrayType& Values() const;
    [[nodiscard]] double                   P() const;
    double&                                P();
    [[nodiscard]] double                   Q() const;
    double&                                Q();

    template <typename VectorType>
    VectorType CopyTo() const
    {
        auto result = VectorType(msArraySize);
        std::ranges::copy(mValues, result.begin());
        return result;
    }

private:
    template <std::forward_iterator Iter>
    PQ(Iter First, Iter Last)
    {
        KRATOS_DEBUG_ERROR_IF(std::distance(First, Last) != msArraySize)
            << "Cannot construct a PQ instance: expected " << msArraySize << " values, but got "
            << std::distance(First, Last) << " value(s)\n";

        std::copy(First, Last, mValues.begin());
    }

    InternalArrayType mValues = {0.0, 0.0};
};

} // namespace Kratos::Geo