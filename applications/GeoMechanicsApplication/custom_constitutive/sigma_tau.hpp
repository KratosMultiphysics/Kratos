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

class KRATOS_API(GEO_MECHANICS_APPLICATION) SigmaTau
{
public:
    static constexpr std::size_t msArraySize = 2;
    using InternalArrayType                  = std::array<double, msArraySize>;

    SigmaTau() = default;

    template <typename VectorType>
    explicit SigmaTau(const VectorType& rValues) : SigmaTau{std::begin(rValues), std::end(rValues)}
    {
    }

    explicit SigmaTau(const std::initializer_list<double>& rValues);

    [[nodiscard]] const InternalArrayType& Values() const;
    [[nodiscard]] double                   Sigma() const;
    double&                                Sigma();
    [[nodiscard]] double                   Tau() const;
    double&                                Tau();

    template <typename VectorType>
    VectorType CopyTo() const
    {
        auto result = VectorType(msArraySize);
        std::ranges::copy(mValues, result.begin());
        return result;
    }

    // See Section "Binary arithmetic operators" at https://en.cppreference.com/w/cpp/language/operators.html
    SigmaTau& operator+=(const SigmaTau& rRhsTraction);
    KRATOS_API(GEO_MECHANICS_APPLICATION)
    friend SigmaTau operator+(SigmaTau LhsTraction, /* passing this one by value helps optimize chained a+b+c */
                              const SigmaTau& rRhsTraction);

private:
    template <std::forward_iterator Iter>
    SigmaTau(Iter First, Iter Last)
    {
        KRATOS_DEBUG_ERROR_IF(std::distance(First, Last) != msArraySize)
            << "Cannot construct a SigmaTau instance: expected " << msArraySize
            << " values, but got " << std::distance(First, Last) << " value(s)\n";

        std::copy(First, Last, mValues.begin());
    }

    InternalArrayType mValues = {0.0, 0.0};
};

} // namespace Kratos::Geo
