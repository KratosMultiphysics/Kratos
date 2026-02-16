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
#include "includes/ublas_interface.h"

#include <algorithm>
#include <initializer_list>
#include <iterator>

namespace Kratos::Geo
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) PQ
{
public:
    static constexpr std::size_t msVectorSize = 2;
    using InternalVectorType                  = BoundedVector<double, msVectorSize>;

    PQ() = default;
    PQ(double P, double Q);

    template <typename VectorType>
    explicit PQ(const VectorType& rValues) : PQ{std::begin(rValues), std::end(rValues)}
    {
    }

    [[nodiscard]] const InternalVectorType& Values() const noexcept;
    [[nodiscard]] double                    P() const noexcept;
    double&                                 P() noexcept;
    [[nodiscard]] double                    Q() const noexcept;
    double&                                 Q() noexcept;

    template <typename VectorType>
    VectorType CopyTo() const
    {
        auto result = VectorType(msVectorSize);
        std::ranges::copy(mValues, result.begin());
        return result;
    }

private:
    template <std::forward_iterator Iter>
    PQ(Iter First, Iter Last)
    {
        KRATOS_DEBUG_ERROR_IF(std::distance(First, Last) != msVectorSize)
            << "Cannot construct a PQ instance: expected " << msVectorSize << " values, but got "
            << std::distance(First, Last) << " value(s)\n";

        std::copy(First, Last, mValues.begin());
    }

    InternalVectorType mValues = ZeroVector{msVectorSize};
};

} // namespace Kratos::Geo