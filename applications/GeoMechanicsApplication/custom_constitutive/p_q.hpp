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
    explicit PQ(const VectorType& rValues)
    {
        KRATOS_DEBUG_ERROR_IF(std::ranges::distance(rValues) != msVectorSize)
            << "Cannot construct a PQ instance: expected " << msVectorSize << " values, but got "
            << std::ranges::distance(rValues) << " value(s)\n";

        std::ranges::copy(rValues, mValues.begin());
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
    InternalVectorType mValues = ZeroVector{msVectorSize};
};

} // namespace Kratos::Geo