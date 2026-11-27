// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Wijtze Pieter Kikstra
//

#pragma once

#include "includes/default_interface.h"
#include "includes/exception.h"
#include "includes/global_variables.h"
#include "includes/kratos_export_api.h"

#include <algorithm>
#include <cstddef> // for std::size_t
#include <iterator>

namespace Kratos::Geo
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) PQTheta
{
public:
    static constexpr std::size_t msVectorSize = 3;
    using InternalVectorType                  = BoundedVector<double, msVectorSize>;

    PQTheta() = default;
    PQTheta(double P, double Q, double Theta); // Pressure, deviatoric stress, Lode angle in radians

    explicit PQTheta(const InternalVectorType& rValues);

    template <typename VectorType>
    explicit PQTheta(const VectorType& rValues) : PQTheta{SizeCheckedCopy(rValues)}
    {
    }

    [[nodiscard]] const InternalVectorType& Values() const noexcept;
    [[nodiscard]] double                    P() const noexcept;
    [[nodiscard]] double                    Q() const noexcept;
    [[nodiscard]] double                    Theta() const noexcept;

    template <typename VectorType>
    [[nodiscard]] VectorType CopyTo() const
    {
        auto result = VectorType(msVectorSize);
        std::ranges::copy(mValues, result.begin());
        return result;
    }

    PQTheta& operator+=(const PQTheta& rRhs);
    KRATOS_API(GEO_MECHANICS_APPLICATION)
    friend PQTheta operator+(PQTheta Lhs, const PQTheta& rRhs);

private:
    InternalVectorType mValues = ZeroVector{msVectorSize};

    static InternalVectorType MakeVector(double P, double Q, double Theta);
    void                      CheckValues() const;

    template <typename VectorType>
    static InternalVectorType SizeCheckedCopy(const VectorType& rValues)
    {
        // For some reason, the `std::ranges` versions of the below algorithms don't play nicely
        // with UBlas vector expressions. Therefore, we're using the iterator-style algorithms.
        auto       first = std::begin(rValues);
        const auto size  = std::distance(first, std::end(rValues));
        KRATOS_ERROR_IF(size != static_cast<std::ptrdiff_t>(msVectorSize))
            << "Cannot construct a PQTheta instance: expected " << msVectorSize
            << " values, but got " << size << " value(s)\n";

        return MakeVector(*first, *std::next(first, 1), *std::next(first, 2));
    }
};

KRATOS_API(GEO_MECHANICS_APPLICATION)
PQTheta operator+(PQTheta Lhs, const PQTheta& rRhs);

} // namespace Kratos::Geo