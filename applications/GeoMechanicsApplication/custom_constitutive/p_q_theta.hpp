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

    template <typename VectorType>
    explicit PQTheta(const VectorType& rValues)
    {
        // For some reason, the `std::ranges` versions of the below algorithms don't play nicely
        // with UBlas vector expressions. Therefore, we're using the iterator-style algorithms.
        auto first = std::begin(rValues);
        auto last  = std::end(rValues);
        KRATOS_DEBUG_ERROR_IF(std::distance(first, last) != msVectorSize)
            << "Cannot construct a PQTheta instance: expected " << msVectorSize
            << " values, but got " << std::distance(first, last) << " value(s)\n";
        CheckQ(*(first + 1));
        CheckTheta(*(first + 2));
        std::copy(first, last, mValues.begin());
    }

    [[nodiscard]] const InternalVectorType& Values() const noexcept;
    [[nodiscard]] double                    P() const noexcept;
    [[nodiscard]] double                    Q() const noexcept;
    [[nodiscard]] double                    Theta() const noexcept;

    template <typename VectorType>
    VectorType CopyTo() const
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

    static void CheckQ(double Q)
    {
        KRATOS_ERROR_IF(Q < 0.0) << "Invalid Q value: " << Q << ". Q must be non-negative." << std::endl;
    }

    static void CheckTheta(double Theta)
    {
        constexpr double max_abs_theta = Globals::Pi / 6.0;
        KRATOS_ERROR_IF(std::abs(Theta) > max_abs_theta)
            << "Invalid Theta value: " << Theta
            << ". Theta must be in the range [-pi/6, pi/6] (radians)." << std::endl;
    }
};

} // namespace Kratos::Geo