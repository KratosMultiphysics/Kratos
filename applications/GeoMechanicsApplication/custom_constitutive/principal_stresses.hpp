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
#include "includes/ublas_interface.h"

#include <algorithm>
#include <initializer_list>
#include <iterator>

namespace Kratos::Geo
{
class KRATOS_API(GEO_MECHANICS_APPLICATION) PrincipalStresses
{
public:
    enum class PrincipalStressesAveragingType {
        NO_AVERAGING,
        LOWEST_PRINCIPAL_STRESSES,
        HIGHEST_PRINCIPAL_STRESSES
    };

    static constexpr std::size_t msVectorSize = 3;
    using InternalVectorType                  = BoundedVector<double, msVectorSize>;

    PrincipalStresses() = default;

    template <typename VectorType>
    explicit PrincipalStresses(const VectorType& rStressVector)
        : PrincipalStresses{std::begin(rStressVector), std::end(rStressVector)}
    {
    }

    template <std::forward_iterator InputIt>
    PrincipalStresses(InputIt First, InputIt Last)
    {
        KRATOS_DEBUG_ERROR_IF(std::distance(First, Last) != msVectorSize)
            << "Cannot construct a PrincipalStresses instance: expected " << msVectorSize
            << " values, but got " << std::distance(First, Last) << " value(s)\n";

        std::copy(First, Last, mValues.begin());
    }

    explicit PrincipalStresses(const std::initializer_list<double>& rValues);

    template <typename VectorType>
    VectorType CopyTo()
    {
        VectorType result(msVectorSize);
        std::ranges::copy(mValues, result.begin());
        return result;
    }

    [[nodiscard]] const InternalVectorType& Values() const;
    [[nodiscard]] InternalVectorType&       Values();

    PrincipalStresses& operator+=(const PrincipalStresses& rRhs);
    KRATOS_API(GEO_MECHANICS_APPLICATION)
    friend PrincipalStresses operator+(PrincipalStresses Lhs, const PrincipalStresses& rRhs);

private:
    InternalVectorType mValues = ZeroVector{msVectorSize};
};
} // namespace Kratos::Geo