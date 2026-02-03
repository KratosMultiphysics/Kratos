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

#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

#include <algorithm>
#include <initializer_list>

namespace Kratos::Geo
{
class KRATOS_API(GEO_MECHANICS_APPLICATION) PrincipalStresses
{
public:
    static constexpr std::size_t msVectorSize = 3;
    using InternalVectorType                  = BoundedVector<double, msVectorSize>;

    PrincipalStresses() = default;

    template <typename VectorType>
    explicit PrincipalStresses(const VectorType& rStressVector)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rStressVector.size() == msVectorSize)
            << "PrincipalStresses can only be initialized with a vector of size " << msVectorSize
            << ", got " << rStressVector.size() << std::endl;
        std::ranges::copy(rStressVector, mValues.begin());
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

private:
    InternalVectorType mValues = ZeroVector{msVectorSize};
};
} // namespace Kratos::Geo