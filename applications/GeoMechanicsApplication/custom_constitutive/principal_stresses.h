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

namespace Kratos::Geo
{
class KRATOS_API(GEO_MECHANICS_APPLICATION) PrincipalStresses
{
public:
    PrincipalStresses() = default;

    template <typename VectorType>
    explicit PrincipalStresses(const VectorType& rStressVector)
    {
        std::ranges::copy(rStressVector, mValues.begin());
    }

    static constexpr std::size_t msVectorSize = 3;
    using InternalVectorType                  = BoundedVector<double, msVectorSize>;
    InternalVectorType Values() const const;

private:
    InternalVectorType mValues = ZeroVector{msVectorSize};
};
} // namespace Kratos::Geo