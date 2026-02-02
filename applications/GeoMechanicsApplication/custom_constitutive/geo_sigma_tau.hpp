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

#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

#include <algorithm>

namespace Kratos::Geo
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SigmaTau
{
public:
    static constexpr std::size_t msVectorSize = 2;
    using InternalVectorType                  = BoundedVector<double, msVectorSize>;

    SigmaTau() = default;

    template <typename VectorType>
    explicit SigmaTau(const VectorType& rValues)
    {
        std::ranges::copy(rValues, mValues.begin());
    }

    const InternalVectorType& Values() const;
    double                    Sigma() const;
    double                    Tau() const;

private:
    InternalVectorType mValues = ZeroVector{msVectorSize};
};

} // namespace Kratos::Geo