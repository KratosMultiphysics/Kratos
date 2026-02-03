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
#include <initializer_list>

namespace Kratos::Geo
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SigmaTau
{
public:
    static constexpr std::size_t msVectorSize = 2;
    using InternalVectorType                  = BoundedVector<double, msVectorSize>;

    SigmaTau() = default;

    template <typename VectorType>
    explicit SigmaTau(const VectorType& rValues) : SigmaTau{std::begin(rValues), std::end(rValues)}
    {
    }

    explicit SigmaTau(const std::initializer_list<double>& rValues);

    [[nodiscard]] const InternalVectorType& Values() const;
    [[nodiscard]] double                    Sigma() const;
    double&                                 Sigma();
    [[nodiscard]] double                    Tau() const;
    double&                                 Tau();

    template <typename VectorType>
    VectorType CopyTo() const
    {
        auto result = VectorType(msVectorSize);
        std::ranges::copy(mValues, result.begin());
        return result;
    }

private:
    template <typename InputIt>
    SigmaTau(InputIt First, InputIt Last)
    {
        KRATOS_DEBUG_ERROR_IF(std::distance(First, Last) != msVectorSize)
            << "Cannot construct a SigmaTau instance: expected " << msVectorSize
            << " values, but got " << std::distance(First, Last) << " value(s)\n";

        std::copy(First, Last, mValues.begin());
    }

    InternalVectorType mValues = ZeroVector{msVectorSize};
};

} // namespace Kratos::Geo