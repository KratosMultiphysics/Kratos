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

#include "includes/ublas_interface.h"

#include <optional>

namespace Kratos
{

template <unsigned int NumberOfRows, unsigned int NumberOfColumns = NumberOfRows>
class ContributionCalculator
{
public:
    virtual ~ContributionCalculator() = default;

    using LHSMatrixType = BoundedMatrix<double, NumberOfRows, NumberOfColumns>;
    using RHSVectorType = BoundedVector<double, NumberOfRows>;

    virtual std::optional<LHSMatrixType>                           LHSContribution()         = 0;
    virtual RHSVectorType                                          RHSContribution()         = 0;
    virtual std::pair<std::optional<LHSMatrixType>, RHSVectorType> LocalSystemContribution() = 0;
};

} // namespace Kratos
