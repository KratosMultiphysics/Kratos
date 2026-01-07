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

template <unsigned int TNumNodes>
class ContributionCalculator
{
public:
    virtual ~ContributionCalculator() = default;

    virtual std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> LHSContribution() = 0;
    virtual BoundedVector<double, TNumNodes>                           RHSContribution() = 0;
    virtual std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> LocalSystemContribution() = 0;
};

} // namespace Kratos
