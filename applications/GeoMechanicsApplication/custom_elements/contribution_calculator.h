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

class ContributionCalculator
{
public:
    virtual ~ContributionCalculator() = default;

    virtual std::optional<Matrix>                    LHSContribution()         = 0;
    virtual std::optional<Vector>                    RHSContribution()         = 0;
    virtual std::pair<std::optional<Matrix>, std::optional<Vector>> LocalSystemContribution() = 0;
};

} // namespace Kratos
