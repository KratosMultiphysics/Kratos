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

namespace Kratos
{

class ContributionCalculator
{
public:
    virtual ~ContributionCalculator() = default;

    virtual Matrix                    LHSContribution()         = 0;
    virtual Vector                    RHSContribution()         = 0;
    virtual std::pair<Matrix, Vector> LocalSystemContribution() = 0;
};

} // namespace Kratos
