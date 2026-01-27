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

#include "contribution_calculator.h"

namespace Kratos
{
    template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
    class UPCouplingCalculator : public ContributionCalculator<NumberOfRows, NumberOfColumns>
    {
    public:
        using BaseType = ContributionCalculator<NumberOfRows, NumberOfColumns>;

        std::optional<typename BaseType::LHSMatrixType> LHSContribution() override
        {
            return {};
        }

        typename BaseType::RHSVectorType RHSContribution() override
        {
            return {};
        }

        std::pair<std::optional<typename BaseType::LHSMatrixType>,
                  typename BaseType::RHSVectorType> LocalSystemContribution() override
        {
            return {LHSContribution(), RHSContribution()};
        };
    };
}
