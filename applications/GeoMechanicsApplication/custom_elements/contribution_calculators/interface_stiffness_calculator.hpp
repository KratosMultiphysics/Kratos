// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Richard Faasse

#pragma once

#include "contribution_calculator.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

#include <utility>

namespace Kratos
{
template <unsigned int TNumNodes>
class InterfaceStiffnessCalculator : public ContributionCalculator<TNumNodes>
{
public:
    struct InputProvider {
    };

    explicit InterfaceStiffnessCalculator(InputProvider AnInputProvider)
        : mInputProvider(std::move(AnInputProvider))
    {
    }

    std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> LHSContribution() override
    {
        return BoundedMatrix<double, TNumNodes, TNumNodes>{};
    }

    BoundedVector<double, TNumNodes> RHSContribution() override
    {
        return BoundedVector<double, TNumNodes>{};
    }

    std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> LocalSystemContribution() override
    {
        return {std::make_optional(LHSContribution()), RHSContribution()};
    }

private:
    InputProvider mInputProvider;
};

} // namespace Kratos
