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
class StiffnessCalculator : public ContributionCalculator<TNumNodes>
{
public:
    struct InputProvider {
        InputProvider(std::function<std::vector<Matrix>()> GetBMatrices,
                      std::function<std::vector<Vector>()> GetRelativeDisplacements,
                      std::function<Vector()>              GetIntegrationCoefficients)
            : GetBMatrices(std::move(GetBMatrices)),
              GetRelativeDisplacements(std::move(GetRelativeDisplacements)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients))
        {
        }

        std::function<std::vector<Vector>()> GetRelativeDisplacements;
        std::function<std::vector<Matrix>()> GetBMatrices;
        std::function<Vector()>              GetIntegrationCoefficients;
    };

    explicit StiffnessCalculator(InputProvider AnInputProvider)
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
