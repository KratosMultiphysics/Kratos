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
#include "geo_aliases.h"

namespace Kratos
{
template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
class UPCouplingCalculator : public ContributionCalculator<NumberOfRows, NumberOfColumns>
{
public:
    struct InputProvider {
        InputProvider(std::function<const Matrix&()>       GetNpContainer,
                      Geo::BMatricesGetter                 GetBMatrices,
                      std::function<Vector()>              GetVoigtVector,
                      Geo::IntegrationCoefficientsGetter   GetIntegrationCoefficients,
                      std::function<std::vector<double>()> GetBiotCoefficients,
                      std::function<std::vector<double>()> GetBishopCoefficients)
            : GetNpContainer(std::move(GetNpContainer)),
              GetBMatrices(std::move(GetBMatrices)),
              GetVoigtVector(std::move(GetVoigtVector)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetBiotCoefficients(std::move(GetBiotCoefficients)),
              GetBishopCoefficients(std::move(GetBishopCoefficients))
        {
        }

        std::function<const Matrix&()>       GetNpContainer;
        Geo::BMatricesGetter                 GetBMatrices;
        std::function<Vector()>              GetVoigtVector;
        Geo::IntegrationCoefficientsGetter   GetIntegrationCoefficients;
        std::function<std::vector<double>()> GetBiotCoefficients;
        std::function<std::vector<double>()> GetBishopCoefficients;
    };

    explicit UPCouplingCalculator(InputProvider CouplingInputProvider)
        : mInputProvider{std::move(CouplingInputProvider)}
    {
    }

    using BaseType = ContributionCalculator<NumberOfRows, NumberOfColumns>;

    std::optional<typename BaseType::LHSMatrixType> LHSContribution() override { return {}; }

    typename BaseType::RHSVectorType RHSContribution() override { return {}; }

    std::pair<std::optional<typename BaseType::LHSMatrixType>, typename BaseType::RHSVectorType> LocalSystemContribution() override
    {
        return {LHSContribution(), RHSContribution()};
    };

    InputProvider mInputProvider;
};
} // namespace Kratos
