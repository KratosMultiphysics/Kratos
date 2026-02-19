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
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_aliases.h"

namespace Kratos
{
template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
class UPCouplingCalculator : public ContributionCalculator<NumberOfRows, NumberOfColumns>
{
public:
    struct InputProvider {
        InputProvider(std::function<const Matrix()>        GetNpContainer,
                      Geo::BMatricesGetter                 GetBMatrices,
                      std::function<Vector()>              GetVoigtVector,
                      Geo::IntegrationCoefficientsGetter   GetIntegrationCoefficients,
                      std::function<std::vector<double>()> GetBiotCoefficients,
                      std::function<std::vector<double>()> GetBishopCoefficients,
                      std::function<Vector()>              GetNodalWaterPressures)
            : GetNpContainer(std::move(GetNpContainer)),
              GetBMatrices(std::move(GetBMatrices)),
              GetVoigtVector(std::move(GetVoigtVector)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetBiotCoefficients(std::move(GetBiotCoefficients)),
              GetBishopCoefficients(std::move(GetBishopCoefficients)),
              GetNodalWaterPressures(std::move(GetNodalWaterPressures))
        {
        }

        std::function<const Matrix()>        GetNpContainer;
        Geo::BMatricesGetter                 GetBMatrices;
        std::function<Vector()>              GetVoigtVector;
        Geo::IntegrationCoefficientsGetter   GetIntegrationCoefficients;
        std::function<std::vector<double>()> GetBiotCoefficients;
        std::function<std::vector<double>()> GetBishopCoefficients;
        std::function<Vector()>              GetNodalWaterPressures;
    };

    explicit UPCouplingCalculator(InputProvider CouplingInputProvider)
        : mInputProvider{std::move(CouplingInputProvider)}
    {
    }

    using BaseType = ContributionCalculator<NumberOfRows, NumberOfColumns>;

    std::optional<typename BaseType::LHSMatrixType> LHSContribution() override
    {
        return CalculateCouplingMatrix();
    }

    typename BaseType::RHSVectorType RHSContribution() override
    {
        auto result = typename BaseType::RHSVectorType{
            prod(CalculateCouplingMatrix(), mInputProvider.GetNodalWaterPressures())};
        KRATOS_INFO("CouplingCalculator::RHSContribution") << result << std::endl;
        return result;
    }

    std::pair<std::optional<typename BaseType::LHSMatrixType>, typename BaseType::RHSVectorType> LocalSystemContribution() override
    {
        return {LHSContribution(), RHSContribution()};
    };

private:
    typename BaseType::LHSMatrixType CalculateCouplingMatrix()
    {
        typename BaseType::LHSMatrixType result       = ZeroMatrix(NumberOfRows, NumberOfColumns);
        const auto                       b_matrices   = mInputProvider.GetBMatrices();
        const auto                       voigt_vector = mInputProvider.GetVoigtVector();
        const auto integration_coefficients           = mInputProvider.GetIntegrationCoefficients();
        const auto biot_coefficients                  = mInputProvider.GetBiotCoefficients();
        const auto bishop_coefficients                = mInputProvider.GetBishopCoefficients();
        const auto np_container                       = mInputProvider.GetNpContainer();
        for (auto i = std::size_t{0}; i < integration_coefficients.size(); ++i) {
            typename BaseType::LHSMatrixType coupling_contribution;
            GeoTransportEquationUtilities::CalculateCouplingMatrix(
                coupling_contribution, b_matrices[i], voigt_vector, row(np_container, i),
                biot_coefficients[i], bishop_coefficients[i], integration_coefficients[i]);
            result += coupling_contribution;
        }
        return result;
    }

    InputProvider mInputProvider;
};
} // namespace Kratos
