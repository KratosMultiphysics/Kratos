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
class PUCouplingCalculator : public ContributionCalculator<NumberOfRows, NumberOfColumns>
{
public:
    struct InputProvider {
        InputProvider(std::function<const Matrix&()>       GetNpContainer,
                      Geo::BMatricesGetter                 GetBMatrices,
                      std::function<Vector()>              GetVoigtVector,
                      Geo::IntegrationCoefficientsGetter   GetIntegrationCoefficients,
                      std::function<std::vector<double>()> GetBiotCoefficients,
                      std::function<std::vector<double>()> GetDegreesOfSaturation,
                      std::function<Vector()>              GetNodalVelocities,
                      std::function<double()>              GetVelocityCoefficient)
            : GetNpContainer(std::move(GetNpContainer)),
              GetBMatrices(std::move(GetBMatrices)),
              GetVoigtVector(std::move(GetVoigtVector)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetBiotCoefficients(std::move(GetBiotCoefficients)),
              GetDegreesOfSaturation(std::move(GetDegreesOfSaturation)),
              GetNodalVelocities(std::move(GetNodalVelocities)),
              GetVelocityCoefficient(std::move(GetVelocityCoefficient))
        {
        }

        std::function<const Matrix&()>       GetNpContainer;
        Geo::BMatricesGetter                 GetBMatrices;
        std::function<Vector()>              GetVoigtVector;
        Geo::IntegrationCoefficientsGetter   GetIntegrationCoefficients;
        std::function<std::vector<double>()> GetBiotCoefficients;
        std::function<std::vector<double>()> GetDegreesOfSaturation;
        std::function<Vector()>              GetNodalVelocities;
        std::function<double()>              GetVelocityCoefficient;
    };

    explicit PUCouplingCalculator(InputProvider CouplingInputProvider)
        : mInputProvider{std::move(CouplingInputProvider)}
    {
    }

    using BaseType = ContributionCalculator<NumberOfRows, NumberOfColumns>;

    typename BaseType::LHSMatrixType CalculateCouplingMatrix()
    {
        // The NumberOfColumns and NumberOfRows are swapped because the coupling matrix is
        // transposed before returning it.
        Matrix     result                   = ZeroMatrix(NumberOfColumns, NumberOfRows);
        const auto b_matrices               = mInputProvider.GetBMatrices();
        const auto voigt_vector             = mInputProvider.GetVoigtVector();
        const auto integration_coefficients = mInputProvider.GetIntegrationCoefficients();
        const auto biot_coefficients        = mInputProvider.GetBiotCoefficients();
        const auto bishop_coefficients      = mInputProvider.GetDegreesOfSaturation();
        const auto np_container             = mInputProvider.GetNpContainer();
        for (int i = 0; i < mInputProvider.GetBMatrices().size(); ++i) {
            Matrix coupling_contribution(NumberOfColumns, NumberOfRows);
            GeoTransportEquationUtilities::CalculateCouplingMatrix(
                coupling_contribution, b_matrices[i], voigt_vector, row(np_container, i),
                biot_coefficients[i], bishop_coefficients[i], integration_coefficients[i]);
            result += coupling_contribution;
        }

        return trans(result) * mInputProvider.GetVelocityCoefficient() * PORE_PRESSURE_SIGN_FACTOR;
    }

    std::optional<typename BaseType::LHSMatrixType> LHSContribution() override
    {
        return CalculateCouplingMatrix();
    }

    typename BaseType::RHSVectorType RHSContribution() override
    {
        return prod(CalculateCouplingMatrix(), mInputProvider.GetNodalVelocities());
    }

    std::pair<std::optional<typename BaseType::LHSMatrixType>, typename BaseType::RHSVectorType> LocalSystemContribution() override
    {
        return {LHSContribution(), RHSContribution()};
    };

    InputProvider mInputProvider;
};
} // namespace Kratos
