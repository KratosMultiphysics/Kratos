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
#include "up_coupling_calculator.hpp"

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
                      std::function<Vector()>              GetNodalVelocities)
            : GetNpContainer(std::move(GetNpContainer)),
              GetBMatrices(std::move(GetBMatrices)),
              GetVoigtVector(std::move(GetVoigtVector)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetBiotCoefficients(std::move(GetBiotCoefficients)),
              GetDegreesOfSaturation(std::move(GetDegreesOfSaturation)),
              GetNodalVelocities(std::move(GetNodalVelocities))
        {
        }

        std::function<const Matrix&()>       GetNpContainer;
        Geo::BMatricesGetter                 GetBMatrices;
        std::function<Vector()>              GetVoigtVector;
        Geo::IntegrationCoefficientsGetter   GetIntegrationCoefficients;
        std::function<std::vector<double>()> GetBiotCoefficients;
        std::function<std::vector<double>()> GetDegreesOfSaturation;
        std::function<Vector()>              GetNodalVelocities;
    };

    explicit PUCouplingCalculator(InputProvider CouplingInputProvider)
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
        return prod(CalculateCouplingMatrix(), mInputProvider.GetNodalVelocities());
    }

    std::pair<std::optional<typename BaseType::LHSMatrixType>, typename BaseType::RHSVectorType> LocalSystemContribution() override
    {
        return {LHSContribution(), RHSContribution()};
    };

private:
    typename BaseType::LHSMatrixType CalculateCouplingMatrix()
    {
        // For the calculation, we can re-use the calculation of the UP coupling matrix.
        // To do this, the NumberOfColumns and NumberOfRows are swapped because the
        // PU coupling matrix is the transpose of the UP coupling matrix.
        // Also note that instead of the Bishop coefficients, the degrees of saturation are passed.
        // Finally, the nodal water pressure getter is not needed for the matrix calculation, so
        // we pass a dummy function.
        auto dummy_nodal_water_pressure_function = []() { return Vector(); };
        typename UPCouplingCalculator<NumberOfColumns, NumberOfRows>::InputProvider input_provider(
            mInputProvider.GetNpContainer, mInputProvider.GetBMatrices, mInputProvider.GetVoigtVector,
            mInputProvider.GetIntegrationCoefficients, mInputProvider.GetBiotCoefficients,
            mInputProvider.GetDegreesOfSaturation, dummy_nodal_water_pressure_function);
        UPCouplingCalculator<NumberOfColumns, NumberOfRows> up_coupling_calculator(input_provider);
        const auto up_coupling_matrix = up_coupling_calculator.LHSContribution().value();

        return trans(up_coupling_matrix) * PORE_PRESSURE_SIGN_FACTOR;
    }

    InputProvider mInputProvider;
};
} // namespace Kratos
