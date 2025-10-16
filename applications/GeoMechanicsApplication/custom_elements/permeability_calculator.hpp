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
#include "custom_retention/retention_law.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "includes/cfd_variables.h"

namespace Kratos
{

template <unsigned int TNumNodes>
class PermeabilityCalculator : public ContributionCalculator<TNumNodes> ///
{
public:
    struct InputProvider {
        InputProvider(std::function<const Properties&()> GetElementProperties,
                      std::function<const std::vector<RetentionLaw::Pointer>&()> GetRetentionLaws,
                      std::function<Vector()>                        GetIntegrationCoefficients,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients,
                      std::function<std::vector<double>()> GetFluidPressures)
            : GetElementProperties(std::move(GetElementProperties)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetNodalValues(std::move(GetNodalValuesOf)),
              GetShapeFunctionGradients(std::move(GetShapeFunctionGradients)),
              GetFluidPressures(std::move(GetFluidPressures))
        {
        }

        std::function<const Properties&()>                           GetElementProperties;
        std::function<const std::vector<RetentionLaw::Pointer>&()>   GetRetentionLaws;
        std::function<Vector()>                                      GetIntegrationCoefficients;
        std::function<Vector(const Variable<double>&)>               GetNodalValues;
        std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients;
        std::function<std::vector<double>()>                         GetFluidPressures;
    };

    explicit PermeabilityCalculator(InputProvider InputProvider)
        : mInputProvider{std::move(InputProvider)}
    {
    }

    std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> LHSContribution() override
    {
        return std::make_optional(CalculatePermeabilityMatrix());
    }

    BoundedVector<double, TNumNodes> RHSContribution() override
    {
        return RHSContribution(CalculatePermeabilityMatrix());
    }

    std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> LocalSystemContribution() override
    {
        const auto permeability_matrix = CalculatePermeabilityMatrix();
        return {std::make_optional(permeability_matrix), RHSContribution(permeability_matrix)};
    }

private:
    InputProvider mInputProvider;

    [[nodiscard]] BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix() const
    {
        RetentionLaw::Parameters retention_parameters(mInputProvider.GetElementProperties());
        const auto&              r_properties      = mInputProvider.GetElementProperties();
        const auto&       integration_coefficients = mInputProvider.GetIntegrationCoefficients();
        const auto&       shape_function_gradients = mInputProvider.GetShapeFunctionGradients();
        const std::size_t local_dimension          = shape_function_gradients[0].size2();
        const Matrix      constitutive_matrix =
            GeoElementUtilities::FillPermeabilityMatrix(r_properties, local_dimension);

        BoundedMatrix<double, TNumNodes, TNumNodes> result = ZeroMatrix(TNumNodes, TNumNodes);

        const double dynamic_viscosity_inverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
        const auto&  fluid_pressures           = mInputProvider.GetFluidPressures();

        for (std::size_t integration_point_index = 0;
             integration_point_index < integration_coefficients.size(); ++integration_point_index) {
            retention_parameters.SetFluidPressure(fluid_pressures[integration_point_index]);
            const double relative_permeability =
                mInputProvider.GetRetentionLaws()[integration_point_index]->CalculateRelativePermeability(
                    retention_parameters);

            noalias(result) += GeoTransportEquationUtilities::CalculatePermeabilityMatrix(
                shape_function_gradients[integration_point_index], dynamic_viscosity_inverse, constitutive_matrix,
                relative_permeability, integration_coefficients[integration_point_index]);
        }
        return result;
    }

    [[nodiscard]] BoundedVector<double, TNumNodes> RHSContribution(
        const BoundedMatrix<double, TNumNodes, TNumNodes>& rPermeabilityMatrix) const
    {
        return -prod(rPermeabilityMatrix, mInputProvider.GetNodalValues(WATER_PRESSURE));
    }
};

} // namespace Kratos
