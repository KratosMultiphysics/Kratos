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
                      std::function<Matrix()>                        GetMaterialPermeability,
                      std::function<Vector()>                        GetIntegrationCoefficients,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients,
                      std::function<std::vector<double>()> GetFluidPressures)
            : GetElementProperties(std::move(GetElementProperties)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetMaterialPermeability(std::move(GetMaterialPermeability)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetNodalValues(std::move(GetNodalValuesOf)),
              GetShapeFunctionGradients(std::move(GetShapeFunctionGradients)),
              GetFluidPressures(std::move(GetFluidPressures))
        {
        }

        std::function<const Properties&()>                           GetElementProperties;
        std::function<const std::vector<RetentionLaw::Pointer>&()>   GetRetentionLaws;
        std::function<Matrix()>                                      GetMaterialPermeability;
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
        const auto&              r_properties  = mInputProvider.GetElementProperties();
        const auto& r_integration_coefficients = mInputProvider.GetIntegrationCoefficients();
        const auto& r_shape_function_gradients = mInputProvider.GetShapeFunctionGradients();
        const auto  material_permeability      = mInputProvider.GetMaterialPermeability();

        BoundedMatrix<double, TNumNodes, TNumNodes> result = ZeroMatrix(TNumNodes, TNumNodes);

        const auto  dynamic_viscosity_inverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
        const auto& r_fluid_pressures         = mInputProvider.GetFluidPressures();

        for (std::size_t integration_point_index = 0;
             integration_point_index < r_integration_coefficients.size(); ++integration_point_index) {
            retention_parameters.SetFluidPressure(r_fluid_pressures[integration_point_index]);
            const auto relative_permeability =
                mInputProvider.GetRetentionLaws()[integration_point_index]->CalculateRelativePermeability(
                    retention_parameters);

            noalias(result) += GeoTransportEquationUtilities::CalculatePermeabilityMatrix(
                r_shape_function_gradients[integration_point_index], dynamic_viscosity_inverse, material_permeability,
                relative_permeability, r_integration_coefficients[integration_point_index]);
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
