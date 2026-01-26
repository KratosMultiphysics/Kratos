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
#include "geo_aliases.h"
#include "includes/cfd_variables.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

#include <utility>
#include <vector>

namespace Kratos
{

template <unsigned int TNumNodes>
class FluidBodyFlowCalculator : public ContributionCalculator<TNumNodes>
{
public:
    struct InputProvider {
        InputProvider(Geo::PropertiesGetter                GetElementProperties,
                      Geo::RetentionLawsGetter             GetRetentionLaws,
                      std::function<Matrix()>              GetMaterialPermeability,
                      std::function<Vector()>              GetIntegrationCoefficients,
                      std::function<std::vector<Vector>()> GetProjectedGravityAtIntegrationPoints,
                      std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients,
                      std::function<std::vector<double>()> GetFluidPressures)

            : GetElementProperties(std::move(GetElementProperties)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetMaterialPermeability(std::move(GetMaterialPermeability)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetProjectedGravityAtIntegrationPoints(std::move(GetProjectedGravityAtIntegrationPoints)),
              GetShapeFunctionGradients(std::move(GetShapeFunctionGradients)),
              GetFluidPressures(std::move(GetFluidPressures))
        {
        }

        Geo::PropertiesGetter                GetElementProperties;
        Geo::RetentionLawsGetter             GetRetentionLaws;
        std::function<Matrix()>              GetMaterialPermeability;
        std::function<Vector()>              GetIntegrationCoefficients;
        std::function<std::vector<Vector>()> GetProjectedGravityAtIntegrationPoints;
        std::function<Geometry<Node>::ShapeFunctionsGradientsType()> GetShapeFunctionGradients;
        std::function<std::vector<double>()>                         GetFluidPressures;
    };

    explicit FluidBodyFlowCalculator(InputProvider AnInputProvider)
        : mInputProvider(std::move(AnInputProvider))
    {
    }

    std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> LHSContribution() override
    {
        return std::nullopt;
    }

    BoundedVector<double, TNumNodes> RHSContribution() override
    {
        const auto& r_shape_function_gradients = mInputProvider.GetShapeFunctionGradients();

        const auto& r_integration_coefficients = mInputProvider.GetIntegrationCoefficients();

        const auto& r_properties          = mInputProvider.GetElementProperties();
        const auto  material_permeability = mInputProvider.GetMaterialPermeability();

        RetentionLaw::Parameters retention_law_parameters(r_properties);
        const auto&              r_projected_gravity_on_integration_points =
            mInputProvider.GetProjectedGravityAtIntegrationPoints();
        const auto&                      r_fluid_pressures = mInputProvider.GetFluidPressures();
        BoundedVector<double, TNumNodes> result            = ZeroVector(TNumNodes);
        const auto bishop_coefficients = CalculateBishopCoefficients(r_fluid_pressures);

        for (unsigned int integration_point_index = 0;
             integration_point_index < r_integration_coefficients.size(); ++integration_point_index) {
            retention_law_parameters.SetFluidPressure(r_fluid_pressures[integration_point_index]);
            const auto relative_permeability =
                mInputProvider.GetRetentionLaws()[integration_point_index]->CalculateRelativePermeability(
                    retention_law_parameters);
            noalias(result) +=
                r_properties[DENSITY_WATER] * bishop_coefficients[integration_point_index] * relative_permeability *
                prod(prod(r_shape_function_gradients[integration_point_index], material_permeability),
                     r_projected_gravity_on_integration_points[integration_point_index]) *
                r_integration_coefficients[integration_point_index] / r_properties[DYNAMIC_VISCOSITY];
        }
        return result;
    }

    std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> LocalSystemContribution() override
    {
        return {LHSContribution(), RHSContribution()};
    }

private:
    InputProvider mInputProvider;

    [[nodiscard]] std::vector<double> CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const
    {
        KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mInputProvider.GetRetentionLaws().size());

        auto retention_law_params = RetentionLaw::Parameters{mInputProvider.GetElementProperties()};

        auto result = std::vector<double>{};
        result.reserve(rFluidPressures.size());
        std::ranges::transform(mInputProvider.GetRetentionLaws(), rFluidPressures, std::back_inserter(result),
                               [&retention_law_params](const auto& pRetentionLaw, auto FluidPressure) {
            retention_law_params.SetFluidPressure(FluidPressure);
            return pRetentionLaw->CalculateBishopCoefficient(retention_law_params);
        });
        return result;
    }
};

} // namespace Kratos
