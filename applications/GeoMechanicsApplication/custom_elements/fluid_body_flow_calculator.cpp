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

#include "fluid_body_flow_calculator.h"
#include "custom_utilities/element_utilities.hpp"
#include "includes/cfd_variables.h"

namespace Kratos
{

FluidBodyFlowCalculator::FluidBodyFlowCalculator(InputProvider AnInputProvider)
    : mInputProvider(std::move(AnInputProvider))
{
}

std::optional<Matrix> FluidBodyFlowCalculator::LHSContribution() { return std::nullopt; }

std::optional<Vector> FluidBodyFlowCalculator::RHSContribution()
{
    const auto    shape_function_gradients = mInputProvider.GetShapeFunctionGradients();
    const Matrix& r_N_container              = mInputProvider.GetNContainer();

    const auto integration_coefficients = mInputProvider.GetIntegrationCoefficients();

    const auto& r_properties        = mInputProvider.GetElementProperties();
    const auto  constitutive_matrix = GeoElementUtilities::FillPermeabilityMatrix(
        r_properties, mInputProvider.GetLocalSpaceDimension());

    RetentionLaw::Parameters RetentionParameters(r_properties);
    const auto projected_gravity_on_integration_points = mInputProvider.GetProjectedGravityForIntegrationPoints();
    Vector fluid_body_vector = ZeroVector(shape_function_gradients[0].size1());
    for (unsigned int integration_point_index = 0;
         integration_point_index < integration_coefficients.size(); ++integration_point_index) {
        const auto N = Vector{row(r_N_container, integration_point_index)};
        const auto relative_permeability =
            mInputProvider.GetRetentionLaws()[integration_point_index]->CalculateRelativePermeability(RetentionParameters);
        noalias(fluid_body_vector) +=
            r_properties[DENSITY_WATER] * relative_permeability *
            prod(prod(shape_function_gradients[integration_point_index], constitutive_matrix),
                 ScalarVector(1, projected_gravity_on_integration_points[integration_point_index])) *
            integration_coefficients[integration_point_index] / r_properties[DYNAMIC_VISCOSITY];
    }
    return fluid_body_vector;
}

std::pair<std::optional<Matrix>, std::optional<Vector>> FluidBodyFlowCalculator::LocalSystemContribution()
{
    return {std::nullopt, RHSContribution()};
}

} // namespace Kratos
