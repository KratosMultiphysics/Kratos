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

Vector FluidBodyFlowCalculator::RHSContribution()
{
    const auto shape_function_gradients = mInputProvider.GetShapeFunctionGradients();

    const auto integration_coefficients = mInputProvider.GetIntegrationCoefficients();

    const auto& r_properties        = mInputProvider.GetElementProperties();
    const auto  constitutive_matrix = GeoElementUtilities::FillPermeabilityMatrix(
        r_properties, mInputProvider.GetLocalSpaceDimension());

    RetentionLaw::Parameters retention_law_parameters(r_properties);
    const auto               projected_gravity_on_integration_points =
        mInputProvider.GetProjectedGravityAtIntegrationPoints();
    std::cout << " projected_gravity_on_integration_points " << projected_gravity_on_integration_points << std::endl;

    const auto fluid_body_vector_length = shape_function_gradients[0].size1();
    auto       result                   = Vector{ZeroVector(fluid_body_vector_length)};
    for (unsigned int integration_point_index = 0;
         integration_point_index < integration_coefficients.size(); ++integration_point_index) {
        const auto relative_permeability =
            mInputProvider.GetRetentionLaws()[integration_point_index]->CalculateRelativePermeability(
                retention_law_parameters);
        noalias(result) +=
            r_properties[DENSITY_WATER] * relative_permeability *
            prod(prod(shape_function_gradients[integration_point_index], constitutive_matrix),
                 projected_gravity_on_integration_points[integration_point_index]) *
            integration_coefficients[integration_point_index] / r_properties[DYNAMIC_VISCOSITY];
    }
    return result;
}

std::pair<std::optional<Matrix>, Vector> FluidBodyFlowCalculator::LocalSystemContribution()
{
    return {LHSContribution(), RHSContribution()};
}

} // namespace Kratos
