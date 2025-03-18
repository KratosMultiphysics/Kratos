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
//
#include "fluid_body_flow_calculator.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"

#include <custom_utilities/element_utilities.hpp>

namespace Kratos
{

FluidBodyFlowCalculator::FluidBodyFlowCalculator(InputProvider AnInputProvider)
    : mInputProvider(std::move(AnInputProvider))
{
}

Matrix FluidBodyFlowCalculator::LHSContribution() { return {}; }

Vector FluidBodyFlowCalculator::RHSContribution()
{
    const auto shape_function_gradients = mInputProvider.GetShapeFunctionGradients();
    const Matrix& N_container = mInputProvider.GetNContainer();

    const auto integration_coefficients = mInputProvider.GetIntegrationCoefficients();

    const auto&                 r_properties = mInputProvider.GetElementProperties();
    const auto constitutive_matrix = GeoElementUtilities::FillPermeabilityMatrix(r_properties, 1);

    RetentionLaw::Parameters RetentionParameters(r_properties);

    const auto projected_gravity = mInputProvider.GetProjectedGravityForIntegrationPoints();

    Vector fluid_body_vector = ZeroVector(shape_function_gradients[0].size1());
    for (unsigned int integration_point_index = 0;
         integration_point_index < integration_coefficients.size(); ++integration_point_index) {
        const auto N = Vector{row(N_container, integration_point_index)};
        double     RelativePermeability =
            mInputProvider.GetRetentionLaws()[integration_point_index]->CalculateRelativePermeability(RetentionParameters);
        fluid_body_vector +=
            r_properties[DENSITY_WATER] * RelativePermeability *
            prod(prod(shape_function_gradients[integration_point_index], constitutive_matrix),
                 ScalarVector(1, projected_gravity[integration_point_index])) *
            integration_coefficients[integration_point_index] / r_properties[DYNAMIC_VISCOSITY];
    }
    return fluid_body_vector;
}

std::pair<Matrix, Vector> FluidBodyFlowCalculator::LocalSystemContribution() { return {{}, RHSContribution()}; }

} // namespace Kratos
