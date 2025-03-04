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
#include "permeability_calculator.h"
#include "custom_retention/retention_law.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "includes/cfd_variables.h"

namespace Kratos
{

PermeabilityCalculator::PermeabilityCalculator(InputProvider InputProvider)
    : mInputProvider{std::move(InputProvider)}
{
}

Matrix PermeabilityCalculator::LHSContribution() { return CalculatePermeabilityMatrix(); }

Vector PermeabilityCalculator::RHSContribution()
{
    return RHSContribution(CalculatePermeabilityMatrix());
}

std::pair<Matrix, Vector> PermeabilityCalculator::LocalSystemContribution()
{
    const auto permeability_matrix = CalculatePermeabilityMatrix();
    return {permeability_matrix, RHSContribution(permeability_matrix)};
}

Vector PermeabilityCalculator::RHSContribution(const Matrix& rPermeabilityMatrix) const
{
    return -prod(rPermeabilityMatrix, mInputProvider.GetNodalValues(WATER_PRESSURE));
}

Matrix PermeabilityCalculator::CalculatePermeabilityMatrix() const
{
    RetentionLaw::Parameters retention_parameters(mInputProvider.GetElementProperties());
    const auto&              r_properties             = mInputProvider.GetElementProperties();
    auto                     integration_coefficients = mInputProvider.GetIntegrationCoefficients();
    const auto               shape_function_gradients = mInputProvider.GetShapeFunctionGradients();
    const auto               local_dimension          = shape_function_gradients[0].size2();
    const Matrix             constitutive_matrix =
        GeoElementUtilities::FillPermeabilityMatrix(r_properties, local_dimension);

    const auto   number_of_nodes           = shape_function_gradients[0].size1();
    auto         result                    = Matrix{ZeroMatrix{number_of_nodes, number_of_nodes}};
    const double dynamic_viscosity_inverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
    for (unsigned int integration_point_index = 0;
         integration_point_index < integration_coefficients.size(); ++integration_point_index) {
        const double relative_permeability =
            mInputProvider.GetRetentionLaws()[integration_point_index]->CalculateRelativePermeability(retention_parameters);
        result += GeoTransportEquationUtilities::CalculatePermeabilityMatrix(
            shape_function_gradients[integration_point_index], dynamic_viscosity_inverse, constitutive_matrix,
            relative_permeability, integration_coefficients[integration_point_index]);
    }
    return result;
}

} // namespace Kratos
