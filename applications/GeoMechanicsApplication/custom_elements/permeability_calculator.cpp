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

template <unsigned int TNumNodes>
PermeabilityCalculator<TNumNodes>::PermeabilityCalculator(InputProvider InputProvider)
    : mInputProvider{std::move(InputProvider)}
{
}

template <unsigned int TNumNodes>
std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> PermeabilityCalculator<TNumNodes>::LHSContribution()
{
    return std::make_optional(CalculatePermeabilityMatrix());
}

template <unsigned int TNumNodes>
BoundedVector<double, TNumNodes> PermeabilityCalculator<TNumNodes>::RHSContribution()
{
    return RHSContribution(CalculatePermeabilityMatrix());
}

template <unsigned int TNumNodes>
std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> PermeabilityCalculator<TNumNodes>::LocalSystemContribution()
{
    const auto permeability_matrix = CalculatePermeabilityMatrix();
    return {std::make_optional(permeability_matrix), RHSContribution(permeability_matrix)};
}

template <unsigned int TNumNodes>
BoundedVector<double, TNumNodes> PermeabilityCalculator<TNumNodes>::RHSContribution(const Matrix& rPermeabilityMatrix) const
{
    return -prod(rPermeabilityMatrix, mInputProvider.GetNodalValues(WATER_PRESSURE));
}

template <unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes> PermeabilityCalculator<TNumNodes>::CalculatePermeabilityMatrix() const
{
    RetentionLaw::Parameters retention_parameters(mInputProvider.GetElementProperties());
    const auto&              r_properties             = mInputProvider.GetElementProperties();
    const auto&              integration_coefficients = mInputProvider.GetIntegrationCoefficients();
    const auto&              shape_function_gradients = mInputProvider.GetShapeFunctionGradients();
    const std::size_t        local_dimension          = shape_function_gradients[0].size2();
    const Matrix             constitutive_matrix =
        GeoElementUtilities::FillPermeabilityMatrix(r_properties, local_dimension);

    const std::size_t number_of_nodes = shape_function_gradients[0].size1();
    BoundedMatrix<double, TNumNodes, TNumNodes> result;

    const double dynamic_viscosity_inverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
    const auto&  fluid_pressures           = mInputProvider.GetFluidPressures();

    // Precompute scaling factors outside the loop if possible

    for (std::size_t integration_point_index = 0;
         integration_point_index < integration_coefficients.size(); ++integration_point_index) {
        retention_parameters.SetFluidPressure(fluid_pressures[integration_point_index]);
        const double relative_permeability =
            mInputProvider.GetRetentionLaws()[integration_point_index]->CalculateRelativePermeability(retention_parameters);

        // Use noalias to avoid temporaries in uBLAS
        BoundedMatrix<double, TNumNodes, TNumNodes> permeability_matrix;
        noalias(permeability_matrix) = GeoTransportEquationUtilities::CalculatePermeabilityMatrix(
            shape_function_gradients[integration_point_index], dynamic_viscosity_inverse, constitutive_matrix,
            relative_permeability, integration_coefficients[integration_point_index]);

        noalias(result) += permeability_matrix;
    }
    return result;
}

} // namespace Kratos
