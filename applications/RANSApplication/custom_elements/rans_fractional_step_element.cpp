//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//  Extended by :    Suneth Warnakulasuriya
//

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "rans_fractional_step_element.h"

namespace Kratos
{
template <unsigned int TDim>
void RansFractionalStepElement<TDim>::CalculateLocalFractionalVelocitySystem(
    Matrix& rLeftHandSideMatrix,
    Vector& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();
    const SizeType local_size = TDim * number_of_nodes;

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size);
    }

    rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);

    if (rRightHandSideVector.size() != local_size) {
        rRightHandSideVector.resize(local_size);
    }

    rRightHandSideVector = ZeroVector(local_size);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    Matrix shape_functions;
    Vector gauss_weights;
    this->CalculateGeometryData(shape_function_derivatives, shape_functions, gauss_weights);
    const unsigned int number_of_gauss_points = gauss_weights.size();

    Matrix mass_matrix = ZeroMatrix(local_size, local_size);

    const double eta = rCurrentProcessInfo[PRESSURE_COEFFICIENT];

    // Stabilization parameters
    const double element_size = this->ElementSize();
    double tau_one, tau_two, density, mass_projection, old_pressure;

    array_1d<double, 3> body_force, momentum_projection, convective_velocity;

    Vector convection_operator(number_of_nodes);

    // Loop on integration points
    for (unsigned int g = 0; g < number_of_gauss_points; g++) {
        const double weight = gauss_weights[g];
        const Vector& r_shape_functions = row(shape_functions, g);
        const Matrix& r_shape_function_derivatives = shape_function_derivatives[g];

        EvaluateInPoint(r_geometry, r_shape_functions,
                        std::tie(density, DENSITY),
                        std::tie(mass_projection, DIVPROJ),
                        std::tie(body_force, BODY_FORCE),
                        std::tie(momentum_projection, CONV_PROJ),
                        std::tie(old_pressure, PRESSURE),
                        std::tie(convective_velocity, VELOCITY));

        const double viscosity = this->EffectiveViscosity(
            density, r_shape_functions, r_shape_function_derivatives,
            element_size, rCurrentProcessInfo);
        this->CalculateTau(tau_one, tau_two, element_size, convective_velocity,
                           density, viscosity, rCurrentProcessInfo);

        // Evaluate convection operator Velocity * Grad(N)
        this->ConvectionOperator(convection_operator, convective_velocity,
                                 r_shape_function_derivatives);

        // Add integration point contribution to the local mass matrix
        this->AddMomentumMassTerm(mass_matrix, r_shape_functions, weight * density);

        // Add convection, stabilization and RHS contributions to the local system equation
        this->AddMomentumSystemTerms(
            rLeftHandSideMatrix, rRightHandSideVector, density, convection_operator,
            body_force, old_pressure * eta, tau_one, tau_two, momentum_projection,
            mass_projection, r_shape_functions, r_shape_function_derivatives, weight);

        // Add viscous term
        const double effective_viscosity = viscosity * weight;
        this->AddViscousTerm(rLeftHandSideMatrix, r_shape_function_derivatives,
                             effective_viscosity);
    }

    // Add residual of previous iteration to RHS
    Vector last_velocity_values;
    this->GetVelocityValues(last_velocity_values, 0);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, last_velocity_values);

    // Add dynamic term
    const Vector& r_bdf_coefficients = rCurrentProcessInfo[BDF_COEFFICIENTS];
    noalias(rLeftHandSideMatrix) += r_bdf_coefficients[0] * mass_matrix;

    Vector time_term = r_bdf_coefficients[0] * last_velocity_values;
    for (SizeType i = 1; i < r_bdf_coefficients.size(); i++) {
        this->GetVelocityValues(last_velocity_values, i);
        noalias(time_term) += r_bdf_coefficients[i] * last_velocity_values;
    }

    noalias(rRightHandSideVector) -= prod(mass_matrix, time_term);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void RansFractionalStepElement<TDim>::CalculateLocalPressureSystem(
    Matrix& rLeftHandSideMatrix,
    Vector& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    using namespace RansCalculationUtilities;

    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != number_of_nodes) {
        rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes);
    }

    rLeftHandSideMatrix = ZeroMatrix(number_of_nodes, number_of_nodes);

    if (rRightHandSideVector.size() != number_of_nodes) {
        rRightHandSideVector.resize(number_of_nodes);
    }

    rRightHandSideVector = ZeroVector(number_of_nodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    Matrix shape_functions;
    Vector gauss_weights;
    this->CalculateGeometryData(shape_function_derivatives, shape_functions, gauss_weights);
    const unsigned int number_of_gauss_points = gauss_weights.size();

    // Stabilization parameters
    const double element_size = this->ElementSize();
    double tau_one, tau_two, density, velocity_divergence;

    array_1d<double, 3> body_force, momentum_projection, convective_velocity;
    array_1d<double, TDim> old_pressure_gradient;

    const double eta = rCurrentProcessInfo[PRESSURE_COEFFICIENT];

    // Loop on integration points
    for (unsigned int g = 0; g < number_of_gauss_points; g++) {
        const double weight = gauss_weights[g];
        const Vector& r_shape_functions = row(shape_functions, g);
        const Matrix& r_shape_function_derivatives = shape_function_derivatives[g];

        // Evaluate required variables at the integration point
        EvaluateInPoint(r_geometry, r_shape_functions,
                        std::tie(density, DENSITY),
                        std::tie(body_force, BODY_FORCE),
                        std::tie(momentum_projection, PRESS_PROJ),
                        std::tie(convective_velocity, VELOCITY));

        this->EvaluateGradientInPoint(old_pressure_gradient, PRESSURE,
                                      r_shape_function_derivatives);

        const double viscosity = this->EffectiveViscosity(
            density, r_shape_functions, r_shape_function_derivatives,
            element_size, rCurrentProcessInfo);
        this->CalculateTau(tau_one, tau_two, element_size, convective_velocity,
                           density, viscosity, rCurrentProcessInfo);

        this->EvaluateDivergenceInPoint(velocity_divergence, VELOCITY,
                                        r_shape_function_derivatives);

        // constant coefficient multiplying the pressure Laplacian (See Codina, Badia 2006 paper for details in case of a BDF2 time scheme)
        const double laplacian_coefficient =
            1.0 / (density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]);

        // Add convection, stabilization and RHS contributions to the local system equation
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            // LHS contribution
            for (SizeType j = 0; j < number_of_nodes; ++j) {
                double l_ij = 0.0;
                for (SizeType d = 0; d < TDim; ++d)
                    l_ij += r_shape_function_derivatives(i, d) *
                            r_shape_function_derivatives(j, d);
                l_ij *= (laplacian_coefficient + tau_one);

                rLeftHandSideMatrix(i, j) += weight * l_ij;
            }

            // RHS contribution

            // Velocity divergence
            double rhs_i = -r_shape_functions[i] * velocity_divergence;

            for (SizeType d = 0; d < TDim; ++d) {
                // Momentum stabilization
                rhs_i += r_shape_function_derivatives(i, d) * tau_one *
                         (density * body_force[d] - old_pressure_gradient[d] -
                          momentum_projection[d]);
                rhs_i += (eta - 1) * laplacian_coefficient *
                         r_shape_function_derivatives(i, d) * old_pressure_gradient[d];
            }

            rRightHandSideVector[i] += weight * rhs_i;
        }
    }
}

/*
 * Template class definition (this should allow us to compile the desired template instantiations)
 */

template class RansFractionalStepElement<2>;
template class RansFractionalStepElement<3>;

} // namespace Kratos
