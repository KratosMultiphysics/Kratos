//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/checks.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_stabilization_adjoint_utilities.h"
#include "custom_elements/convection_diffusion_reaction_stabilization_utilities.h"

namespace Kratos
{
namespace Testing
{
namespace ConvectionDiffusionReactionStabilizationAdjointUtilities
{
BoundedMatrix<double, 3, 3> CalculateInputMatrix(
    const double X,
    const double Y)
{
    BoundedMatrix<double, 3, 3> output;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            output(i, j) = (i + 1) * std::pow(X, j) * Y + (j + 1) * std::pow(Y, i);
        }
    }

    for (int i = 0; i < 3; ++i) {
        output(1, i) *= -1.0;
    }

    return output;
}

BoundedVector<BoundedMatrix<double, 3, 3>, 2> CalculateInputMatrixDerivatives(
    const double X,
    const double Y)
{
    BoundedVector<BoundedMatrix<double, 3, 3>, 2> output;
    auto& derivatives_0 = output[0];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            derivatives_0(i, j) = (i + 1) * j * std::pow(X, j - 1) * Y;
        }
    }

    for (int i = 0; i < 3; ++i) {
        derivatives_0(1, i) *= -1.0;
    }

    auto& derivatives_1 = output[1];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            derivatives_1(i, j) =
                (i + 1) * std::pow(X, j) + (j + 1) * i * std::pow(Y, i - 1);
        }
    }

    for (int i = 0; i < 3; ++i) {
        derivatives_1(1, i) *= -1.0;
    }

    return output;
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansAdjointStabilizedRFCCalculateDiscreteUpwindOperatorResidualContributionDerivatives, KratosRansFastSuite)
{
    using namespace ConvectionDiffusionReactionStabilizationAdjointUtilities;
    using adjoint_utilities = ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<2, 3>;

    BoundedVector<double, 3> input_values;
    for (int i = 0; i < 3; ++i) {
        input_values[i] = i+1;
    }

    const auto& calculate_residual = [&](const double X, const double Y) -> BoundedVector<double, 3> {
        const auto& matrix = CalculateInputMatrix(X, Y);

        BoundedMatrix<double, 3, 3> discrete_diffusion_matrix;
        double scalar_multiplier;
        ConvectionDiffusionReactionStabilizationUtilities::CalculateDiscreteUpwindOperator(
            scalar_multiplier, discrete_diffusion_matrix, matrix);

        return prod(discrete_diffusion_matrix, input_values);
    };

    double x_ref = 1.3;
    double y_ref = 1.5;

    const auto& values_ref = calculate_residual(x_ref, y_ref);
    const auto& input_matrix_ref = CalculateInputMatrix(x_ref, y_ref);
    const auto& input_matrix_derivatives = CalculateInputMatrixDerivatives(x_ref, y_ref);

    BoundedMatrix<double, 2, 3> adjoint_residual_derivatives;
    adjoint_utilities::CalculateDiscreteUpwindOperatorResidualContributionDerivatives(
        adjoint_residual_derivatives, input_values, input_matrix_ref, input_matrix_derivatives);

    // finite difference calculation
    const double delta = 1e-8;

    x_ref += delta;
    const BoundedVector<double, 3>& fd_x = (calculate_residual(x_ref, y_ref) - values_ref) / delta;
    x_ref -= delta;

    y_ref += delta;
    const BoundedVector<double, 3>& fd_y = (calculate_residual(x_ref, y_ref) - values_ref) / delta;
    y_ref -= delta;

    KRATOS_CHECK_VECTOR_NEAR(fd_x, row(adjoint_residual_derivatives, 0), 1e-6);
    KRATOS_CHECK_VECTOR_NEAR(fd_y, row(adjoint_residual_derivatives, 1), 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansAdjointStabilizedRFCCalculatePositivityPreservingMatrixDerivatives, KratosRansFastSuite)
{
    using namespace ConvectionDiffusionReactionStabilizationAdjointUtilities;
    using adjoint_utilities = ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<2, 3>;

    const auto& calculate_coefficient = [&](const double X, const double Y) -> double {
        const auto& matrix = CalculateInputMatrix(X, Y);
        BoundedMatrix<double, 3, 3> discrete_diffusion_matrix;
        return ConvectionDiffusionReactionStabilizationUtilities::CalculatePositivityPreservingMatrix(matrix);
    };

    double x_ref = 1.3;
    double y_ref = 1.5;

    const auto& values_ref = calculate_coefficient(x_ref, y_ref);
    const auto& input_matrix_ref = CalculateInputMatrix(x_ref, y_ref);
    const auto& input_matrix_derivatives = CalculateInputMatrixDerivatives(x_ref, y_ref);

    BoundedVector<double, 2> coefficient_derivatives;
    adjoint_utilities::CalculatePositivityPreservingCoefficientDerivatives(
        coefficient_derivatives, values_ref, input_matrix_ref, input_matrix_derivatives);

    // finite difference calculation
    const double delta = 1e-8;

    x_ref += delta;
    const double& fd_x = (calculate_coefficient(x_ref, y_ref) - values_ref) / delta;
    x_ref -= delta;

    y_ref += delta;
    const double& fd_y = (calculate_coefficient(x_ref, y_ref) - values_ref) / delta;
    y_ref -= delta;

    KRATOS_CHECK_NEAR(fd_x, coefficient_derivatives[0], 1e-6);
    KRATOS_CHECK_NEAR(fd_y, coefficient_derivatives[1], 1e-6);
}
} // namespace Testing
} // namespace Kratos