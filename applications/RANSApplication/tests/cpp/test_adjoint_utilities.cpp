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

// External includes

// Project includes
#include "includes/checks.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/rans_adjoint_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(RansAdjointUtilitiesCalculateYPlusAndUtauDerivative, KratosRansFastSuite)
{
    const double kappa = 0.41;
    const double beta = 5.2;
    const double nu = 2.3;

    const double y_plus_limit =
        RansCalculationUtilities::CalculateLogarithmicYPlusLimit(kappa, beta);

    const auto& u_method = [](const double x) -> double { return x * x; };
    const auto& u_derivative_method = [](const double x) -> double { return 2.0 * x; };

    const auto& y_method = [](const double x) -> double { return 3 * x + 4 * x * x; };
    const auto& y_derivative_method = [](const double x) -> double { return 3.0 + 8.0 * x; };

    const auto& y_plus_u_tau_method = [&](const double x) {
        const double u = u_method(x);
        const double y = y_method(x);

        double u_tau, y_plus;
        RansCalculationUtilities::CalculateYPlusAndUtau(y_plus, u_tau, u, y, nu, kappa, beta);

        return std::make_pair(y_plus, u_tau);
    };

    const auto& y_plus_u_tau_derivative_method = [&](const double x) {
        double y_plus, u_tau;
        std::tie(y_plus, u_tau) = y_plus_u_tau_method(x);

        double y_plus_derivative, u_tau_derivative;
        RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
            y_plus_derivative, u_tau_derivative, y_plus, u_tau, u_method(x),
            u_derivative_method(x), y_method(x), y_derivative_method(x), nu,
            kappa, beta, y_plus_limit);

        return std::make_pair(y_plus_derivative, u_tau_derivative);
    };

    double delta, x, y_plus_ref, u_tau_ref, y_plus, u_tau, y_plus_derivative,
        u_tau_derivative, fd_y_plus_derivative, fd_u_tau_derivative;

    // checking in the logarithmic region
    x = 4.3;
    std::tie(y_plus_ref, u_tau_ref) = y_plus_u_tau_method(x);
    std::tie(y_plus_derivative, u_tau_derivative) = y_plus_u_tau_derivative_method(x);

    delta = 1e-9;
    x += delta;
    std::tie(y_plus, u_tau) = y_plus_u_tau_method(x);

    // compute fd sensitivities
    fd_y_plus_derivative = (y_plus - y_plus_ref) / delta;
    fd_u_tau_derivative = (u_tau - u_tau_ref) / delta;

    KRATOS_CHECK(y_plus_ref > y_plus_limit); // checks whether it is in logarithmic region
    KRATOS_CHECK_RELATIVE_NEAR(fd_y_plus_derivative, y_plus_derivative, 1e-6);
    KRATOS_CHECK_RELATIVE_NEAR(fd_u_tau_derivative, u_tau_derivative, 1e-6);

    // checking for linear region
    x = 1.1;
    std::tie(y_plus_ref, u_tau_ref) = y_plus_u_tau_method(x);
    std::tie(y_plus_derivative, u_tau_derivative) = y_plus_u_tau_derivative_method(x);

    delta = 1e-9;
    x += delta;
    std::tie(y_plus, u_tau) = y_plus_u_tau_method(x);

    // compute fd sensitivities
    fd_y_plus_derivative = (y_plus - y_plus_ref) / delta;
    fd_u_tau_derivative = (u_tau - u_tau_ref) / delta;

    KRATOS_CHECK(y_plus_ref < y_plus_limit); // checks whether it is in linear region
    KRATOS_CHECK_RELATIVE_NEAR(fd_y_plus_derivative, y_plus_derivative, 1e-6);
    KRATOS_CHECK_RELATIVE_NEAR(fd_u_tau_derivative, u_tau_derivative, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansAdjointUtilitiesCalculateUnitVectorDerivative, KratosRansFastSuite)
{
    const auto& vector_method = [](const double x) -> array_1d<double, 3> {
        array_1d<double, 3> result;
        result[0] = x * x;
        result[1] = 2.0 * x;
        result[2] = 3.0 * x + x * x * x;
        return result;
    };

    const auto& vector_derivative = [](const double x) -> array_1d<double, 3> {
        array_1d<double, 3> result;
        result[0] = 2.0 * x;
        result[1] = 2.0;
        result[2] = 3.0 + 3.0 * x * x;
        return result;
    };

    double x = 2.0;
    array_1d<double, 3> u, u_ref, unit_u_fd_sensitivity, unit_u_analytical_sensitivity;
    u_ref = vector_method(x);
    double u_ref_magnitude = norm_2(u_ref);
    u_ref /= u_ref_magnitude;

    unit_u_analytical_sensitivity = RansAdjointUtilities::CalculateUnitVectorDerivative(
        u_ref_magnitude, u_ref, vector_derivative(x));

    double delta = 1e-8;
    x += delta;
    u = vector_method(x);
    u /= norm_2(u);

    unit_u_fd_sensitivity = (u - u_ref) / delta;

    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(unit_u_fd_sensitivity, unit_u_analytical_sensitivity, 1e-6);
}

} // namespace Testing

} // namespace Kratos