//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes
#include <functional>
#include <cmath>

// External includes

// Project includes
#include "includes/checks.h"
#include "testing/testing.h"

// Application includes
#include "custom_elements/data_containers/k_omega_sst/element_data_utilities.h"
#include "custom_elements/data_containers/k_omega_sst/element_data_derivative_utilities.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTUtilitiesCalculateTanhDerivative, KratosRansFastSuite)
{
    KRATOS_TRY

    const auto& f_h = [](const double x) { return x * x; };
    const auto& f_h_derivative = [](const double x) { return 2.0 * x; };

    double x = 1.2;

    // calculate reference values
    const double h = f_h(x);
    const double tanh_ref = std::tanh(h);

    // calculate analytical derivatives
    const double h_derivative = f_h_derivative(x);

    const double tanh_analytical_derivative = KOmegaSSTElementData::AdjointUtilities<3>::CalculateTanhDerivative(h, h_derivative);

    const double delta = 1e-9;

    x += delta;

    const double tanh_finite_difference_derivative = (std::tanh(f_h(x)) - tanh_ref) / delta;

    KRATOS_CHECK_RELATIVE_NEAR(tanh_analytical_derivative, tanh_finite_difference_derivative, 1e-7);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTUtilitiesCalculateArg1Derivative, KratosRansFastSuite)
{
    KRATOS_TRY

    double beta_star = 2.1;

    const auto& f_tke = [](const double x) { return x * x * x; };
    const auto& f_tke_derivative = [](const double x) { return 3.0 * x * x; };

    const auto& f_omega = [](const double x) { return x * x; };
    const auto& f_omega_derivative = [](const double x) { return 2.0 * x; };

    const auto& f_y = [](const double x) { return 3.0 * x; };
    const auto& f_y_derivative = [](const double x) { return 3.0; };

    const auto& f_arg1 = [&](const double x) {
        const double tke = f_tke(x);
        const double omega = f_omega(x);
        const double y = f_y(x);

        return KOmegaSSTElementData::CalculateArg1(beta_star, tke, omega, y);
    };

    const auto& f_arg1_analytical_derivative = [&](const double x) {
        const double tke = f_tke(x);
        const double omega = f_omega(x);
        const double y = f_y(x);

        const double tke_derivative = f_tke_derivative(x);
        const double omega_derivative = f_omega_derivative(x);
        const double y_derivative = f_y_derivative(x);

        return KOmegaSSTElementData::AdjointUtilities<3>::CalculateArg1Derivative(
            beta_star, tke, tke_derivative, omega, omega_derivative, y, y_derivative);
    };

    const double x = 1.5;
    const double delta = 1e-9;

    const double finite_difference_derivative = (f_arg1(x + delta) - f_arg1(x)) / delta;
    const double analytical_derivative = f_arg1_analytical_derivative(x);

    KRATOS_CHECK_RELATIVE_NEAR(analytical_derivative, finite_difference_derivative, 1e-7);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTUtilitiesCalculateArg2Derivative, KratosRansFastSuite)
{
    KRATOS_TRY

    double nu = 2.1;

    const auto& f_omega = [](const double x) { return x * x; };
    const auto& f_omega_derivative = [](const double x) { return 2.0 * x; };

    const auto& f_y_2 = [](const double x) { return 3.0 * x; };
    const auto& f_y_2_derivative = [](const double x) { return 3.0; };

    const auto& f_arg2 = [&](const double x) {
        const double omega = f_omega(x);
        const double y_2 = f_y_2(x);

        return KOmegaSSTElementData::CalculateArg2(nu, omega, y_2);
    };

    const auto& f_arg2_analytical_derivative = [&](const double x) {
        const double omega = f_omega(x);
        const double y_2 = f_y_2(x);

        const double omega_derivative = f_omega_derivative(x);
        const double y_2_derivative = f_y_2_derivative(x);

        return KOmegaSSTElementData::AdjointUtilities<3>::CalculateArg2Derivative(nu, omega, omega_derivative, y_2, y_2_derivative);
    };

    const double x = 5.5;
    const double delta = 1e-8;

    const double finite_difference_derivative = (f_arg2(x + delta) - f_arg2(x)) / delta;
    const double analytical_derivative = f_arg2_analytical_derivative(x);

    KRATOS_CHECK_RELATIVE_NEAR(analytical_derivative, finite_difference_derivative, 1e-7);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTUtilitiesCalculateArg3Derivative, KratosRansFastSuite)
{
    KRATOS_TRY

    double sigma_omega_2 = 1.8;

    const auto& f_tke = [](const double x) { return x * x * x; };
    const auto& f_tke_derivative = [](const double x) { return 3.0 * x * x; };

    const auto& f_y_2 = [](const double x) { return 3.0 * x; };
    const auto& f_y_2_derivative = [](const double x) { return 3.0; };

    const auto& f_cross_diffusion = [](const double x) { return 3.0 * x + x * x; };
    const auto& f_cross_diffusion_derivative = [](const double x) { return 3.0 + 2.0 * x; };

    const auto& f_arg3 = [&](const double x) {
        const double tke = f_tke(x);
        const double cross_diffusion = f_cross_diffusion(x);
        const double y_2 = f_y_2(x);

        return KOmegaSSTElementData::CalculateArg3(sigma_omega_2, tke, cross_diffusion, y_2);
    };

    const auto& f_arg3_analytical_derivative = [&](const double x) {
        const double tke = f_tke(x);
        const double cross_diffusion = f_cross_diffusion(x);
        const double y_2 = f_y_2(x);

        const double tke_derivative = f_tke_derivative(x);
        const double cross_diffusion_derivative = f_cross_diffusion_derivative(x);
        const double y_2_derivative = f_y_2_derivative(x);

        return KOmegaSSTElementData::AdjointUtilities<3>::CalculateArg3Derivative(sigma_omega_2, tke, tke_derivative, cross_diffusion, cross_diffusion_derivative, y_2, y_2_derivative);
    };

    const double x = 5.5;
    const double delta = 1e-8;

    const double finite_difference_derivative = (f_arg3(x + delta) - f_arg3(x)) / delta;
    const double analytical_derivative = f_arg3_analytical_derivative(x);

    KRATOS_CHECK_RELATIVE_NEAR(analytical_derivative, finite_difference_derivative, 1e-7);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTUtilitiesCalculateBlendedPhiDerivative, KratosRansFastSuite)
{
    KRATOS_TRY

    const auto& f_phi1 = [](const double x) { return x * x; };
    const auto& f_phi1_derivative = [](const double x) { return 2.0 * x; };

    const auto& f_phi2 = [](const double x) { return 2.0 * x; };
    const auto& f_phi2_derivative = [](const double x) { return 2.0; };

    const auto& f_f1 = [](const double x) { return 2.0 * x + std::pow(x, 3); };
    const auto& f_f1_derivative = [](const double x) { return 2.0 + 3.0 * std::pow(x, 2); };

    double x = 1.3;

    // calculating reference values
    const double phi1_ref = f_phi1(x);
    const double phi2_ref = f_phi2(x);
    const double f1_ref = f_f1(x);

    // calculating analytical derivatives
    const double phi1_derivative_analytical = f_phi1_derivative(x);
    const double phi2_derivative_analytical = f_phi2_derivative(x);
    const double f1_derivative_analytical = f_f1_derivative(x);
    const double blended_phi_derivative_analytical = KOmegaSSTElementData::AdjointUtilities<2>::CalculateBlendedPhiDerivative(phi1_ref, phi1_derivative_analytical, phi2_ref, phi2_derivative_analytical, f1_ref, f1_derivative_analytical);

    // calculating finite difference derivatives

    const double blended_phi_ref = KOmegaSSTElementData::CalculateBlendedPhi(phi1_ref, phi2_ref, f1_ref);

    const double delta = 1e-8;
    x += delta;

    const double phi1 = f_phi1(x);
    const double phi2 = f_phi2(x);
    const double f1 = f_f1(x);

    const double blended_phi = KOmegaSSTElementData::CalculateBlendedPhi(phi1, phi2, f1);

    x -= delta;

    const double phi1_derivative_finite_difference = (phi1 - phi1_ref) / delta;
    const double phi2_derivative_finite_difference = (phi2 - phi2_ref) / delta;
    const double f1_derivative_finite_difference = (f1 - f1_ref) / delta;
    const double blended_phi_derivative_finite_difference = (blended_phi - blended_phi_ref) / delta;

    KRATOS_CHECK_RELATIVE_NEAR(phi1_derivative_finite_difference, phi1_derivative_analytical, 1e-7);
    KRATOS_CHECK_RELATIVE_NEAR(phi2_derivative_finite_difference, phi2_derivative_analytical, 1e-7);
    KRATOS_CHECK_RELATIVE_NEAR(f1_derivative_analytical, f1_derivative_finite_difference, 1e-7);
    KRATOS_CHECK_RELATIVE_NEAR(blended_phi_derivative_analytical, blended_phi_derivative_finite_difference, 1e-6);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTUtilitiesCalculateCrossDiffusionTermDerivative, KratosRansFastSuite)
{
    KRATOS_TRY

    const double sigma_omega_2 = 1.8;

    const auto& f_omega = [](const double x) { return x * x; };
    const auto& f_omega_derivative = [](const double x) { return 2.0 * x; };

    const auto& f_tke_gradient = [](const double x) { return array_1d<double, 3>{x, 2.0 * x, std::pow(x, 3)}; };
    const auto& f_tke_gradient_derivative = [](const double x) { return array_1d<double, 3>{1.0, 2.0, 3.0 * std::pow(x, 2)}; };

    const auto& f_omega_gradient = [](const double x) { return array_1d<double, 3>{2.0 * x, 2.0 * x * x, std::pow(x, 2)}; };
    const auto& f_omega_gradient_derivative = [](const double x) { return array_1d<double, 3>{2.0, 4.0 * x, 2.0 * x}; };

    double x = 1.6;

    // calculating reference values
    const double omega_ref = f_omega(x);
    const array_1d<double, 3>& tke_gradient_ref = f_tke_gradient(x);
    const array_1d<double, 3>& omega_gradient_ref = f_omega_gradient(x);

    const double cross_wind_diffusion_ref =
        KOmegaSSTElementData::CalculateCrossDiffusionTerm<3>(
            sigma_omega_2, omega_ref, tke_gradient_ref, omega_gradient_ref);

    // calculate analytical derivatives
    const double omega_analytical_derivative = f_omega_derivative(x);
    const array_1d<double, 3> tke_gradient_analytical_derivative = f_tke_gradient_derivative(x);
    const array_1d<double, 3> omega_gradient_analytical_derivative = f_omega_gradient_derivative(x);

    const double cross_wind_analytical_derivative =
        KOmegaSSTElementData::AdjointUtilities<3>::CalculateCrossDiffusionTermDerivative(
            sigma_omega_2, omega_ref, omega_analytical_derivative,
            tke_gradient_ref, tke_gradient_analytical_derivative,
            omega_gradient_ref, omega_gradient_analytical_derivative);

    const double delta = 1e-8;

    // calculate finite difference derivatives
    x += delta;

    const double omega = f_omega(x);
    const array_1d<double, 3> tke_gradient = f_tke_gradient(x);
    const array_1d<double, 3> omega_gradient = f_omega_gradient(x);

    const double cross_wind_diffusion =
        KOmegaSSTElementData::CalculateCrossDiffusionTerm<3>(
            sigma_omega_2, omega, tke_gradient, omega_gradient);

    const double cross_wind_diffusion_finite_difference_derivative = (cross_wind_diffusion - cross_wind_diffusion_ref) / delta;

    KRATOS_CHECK_RELATIVE_NEAR(cross_wind_analytical_derivative, cross_wind_diffusion_finite_difference_derivative, 1e-7);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTUtilitiesCalculateF1Derivative, KratosRansFastSuite)
{
    KRATOS_TRY

    double beta_star = 2.1;
    double sigma_omega_2 = 1.8;
    double nu = 2.1;

    const auto& f_tke = [](const double x) { return x * x * x; };
    const auto& f_tke_derivative = [](const double x) { return 3.0 * x * x; };

    const auto& f_omega = [](const double x) { return x * x; };
    const auto& f_omega_derivative = [](const double x) { return 2.0 * x; };

    const auto& f_y = [](const double x) { return 3.0 * x; };
    const auto& f_y_derivative = [](const double x) { return 3.0; };

    const auto& f_cross_diffusion = [](const double x) { return 3.0 * x + x * x; };
    const auto& f_cross_diffusion_derivative = [](const double x) { return 3.0 + 2.0 * x; };

    const auto& f_f1 = [&](const double x) {
        const double tke_ref = f_tke(x);
        const double omega_ref = f_omega(x);
        const double y_ref = f_y(x);
        const double cross_diffusion_ref = f_cross_diffusion(x);

        return KOmegaSSTElementData::CalculateF1(tke_ref, omega_ref, nu, y_ref, beta_star, cross_diffusion_ref, sigma_omega_2);
    };

    const auto& f_f1_analytical_derivative = [&](const double x) {
        const double tke = f_tke(x);
        const double omega = f_omega(x);
        const double y = f_y(x);
        const double cross_diffusion = f_cross_diffusion(x);

        const double tke_derivative = f_tke_derivative(x);
        const double omega_derivative = f_omega_derivative(x);
        const double y_derivative = f_y_derivative(x);
        const double cross_diffusion_derivative = f_cross_diffusion_derivative(x);

        return KOmegaSSTElementData::AdjointUtilities<3>::CalculateF1Derivative(
            tke, tke_derivative, omega, omega_derivative, nu, y, y_derivative,
            beta_star, cross_diffusion, cross_diffusion_derivative, sigma_omega_2);
    };

    const auto& check_analytical_derivative = [&](const double x, const double delta, const double tolerance) {
        const double f1_finite_difference_derivative = (
            f_f1(x + delta) - f_f1(x)) / delta;
        const double f1_analytical_derivative = f_f1_analytical_derivative(x);

        KRATOS_CHECK_RELATIVE_NEAR(f1_analytical_derivative, f1_finite_difference_derivative, tolerance);
    };

    const double delta = 1e-7;
    const double tolerance = 10.0;

    check_analytical_derivative(1.0, delta, tolerance);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTUtilitiesCalculateF2Derivative, KratosRansFastSuite)
{
    KRATOS_TRY

    double beta_star = 2.1;
    double nu = 2.1;

    const auto& f_tke = [](const double x) { return x * x * x; };
    const auto& f_tke_derivative = [](const double x) { return 3.0 * x * x; };

    const auto& f_omega = [](const double x) { return x * x; };
    const auto& f_omega_derivative = [](const double x) { return 2.0 * x; };

    const auto& f_y = [](const double x) { return 3.0 * x; };
    const auto& f_y_derivative = [](const double x) { return 3.0; };

    const auto& f_f1 = [&](const double x) {
        const double tke_ref = f_tke(x);
        const double omega_ref = f_omega(x);
        const double y_ref = f_y(x);

        return KOmegaSSTElementData::CalculateF2(tke_ref, omega_ref, nu, y_ref, beta_star);
    };

    const auto& f_f2_analytical_derivative = [&](const double x) {
        const double tke = f_tke(x);
        const double omega = f_omega(x);
        const double y = f_y(x);

        const double tke_derivative = f_tke_derivative(x);
        const double omega_derivative = f_omega_derivative(x);
        const double y_derivative = f_y_derivative(x);

        return KOmegaSSTElementData::AdjointUtilities<3>::CalculateF2Derivative(
            tke, tke_derivative, omega, omega_derivative, nu, y, y_derivative,
            beta_star);
    };

    const auto& check_analytical_derivative = [&](const double x, const double delta, const double tolerance) {
        const double f2_finite_difference_derivative = (
            f_f1(x + delta) - f_f1(x)) / delta;
        const double f2_analytical_derivative = f_f2_analytical_derivative(x);


        KRATOS_CHECK_RELATIVE_NEAR(f2_analytical_derivative, f2_finite_difference_derivative, tolerance);
    };

    const double delta = 1e-7;
    const double tolerance = 10.0;

    check_analytical_derivative(1e+2, delta, tolerance);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTUtilitiesCalculateTurbulentKinematicViscosityDerivative, KratosRansFastSuite)
{
    KRATOS_TRY

    double a1 = 2.1;

    const auto& f_tke = [](const double x) { return x * x * x; };
    const auto& f_tke_derivative = [](const double x) { return 3.0 * x * x; };

    const auto& f_omega = [](const double x) { return x * x; };
    const auto& f_omega_derivative = [](const double x) { return 2.0 * x; };

    const auto& f_u_norm = [](const double x) { return 4.0 * x; };
    const auto& f_u_norm_derivative = [](const double x) { return 4.0; };

    const auto& f_f2 = [](const double x) { return 3.0 * x; };
    const auto& f_f2_derivative = [](const double x) { return 3.0; };

    const auto& f_nut = [&](const double x) {
        const double tke = f_tke(x);
        const double omega = f_omega(x);
        const double u_norm = f_u_norm(x);
        const double f2 = f_f2(x);

        return KOmegaSSTElementData::CalculateTurbulentKinematicViscosity(tke, omega, u_norm, f2, a1);
    };

    const auto& f_nut_analytical_derivative = [&](const double x) {
        const double tke = f_tke(x);
        const double omega = f_omega(x);
        const double u_norm = f_u_norm(x);
        const double f2 = f_f2(x);

        const double tke_derivative = f_tke_derivative(x);
        const double omega_derivative = f_omega_derivative(x);
        const double u_norm_derivative = f_u_norm_derivative(x);
        const double f2_derivative = f_f2_derivative(x);

        return KOmegaSSTElementData::AdjointUtilities<3>::CalculateTurbulentKinematicViscosityDerivative(
            tke, tke_derivative, omega, omega_derivative, u_norm, u_norm_derivative, f2, f2_derivative, a1);
    };

    const auto& check_analytical_derivative = [&](const double x, const double delta, const double tolerance) {
        const double nut_finite_difference_derivative = (
            f_nut(x + delta) - f_nut(x)) / delta;
        const double nut_analytical_derivative = f_nut_analytical_derivative(x);

        KRATOS_CHECK_RELATIVE_NEAR(nut_analytical_derivative, nut_finite_difference_derivative, tolerance);
    };

    const double delta = 1e-7;
    const double tolerance = 10.0;

    check_analytical_derivative(1e+2, delta, tolerance);

    KRATOS_CATCH("");
}
} // namespace Testing
} // namespace Kratos