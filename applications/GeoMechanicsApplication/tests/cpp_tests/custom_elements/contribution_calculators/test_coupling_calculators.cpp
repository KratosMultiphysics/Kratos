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

#include "custom_elements/contribution_calculators/pu_coupling_calculator.hpp"
#include "custom_elements/contribution_calculators/up_coupling_calculator.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingMatrixContributionSmaller)
{
    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto   voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const double integration_coefficient = 0.5;
    const double biot_coefficient        = 2.0;
    const double bishop_coefficient      = 0.1;

    // clang-format off
    // Checked by hand
    const auto expected_coupling_matrix = UblasUtilities::CreateMatrix(
    {{1.5, 3.0},
     {1.8, 3.6},
     {2.1, 4.2},
     {2.4, 4.8}});
    // clang-format on

    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    auto get_b_matrices               = [b_matrix]() { return std::vector{b_matrix}; };
    auto get_integration_coefficients = [integration_coefficient]() {
        return std::vector{integration_coefficient};
    };
    const auto np_container     = UblasUtilities::CreateMatrix({{1.0, 2.0}});
    auto       get_np_container = [&np_container]() -> const Matrix& { return np_container; };
    auto get_biot_coefficients  = [biot_coefficient]() { return std::vector{biot_coefficient}; };
    auto get_bishop_coefficients = [bishop_coefficient]() { return std::vector{bishop_coefficient}; };
    auto get_voigt_vector    = [voigt_vector]() { return voigt_vector; };
    auto get_fluid_pressures = []() { return UblasUtilities::CreateVector({1.0, 2.0}); };

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof>::InputProvider input_provider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_bishop_coefficients, get_fluid_pressures);

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator(input_provider);

    const auto calculated_coupling_matrix = coupling_calculator.LHSContribution().value();
    KRATOS_CHECK_MATRIX_NEAR(calculated_coupling_matrix, expected_coupling_matrix, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingVectorContribution)
{
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto   voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const double integration_coefficient = 0.5;
    const double biot_coefficient        = 2.0;
    const double bishop_coefficient      = 0.1;

    // Checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({7.5, 9.0, 10.5, 12.0});

    auto get_b_matrices               = [b_matrix]() { return std::vector{b_matrix}; };
    auto get_integration_coefficients = [integration_coefficient]() {
        return std::vector{integration_coefficient};
    };
    const auto np_container     = UblasUtilities::CreateMatrix({{1.0, 2.0}});
    auto       get_np_container = [&np_container]() -> const Matrix& { return np_container; };
    auto get_biot_coefficients  = [biot_coefficient]() { return std::vector{biot_coefficient}; };
    auto get_bishop_coefficients = [bishop_coefficient]() { return std::vector{bishop_coefficient}; };
    auto get_voigt_vector    = [voigt_vector]() { return voigt_vector; };
    auto get_fluid_pressures = []() { return UblasUtilities::CreateVector({1.0, 2.0}); };

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof>::InputProvider input_provider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_bishop_coefficients, get_fluid_pressures);

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator(input_provider);

    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();
    KRATOS_CHECK_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingVectorContributionMultipleIntegrationPoints)
{
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto   voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const double integration_coefficient = 0.5;
    const double biot_coefficient        = 2.0;
    const double bishop_coefficient      = 0.1;

    // Checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({15.0, 18.0, 21.0, 24.0});

    auto get_b_matrices               = [b_matrix]() { return std::vector{b_matrix, b_matrix}; };
    auto get_integration_coefficients = [integration_coefficient]() {
        return std::vector{integration_coefficient, integration_coefficient};
    };
    const auto np_container          = UblasUtilities::CreateMatrix({{1.0, 2.0}, {1.0, 2.0}});
    auto       get_np_container      = [&np_container]() -> const Matrix& { return np_container; };
    auto       get_biot_coefficients = [biot_coefficient]() {
        return std::vector{biot_coefficient, biot_coefficient};
    };
    auto get_bishop_coefficients = [bishop_coefficient]() {
        return std::vector{bishop_coefficient, bishop_coefficient};
    };
    auto get_voigt_vector    = [voigt_vector]() { return voigt_vector; };
    auto get_fluid_pressures = []() { return UblasUtilities::CreateVector({1.0, 2.0}); };

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof>::InputProvider input_provider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_bishop_coefficients, get_fluid_pressures);

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator(input_provider);

    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();
    KRATOS_CHECK_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestPUCouplingMatrixContributionSmaller)
{
    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto   voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const double integration_coefficient = 0.5;
    const double biot_coefficient        = 2.0;
    const double degree_of_saturation    = 0.1;

    // clang-format off
    // Checked by hand
    const auto expected_coupling_matrix = UblasUtilities::CreateMatrix(
    {{3.0, 3.6, 4.2, 4.8},
     {6.0, 7.2, 8.4, 9.6}});
    // clang-format on

    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    auto get_b_matrices               = [b_matrix]() { return std::vector{b_matrix}; };
    auto get_integration_coefficients = [integration_coefficient]() {
        return std::vector{integration_coefficient};
    };
    const auto np_container        = UblasUtilities::CreateMatrix({{1.0, 2.0}});
    auto       get_np_container    = [&np_container]() -> const Matrix& { return np_container; };
    auto get_biot_coefficients     = [biot_coefficient]() { return std::vector{biot_coefficient}; };
    auto get_degrees_of_saturation = [degree_of_saturation]() {
        return std::vector{degree_of_saturation};
    };
    auto get_voigt_vector         = [voigt_vector]() { return voigt_vector; };
    auto get_fluid_pressures      = []() { return UblasUtilities::CreateVector({1.0, 2.0}); };
    auto get_velocity_coefficient = []() { return 2.0; };

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof>::InputProvider input_provider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_degrees_of_saturation, get_fluid_pressures, get_velocity_coefficient);

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof> coupling_calculator(input_provider);

    const auto calculated_coupling_matrix = coupling_calculator.LHSContribution().value();
    KRATOS_CHECK_MATRIX_NEAR(calculated_coupling_matrix, expected_coupling_matrix, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestPUCouplingVectorContribution)
{
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto   voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const double integration_coefficient = 0.5;
    const double biot_coefficient        = 2.0;
    const double degree_of_saturation    = 0.1;

    // Checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({42.0, 84.0});

    auto get_b_matrices               = [b_matrix]() { return std::vector{b_matrix}; };
    auto get_integration_coefficients = [integration_coefficient]() {
        return std::vector{integration_coefficient};
    };
    const auto np_container        = UblasUtilities::CreateMatrix({{1.0, 2.0}});
    auto       get_np_container    = [&np_container]() -> const Matrix& { return np_container; };
    auto get_biot_coefficients     = [biot_coefficient]() { return std::vector{biot_coefficient}; };
    auto get_degrees_of_saturation = [degree_of_saturation]() {
        return std::vector{degree_of_saturation};
    };
    auto get_voigt_vector = [voigt_vector]() { return voigt_vector; };
    auto get_velocities   = []() { return UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0}); };
    auto get_velocity_coefficient = []() { return 2.0; };

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof>::InputProvider input_provider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_degrees_of_saturation, get_velocities, get_velocity_coefficient);

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof> coupling_calculator(input_provider);

    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();
    KRATOS_CHECK_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestPUCouplingVectorContributionMultipleIntegrationPoints)
{
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto   voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    const double integration_coefficient = 0.5;
    const double biot_coefficient        = 2.0;
    const double degree_of_saturation    = 0.1;

    // Checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({84.0, 168.0});

    auto get_b_matrices               = [b_matrix]() { return std::vector{b_matrix, b_matrix}; };
    auto get_integration_coefficients = [integration_coefficient]() {
        return std::vector{integration_coefficient, integration_coefficient};
    };
    const auto np_container          = UblasUtilities::CreateMatrix({{1.0, 2.0}, {1.0, 2.0}});
    auto       get_np_container      = [&np_container]() -> const Matrix& { return np_container; };
    auto       get_biot_coefficients = [biot_coefficient]() {
        return std::vector{biot_coefficient, biot_coefficient};
    };
    auto get_degrees_of_saturation = [degree_of_saturation]() {
        return std::vector{degree_of_saturation, degree_of_saturation};
    };
    auto get_voigt_vector = [voigt_vector]() { return voigt_vector; };
    auto get_velocities   = []() { return UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0}); };
    auto get_velocity_coefficient = []() { return 2.0; };

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof>::InputProvider input_provider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_degrees_of_saturation, get_velocities, get_velocity_coefficient);

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof> coupling_calculator(input_provider);

    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();
    KRATOS_CHECK_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
