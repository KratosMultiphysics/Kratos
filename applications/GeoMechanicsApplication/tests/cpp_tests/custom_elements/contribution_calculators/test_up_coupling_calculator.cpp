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

#include "custom_elements/contribution_calculators/up_coupling_calculator.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingMatrixContribution)
{
    constexpr std::size_t size_D = 3;
    constexpr std::size_t size_N = 4;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{ 1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 11.0, 12.0},
                 {11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0},
                 {21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0},
                 {31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0},
                 {41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0},
                 {51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0}});

    // clang-format on

    const auto   voigt_vector = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0, 0.0, 0.0});
    Vector       n_p          = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0});
    const double integration_coefficient = 1.0;
    const double biot_coefficient        = 0.02;
    const double bishop_coefficient      = 0.1;

    BoundedMatrix<double, size_D * size_N, size_N> expected_coupling_matrix;
    // clang-format off
    expected_coupling_matrix = UblasUtilities::CreateMatrix(
    {{0.066,0.132,0.198,0.264},
     {0.072,0.144,0.216,0.288},
     {0.078,0.156,0.234,0.312},
     {0.084,0.168,0.252,0.336},
     {0.090,0.180,0.270,0.360},
     {0.096,0.192,0.288,0.384},
     {0.102,0.204,0.306,0.408},
     {0.108,0.216,0.324,0.432},
     {0.114,0.228,0.342,0.456},
     {0.120,0.240,0.360,0.480},
     {0.126,0.252,0.378,0.504},
     {0.132,0.264,0.396,0.528}});
    // clang-format on

    constexpr auto number_of_u_dof  = 12;
    constexpr auto number_of_pw_dof = 4;

    auto get_b_matrices               = [b_matrix]() { return std::vector{b_matrix}; };
    auto get_integration_coefficients = [integration_coefficient]() {
        return std::vector{integration_coefficient};
    };
    const auto np_container     = UblasUtilities::CreateMatrix({{1.0, 2.0, 3.0, 4.0}});
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
    Vector       n_p                     = UblasUtilities::CreateVector({1.0, 2.0});
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
    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto   voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    Vector       n_p                     = UblasUtilities::CreateVector({1.0, 2.0});
    const double integration_coefficient = 0.5;
    const double biot_coefficient        = 2.0;
    const double bishop_coefficient      = 0.1;

    // Checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({7.5, 9.0, 10.5, 12.0});

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

    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();
    KRATOS_CHECK_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
