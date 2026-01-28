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

namespace
{

using namespace Kratos;

template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
UPCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider CreateUPInputProvider(
    const Matrix& b_matrix,
    const Vector& voigt_vector,
    double        integration_coefficient,
    double        biot_coefficient,
    double        bishop_coefficient,
    const Vector& fluid_pressures,
    const Matrix& rNpContainer,
    std::size_t   NumberOfIntegrationPoints = 1)
{
    auto get_b_matrices = [b_matrix, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, b_matrix);
    };
    auto get_integration_coefficients = [integration_coefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, integration_coefficient);
    };
    auto get_np_container      = [&rNpContainer]() -> const Matrix& { return rNpContainer; };
    auto get_biot_coefficients = [biot_coefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, biot_coefficient);
    };
    auto get_bishop_coefficients = [bishop_coefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, bishop_coefficient);
    };
    auto get_voigt_vector    = [voigt_vector]() { return voigt_vector; };
    auto get_fluid_pressures = [fluid_pressures]() { return fluid_pressures; };

    return typename UPCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_bishop_coefficients, get_fluid_pressures);
}

template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
PUCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider CreatePUInputProvider(
    const Matrix& b_matrix,
    const Vector& voigt_vector,
    double        integration_coefficient,
    double        biot_coefficient,
    double        degree_of_saturation,
    const Vector& velocities,
    double        velocity_coefficient,
    const Matrix& rNpContainer,
    std::size_t   NumberOfIntegrationPoints = 1)
{
    auto get_b_matrices               = [b_matrix, NumberOfIntegrationPoints]() { return std::vector(NumberOfIntegrationPoints, b_matrix); };
    auto get_integration_coefficients = [integration_coefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, integration_coefficient);
    };
    auto get_np_container          = [&rNpContainer]() -> const Matrix& { return rNpContainer; };
    auto get_biot_coefficients     = [biot_coefficient, NumberOfIntegrationPoints]() { return std::vector(NumberOfIntegrationPoints, biot_coefficient); };
    auto get_degrees_of_saturation = [degree_of_saturation, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, degree_of_saturation);
    };
    auto get_voigt_vector         = [voigt_vector]() { return voigt_vector; };
    auto get_velocities           = [velocities]() { return velocities; };
    auto get_velocity_coefficient = [velocity_coefficient]() { return velocity_coefficient; };

    return typename PUCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_degrees_of_saturation, get_velocities, get_velocity_coefficient);
}

} // namespace

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingMatrixContribution)
{
    // Given
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
    const auto   np_container            = UblasUtilities::CreateMatrix({{1.0, 2.0}});
    auto         fluid_pressures         = UblasUtilities::CreateVector({1.0, 2.0});

    const auto input_provider = CreateUPInputProvider<number_of_u_dof, number_of_pw_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, bishop_coefficient,
        fluid_pressures, np_container);

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator(input_provider);

    // When
    const auto calculated_coupling_matrix = coupling_calculator.LHSContribution().value();

    // Then
    // clang-format off
    // Expected result is checked by hand
    const auto expected_coupling_matrix = UblasUtilities::CreateMatrix(
    {{1.5, 3.0},
     {1.8, 3.6},
     {2.1, 4.2},
     {2.4, 4.8}});
    // clang-format on
    KRATOS_CHECK_MATRIX_NEAR(calculated_coupling_matrix, expected_coupling_matrix, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingVectorContribution)
{
    // Given
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
    const auto   np_container            = UblasUtilities::CreateMatrix({{1.0, 2.0}});
    auto         fluid_pressures         = UblasUtilities::CreateVector({1.0, 2.0});

    const auto input_provider = CreateUPInputProvider<number_of_u_dof, number_of_pw_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, bishop_coefficient,
        fluid_pressures, np_container);

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator(input_provider);

    // When
    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();

    // Then
    // Expected result is checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({7.5, 9.0, 10.5, 12.0});
    KRATOS_CHECK_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingVectorContributionMultipleIntegrationPoints)
{
    // Given
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
    // The Np container now has two rows since there are two integration points
    const auto     np_container    = UblasUtilities::CreateMatrix({{1.0, 2.0}, {1.0, 2.0}});
    const auto     fluid_pressures = UblasUtilities::CreateVector({1.0, 2.0});
    constexpr auto number_of_integration_points = std::size_t{2};

    const auto input_provider = CreateUPInputProvider<number_of_u_dof, number_of_pw_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, bishop_coefficient,
        fluid_pressures, np_container, number_of_integration_points);

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator(input_provider);

    // When
    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();

    // Then
    // Expected result is checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({15.0, 18.0, 21.0, 24.0});
    KRATOS_CHECK_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestPUCouplingMatrixContribution)
{
    // Given
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
    const auto   velocity_coefficient    = 2.0;
    const auto   velocities              = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0});
    const auto   np_container            = UblasUtilities::CreateMatrix({{1.0, 2.0}});

    const auto input_provider = CreatePUInputProvider<number_of_pw_dof, number_of_u_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, degree_of_saturation,
        velocities, velocity_coefficient, np_container);

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof> coupling_calculator(input_provider);

    // When
    const auto calculated_coupling_matrix = coupling_calculator.LHSContribution().value();

    // Then
    // clang-format off
    // Expected result is checked by hand
    const auto expected_coupling_matrix = UblasUtilities::CreateMatrix(
    {{3.0, 3.6, 4.2, 4.8},
     {6.0, 7.2, 8.4, 9.6}});
    // clang-format on
    KRATOS_CHECK_MATRIX_NEAR(calculated_coupling_matrix, expected_coupling_matrix, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestPUCouplingVectorContribution)
{
    // Given
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
    const auto   velocity_coefficient    = 2.0;
    const auto   velocities              = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0});
    const auto   np_container            = UblasUtilities::CreateMatrix({{1.0, 2.0}});

    const auto input_provider = CreatePUInputProvider<number_of_pw_dof, number_of_u_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, degree_of_saturation,
        velocities, velocity_coefficient, np_container);

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof> coupling_calculator(input_provider);

    // When
    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();

    // Then
    // Expected result is checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({42.0, 84.0});
    KRATOS_CHECK_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestPUCouplingVectorContributionMultipleIntegrationPoints)
{
    // Given
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
    const auto   velocity_coefficient    = 2.0;
    const auto   velocities              = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0});
    const auto np_container          = UblasUtilities::CreateMatrix({{1.0, 2.0}, {1.0, 2.0}});
    const auto number_of_integration_points = std::size_t{2};
    const auto input_provider = CreatePUInputProvider<number_of_pw_dof, number_of_u_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, degree_of_saturation,
        velocities, velocity_coefficient, np_container, number_of_integration_points);

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof> coupling_calculator(input_provider);

    // When
    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();

    // Then
    // Expected result is checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({84.0, 168.0});
    KRATOS_CHECK_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
