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
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

namespace
{

using namespace Kratos;

template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
typename UPCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider CreateUPInputProvider(
    const Matrix& rBMatrix,
    const Vector& rVoigtVector,
    double        IntegrationCoefficient,
    double        BiotCoefficient,
    double        BishopCoefficient,
    const Vector& rFluidPressures,
    const Matrix& rNpContainer,
    std::size_t   NumberOfIntegrationPoints = 1)
{
    auto get_b_matrices = [&rBMatrix, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, rBMatrix);
    };
    auto get_integration_coefficients = [IntegrationCoefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, IntegrationCoefficient);
    };
    auto get_np_container      = [&rNpContainer]() -> const auto& { return rNpContainer; };
    auto get_biot_coefficients = [BiotCoefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, BiotCoefficient);
    };
    auto get_bishop_coefficients = [BishopCoefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, BishopCoefficient);
    };
    auto get_voigt_vector    = [rVoigtVector]() { return rVoigtVector; };
    auto get_fluid_pressures = [rFluidPressures]() { return rFluidPressures; };

    return typename UPCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_bishop_coefficients, get_fluid_pressures);
}

template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
typename PUCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider CreatePUInputProvider(
    const Matrix& rBMatrix,
    const Vector& rVoigtVector,
    double        IntegrationCoefficient,
    double        BiotCoefficient,
    double        DegreeOfSaturation,
    const Vector& rVelocities,
    const Matrix& rNpContainer,
    std::size_t   NumberOfIntegrationPoints = 1)
{
    auto get_b_matrices = [&rBMatrix, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, rBMatrix);
    };
    auto get_integration_coefficients = [IntegrationCoefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, IntegrationCoefficient);
    };
    auto get_np_container      = [&rNpContainer]() -> const auto& { return rNpContainer; };
    auto get_biot_coefficients = [BiotCoefficient, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, BiotCoefficient);
    };
    auto get_degrees_of_saturation = [DegreeOfSaturation, NumberOfIntegrationPoints]() {
        return std::vector(NumberOfIntegrationPoints, DegreeOfSaturation);
    };
    auto get_voigt_vector = [rVoigtVector]() { return rVoigtVector; };
    auto get_velocities   = [rVelocities]() { return rVelocities; };

    return typename PUCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider(
        get_np_container, get_b_matrices, get_voigt_vector, get_integration_coefficients,
        get_biot_coefficients, get_degrees_of_saturation, get_velocities);
}

} // namespace

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingMatrixContribution)
{
    // Arrange
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto       voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    constexpr double integration_coefficient = 0.5;
    constexpr double biot_coefficient        = 2.0;
    constexpr double bishop_coefficient      = 0.1;
    const auto       np_container            = UblasUtilities::CreateMatrix({{1.0, 2.0}});
    auto             fluid_pressures         = UblasUtilities::CreateVector({1.0, 2.0});

    const auto input_provider = CreateUPInputProvider<number_of_u_dof, number_of_pw_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, bishop_coefficient,
        fluid_pressures, np_container);

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator(input_provider);

    // Act
    const auto calculated_coupling_matrix = coupling_calculator.LHSContribution().value();

    // Assert
    // clang-format off
    // Expected result is checked by hand
    const auto expected_coupling_matrix = UblasUtilities::CreateMatrix(
    {{1.5, 3.0},
     {1.8, 3.6},
     {2.1, 4.2},
     {2.4, 4.8}});
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(calculated_coupling_matrix, expected_coupling_matrix, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingVectorContribution)
{
    // Arrange
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto       voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    constexpr double integration_coefficient = 0.5;
    constexpr double biot_coefficient        = 2.0;
    constexpr double bishop_coefficient      = 0.1;
    const auto       np_container            = UblasUtilities::CreateMatrix({{1.0, 2.0}});
    auto             fluid_pressures         = UblasUtilities::CreateVector({1.0, 2.0});

    const auto input_provider = CreateUPInputProvider<number_of_u_dof, number_of_pw_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, bishop_coefficient,
        fluid_pressures, np_container);

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator(input_provider);

    // Act
    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();

    // Assert
    // Expected result is checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({7.5, 9.0, 10.5, 12.0});
    KRATOS_EXPECT_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestUPCouplingVectorContributionMultipleIntegrationPoints)
{
    // Arrange
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto       voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    constexpr double integration_coefficient = 0.5;
    constexpr double biot_coefficient        = 2.0;
    constexpr double bishop_coefficient      = 0.1;
    // The Np container now has two rows since there are two integration points
    const auto     np_container    = UblasUtilities::CreateMatrix({{1.0, 2.0}, {1.0, 2.0}});
    const auto     fluid_pressures = UblasUtilities::CreateVector({1.0, 2.0});
    constexpr auto number_of_integration_points = std::size_t{2};

    const auto input_provider = CreateUPInputProvider<number_of_u_dof, number_of_pw_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, bishop_coefficient,
        fluid_pressures, np_container, number_of_integration_points);

    UPCouplingCalculator<number_of_u_dof, number_of_pw_dof> coupling_calculator(input_provider);

    // Act
    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();

    // Assert
    // Expected result is checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({15.0, 18.0, 21.0, 24.0});
    KRATOS_EXPECT_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestPUCouplingMatrixContribution)
{
    // Arrange
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto       voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    constexpr double integration_coefficient = 0.5;
    constexpr double biot_coefficient        = 2.0;
    constexpr double degree_of_saturation    = 0.1;
    const auto       velocities              = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0});
    const auto       np_container            = UblasUtilities::CreateMatrix({{1.0, 2.0}});

    const auto input_provider = CreatePUInputProvider<number_of_pw_dof, number_of_u_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, degree_of_saturation,
        velocities, np_container);

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof> coupling_calculator(input_provider);

    // Act
    const auto calculated_coupling_matrix = coupling_calculator.LHSContribution().value();

    // Assert
    // clang-format off
    // Expected result is checked by hand
    const auto expected_coupling_matrix = UblasUtilities::CreateMatrix(
    {{1.5, 1.8, 2.1, 2.4},
     {3.0, 3.6, 4.2, 4.8}});
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(calculated_coupling_matrix, expected_coupling_matrix, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestPUCouplingVectorContribution)
{
    // Arrange
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto       voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    constexpr double integration_coefficient = 0.5;
    constexpr double biot_coefficient        = 2.0;
    constexpr double degree_of_saturation    = 0.1;
    const auto       velocities              = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0});
    const auto       np_container            = UblasUtilities::CreateMatrix({{1.0, 2.0}});

    const auto input_provider = CreatePUInputProvider<number_of_pw_dof, number_of_u_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, degree_of_saturation,
        velocities, np_container);

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof> coupling_calculator(input_provider);

    // Act
    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();

    // Assert
    // Expected result is checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({21.0, 42.0});
    KRATOS_EXPECT_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestPUCouplingVectorContributionMultipleIntegrationPoints)
{
    // Arrange
    constexpr auto number_of_u_dof  = 4;
    constexpr auto number_of_pw_dof = 2;

    // clang-format off
    const auto b_matrix = UblasUtilities::CreateMatrix(
                {{1.0,  2.0, 3.0, 4.0},
                 {5.0,  6.0, 7.0, 8.0},
                 {9.0, 10.0,11.0,12.0},
                 {13.0,14.0,15.0,16.0}});
    // clang-format on

    const auto       voigt_vector            = UblasUtilities::CreateVector({1.0, 1.0, 1.0, 0.0});
    constexpr double integration_coefficient = 0.5;
    constexpr double biot_coefficient        = 2.0;
    constexpr double degree_of_saturation    = 0.1;
    const auto       velocities              = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0});
    const auto       np_container = UblasUtilities::CreateMatrix({{1.0, 2.0}, {1.0, 2.0}});
    const auto       number_of_integration_points = std::size_t{2};
    const auto       input_provider = CreatePUInputProvider<number_of_pw_dof, number_of_u_dof>(
        b_matrix, voigt_vector, integration_coefficient, biot_coefficient, degree_of_saturation,
        velocities, np_container, number_of_integration_points);

    PUCouplingCalculator<number_of_pw_dof, number_of_u_dof> coupling_calculator(input_provider);

    // Act
    const auto calculated_coupling_vector = coupling_calculator.RHSContribution();

    // Assert
    // Expected result is checked by hand
    const auto expected_coupling_vector = UblasUtilities::CreateVector({42.0, 84.0});
    KRATOS_EXPECT_VECTOR_NEAR(calculated_coupling_vector, expected_coupling_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
