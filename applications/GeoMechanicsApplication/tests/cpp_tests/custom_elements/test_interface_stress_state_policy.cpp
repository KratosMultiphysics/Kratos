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

#include "custom_elements/interface_stress_state.h"
#include "custom_geometries/line_interface_geometry.h"
#include "custom_utilities/registration_utilities.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <string>
#include <type_traits>

using namespace Kratos;
using namespace std::string_literals;

namespace
{

auto CreateThreePlusThree2DLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 2.5, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 2.5, 0.0, 0.0));
    return LineInterfaceGeometry<Line2D3<Node>>{1, nodes};
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_CannotBeCopiedButItCanBeMoved, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    EXPECT_FALSE(std::is_copy_constructible_v<InterfaceStressState>);
    EXPECT_FALSE(std::is_copy_assignable_v<InterfaceStressState>);
    EXPECT_TRUE(std::is_move_constructible_v<InterfaceStressState>);
    EXPECT_TRUE(std::is_move_assignable_v<InterfaceStressState>);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_CloneCreatesCorrectInstance, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<InterfaceStressState>();

    const auto p_cloned_policy = p_stress_state_policy->Clone();
    KRATOS_EXPECT_NE(dynamic_cast<InterfaceStressState*>(p_cloned_policy.get()), nullptr);
    KRATOS_EXPECT_NE(p_cloned_policy.get(), p_stress_state_policy.get());
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_ThrowsWhenInputtingEmptyShapeFunctionValues,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = InterfaceStressState{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto b_matrix = stress_state_policy.CalculateBMatrix({}, {}, {}),
        "Shape function values are empty. Therefore, the B matrix can not be computed.\n");
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_ThrowsWhenNumberOfShapeFunctionsIsNotEqualToNumberOfNodePairs,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = InterfaceStressState{};
    const auto geometry            = CreateThreePlusThree2DLineInterfaceGeometry();

    Vector shape_function_values(2);
    shape_function_values <<= 1.0, 2.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto b_matrix =
            stress_state_policy.CalculateBMatrix({}, shape_function_values, geometry),
        "The number of shape functions should be equal to the number of node pairs. Therefore, "
        "the B matrix can not be computed.\n");
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_ReturnsExpectedVoigtSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = InterfaceStressState{};

    KRATOS_EXPECT_EQ(stress_state_policy.GetVoigtSize(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_ReturnsExpectedVoigtVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = InterfaceStressState{};

    const auto& voigt_vector = stress_state_policy.GetVoigtVector();

    Vector expected_voigt_vector(2);
    expected_voigt_vector <<= 1.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(expected_voigt_vector, voigt_vector, 1.0e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_ReturnsCorrectBMatrixForThreePlusThreeNodesGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = InterfaceStressState{};
    const auto geometry            = CreateThreePlusThree2DLineInterfaceGeometry();

    Vector shape_function_values(3);
    shape_function_values <<= -0.125, 0.375, 0.75; // Shape function values for xi = 0.5

    const auto b_matrix = stress_state_policy.CalculateBMatrix({}, shape_function_values, geometry);

    // clang-format off
    Matrix expected_b_matrix(2, 12);
    expected_b_matrix <<= 0,   0.125, 0, -0.375, 0, -0.75, 0, -0.125, 0, 0.375, 0, 0.75,
                          0.125, 0, -0.375, 0, -0.75, 0, -0.125, 0, 0.375, 0, 0.75, 0;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(b_matrix, expected_b_matrix, 1.0e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_ReturnsCorrectIntegrationCoefficient, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto interface_stress_state = InterfaceStressState{};

    const Geometry<Node>::IntegrationPointType integration_point(0.5, 0.3, 0.0, 0.5);

    constexpr auto detJ                   = 2.0;
    const auto     calculated_coefficient = interface_stress_state.CalculateIntegrationCoefficient(
        integration_point, detJ, LineInterfaceGeometry<Line2D3<Node>>());

    // The expected number is calculated as follows:
    // weight * detJ = 0.5 * 2.0 = 1.0
    KRATOS_EXPECT_NEAR(calculated_coefficient, 1.0, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_Throws_WhenAskingForStrain, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = InterfaceStressState{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto strain = stress_state_policy.CalculateGreenLagrangeStrain({}),
        "For interfaces, it is not possible to calculate the Green-Lagrange "
        "strain based on a deformation gradient.")
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_Throws_WhenAskingForStressTensorSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = InterfaceStressState{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto strain = stress_state_policy.GetStressTensorSize(),
        "For interfaces, the stress tensor size is not implemented.")
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceStressState_CanBeSavedAndLoadedThroughInterface, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("InterfaceStressState"s, InterfaceStressState{})};
    const auto p_policy = std::unique_ptr<StressStatePolicy>{std::make_unique<InterfaceStressState>()};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_policy);
    auto p_loaded_policy = std::unique_ptr<StressStatePolicy>{};
    serializer.load("test_tag"s, p_loaded_policy);

    // Assert
    ASSERT_NE(p_loaded_policy, nullptr);
    KRATOS_EXPECT_EQ(p_loaded_policy->GetVoigtSize(), VOIGT_SIZE_2D_INTERFACE);
    auto expected_voigt_vector = Vector{2};
    expected_voigt_vector <<= 1.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(p_loaded_policy->GetVoigtVector(), expected_voigt_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
