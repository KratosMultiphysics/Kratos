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
#include "custom_geometries/interface_geometry.hpp"
#include "custom_utilities/registration_utilities.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <string>

using namespace Kratos;
using namespace std::string_literals;

namespace
{

auto CreateThreePlusThree2DLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 2.5, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(5, 5.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(6, 2.5, 0.0, 0.0));
    return InterfaceGeometry<Line2D3<Node>>{1, nodes};
}

auto CreateThreePlusThreeSurfaceInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 2.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 0.0, 2.0, 0.0));
    nodes.push_back(make_intrusive<Node>(4, 0.0, 0.0, 0.2));
    nodes.push_back(make_intrusive<Node>(5, 2.0, 0.0, 0.2));
    nodes.push_back(make_intrusive<Node>(6, 0.0, 2.0, 0.2));
    return InterfaceGeometry<Triangle3D3<Node>>{1, nodes};
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(Line2DInterfaceStressState_CloneCreatesCorrectInstance, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::unique_ptr<StressStatePolicy> p_stress_state_policy =
        std::make_unique<Line2DInterfaceStressState>();

    const auto p_cloned_policy = p_stress_state_policy->Clone();
    KRATOS_EXPECT_NE(dynamic_cast<Line2DInterfaceStressState*>(p_cloned_policy.get()), nullptr);
    KRATOS_EXPECT_NE(p_cloned_policy.get(), p_stress_state_policy.get());
}

KRATOS_TEST_CASE_IN_SUITE(Line2DInterfaceStressState_ThrowsWhenInputtingEmptyShapeFunctionValues,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = Line2DInterfaceStressState{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto b_matrix = stress_state_policy.CalculateBMatrix({}, {}, {}),
        "Shape function values are empty. Therefore, the B matrix can not be computed.\n");
}

KRATOS_TEST_CASE_IN_SUITE(Line2DInterfaceStressState_ThrowsWhenNumberOfShapeFunctionsIsNotEqualToNumberOfNodePairs,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = Line2DInterfaceStressState{};
    const auto geometry            = CreateThreePlusThree2DLineInterfaceGeometry();

    Vector shape_function_values(2);
    shape_function_values <<= 1.0, 0.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto b_matrix =
            stress_state_policy.CalculateBMatrix({}, shape_function_values, geometry),
        "The number of shape functions should be equal to the number of node pairs. Therefore, "
        "the B matrix can not be computed.\n");
}

KRATOS_TEST_CASE_IN_SUITE(Line2DInterfaceStressState_ReturnsExpectedVoigtSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = Line2DInterfaceStressState{};

    KRATOS_EXPECT_EQ(stress_state_policy.GetVoigtSize(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(Line2DInterfaceStressState_ReturnsExpectedVoigtVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = Line2DInterfaceStressState{};

    const auto& r_voigt_vector = stress_state_policy.GetVoigtVector();

    Vector expected_voigt_vector(2);
    expected_voigt_vector <<= 1.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(expected_voigt_vector, r_voigt_vector, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(Line2DInterfaceStressState_ReturnsCorrectBMatrixForThreePlusThreeNodesGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = Line2DInterfaceStressState{};
    const auto geometry            = CreateThreePlusThree2DLineInterfaceGeometry();

    Vector shape_function_values(3);
    shape_function_values <<= -0.125, 0.375, 0.75; // Shape function values for xi = 0.5

    const auto b_matrix = stress_state_policy.CalculateBMatrix({}, shape_function_values, geometry);

    // clang-format off
    Matrix expected_b_matrix(2, 12);
    expected_b_matrix <<= 0,   0.125, 0, -0.375, 0, -0.75, 0, -0.125, 0, 0.375, 0, 0.75,
                          0.125, 0, -0.375, 0, -0.75, 0, -0.125, 0, 0.375, 0, 0.75, 0;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(b_matrix, expected_b_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(Line2DInterfaceStressState_Throws_WhenAskingForStrain, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = Line2DInterfaceStressState{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto strain = stress_state_policy.CalculateGreenLagrangeStrain({}),
        "For line interfaces, it is not possible to calculate the Green-Lagrange "
        "strain.")
}

KRATOS_TEST_CASE_IN_SUITE(Line2DInterfaceStressState_Throws_WhenAskingForStressTensorSize,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = Line2DInterfaceStressState{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto size = stress_state_policy.GetStressTensorSize(),
        "For line interfaces, the stress tensor size is not implemented.")
}

KRATOS_TEST_CASE_IN_SUITE(Line2DInterfaceStressState_CanBeSavedAndLoadedThroughInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration = ScopedSerializerRegistration{
        std::make_pair("Line2DInterfaceStressState"s, Line2DInterfaceStressState{})};
    const auto p_policy =
        std::unique_ptr<StressStatePolicy>{std::make_unique<Line2DInterfaceStressState>()};
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

KRATOS_TEST_CASE_IN_SUITE(SurfaceInterfaceStressState_CloneCreatesCorrectInstance, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto p_stress_state_policy =
        std::unique_ptr<StressStatePolicy>{std::make_unique<SurfaceInterfaceStressState>()};

    const auto p_cloned_policy = p_stress_state_policy->Clone();
    KRATOS_EXPECT_NE(dynamic_cast<SurfaceInterfaceStressState*>(p_cloned_policy.get()), nullptr);
    KRATOS_EXPECT_NE(p_cloned_policy.get(), p_stress_state_policy.get());
}

KRATOS_TEST_CASE_IN_SUITE(SurfaceInterfaceStressState_ThrowsWhenInputtingEmptyShapeFunctionValues,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = SurfaceInterfaceStressState{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto b_matrix = stress_state_policy.CalculateBMatrix({}, {}, {}),
        "Shape function values are empty. Therefore, the B matrix can not be computed.\n");
}

KRATOS_TEST_CASE_IN_SUITE(SurfaceInterfaceStressState_ThrowsWhenNumberOfShapeFunctionsIsNotEqualToNumberOfNodePairs,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = SurfaceInterfaceStressState{};
    const auto geometry            = CreateThreePlusThreeSurfaceInterfaceGeometry();

    Vector shape_function_values(2);
    shape_function_values <<= 1.0, 0.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto b_matrix =
            stress_state_policy.CalculateBMatrix({}, shape_function_values, geometry),
        "The number of shape functions should be equal to the number of node pairs. Therefore, "
        "the B matrix can not be computed.\n");
}

KRATOS_TEST_CASE_IN_SUITE(SurfaceInterfaceStressState_ReturnsExpectedVoigtSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = SurfaceInterfaceStressState{};

    KRATOS_EXPECT_EQ(stress_state_policy.GetVoigtSize(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(SurfaceInterfaceStressState_ReturnsExpectedVoigtVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = SurfaceInterfaceStressState{};

    const auto& r_voigt_vector = stress_state_policy.GetVoigtVector();

    Vector expected_voigt_vector(3);
    expected_voigt_vector <<= 1.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(expected_voigt_vector, r_voigt_vector, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(SurfaceInterfaceStressState_ReturnsCorrectBMatrixForThreePlusThreeNodesGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = SurfaceInterfaceStressState{};
    const auto geometry            = CreateThreePlusThreeSurfaceInterfaceGeometry();

    Vector shape_function_values(3);
    shape_function_values <<= -0.125, 0.375, 0.75;

    const auto b_matrix = stress_state_policy.CalculateBMatrix({}, shape_function_values, geometry);

    // clang-format off
    Matrix expected_b_matrix(3, 18);
    expected_b_matrix <<= 0.0,   0.0,   0.125,  0.0,    0.0,   -0.375,  0.0,   0.0,  -0.75,  0.0,    0.0,   -0.125, 0.0,   0.0,   0.375, 0.0,  0.0,  0.75,
                          0.125, 0.0,   0.0,   -0.375,  0.0,    0.0,   -0.75,  0.0,   0.0,  -0.125,  0.0,    0.0,   0.375, 0.0,   0.0,   0.75, 0.0,  0.0,
                          0.0,   0.125, 0.0,    0.0,   -0.375,  0.0,    0.0,  -0.75,  0.0,   0.0,   -0.125,  0.0,   0.0,   0.375, 0.0,   0.0,  0.75, 0.0;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(b_matrix, expected_b_matrix, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(SurfaceInterfaceStressState_Throws_WhenAskingForStrain, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = SurfaceInterfaceStressState{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto strain = stress_state_policy.CalculateGreenLagrangeStrain({}),
        "For surface interfaces, it is not possible to calculate the Green-Lagrange "
        "strain.")
}

KRATOS_TEST_CASE_IN_SUITE(SurfaceInterfaceStressState_Throws_WhenAskingForStressTensorSize,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_state_policy = SurfaceInterfaceStressState{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto size = stress_state_policy.GetStressTensorSize(),
        "For surface interfaces, the stress tensor size is not implemented.")
}

KRATOS_TEST_CASE_IN_SUITE(SurfaceInterfaceStressState_CanBeSavedAndLoadedThroughInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration = ScopedSerializerRegistration{
        std::make_pair("SurfaceInterfaceStressState"s, SurfaceInterfaceStressState{})};
    const auto p_policy =
        std::unique_ptr<StressStatePolicy>{std::make_unique<SurfaceInterfaceStressState>()};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_policy);
    auto p_loaded_policy = std::unique_ptr<StressStatePolicy>{};
    serializer.load("test_tag"s, p_loaded_policy);

    // Assert
    ASSERT_NE(p_loaded_policy, nullptr);
    KRATOS_EXPECT_EQ(p_loaded_policy->GetVoigtSize(), VOIGT_SIZE_3D_INTERFACE);
    auto expected_voigt_vector = Vector{3};
    expected_voigt_vector <<= 1.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(p_loaded_policy->GetVoigtVector(), expected_voigt_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
