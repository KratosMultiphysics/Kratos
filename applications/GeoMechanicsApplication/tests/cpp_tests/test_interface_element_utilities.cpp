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

#include "geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <custom_geometries/line_interface_geometry.h>
#include <custom_utilities/interface_element_utilities.h>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_RotationMatrixForHorizontalInterfaceIsIdentityMatrix,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 5.0, 0.0, 0.0));
    const LineInterfaceGeometry geometry(1, nodes);

    // Act
    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry, local_coordinate);

    // Assert
    KRATOS_EXPECT_MATRIX_NEAR(Matrix{IdentityMatrix{2}}, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationMatrixForInclinedInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.5 * std::sqrt(3), -0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.5 * std::sqrt(3), -0.5, 0.0));
    const LineInterfaceGeometry geometry(1, nodes);

    // Act
    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry, local_coordinate);

    // Assert
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5 * sqrt(3), 0.5, -0.5, 0.5 * sqrt(3); // Rotation of 30 degrees clockwise
    KRATOS_INFO("Expected rotation matrix") << expected_rotation_matrix << std::endl;
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
    KRATOS_INFO("End test") << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationMatrixForInclinedInterface2,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, -0.5, 0.5 * std::sqrt(3), 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -0.5, 0.5 * std::sqrt(3), 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    const LineInterfaceGeometry geometry(1, nodes);

    // Act
    // Since the gradient of the shape functions is constant, the rotation matrix is the same at
    // point, meaning the local_coordinate should not have an effect
    const auto local_coordinate = array_1d<double, 3>{0.5, 0.0, 0.0};
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry, local_coordinate);

    // Assert
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5, 0.5 * sqrt(3), -0.5 * sqrt(3), 0.5; // rotation of 60 degrees clockwise
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationMatrixForInclinedInterface_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, -0.5 * std::sqrt(3), 0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.5 * std::sqrt(3), -0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -0.5 * std::sqrt(3), 0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.5 * std::sqrt(3), -0.5, 0.0));
    const LineInterfaceGeometry geometry(1, nodes);

    // Act
    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry, local_coordinate);

    // Assert
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5 * sqrt(3), 0.5, -0.5, 0.5 * sqrt(3); // Rotation of 30 degrees clockwise
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationMatrixForVerticalElement_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, -7.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 10.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, -7.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 10.0, 0.0));
    const LineInterfaceGeometry geometry(1, nodes);

    // Act
    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry, local_coordinate);

    // Assert
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.0, -1.0, 1.0, 0.0; // Rotation of 90 degrees counterclockwise
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_UnityRotationFor3Plus3NodedInterface_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 2.0, 0.0, 0.0));
    const LineInterfaceGeometry geometry(1, nodes);

    // Act
    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry, local_coordinate);

    // Assert
    KRATOS_EXPECT_MATRIX_NEAR(Matrix{IdentityMatrix{2}}, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_UnityRotationForCurved3Plus3NodedInterface_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 1.0, -1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 1.0, -1.0, 0.0));
    const LineInterfaceGeometry geometry(1, nodes);

    // Act
    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry, local_coordinate);

    // Assert
    KRATOS_EXPECT_MATRIX_NEAR(Matrix{IdentityMatrix{2}}, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationForInclinedCurved3Plus3NodedInterface_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 2.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.8, 0.8, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 2.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 0.8, 0.8, 0.0));
    const LineInterfaceGeometry geometry(1, nodes);

    // Act
    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry, local_coordinate);

    // Assert
    // clang-format off
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5 * sqrt(2), 0.5 * sqrt(2),
                                -0.5 * sqrt(2), 0.5 * sqrt(2); // Rotation of 45 degrees clockwise
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationAtEdgeOfQuadraticElement_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, -1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, -1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 0.0, 0.0, 0.0));
    const LineInterfaceGeometry geometry(1, nodes);

    // Act
    const auto local_coordinate = array_1d<double, 3>{1.0, 0.0, 0.0};
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry, local_coordinate);

    // Assert
    // clang-format off
    Matrix expected_rotation_matrix(2, 2);
    // Rotation of 63.43 (atan(2)) degrees counterclockwise
    expected_rotation_matrix <<= std::cos(1.1071), -std::sin(1.1071),
                                 std::sin(1.1071), std::cos(1.1071);
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-3)
}

} // namespace Kratos::Testing