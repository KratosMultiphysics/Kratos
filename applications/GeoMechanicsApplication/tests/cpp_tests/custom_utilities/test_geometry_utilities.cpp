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

#include "custom_geometries/line_interface_geometry.h"
#include "custom_utilities/geometry_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_RotationMatrixForHorizontal2Plus2LineInterfaceIsIdentityMatrix,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 5.0, 0.0, 0.0));
    const LineInterfaceGeometry<Line2D2<Node>> geometry(1, nodes);

    // Note that only the first component of the local coordinate is used
    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    KRATOS_EXPECT_MATRIX_NEAR(Matrix{IdentityMatrix{2}}, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_ReturnsCorrectRotationMatrixForInclinedInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.5 * std::sqrt(3), -0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.5 * std::sqrt(3), -0.5, 0.0));
    const LineInterfaceGeometry<Line2D2<Node>> geometry(1, nodes);

    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    // clang-format off
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5 * sqrt(3), 0.5,
                                -0.5,           0.5 * sqrt(3); // Rotation of 30 degrees clockwise
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_ReturnsCorrectRotationMatrixForInclinedInterface2,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, -0.5, 0.5 * std::sqrt(3), 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -0.5, 0.5 * std::sqrt(3), 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    const LineInterfaceGeometry<Line2D2<Node>> geometry(1, nodes);

    // Since the gradient of the shape functions is constant, the rotation matrix is the same at
    // each point, meaning the local_coordinate should not have an effect
    const auto local_coordinate = array_1d<double, 3>{0.5, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    // clang-format off
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5,           0.5 * sqrt(3),
                                -0.5 * sqrt(3), 0.5; // rotation of 60 degrees clockwise
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_ReturnsCorrectRotationMatrixForInclinedInterface_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, -0.5 * std::sqrt(3), 0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.5 * std::sqrt(3), -0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -0.5 * std::sqrt(3), 0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.5 * std::sqrt(3), -0.5, 0.0));
    const LineInterfaceGeometry<Line2D2<Node>> geometry(1, nodes);

    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    // clang-format off
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5 * sqrt(3), 0.5,
                                -0.5,           0.5 * sqrt(3); // Rotation of 30 degrees clockwise
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_ReturnsCorrectRotationMatrixForVerticalElement_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, -7.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 10.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, -7.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 10.0, 0.0));
    const LineInterfaceGeometry<Line2D2<Node>> geometry(1, nodes);

    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    // clang-format off
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.0, -1.0,
                                 1.0,  0.0; // Rotation of 90 degrees counterclockwise
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_RotationMatrixForHorizontal3Plus3LineInterfaceIsIdentityMatrix,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 2.0, 0.0, 0.0));
    const LineInterfaceGeometry<Line2D3<Node>> geometry(1, nodes);

    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    KRATOS_EXPECT_MATRIX_NEAR(Matrix{IdentityMatrix{2}}, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_UnityRotationForCenterOfCurved3Plus3NodedInterface_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 1.0, -1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 1.0, -1.0, 0.0));
    const LineInterfaceGeometry<Line2D3<Node>> geometry(1, nodes);

    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    KRATOS_EXPECT_MATRIX_NEAR(Matrix{IdentityMatrix{2}}, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_ReturnsCorrectRotationForInclinedCurved3Plus3NodedInterface_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 2.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.8, 0.8, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 2.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 2.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 0.8, 0.8, 0.0));
    const LineInterfaceGeometry<Line2D3<Node>> geometry(1, nodes);

    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    // clang-format off
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5 * sqrt(2), 0.5 * sqrt(2),
                                -0.5 * sqrt(2), 0.5 * sqrt(2); // Rotation of 45 degrees clockwise
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_ReturnsCorrectRotationAtEdgeOfQuadraticElement_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, -1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, -1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 0.0, 0.0, 0.0));
    const LineInterfaceGeometry<Line2D3<Node>> geometry(1, nodes);

    const auto local_coordinate = array_1d<double, 3>{1.0, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    // clang-format off
    Matrix expected_rotation_matrix(2, 2);
    // Rotation of 63.43 (atan(2)) degrees counterclockwise
    expected_rotation_matrix <<= std::cos(1.1071), -std::sin(1.1071),
                                 std::sin(1.1071), std::cos(1.1071);
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-3)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_ReturnsCorrectRotationAtArbitraryXiForQuadraticElement_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, -1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, -1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 0.0, 0.0, 0.0));
    const LineInterfaceGeometry<Line2D3<Node>> geometry(1, nodes);

    const auto local_coordinate = array_1d<double, 3>{-0.5, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    // clang-format off
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5 * sqrt(2), 0.5 * sqrt(2),
                                -0.5 * sqrt(2), 0.5 * sqrt(2); // Rotation of 45 degrees clockwise
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-3)
}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_RotationMatrixForOpenHorizontalInterfaceIsIdentityMatrix,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    // Make sure the two sides of the interface are separated by a distance of 2
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, -1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 5.0, -1.0, 0.0));
    const LineInterfaceGeometry<Line2D2<Node>> geometry(1, nodes);

    const auto local_coordinate = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto rotation_matrix =
        GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(geometry, local_coordinate);

    KRATOS_EXPECT_MATRIX_NEAR(Matrix{IdentityMatrix{2}}, rotation_matrix, 1e-6)
}

} // namespace Kratos::Testing