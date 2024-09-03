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
#include <custom_utilities/interface_element_utilities.hpp>

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
    LineInterfaceGeometry geometry(1, nodes);

    // Act
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry);

    // Assert
    KRATOS_EXPECT_MATRIX_NEAR(Matrix{IdentityMatrix{2}}, rotation_matrix, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationMatrixForInclinedInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.5 * std::sqrt(3), -0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.5 * std::sqrt(3), -0.5, 0.0));
    LineInterfaceGeometry geometry(1, nodes);

    // Act
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry);

    // Assert
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5 * sqrt(3), -0.5, 0.5, 0.5 * sqrt(3); // Rotation of 30 degrees counterclockwise
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationMatrixForInclinedInterface2,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, -0.5, 0.5 * std::sqrt(3), 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -0.5, 0.5 * std::sqrt(3), 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    LineInterfaceGeometry geometry(1, nodes);

    // Act
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry);

    // Assert
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5, -0.5 * sqrt(3), 0.5 * sqrt(3), 0.5; // rotation of 60 degrees counterclockwise
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationMatrixForInclinedInterface_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, -0.5 * std::sqrt(3), 0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.5 * std::sqrt(3), -0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -0.5 * std::sqrt(3), 0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.5 * std::sqrt(3), -0.5, 0.0));
    LineInterfaceGeometry geometry(1, nodes);

    // Act
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry);

    // Assert
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.5 * sqrt(3), -0.5, 0.5, 0.5 * sqrt(3); // Rotation of 30 degrees counterclockwise
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_ReturnsCorrectRotationMatrixForVerticalElement_WithNonUnitLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, -7.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 10.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, -7.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 10.0, 0.0));
    LineInterfaceGeometry geometry(1, nodes);

    // Act
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix(geometry);

    // Assert
    Matrix expected_rotation_matrix(2, 2);
    expected_rotation_matrix <<= 0.0, 1.0, -1.0, 0.0; // Rotation of 270 degrees counterclockwise
    KRATOS_EXPECT_MATRIX_NEAR(expected_rotation_matrix, rotation_matrix, 1e-6);
}

} // namespace Kratos::Testing