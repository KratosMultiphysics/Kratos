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

#include "custom_utilities/nodal_extrapolator.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/triangle_2d_3.h"
#include "testing/testing.h"
#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D3NTriangle, KratosGeoMechanicsFastSuite)
{
    // Create a triangle
    Kratos::Triangle2D3<Node> geometry(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
                                       Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
                                       Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0));

    NodalExtrapolator nodal_extrapolator;

    auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto integration_points = geometry.IntegrationPoints(integration_method);
    // Calculate the extrapolation matrix
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method).size(), integration_points,
        integration_method, geometry.size());

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(3, 3);
    expected_extrapolation_matrix <<= 5.0/3.0, -1.0/3.0, -1.0/3.0,
                                     -1.0/3.0,  5.0/3.0, -1.0/3.0,
                                     -1.0/3.0, -1.0/3.0,  5.0/3.0;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolator_GivesCorrectExtrapolationMatrix_For2D4NQuadrilateral, KratosGeoMechanicsFastSuite)
{
    Kratos::Quadrilateral2D4<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0), Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0), Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0));

    NodalExtrapolator nodal_extrapolator;

    auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    auto integration_points = geometry.IntegrationPoints(integration_method);
    // Calculate the extrapolation matrix
    auto extrapolation_matrix = nodal_extrapolator.CalculateElementExtrapolationMatrix(
        geometry, geometry.IntegrationPoints(integration_method).size(), integration_points,
        integration_method, geometry.size());

    // clang-format off
    Matrix expected_extrapolation_matrix = ZeroMatrix(4, 4);
    expected_extrapolation_matrix <<= 1.8660254037844386, -0.5,  0.13397459621556132, -0.5,
                                      -0.5, 1.8660254037844386, -0.5, 0.13397459621556132,
                                      0.13397459621556132, -0.5, 1.8660254037844386, -0.5,
                                      -0.5, 0.13397459621556132, -0.5, 1.8660254037844386;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(extrapolation_matrix, expected_extrapolation_matrix, 1e-12);
}

} // namespace Kratos::Testing