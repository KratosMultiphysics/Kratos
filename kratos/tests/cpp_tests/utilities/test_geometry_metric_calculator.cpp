//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/geometry_metric_calculator.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(MetricTensorDataEquilateralTriangle, KratosCoreFastSuite)
{
    // Set triangle geometry
    Geometry<Node>::PointsArrayType nodes;
    nodes.push_back(Node::Pointer(new Node(1, 0.0, 0.0, 0.0)));
    nodes.push_back(Node::Pointer(new Node(2, 1.0, 0.0, 0.0)));
    nodes.push_back(Node::Pointer(new Node(3, 0.5, 0.866025404, 0.0)));
    const auto p_triangle = Geometry<Node>::Pointer(new Triangle2D3<Node>(nodes));

    // Call the triangle metric calculator utility
    double h_ref, met_inf, met_sup;
    BoundedMatrix<double,2,2> metric_tensor;
    GeometryMetricCalculator<2,3>::CalculateMetricTensorDimensionless(*p_triangle, metric_tensor, h_ref, met_inf, met_sup);

    // Check results
    const double tolerance = 1.0e-8;
    KRATOS_CHECK_NEAR(h_ref, 1.0, tolerance);
    KRATOS_CHECK_NEAR(met_inf, 1.0, tolerance);
    KRATOS_CHECK_NEAR(met_sup, 1.0, tolerance);
    BoundedMatrix<double,2,2> expected_metric_tensor;
    expected_metric_tensor(0,0) = 1.0; expected_metric_tensor(0,1) = 0.0;
    expected_metric_tensor(1,0) = 0.0; expected_metric_tensor(1,1) = 1.0;
    KRATOS_CHECK_MATRIX_NEAR(metric_tensor, expected_metric_tensor, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MetricTensorDataUnitTriangle2D3N, KratosCoreFastSuite)
{
    // Set triangle geometry
    Geometry<Node>::PointsArrayType nodes;
    nodes.push_back(Node::Pointer(new Node(1, 0.0, 0.0, 0.0)));
    nodes.push_back(Node::Pointer(new Node(2, 1.0, 0.0, 0.0)));
    nodes.push_back(Node::Pointer(new Node(3, 0.0, 1.0, 0.0)));
    const auto p_triangle = Geometry<Node>::Pointer(new Triangle2D3<Node>(nodes));

    // Call the triangle metric calculator utility
    double h_ref, met_inf, met_sup;
    BoundedMatrix<double,2,2> metric_tensor;
    GeometryMetricCalculator<2,3>::CalculateMetricTensorDimensionless(*p_triangle, metric_tensor, h_ref, met_inf, met_sup);

    // Check results
    const double tolerance = 1.0e-5;
    KRATOS_CHECK_NEAR(h_ref, 1.115360, tolerance);
    KRATOS_CHECK_NEAR(met_inf, 0.5, tolerance);
    KRATOS_CHECK_NEAR(met_sup, 1.5, tolerance);
    BoundedMatrix<double,2,2> expected_metric_tensor;
    expected_metric_tensor(0,0) = 1.244020; expected_metric_tensor(0,1) = 0.622008;
    expected_metric_tensor(1,0) = 0.622008; expected_metric_tensor(1,1) = 1.244020;
    KRATOS_CHECK_MATRIX_NEAR(metric_tensor, expected_metric_tensor, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MetricTensorDataEquilateralTetrahedra3D4N, KratosCoreFastSuite)
{
    // Set triangle geometry
    Geometry<Node>::PointsArrayType nodes;
    nodes.push_back(Node::Pointer(new Node(1, 0.0, 0.0, 0.0)));
    nodes.push_back(Node::Pointer(new Node(2, 0.5*std::sqrt(3), 0.5, 0.0)));
    nodes.push_back(Node::Pointer(new Node(3, 0.5*std::sqrt(3), -0.5, 0.0)));
    nodes.push_back(Node::Pointer(new Node(3, std::sqrt(3) / 3.0, 0.0, std::sqrt(6) / 3.0)));
    const auto p_tetrahedra = Geometry<Node>::Pointer(new Tetrahedra3D4<Node>(nodes));

    // Call the triangle metric calculator utility
    double h_ref, met_inf, met_sup;
    BoundedMatrix<double,3,3> metric_tensor;
    GeometryMetricCalculator<3,4>::CalculateMetricTensorDimensionless(*p_tetrahedra, metric_tensor, h_ref, met_inf, met_sup);

    // Check results
    const double tolerance = 1.0e-5;
    KRATOS_CHECK_NEAR(h_ref, 1.0, tolerance);
    KRATOS_CHECK_NEAR(met_inf, 1.0, tolerance);
    KRATOS_CHECK_NEAR(met_sup, 1.0, tolerance);
    BoundedMatrix<double,3,3> expected_metric_tensor;
    expected_metric_tensor(0,0) = 1.0; expected_metric_tensor(0,1) = 0.0; expected_metric_tensor(0,2) = 0.0;
    expected_metric_tensor(1,0) = 0.0; expected_metric_tensor(1,1) = 1.0; expected_metric_tensor(1,2) = 0.0;
    expected_metric_tensor(2,0) = 0.0; expected_metric_tensor(2,1) = 0.0; expected_metric_tensor(2,2) = 1.0;
    KRATOS_CHECK_MATRIX_NEAR(metric_tensor, expected_metric_tensor, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MetricTensorDataUnitTetrahedra3D4N, KratosCoreFastSuite)
{
    // Set triangle geometry
    Geometry<Node>::PointsArrayType nodes;
    nodes.push_back(Node::Pointer(new Node(1, 0.0, 0.0, 0.0)));
    nodes.push_back(Node::Pointer(new Node(2, 1.0, 0.0, 0.0)));
    nodes.push_back(Node::Pointer(new Node(3, 0.0, 1.0, 0.0)));
    nodes.push_back(Node::Pointer(new Node(3, 0.0, 0.0, 1.0)));
    const auto p_tetrahedra = Geometry<Node>::Pointer(new Tetrahedra3D4<Node>(nodes));

    // Call the triangle metric calculator utility
    double h_ref, met_inf, met_sup;
    BoundedMatrix<double,3,3> metric_tensor;
    GeometryMetricCalculator<3,4>::CalculateMetricTensorDimensionless(*p_tetrahedra, metric_tensor, h_ref, met_inf, met_sup);

    // Check results
    const double tolerance = 1.0e-5;
    KRATOS_CHECK_NEAR(h_ref, 1.17851, tolerance);
    KRATOS_CHECK_NEAR(met_inf, 0.5, tolerance);
    KRATOS_CHECK_NEAR(met_sup, 2.0, tolerance);
    BoundedMatrix<double,3,3> expected_metric_tensor;
    expected_metric_tensor(0,0) = 1.38889; expected_metric_tensor(0,1) = 0.694444; expected_metric_tensor(0,2) = 0.694444;
    expected_metric_tensor(1,0) = 0.694444; expected_metric_tensor(1,1) = 1.38889; expected_metric_tensor(1,2) = 0.694444;
    expected_metric_tensor(2,0) = 0.694444; expected_metric_tensor(2,1) = 0.694444; expected_metric_tensor(2,2) = 1.38889;
    KRATOS_CHECK_MATRIX_NEAR(metric_tensor, expected_metric_tensor, tolerance);
}

} // namespace Testing

} // namespace Kratos
