//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla Martinez
//                   Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "includes/checks.h"
#include "utilities/geometry_utilities.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
// #include "includes/gid_io.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    Line2D2 <NodeType> GenerateExampleLine()
    {
        // First we create the nodes
        NodeType::Pointer p_node_1 = Kratos::make_intrusive<Node<3>>(1, 0.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_2 = Kratos::make_intrusive<Node<3>>(2, 1.0 , 0.0 , 0.0);

        // Now we create the geometry
        std::vector<NodeType::Pointer> line_nodes (2);
        line_nodes[0] = p_node_1;
        line_nodes[1] = p_node_2;
        Line2D2 <NodeType> line( PointerVector<NodeType>{line_nodes} );

        return line;
    }

    Triangle2D3 <NodeType> GenerateExampleTriangle()
    {
        // First we create the nodes
        NodeType::Pointer p_node_1 = Kratos::make_intrusive<Node<3>>(1, 0.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_2 = Kratos::make_intrusive<Node<3>>(2, 1.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_3 = Kratos::make_intrusive<Node<3>>(3, 1.0 , 1.0 , 0.0);

        // Now we create the geometry
        std::vector<NodeType::Pointer> triangle_nodes (3);
        triangle_nodes[0] = p_node_1;
        triangle_nodes[1] = p_node_2;
        triangle_nodes[2] = p_node_3;
        Triangle2D3 <NodeType> triangle( PointerVector<NodeType>{triangle_nodes} );

        return triangle;
    }

    Triangle3D3 <NodeType> GenerateExampleTriangle3D()
    {
        // First we create the nodes
        NodeType::Pointer p_node_1 = Kratos::make_intrusive<Node<3>>(1, 0.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_2 = Kratos::make_intrusive<Node<3>>(2, 1.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_3 = Kratos::make_intrusive<Node<3>>(3, 1.0 , 1.0 , 0.0);

        // Now we create the geometry
        std::vector<NodeType::Pointer> triangle_nodes (3);
        triangle_nodes[0] = p_node_1;
        triangle_nodes[1] = p_node_2;
        triangle_nodes[2] = p_node_3;
        Triangle3D3 <NodeType> triangle( PointerVector<NodeType>{triangle_nodes} );

        return triangle;
    }

    Tetrahedra3D4 <NodeType> GenerateExampleTetrahedra()
    {
        // First we create the nodes
        NodeType::Pointer p_node_1 = Kratos::make_intrusive<Node<3>>(1, 0.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_2 = Kratos::make_intrusive<Node<3>>(2, 1.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_3 = Kratos::make_intrusive<Node<3>>(3, 1.0 , 1.0 , 0.0);
        NodeType::Pointer p_node_4 = Kratos::make_intrusive<Node<3>>(4, 1.0 , 1.0 , 1.0);

        // Now we create the geometry
        std::vector<NodeType::Pointer> tetrahedra_nodes (4);
        tetrahedra_nodes[0] = p_node_1;
        tetrahedra_nodes[1] = p_node_2;
        tetrahedra_nodes[2] = p_node_3;
        tetrahedra_nodes[3] = p_node_4;
        Tetrahedra3D4 <NodeType> tetrahedra( PointerVector<NodeType>{tetrahedra_nodes} );

        return tetrahedra;
    }

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsPointDistanceToTriangleOutOfPlane, KratosCoreFastSuite)
    {
        // Triangle plane points
        Point triangle_point_1(-1.0,-1.0, 0.0);
        Point triangle_point_2( 1.0,-1.0, 0.0);
        Point triangle_point_3(-1.0, 1.0, 0.0);

        // Point out of plane to compute the distance to
        Point distance_point(0.357143, -0.214286, 0.0714286);


        const double dist = GeometryUtils::PointDistanceToTriangle3D(triangle_point_1, triangle_point_2, triangle_point_3, distance_point);

        KRATOS_CHECK_NEAR(dist, 0.123718, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsPointDistanceToTriangleInPlane, KratosCoreFastSuite)
    {
        // Triangle plane points
        Point triangle_point_1( 1.0,-1.0, 0.0);
        Point triangle_point_2( 1.0, 1.0, 0.0);
        Point triangle_point_3(-1.0, 1.0, 0.0);

        // Point over the plane to compute the distance to
        Point distance_point(0.357143, -0.214286, 0.0714286);

        const double dist = GeometryUtils::PointDistanceToTriangle3D(triangle_point_1, triangle_point_2, triangle_point_3, distance_point);

        KRATOS_CHECK_NEAR(dist, distance_point.Z(), 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsCalculateGeometryDataTriangle, KratosCoreFastSuite)
    {
        Triangle2D3 <NodeType> triangle = GenerateExampleTriangle();

        // Computing the info
        BoundedMatrix<double,3,2> DN_DX;
        array_1d<double,3> N;
        double area;

        GeometryUtils::CalculateGeometryData(triangle, DN_DX, N, area);

        const double tolerance = 1.0e-8;

        for (int i = 0; i < 3; ++i)
            KRATOS_CHECK_NEAR(N[i], 1.0/3.0, tolerance);

        KRATOS_CHECK_NEAR(area, 0.5, tolerance);

        KRATOS_CHECK_NEAR(DN_DX(0,0), -1.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(0,1), 0.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(1,0), 1.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(1,1), -1.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(0,1), 0.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(2,0), 0.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(2,1), 1.0,tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsCalculateGeometryDataTetrahedra, KratosCoreFastSuite)
    {
        Tetrahedra3D4 <NodeType> tetrahedra = GenerateExampleTetrahedra();

        // Computing the info
        BoundedMatrix<double,4,3> DN_DX;
        array_1d<double,4> N;
        double volume;

        GeometryUtils::CalculateGeometryData(tetrahedra, DN_DX, N, volume);

        const double tolerance = 1.0e-8;

        for (int i = 0; i < 4; ++i)
            KRATOS_CHECK_NEAR(N[i], 1.0/4.0, tolerance);

        KRATOS_CHECK_NEAR(volume, 1.0/6.0, tolerance);

        KRATOS_CHECK_NEAR(DN_DX(0,0), -1.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(0,1), 0.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(0,2), 0.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(1,0), 1.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(1,1), -1.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(1,2), 0.0, tolerance);
        KRATOS_CHECK_NEAR(DN_DX(2,0), 0.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(2,1), 1.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(2,2), -1.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(3,0), 0.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(3,1), 0.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(3,2), 1.0, tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsCalculateGeometryDataLine, KratosCoreFastSuite)
    {
        Line2D2 <NodeType> line = GenerateExampleLine();

        // Computing the info
        BoundedMatrix<double,2,1> DN_DX;
        array_1d<double,2> N;
        double length;

        GeometryUtils::CalculateGeometryData(line, DN_DX, N, length);

        const double tolerance = 1.0e-8;

        for (int i = 0; i < 2; ++i)
            KRATOS_CHECK_NEAR(N[i], 1.0/2.0, tolerance);

        KRATOS_CHECK_NEAR(length, 1.0, tolerance);

        KRATOS_CHECK_NEAR(DN_DX(0,0), -1.0,tolerance);
        KRATOS_CHECK_NEAR(DN_DX(1,0), 1.0,tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsSideLenghts2D, KratosCoreFastSuite)
    {
        Triangle2D3 <NodeType> triangle = GenerateExampleTriangle();

        // Computing the info
        double hmin, hmax;
        GeometryUtils::SideLenghts2D(triangle, hmin, hmax);

        const double tolerance = 1.0e-8;

        KRATOS_CHECK_NEAR(hmin, 1.0,tolerance);
        KRATOS_CHECK_NEAR(hmax, std::sqrt(2.0),tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsCalculateTriangleDistances, KratosCoreFastSuite)
    {
        Triangle3D3 <NodeType> triangle = GenerateExampleTriangle3D();

        // Computing the info
        array_1d<double, 3> distances;
        distances[0] = 0.1;
        distances[1] = 0.2;
        distances[2] = -0.3;
        GeometryUtils::CalculateTriangleDistances(triangle, distances);

        const double tolerance = 1.0e-6;

        KRATOS_CHECK_NEAR(distances[0], 0.353553, tolerance);
        KRATOS_CHECK_NEAR(distances[1], 0.392232, tolerance);
        KRATOS_CHECK_NEAR(distances[2], 0.6, tolerance);
    }

//     KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsCalculateTetrahedraIntersectionPoints, KratosCoreFastSuite)
//     {
//         Tetrahedra3D4 <NodeType> tetrahedra = GenerateExampleTetrahedra()
//         // Computing the info
//         array_1d<double, 4> distances;
//         distances[0] = 0.1;
//         distances[1] = 0.2;
//         distances[2] = -0.3;
//         distances[2] = 0.2;
//         array_1d<Point, 4> intersection_points;
//         const int intersections = GeometryUtils::CalculateTetrahedraIntersectionPoints(tetrahedra, distances, intersection_points);
// 
//        KRATOS_CHECK_EQUAL(intersections, 1);
//     }

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsSeveralUtilities, KratosCoreFastSuite)
    {
        Triangle2D3 <NodeType> triangle = GenerateExampleTriangle();

        // Getting reference solution
        BoundedMatrix<double,3,2> reference_DN_DX; // NOTE: On reference configuration
        array_1d<double,3> N;
        double area;

        GeometryUtils::CalculateGeometryData(triangle, reference_DN_DX, N, area); // NOTE: COmputing here, this method always computes on current configuration

        // Adding some deformation
        triangle[0].X() += 0.01;

        // Computing jacobian
        const auto integration_method = triangle.GetDefaultIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points = triangle.IntegrationPoints( integration_method );
        Matrix J0(2, 2);
        GeometryUtils::JacobianOnInitialConfiguration(triangle, integration_points[0], J0);

        const double tolerance = 1.0e-6;

        // Checking jacobian
        KRATOS_CHECK_NEAR(J0(0, 0), 1.0, tolerance);
        KRATOS_CHECK_NEAR(J0(1, 0), 0.0, tolerance);
        KRATOS_CHECK_NEAR(J0(0, 1), 1.0, tolerance);
        KRATOS_CHECK_NEAR(J0(1, 1), 1.0, tolerance);

        Matrix auxiliar_J0(2, 2);
        GeometryUtils::DirectJacobianOnInitialConfiguration<Matrix>(triangle, auxiliar_J0, 0, integration_method);

        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                KRATOS_CHECK_NEAR(J0(i, j), auxiliar_J0(i, j), tolerance);
            }
        }

        // Computing
        double detJ0;
        Matrix InvJ0(2, 2);
        MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);
        Matrix DN_De(3, 2);
        triangle.ShapeFunctionsLocalGradients( DN_De, integration_points[0].Coordinates() );
        Matrix DN_DX(3, 2);
        GeometryUtils::ShapeFunctionsGradients(DN_De, InvJ0, DN_DX);

        // Comparing
        for (std::size_t i = 0; i < DN_DX.size1(); ++i) {
            for (std::size_t j = 0; j < DN_DX.size2(); ++j) {
                KRATOS_CHECK_NEAR(DN_DX(i, j), reference_DN_DX(i, j), tolerance);
            }
        }

        // Computing deformation gradient
        Matrix J;
        J = triangle.Jacobian(J, 0, integration_method);
        Matrix F(2, 2);
        GeometryUtils::DeformationGradient(J, InvJ0, F);

        // Checking deformation gradient
        KRATOS_CHECK_NEAR(F(0, 0), 0.99, tolerance);
        KRATOS_CHECK_NEAR(F(1, 0), 0.0, tolerance);
        KRATOS_CHECK_NEAR(F(0, 1), 0.0, tolerance);
        KRATOS_CHECK_NEAR(F(1, 1), 1.0, tolerance);

        Matrix auxiliar_J(2, 2);
        GeometryUtils::DirectJacobianOnCurrentConfiguration<Matrix>(triangle, integration_points[0].Coordinates(), auxiliar_J);

        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                KRATOS_CHECK_NEAR(J(i, j), auxiliar_J(i, j), tolerance);
            }
        }
    }
}  // namespace Testing.
}  // namespace Kratos.
