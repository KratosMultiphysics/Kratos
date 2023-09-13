//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Vicente Mataix Ferrandiz
//
//

// Project includes
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/expect.h"
#include "testing/testing.h"
#include "utilities/intersection_utilities.h"

namespace Kratos {
namespace Testing {

    Line2D2<Point> GenerateSlopedLine2D2(){
        Point::Pointer p_point_1 = Kratos::make_shared<Point>(1.0, 0.0, 0.0);
        Point::Pointer p_point_2 = Kratos::make_shared<Point>(0.0, 1.0, 0.0);
        Line2D2<Point> line_geom(p_point_1, p_point_2);
        return line_geom;
    }

    Line2D2<Point> GenerateVerticalLine2D2(){
        Point::Pointer p_point_1 = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
        Point::Pointer p_point_2 = Kratos::make_shared<Point>(0.0, 1.0, 0.0);
        Line2D2<Point> line_geom(p_point_1, p_point_2);
        return line_geom;
    }

    Triangle3D3<Point> GenerateStraightTriangle3D3(){
        Point::Pointer p_point_1 = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
        Point::Pointer p_point_2 = Kratos::make_shared<Point>(0.0, 1.0, 0.0);
        Point::Pointer p_point_3 = Kratos::make_shared<Point>(0.0, 0.0, 1.0);
        Triangle3D3<Point> triang_geom(p_point_1, p_point_2, p_point_3);
        return triang_geom;
    }

    Triangle3D3<Point> GenerateSkewedTriangle3D3(){
        Point::Pointer p_point_1 = Kratos::make_shared<Point>(1.0, 0.0, 0.0);
        Point::Pointer p_point_2 = Kratos::make_shared<Point>(0.0, 0.0, 1.0);
        Point::Pointer p_point_3 = Kratos::make_shared<Point>(0.0, 1.0, 0.5);
        Triangle3D3<Point> triang_geom(p_point_1, p_point_2, p_point_3);
        return triang_geom;
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineLineIntersectionCentered, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto line_geom = GenerateVerticalLine2D2();

        // Set the points that define the intersection line
        const Point line_pt_1(-1.0,0.5,0.0);
        const Point line_pt_2( 1.0,0.5,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<Point>>(
            line_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        const array_1d<double,3>& int_pt_coords = int_pt.Coordinates();
        KRATOS_EXPECT_EQ(int_id, 1);
        KRATOS_EXPECT_NEAR(int_pt_coords[0], 0.0, 1e-10);
        KRATOS_EXPECT_NEAR(int_pt_coords[1], 0.5, 1e-10);
        KRATOS_EXPECT_NEAR(int_pt_coords[2], 0.0, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineLineIntersection, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto line_geom = GenerateSlopedLine2D2();

        // Set the points that define the intersection line
        const Point line_pt_1(-0.5,0.0,0.0);
        const Point line_pt_2( 1.0,1.0,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<Point>>(
            line_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        const array_1d<double,3>& int_pt_coords = int_pt.Coordinates();
        KRATOS_EXPECT_EQ(int_id, 1);
        KRATOS_EXPECT_NEAR(int_pt_coords[0], 0.4, 1e-10);
        KRATOS_EXPECT_NEAR(int_pt_coords[1], 0.6, 1e-10);
        KRATOS_EXPECT_NEAR(int_pt_coords[2], 0.0, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineLineParallel, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto line_geom = GenerateVerticalLine2D2();

        // Set the points that define the intersection line
        const Point line_pt_1(1.0,2.0,0.0);
        const Point line_pt_2(1.0,3.0,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<Point>>(
            line_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        KRATOS_EXPECT_EQ(int_id, 0);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineLineCollinear, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto line_geom = GenerateVerticalLine2D2();

        // Set the points that define the intersection line
        const Point line_pt_1(0.0,0.25,0.0);
        const Point line_pt_2(0.0,0.75,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<Point>>(
            line_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        KRATOS_EXPECT_EQ(int_id, 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineLineNonParallelNoIntersection, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto line_geom = GenerateSlopedLine2D2();

        // Set the points that define the intersection line
        const Point line_pt_1(0.0,0.0,0.0);
        const Point line_pt_2(0.25,0.25,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<Point>>(
            line_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        KRATOS_EXPECT_EQ(int_id, 0);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineLineThroughPoint, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto line_geom = GenerateVerticalLine2D2();

        // Set the points that define the intersection line
        const Point line_pt_1(-1.0,0.0,0.0);
        const Point line_pt_2(0.0,0.0,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<Point>>(
            line_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        const array_1d<double,3>& int_pt_coords = int_pt.Coordinates();
        KRATOS_EXPECT_EQ(int_id, 3);
        KRATOS_EXPECT_NEAR(int_pt_coords[0], 0.0, 1e-6);
        KRATOS_EXPECT_NEAR(int_pt_coords[1], 0.0, 1e-6);
        KRATOS_EXPECT_NEAR(int_pt_coords[2], 0.0, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesTriangleLineCenteredIntersection, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto triang_geom = GenerateStraightTriangle3D3();

        // Set the points that define the intersection line
        const Point line_pt_1(1.0,0.25,0.25);
        const Point line_pt_2(-1.0,0.25,0.25);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeTriangleLineIntersection<Triangle3D3<Point>>(
            triang_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        const array_1d<double,3>& int_pt_coords = int_pt.Coordinates();
        KRATOS_EXPECT_EQ(int_id, 1);
        KRATOS_EXPECT_NEAR(int_pt_coords[0], 0.0, 1e-10);
        KRATOS_EXPECT_NEAR(int_pt_coords[1], 0.25, 1e-10);
        KRATOS_EXPECT_NEAR(int_pt_coords[2], 0.25, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesTriangleLineSkewedIntersection, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto triang_geom = GenerateSkewedTriangle3D3();

        // Set the points that define the intersection line
        const Point line_pt_1(2.0,0.0,0.0);
        const Point line_pt_2(-1.0,0.25,0.25);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeTriangleLineIntersection<Triangle3D3<Point>>(
            triang_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        const array_1d<double,3>& int_pt_coords = int_pt.Coordinates();
        KRATOS_EXPECT_EQ(int_id, 1);
        KRATOS_EXPECT_NEAR(int_pt_coords[0], 0.857143, 1e-6);
        KRATOS_EXPECT_NEAR(int_pt_coords[1], 0.0952381, 1e-6);
        KRATOS_EXPECT_NEAR(int_pt_coords[2], 0.0952381, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesTriangleLineNoIntersection, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto triang_geom = GenerateStraightTriangle3D3();

        // Set the points that define the intersection line
        const Point line_pt_1(1.0,0.25,0.25);
        const Point line_pt_2(0.1,0.25,0.25);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeTriangleLineIntersection<Triangle3D3<Point>>(
            triang_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Check that there is no intersection
        KRATOS_EXPECT_EQ(int_id, 0);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesTriangleLineThroughPoint, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto triang_geom = GenerateStraightTriangle3D3();

        // Set the points that define the intersection line
        const Point line_pt_1(1.0,0.0,0.0);
        const Point line_pt_2(-1.0,0.0,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeTriangleLineIntersection<Triangle3D3<Point>>(
            triang_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point passes through the node
        KRATOS_EXPECT_EQ(int_id, 1);
        KRATOS_EXPECT_NEAR(int_pt.Coordinates()[0], 0.0, 1e-6);
        KRATOS_EXPECT_NEAR(int_pt.Coordinates()[1], 0.0, 1e-6);
        KRATOS_EXPECT_NEAR(int_pt.Coordinates()[2], 0.0, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesTriangleLineCoplanar, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto triang_geom = GenerateStraightTriangle3D3();

        // Set the points that define the intersection line
        const Point line_pt_1(0.0,0.0,0.0);
        const Point line_pt_2(0.0,1.0,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeTriangleLineIntersection<Triangle3D3<Point>>(
            triang_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Check that there is no intersection
        KRATOS_EXPECT_EQ(int_id, 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesTriangleLineCoplanarIntersection, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto triang_geom = GenerateStraightTriangle3D3();

        // Set the points that define the intersection line
        const Point line_pt_1(0.0,0.0,0.0);
        const Point line_pt_2(0.0,1.0,1.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeTriangleLineIntersection<Triangle3D3<Point>>(
            triang_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Check that are co-planar
        KRATOS_EXPECT_EQ(int_id, 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesTriangleLineCoplanarNoIntersection, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto triang_geom = GenerateStraightTriangle3D3();

        // Set the points that define the intersection line
        const Point line_pt_1(0.0,2.0,0.0);
        const Point line_pt_2(0.0,2.0,1.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeTriangleLineIntersection<Triangle3D3<Point>>(
            triang_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Check that are co-planar
        KRATOS_EXPECT_EQ(int_id, 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesTriangleLineBoundaryIntersection, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        Point::Pointer p_point_1 = Kratos::make_shared<Point>(0.302838, 0.210816, 0.5);
        Point::Pointer p_point_2 = Kratos::make_shared<Point>(0.325, 0.196891, 0.5);
        Point::Pointer p_point_3 = Kratos::make_shared<Point>(0.31342, 0.204924, 0.475141);
        Triangle3D3<Point> triang_geom(p_point_1, p_point_2, p_point_3);

        // Set the points that define the intersection line
        const Point line_pt_1(0.3,0.2,0.5);
        const Point line_pt_2(0.4,0.2,0.5);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeTriangleLineIntersection<Triangle3D3<Point>>(
            triang_geom,
            line_pt_1.Coordinates(),
            line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // The triangle and edge are set such that the intersection occurs close to the triangle boundary
        KRATOS_EXPECT_EQ(int_id, 1);
        KRATOS_EXPECT_NEAR(int_pt[0], 0.320052, 1e-6);
        KRATOS_EXPECT_NEAR(int_pt[1], 0.2, 1e-6);
        KRATOS_EXPECT_NEAR(int_pt[2], 0.5, 1e-6);
    }

        KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesPlaneLineCenteredIntersection, KratosCoreFastSuite)
    {
        // Set the plane to be intersected
        const Point plane_base_pt(0.0,0.0,0.0);
        const Point plane_normal(1.0,0.0,0.0);

        // Set the points that define the intersection line
        const Point line_pt_1(1.0,0.25,0.25);
        const Point line_pt_2(-1.0,0.25,0.25);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputePlaneLineIntersection(
            plane_base_pt.Coordinates(), plane_normal.Coordinates(),
            line_pt_1.Coordinates(),line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        KRATOS_EXPECT_EQ(int_id, 1);
        const std::vector<double> expected_values = {0.0,0.25,0.25};
        KRATOS_EXPECT_VECTOR_NEAR(int_pt.Coordinates(), expected_values, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesPlaneLineSkewedIntersection, KratosCoreFastSuite)
    {
        // Set the plane to be intersected
        const Point plane_base_pt(1.0,0.0,0.0);
        const Point plane_normal(-1.0,0.0-0.5,-1.0);

        // Set the points that define the intersection line
        const Point line_pt_1(2.0,0.0,0.0);
        const Point line_pt_2(-1.0,0.25,0.25);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputePlaneLineIntersection(
            plane_base_pt.Coordinates(), plane_normal.Coordinates(),
            line_pt_1.Coordinates(), line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        KRATOS_EXPECT_EQ(int_id, 1);
        const std::vector<double> expected_values = {0.857143,0.0952381,0.0952381};
        KRATOS_EXPECT_VECTOR_NEAR(int_pt.Coordinates(), expected_values, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesPlaneLineNoIntersection, KratosCoreFastSuite)
    {
        // Set the plane to be intersected
        const Point plane_base_pt(0.0,0.0,0.0);
        const Point plane_normal(1.0,0.0,0.0);

        // Set the points that define the intersection line
        const Point line_pt_1(1.0,0.25,0.25);
        const Point line_pt_2(0.1,0.25,0.25);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputePlaneLineIntersection(
            plane_base_pt.Coordinates(), plane_normal.Coordinates(),
            line_pt_1.Coordinates(), line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Check that there is no intersection
        KRATOS_EXPECT_EQ(int_id, 0);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesPlaneLineThroughPoint, KratosCoreFastSuite)
    {
        // Set the plane to be intersected
        const Point plane_base_pt(0.0,0.0,0.0);
        const Point plane_normal(1.0,0.0,0.0);

        // Set the points that define the intersection line
        const Point line_pt_1(1.0,1.0,0.0);
        const Point line_pt_2(0.0,1.0,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputePlaneLineIntersection(
            plane_base_pt.Coordinates(), plane_normal.Coordinates(),
            line_pt_1.Coordinates(), line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Compute and check the obtained intersection point passes through the node
        KRATOS_EXPECT_EQ(int_id, 1);
        KRATOS_EXPECT_VECTOR_NEAR(int_pt.Coordinates(), line_pt_2.Coordinates(), 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesPlaneLineCoplanar, KratosCoreFastSuite)
    {
        // Set the plane to be intersected
        const Point plane_base_pt(0.0,0.0,0.0);
        const Point plane_normal(1.0,0.0,0.0);

        // Set the points that define the intersection line
        const Point line_pt_1(0.0,0.0,0.0);
        const Point line_pt_2(0.0,1.0,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputePlaneLineIntersection(
            plane_base_pt.Coordinates(), plane_normal.Coordinates(),
            line_pt_1.Coordinates(), line_pt_2.Coordinates(),
            int_pt.Coordinates());

        // Check that there is no intersection
        KRATOS_EXPECT_EQ(int_id, 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineBoxIntersectionNoHitpoint, KratosCoreFastSuite)
    {
        // Set the origin and endpoint of the segment
        const Point line_origin(0.5,0.5,2.0);
        const Point line_endpoint(0.6,0.8,2.5);

        // Set the box minimum and maximum point
        const Point box_min_point(0.0,0.0,0.0);
        const Point box_max_point(1.0,1.0,1.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineBoxIntersection(
            box_min_point.Coordinates(),
            box_max_point.Coordinates(),
            line_origin.Coordinates(),
            line_endpoint.Coordinates());

        KRATOS_EXPECT_EQ(int_id, 0);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineBoxIntersectionSingle, KratosCoreFastSuite)
    {
        // Set the origin and endpoint of the segment
        const Point line_origin(0.5,0.5,0.5);
        const Point line_endpoint(-0.5,0.3,0.3);

        // Set the box minimum and maximum point
        const Point box_min_point(0.0,0.0,0.0);
        const Point box_max_point(1.0,1.0,1.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineBoxIntersection(
            box_min_point.Coordinates(),
            box_max_point.Coordinates(),
            line_origin.Coordinates(),
            line_endpoint.Coordinates());

        KRATOS_EXPECT_EQ(int_id, 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineBoxIntersectionDouble, KratosCoreFastSuite)
    {
        // Set the origin and endpoint of the segment
        const Point line_origin(-0.5,0.5,0.5);
        const Point line_endpoint(1.5,0.5,0.5);

        // Set the box minimum and maximum point
        const Point box_min_point(0.0,0.0,0.0);
        const Point box_max_point(1.0,1.0,1.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineBoxIntersection(
            box_min_point.Coordinates(),
            box_max_point.Coordinates(),
            line_origin.Coordinates(),
            line_endpoint.Coordinates());

        KRATOS_EXPECT_EQ(int_id, 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(ComputeTetrahedraLineIntersection1, KratosCoreFastSuite)
    {
        Point point_1 = Point(0.0, 0.0, 1.0);
        Point point_2 = Point(1.0, 0.0, 1.0);

        Point::Pointer p_point_3 = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
        Point::Pointer p_point_4 = Kratos::make_shared<Point>(1.0, 0.0, 0.0);
        Point::Pointer p_point_5 = Kratos::make_shared<Point>(0.0, 1.0, 0.0);
        Point::Pointer p_point_6 = Kratos::make_shared<Point>(0.0, 0.0, 1.0);
        Tetrahedra3D4<Point> tetrahedra(p_point_3, p_point_4, p_point_5, p_point_6);

        // Intersecting line (corner)
        array_1d<double,3> intersection_point1, intersection_point2;
        auto intersection = IntersectionUtilities::ComputeTetrahedraLineIntersection(tetrahedra, point_1.Coordinates(), point_2.Coordinates(), intersection_point1, intersection_point2);
        KRATOS_EXPECT_EQ(intersection, 8);//KRATOS_EXPECT_EQ(intersection, IntersectionUtilitiesTetrahedraLineIntersectionStatus::FOURTH_CORNER);
        KRATOS_EXPECT_VECTOR_EQ(intersection_point1, point_1.Coordinates());
    }

    KRATOS_TEST_CASE_IN_SUITE(ComputeTetrahedraLineIntersection2, KratosCoreFastSuite)
    {
        Point point_1 = Point(0.0, 0.0, 0.5);
        Point point_2 = Point(1.0, 0.0, 0.5);

        Point::Pointer p_point_3 = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
        Point::Pointer p_point_4 = Kratos::make_shared<Point>(1.0, 0.0, 0.0);
        Point::Pointer p_point_5 = Kratos::make_shared<Point>(0.0, 1.0, 0.0);
        Point::Pointer p_point_6 = Kratos::make_shared<Point>(0.0, 0.0, 1.0);
        Tetrahedra3D4<Point> tetrahedra(p_point_3, p_point_4, p_point_5, p_point_6);

        // Intersecting line (face)
        array_1d<double,3> intersection_point1, intersection_point2;
        auto intersection = IntersectionUtilities::ComputeTetrahedraLineIntersection(tetrahedra, point_1.Coordinates(), point_2.Coordinates(), intersection_point1, intersection_point2);
        KRATOS_EXPECT_EQ(intersection, 1);//KRATOS_EXPECT_EQ(intersection, IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION);
        array_1d<double,3> expected_intersection_point1 = ZeroVector(3);
        expected_intersection_point1[0] = 0.5;
        expected_intersection_point1[2] = 0.5;
        array_1d<double,3>  expected_intersection_point2 = ZeroVector(3);
        expected_intersection_point2[2] = 0.5;
        KRATOS_EXPECT_VECTOR_EQ(intersection_point1, expected_intersection_point1);
        KRATOS_EXPECT_VECTOR_EQ(intersection_point2, expected_intersection_point2);
    }

    KRATOS_TEST_CASE_IN_SUITE(ComputeTetrahedraLineIntersection3, KratosCoreFastSuite)
    {
        Point point_1 = Point(0.25, 0.0, 0.0);
        Point point_2 = Point(1.0, 0.0, 0.0);

        Point::Pointer p_point_3 = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
        Point::Pointer p_point_4 = Kratos::make_shared<Point>(1.0, 0.0, 0.0);
        Point::Pointer p_point_5 = Kratos::make_shared<Point>(0.0, 1.0, 0.0);
        Point::Pointer p_point_6 = Kratos::make_shared<Point>(0.0, 0.0, 1.0);
        Tetrahedra3D4<Point> tetrahedra(p_point_3, p_point_4, p_point_5, p_point_6);

        // Intersecting line (edge first)
        array_1d<double,3> intersection_point1, intersection_point2;
        auto intersection = IntersectionUtilities::ComputeTetrahedraLineIntersection(tetrahedra, point_1.Coordinates(), point_2.Coordinates(), intersection_point1, intersection_point2);
        KRATOS_EXPECT_EQ(intersection, 1);//KRATOS_EXPECT_EQ(intersection, IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION);
        KRATOS_EXPECT_VECTOR_EQ(intersection_point1, Point(1.0, 0.0, 0.0).Coordinates());
        KRATOS_EXPECT_VECTOR_EQ(intersection_point2, point_1.Coordinates());

        // Intersecting line (edge second)
        point_2.X() = 0.75;
        intersection = IntersectionUtilities::ComputeTetrahedraLineIntersection(tetrahedra, point_1.Coordinates(), point_2.Coordinates(), intersection_point1, intersection_point2);
        KRATOS_EXPECT_EQ(intersection, 1);//KRATOS_EXPECT_EQ(intersection, IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION);
        KRATOS_EXPECT_VECTOR_EQ(intersection_point1, point_1.Coordinates());
        KRATOS_EXPECT_VECTOR_EQ(intersection_point2, point_2.Coordinates());
    }

    KRATOS_TEST_CASE_IN_SUITE(ComputeTetrahedraLineIntersection4, KratosCoreFastSuite)
    {
        Point point_1 = Point(0.0, 0.25, 0.5);
        Point point_2 = Point(1.0, 0.25, 0.5);

        Point::Pointer p_point_3 = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
        Point::Pointer p_point_4 = Kratos::make_shared<Point>(1.0, 0.0, 0.0);
        Point::Pointer p_point_5 = Kratos::make_shared<Point>(0.0, 1.0, 0.0);
        Point::Pointer p_point_6 = Kratos::make_shared<Point>(0.0, 0.0, 1.0);
        Tetrahedra3D4<Point> tetrahedra(p_point_3, p_point_4, p_point_5, p_point_6);

        // Intersecting line 
        array_1d<double,3> intersection_point1, intersection_point2;
        auto intersection = IntersectionUtilities::ComputeTetrahedraLineIntersection(tetrahedra, point_1.Coordinates(), point_2.Coordinates(), intersection_point1, intersection_point2);
        KRATOS_EXPECT_EQ(intersection, 1);//KRATOS_EXPECT_EQ(intersection, IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION);
        array_1d<double,3> expected_intersection_point1 = ZeroVector(3);
        expected_intersection_point1[0] = 0.25;
        expected_intersection_point1[1] = 0.25;
        expected_intersection_point1[2] = 0.5;
        array_1d<double,3>  expected_intersection_point2 = ZeroVector(3);
        expected_intersection_point2[1] = 0.25;
        expected_intersection_point2[2] = 0.5;
        KRATOS_EXPECT_VECTOR_EQ(intersection_point1, expected_intersection_point1);
        KRATOS_EXPECT_VECTOR_EQ(intersection_point2, expected_intersection_point2);

        // Not intersecting line 
        point_1.Z() = 1.25; 
        point_2.Z() = 1.25;
        intersection = IntersectionUtilities::ComputeTetrahedraLineIntersection(tetrahedra, point_1.Coordinates(), point_2.Coordinates(), intersection_point1, intersection_point2);
        KRATOS_EXPECT_EQ(intersection, 0);//KRATOS_EXPECT_EQ(intersection, IntersectionUtilitiesTetrahedraLineIntersectionStatus::NO_INTERSECTION);
    }

    KRATOS_TEST_CASE_IN_SUITE(ComputeTetrahedraLineIntersection5, KratosCoreFastSuite)
    {
        Point point_1 = Point(0.0, 0.25, 1.25);
        Point point_2 = Point(1.0, 0.25, 1.25);

        Point::Pointer p_point_3 = Kratos::make_shared<Point>(-10.0, -10.0, -10.0);
        Point::Pointer p_point_4 = Kratos::make_shared<Point>(100.0, 0.0, 0.0);
        Point::Pointer p_point_5 = Kratos::make_shared<Point>(0.0, 100.0, 0.0);
        Point::Pointer p_point_6 = Kratos::make_shared<Point>(0.0, 0.0, 100.0);
        Tetrahedra3D4<Point> tetrahedra(p_point_3, p_point_4, p_point_5, p_point_6);

        // Intersecting line (totally inside)
        array_1d<double,3> intersection_point1, intersection_point2;
        auto intersection = IntersectionUtilities::ComputeTetrahedraLineIntersection(tetrahedra, point_1.Coordinates(), point_2.Coordinates(), intersection_point1, intersection_point2);
        KRATOS_EXPECT_EQ(intersection, 3);//KRATOS_EXPECT_EQ(intersection, IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION_BOTH_INSIDE);
        KRATOS_EXPECT_VECTOR_EQ(intersection_point1, point_1.Coordinates());
        KRATOS_EXPECT_VECTOR_EQ(intersection_point2, point_2.Coordinates());

        // Intersecting line (partially inside)
        point_2.Z() = 125.0;
        intersection = IntersectionUtilities::ComputeTetrahedraLineIntersection(tetrahedra, point_1.Coordinates(), point_2.Coordinates(), intersection_point1, intersection_point2);
        KRATOS_EXPECT_EQ(intersection, 4);//KRATOS_EXPECT_EQ(intersection, IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION_ONE_INSIDE);
        array_1d<double,3> expected_intersection_point1 = ZeroVector(3);
        expected_intersection_point1[0] = 0.789579;
        expected_intersection_point1[1] = 0.25;
        expected_intersection_point1[2] = 98.9604;
        KRATOS_EXPECT_VECTOR_NEAR(intersection_point1, expected_intersection_point1, 1.0e-4);
        KRATOS_EXPECT_VECTOR_EQ(intersection_point2, point_1.Coordinates());
    }

    KRATOS_TEST_CASE_IN_SUITE(ComputeShortestLineBetweenTwoLines, KratosCoreFastSuite)
    {
        Point::Pointer p_point_1 = Kratos::make_shared<Point>(1.0, 0.0, 0.0);
        Point::Pointer p_point_2 = Kratos::make_shared<Point>(0.0, 1.0, 0.0);
        Line3D2<Point> line1(p_point_1, p_point_2);

        Point::Pointer p_point_3 = Kratos::make_shared<Point>(1.0, 0.0, 0.5);
        Point::Pointer p_point_4 = Kratos::make_shared<Point>(0.0, 2.0, 1.0);
        Line3D2<Point> line2(p_point_3, p_point_4);

        Point::Pointer p_point_5 = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
        Point::Pointer p_point_6 = Kratos::make_shared<Point>(2.0, 2.0, 0.0);
        Line3D2<Point> line3(p_point_5, p_point_6);

        Point::Pointer p_point_7 = Kratos::make_shared<Point>(2.0, 0.0, 0.0);
        Point::Pointer p_point_8 = Kratos::make_shared<Point>(0.0, 2.0, 0.0);
        Line3D2<Point> line4(p_point_7, p_point_8);

        Point::Pointer p_point_9 = Kratos::make_shared<Point>(1.75, 0.25, 0.0);
        Point::Pointer p_point_10 = Kratos::make_shared<Point>(0.25, 1.75, 0.0);
        Line3D2<Point> line5(p_point_9, p_point_10);

        // Not intersecting lines
        const auto line1_2 = Line3D2<Point>(IntersectionUtilities::ComputeShortestLineBetweenTwoLines(line1, line2));
        KRATOS_EXPECT_NEAR(line1_2.Length(), 0.408248, 1.0e-6);

        // Intersecting lines
        const auto line1_3 = Line3D2<Point>(IntersectionUtilities::ComputeShortestLineBetweenTwoLines(line1, line3));
        KRATOS_EXPECT_DOUBLE_EQ(line1_3.Length(), 0.0);
        array_1d<double, 3> expected_center = ZeroVector(3);
        expected_center[0] = 0.5;
        expected_center[1] = 0.5;
        KRATOS_EXPECT_VECTOR_NEAR(line1_3.Center().Coordinates(), expected_center, 1e-10);

        // Not intersecting lines (parallel lines)
        const auto line1_4 = Line3D2<Point>(IntersectionUtilities::ComputeShortestLineBetweenTwoLines(line1, line4));
        KRATOS_EXPECT_NEAR(line1_4.Length(), std::sqrt(2.0)/2.0, 1.0e-6);

        const auto line1_5 = Line3D2<Point>(IntersectionUtilities::ComputeShortestLineBetweenTwoLines(line1, line5));
        KRATOS_EXPECT_NEAR(line1_5.Length(), std::sqrt(2.0)/2.0, 1.0e-6);
    }

}  // namespace Testing.
}  // namespace Kratos.