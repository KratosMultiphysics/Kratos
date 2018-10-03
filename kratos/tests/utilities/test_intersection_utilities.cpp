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
//
//

// Project includes
#include "geometries/point.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "includes/checks.h"
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
        const array_1d<double,3> int_pt_coords = int_pt.Coordinates();
        KRATOS_CHECK_EQUAL(int_id, 1);
        KRATOS_CHECK_NEAR(int_pt_coords(0), 0.0, 1e-10);
        KRATOS_CHECK_NEAR(int_pt_coords(1), 0.5, 1e-10);
        KRATOS_CHECK_NEAR(int_pt_coords(2), 0.0, 1e-10);
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
        const array_1d<double,3> int_pt_coords = int_pt.Coordinates();
        KRATOS_CHECK_EQUAL(int_id, 1);
        KRATOS_CHECK_NEAR(int_pt_coords(0), 0.4, 1e-10);
        KRATOS_CHECK_NEAR(int_pt_coords(1), 0.6, 1e-10);
        KRATOS_CHECK_NEAR(int_pt_coords(2), 0.0, 1e-10);
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
        const array_1d<double,3> int_pt_coords = int_pt.Coordinates();
        KRATOS_CHECK_EQUAL(int_id, 0);
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
        const array_1d<double,3> int_pt_coords = int_pt.Coordinates();
        KRATOS_CHECK_EQUAL(int_id, 2);
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
        const array_1d<double,3> int_pt_coords = int_pt.Coordinates();
        KRATOS_CHECK_EQUAL(int_id, 0);
    }

    KRATOS_TEST_CASE_IN_SUITE(IntersectionUtilitiesLineLineThroughPoint, KratosCoreFastSuite)
    {
        // Set the triangle to be intersected
        auto line_geom = GenerateVerticalLine2D2();

        // Set the points that define the intersection line
        const Point line_pt_1(0.0,0.0,0.0);
        const Point line_pt_2(1.0,0.0,0.0);

        // Initialize the intersection point
        Point int_pt(0.0,0.0,0.0);

        // Call the intersection utility
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<Point>>(
            line_geom, 
            line_pt_1.Coordinates(), 
            line_pt_2.Coordinates(), 
            int_pt.Coordinates());

        // Compute and check the obtained intersection point coordinates
        const array_1d<double,3> int_pt_coords = int_pt.Coordinates();
        KRATOS_CHECK_EQUAL(int_id, 1);
        KRATOS_CHECK_NEAR(int_pt.Coordinates()[0], 0.0, 1e-6);
        KRATOS_CHECK_NEAR(int_pt.Coordinates()[1], 0.0, 1e-6);
        KRATOS_CHECK_NEAR(int_pt.Coordinates()[2], 0.0, 1e-6);
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
        const array_1d<double,3> int_pt_coords = int_pt.Coordinates();
        KRATOS_CHECK_EQUAL(int_id, 1);
        KRATOS_CHECK_NEAR(int_pt_coords(0), 0.0, 1e-10);
        KRATOS_CHECK_NEAR(int_pt_coords(1), 0.25, 1e-10);
        KRATOS_CHECK_NEAR(int_pt_coords(2), 0.25, 1e-10);
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
        const array_1d<double,3> int_pt_coords = int_pt.Coordinates();
        KRATOS_CHECK_EQUAL(int_id, 1);
        KRATOS_CHECK_NEAR(int_pt_coords(0), 0.857143, 1e-6);
        KRATOS_CHECK_NEAR(int_pt_coords(1), 0.0952381, 1e-6);
        KRATOS_CHECK_NEAR(int_pt_coords(2), 0.0952381, 1e-6);
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
        KRATOS_CHECK_EQUAL(int_id, 0);
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
        KRATOS_CHECK_EQUAL(int_id, 1);
        KRATOS_CHECK_NEAR(int_pt.Coordinates()[0], 0.0, 1e-6);
        KRATOS_CHECK_NEAR(int_pt.Coordinates()[1], 0.0, 1e-6);
        KRATOS_CHECK_NEAR(int_pt.Coordinates()[2], 0.0, 1e-6);
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
        KRATOS_CHECK_EQUAL(int_id, 2);
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
        KRATOS_CHECK_EQUAL(int_id, 1);
        KRATOS_CHECK_NEAR(int_pt[0], 0.320052, 1e-6);
        KRATOS_CHECK_NEAR(int_pt[1], 0.2, 1e-6);
        KRATOS_CHECK_NEAR(int_pt[2], 0.5, 1e-6);
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

        KRATOS_CHECK_EQUAL(int_id, 0);
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

        KRATOS_CHECK_EQUAL(int_id, 1);
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

        KRATOS_CHECK_EQUAL(int_id, 1);
    }

}  // namespace Testing.
}  // namespace Kratos.
