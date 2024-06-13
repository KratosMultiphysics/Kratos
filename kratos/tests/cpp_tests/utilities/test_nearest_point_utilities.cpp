//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(LineNearestPoint, KratosCoreFastSuite)
{
    constexpr double length = 1.2;

    Point::Pointer p_point_1(make_shared<Point>( 0.00, 0.00, 0.00));
    Point::Pointer p_point_2(make_shared<Point>( length, 0.00, 0.00));
    Point::Pointer p_point_3(make_shared<Point>( 0.00, length, 0.00));

    Point nearest_point( 0.00, 0.00, 0.00);

    Line3D2<Point> line_1(p_point_1, p_point_2);
    Line3D2<Point> line_2(p_point_1, p_point_3);

    Point point_1(0.2 * length, 0.1 * length, 0.00);
    nearest_point = NearestPointUtilities::LineNearestPoint(point_1, line_1);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.00, 0.00), 1e-6);

    nearest_point = NearestPointUtilities::LineNearestPoint(point_1, line_2);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.00, 0.1 * length, 0.00), 1e-6);

    Point point_2(-0.2 * length, -0.1 * length, 0.00);
    nearest_point = NearestPointUtilities::LineNearestPoint(point_2, line_1);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point,Point(0.00, 0.00, 0.00), 1e-6);

    nearest_point = NearestPointUtilities::LineNearestPoint(point_2, line_2);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.00, 0.00, 0.00), 1e-6);

    Point point_3(1.2 * length, 1.1 * length, 0.00);
    nearest_point = NearestPointUtilities::LineNearestPoint(point_3, line_1);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point,Point(length, 0.00, 0.00), 1e-6);

    nearest_point = NearestPointUtilities::LineNearestPoint(point_3, line_2);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.00, length, 0.00), 1e-6);     
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInPlaneNearestPoint, KratosCoreFastSuite)
{
    constexpr double length = 1.2;

    Point::Pointer p_point_1(make_shared<Point>( 0.00, 0.00, 0.00));
    Point::Pointer p_point_2(make_shared<Point>( length, 0.00, 0.00));
    Point::Pointer p_point_3(make_shared<Point>( length, length, 0.00));

    Point nearest_point( 0.00, 0.00, 0.00);

    Triangle3D3<Point> triangle(p_point_1, p_point_3, p_point_2);

    Point inside_point(0.2 * length, 0.1 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(inside_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, inside_point, 1e-6);

    Point border_point(0.2 * length, 0.2 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(border_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, border_point, 1e-6);

    Point lower_corner(-2 * length, -0.1 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(lower_corner, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_1), 1e-6);

    Point below_point(0.2 * length, -0.2 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(below_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.0, 0.00), 1e-6);

    Point below_left_point(-0.1 * length, -0.2 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(below_left_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_1), 1e-6);

    Point below_point_right(1.2 * length, -0.2 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(below_point_right, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_2), 1e-6);

    Point right_side(1.2 * length, 0.2 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(right_side, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(length, 0.2 * length, 0.00), 1e-6);       

    Point upper_corner(1.2 * length, 1.2 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(upper_corner, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_3), 1e-6);

    Point left_point(0.2 * length, 0.4 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(left_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.3 * length, 0.3 * length, 0.00), 1e-6);       

    Point far_left_point(-0.2 * length, 0.2 * length, 0.00);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(far_left_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.00, 0.00, 0.00), 1e-6);       
}

KRATOS_TEST_CASE_IN_SUITE(TriangleOutOfPlaneNearestPoint, KratosCoreFastSuite)
{
    constexpr double length = 1.2;
    constexpr double distance = 2.1;

    Point::Pointer p_point_1(make_shared<Point>( 0.00, 0.00, 0.00));
    Point::Pointer p_point_2(make_shared<Point>( length, 0.00, 0.00));
    Point::Pointer p_point_3(make_shared<Point>( length, length, 0.00));

    Point nearest_point( 0.00, 0.00, 0.00);

    Triangle3D3<Point> triangle(p_point_1, p_point_3, p_point_2);

    Point inside_point(0.2 * length, 0.1 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(inside_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.1 * length, 0.00), 1e-6);

    Point border_point(0.2 * length, 0.2 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(border_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.2 * length, 0.00), 1e-6);

    Point lower_corner(-2 * length, -0.1 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(lower_corner, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_1), 1e-6);

    Point below_point(0.2 * length, -0.2 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(below_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.0, 0.00), 1e-6);

    Point below_left_point(-0.1 * length, -0.2 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(below_left_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_1), 1e-6);

    Point below_point_right(1.2 * length, -0.2 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(below_point_right, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_2), 1e-6);

    Point right_side(1.2 * length, 0.2 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(right_side, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(length, 0.2 * length, 0.00), 1e-6);       

    Point upper_corner(1.2 * length, 1.2 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(upper_corner, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_3), 1e-6);

    Point left_point(0.2 * length, 0.4 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(left_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.3 * length, 0.3 * length, 0.00), 1e-6);       

    Point far_left_point(-0.2 * length, 0.2 * length, distance);
    nearest_point = NearestPointUtilities::TriangleNearestPoint(far_left_point, triangle);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.00, 0.00, 0.00), 1e-6);       
}

} // namespace Kratos::Testing
