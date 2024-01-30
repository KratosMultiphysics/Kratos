//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//  			 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "utilities/geometrical_projection_utilities.h"
#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"

namespace Kratos::Testing
{

using NodeType = Node;
using GeometryNodeType = Geometry<NodeType>;
using GeometryPointType = Geometry<Point>;

namespace
{

GeometryNodeType::Pointer CreateLine3D2NForTestNode2D()
{
    GeometryNodeType::PointsArrayType points;
    points.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    points.push_back(Kratos::make_intrusive<Node>(2, 2.0, 0.0, 0.0));

    return GeometryNodeType::Pointer(new Line3D2<NodeType>(points));
}

GeometryPointType::Pointer CreateLine3D2NForTestPoint2D()
{
    GeometryPointType::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(0.0, 0.0, 0.0));
    points.push_back(Kratos::make_shared<Point>(2.0, 0.0, 0.0));

    return GeometryPointType::Pointer(new Line3D2<Point>(points));
}

GeometryNodeType::Pointer CreateLine3D2NForTestNode3D()
{
    GeometryNodeType::PointsArrayType points;
    points.push_back(Kratos::make_intrusive<Node>(1, 1.0, 3.0, -1.0));
    points.push_back(Kratos::make_intrusive<Node>(2, 3.0, 6.0, 0.0));

    return GeometryNodeType::Pointer(new Line3D2<NodeType>(points));
}

GeometryPointType::Pointer CreateLine3D2NForTestPoint3D()
{
    GeometryPointType::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(1.0, 3.0, -1.0));
    points.push_back(Kratos::make_shared<Point>(3.0, 6.0, 0.0));

    return GeometryPointType::Pointer(new Line3D2<Point>(points));
}

GeometryNodeType::Pointer CreateTriangle3D3NForTestNode()
{
    GeometryNodeType::PointsArrayType points;
    points.push_back(Kratos::make_intrusive<Node>(1,0.04, 0.02, 0.0));
    points.push_back(Kratos::make_intrusive<Node>(2,1.1, 0.03, 0.0));
    points.push_back(Kratos::make_intrusive<Node>(3,1.08, 1.0, 0.0));

    return GeometryNodeType::Pointer(new Triangle3D3<NodeType>(points));
}

GeometryPointType::Pointer CreateTriangle3D3NForTestPoint()
{
    GeometryPointType::PointsArrayType points;
    points.push_back(Kratos::make_shared<Point>(0.04, 0.02, 0.0));
    points.push_back(Kratos::make_shared<Point>(1.1, 0.03, 0.0));
    points.push_back(Kratos::make_shared<Point>(1.08, 1.0, 0.0));

    return GeometryPointType::Pointer(new Triangle3D3<Point>(points));
}

template<class TGeometryType>
void TestFastProjectDirection(const TGeometryType& rGeom)
{
    const double expected_proj_dist = 1.258;

    const double x_coord = 0.325;
    const double y_coord = 0.147;

    const Point point_to_proj(x_coord, y_coord, expected_proj_dist);

    array_1d<double,3> dir_vector;
    array_1d<double,3> normal_vector;

    dir_vector[0] = 0.0;
    dir_vector[1] = 0.0;
    dir_vector[2] = -1.0;

    normal_vector[0] = 0.0;
    normal_vector[1] = 0.0;
    normal_vector[2] = 1.0;

    Point projected_point;

    double proj_distance = GeometricalProjectionUtilities::FastProjectDirection(
        rGeom,
        point_to_proj,
        projected_point,
        normal_vector,
        dir_vector);

    KRATOS_EXPECT_DOUBLE_EQ(expected_proj_dist, proj_distance);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.X(), x_coord);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Y(), y_coord);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Z(), 0.0);
}

template<class TGeometryType>
void TestFastProjectOnGeometry(TGeometryType& rGeom)
{
    const double expected_proj_dist = -1.258;

    const double x_coord = 0.325;
    const double y_coord = 0.147;

    const Point point_to_proj(x_coord, y_coord, -expected_proj_dist);
    Point projected_point;

    const double proj_distance = GeometricalProjectionUtilities::FastProjectOnGeometry(
        rGeom,
        point_to_proj,
        projected_point);

    KRATOS_EXPECT_DOUBLE_EQ(expected_proj_dist, proj_distance);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.X(), x_coord);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Y(), y_coord);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Z(), 0.0);
}

template<class TGeometryType>
void TestFastProjectOnLine2D(TGeometryType& rGeom)
{
    const double expected_proj_dist = 1.258;

    const double x_coord = 0.325;

    const Point point_to_proj(x_coord, expected_proj_dist, 0.0);
    Point projected_point;

    double proj_distance = GeometricalProjectionUtilities::FastProjectOnLine(
        rGeom,
        point_to_proj,
        projected_point);

    KRATOS_EXPECT_DOUBLE_EQ(expected_proj_dist, proj_distance);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.X(), x_coord);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Y(), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Z(), 0.0);

    proj_distance = GeometricalProjectionUtilities::FastProjectOnLine2D(
        rGeom,
        point_to_proj,
        projected_point);

    KRATOS_EXPECT_DOUBLE_EQ(expected_proj_dist, proj_distance);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.X(), x_coord);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Y(), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Z(), 0.0);
}

template<class TGeometryType>
void TestFastProjectOnLine3D(TGeometryType& rGeom)
{
    // Reference: Ref: https://www.qc.edu.hk/math/Advanced%20Level/Point_to_line.htm "Method 3 Using Dot Product"
    const Point point_to_proj(-2.0, 4.0, -3.0);
    Point projected_point;

    const double proj_distance = GeometricalProjectionUtilities::FastProjectOnLine(
        rGeom,
        point_to_proj,
        projected_point);

    const double expected_proj_dist = std::sqrt(171.0/14.0);

    const double expected_coord_x = 4.0/14.0;
    const double expected_coord_y = 27.0/14.0;
    const double expected_coord_z = -19.0/14.0;

    KRATOS_EXPECT_DOUBLE_EQ(expected_proj_dist, proj_distance);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.X(), expected_coord_x);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Y(), expected_coord_y);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Z(), expected_coord_z);
}

template<class TGeometryType>
void TestFastMinimalDistanceOnLine(TGeometryType& rGeom)
{
    const Point point(2.0, 5.0, 0.5);

    const double proj_distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(rGeom, point);

    Point projected_point;
    const double expected_proj_dist = GeometricalProjectionUtilities::FastProjectOnLine(rGeom,point, projected_point);

    KRATOS_EXPECT_DOUBLE_EQ(expected_proj_dist, proj_distance);

    const Point non_inside_projected_point(-2.0, 4.0, -3.0);

    const double non_inside_projected_distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(rGeom, non_inside_projected_point);

    const double expected_non_inside_projected_dist = 3.74166;

    KRATOS_EXPECT_RELATIVE_NEAR(expected_non_inside_projected_dist, non_inside_projected_distance, 1.0e-4);

    const Point far_point(100.0, 100.0, 100.0);

    const double far_distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(rGeom, far_point);

    const double expected_far_dist = 168.062;

    KRATOS_EXPECT_RELATIVE_NEAR(expected_far_dist, far_distance, 1.0e-4);
}

}

KRATOS_TEST_CASE_IN_SUITE(GeometricalProjectionUtilitiesFastProjectDirectionNode, KratosCoreFastSuite)
{
    GeometryNodeType::Pointer p_geom = CreateTriangle3D3NForTestNode();

    TestFastProjectDirection(*p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalProjectionUtilitiesFastProjectDirectionPoint, KratosCoreFastSuite)
{
    GeometryPointType::Pointer p_geom = CreateTriangle3D3NForTestPoint();

    TestFastProjectDirection(*p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalProjectionUtilitiesFastProjectOnGeometryNode, KratosCoreFastSuite)
{
    GeometryNodeType::Pointer p_geom = CreateTriangle3D3NForTestNode();

    TestFastProjectOnGeometry(*p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalProjectionUtilitiesFastProjectOnGeometryPoint, KratosCoreFastSuite)
{
    GeometryPointType::Pointer p_geom = CreateTriangle3D3NForTestPoint();

    TestFastProjectOnGeometry(*p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalProjectionUtilitiesFastProjectOnLineNode, KratosCoreFastSuite)
{
    GeometryNodeType::Pointer p_geom = CreateLine3D2NForTestNode2D();
    TestFastProjectOnLine2D(*p_geom);

    p_geom = CreateLine3D2NForTestNode3D();
    TestFastProjectOnLine3D(*p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalProjectionUtilitiesFastProjectOnLinePoint, KratosCoreFastSuite)
{
    GeometryPointType::Pointer p_geom = CreateLine3D2NForTestPoint2D();
    TestFastProjectOnLine2D(*p_geom);

    p_geom = CreateLine3D2NForTestPoint3D();
    TestFastProjectOnLine3D(*p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalProjectionUtilitiesFastMinimalDistanceOnLine, KratosCoreFastSuite)
{
    GeometryPointType::Pointer p_geom = CreateLine3D2NForTestPoint3D();
    TestFastMinimalDistanceOnLine(*p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalProjectionUtilitiesFastProject, KratosCoreFastSuite)
{
    const double expected_proj_dist = 1.258;

    const double x_coord = 0.325;
    const double y_coord = 0.147;

    const Point point_in_plane(-1.274, 10.478, 0.0);
    const Point point_to_proj(x_coord, y_coord, expected_proj_dist);

    array_1d<double,3> normal_vector;

    normal_vector[0] = 0.0;
    normal_vector[1] = 0.0;
    normal_vector[2] = 1.0;

    double proj_distance;

    Point projected_point = GeometricalProjectionUtilities::FastProject(
        point_in_plane,
        point_to_proj,
        normal_vector,
        proj_distance);

    KRATOS_EXPECT_DOUBLE_EQ(expected_proj_dist, proj_distance);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.X(), x_coord);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Y(), y_coord);
    KRATOS_EXPECT_DOUBLE_EQ(projected_point.Z(), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(FastMinimalDistanceOnLineWithRadius, KratosCoreFastSuite)
{
    constexpr double TOLERANCE_DISTANCE_PATH = 1e-6;

    double distance;
    double radius = 0.0;
    auto line = Kratos::make_shared<Line3D2<Node>>(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0)
    );
    Point point1(0.0,0.0,0.1);
    distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*line, point1);
    auto distance_computed_type = GeometricalProjectionUtilities::GeometricalProjectionUtilities::FastMinimalDistanceOnLineWithRadius(distance, *line, point1, radius);
    KRATOS_EXPECT_NEAR(distance, 0.1, TOLERANCE_DISTANCE_PATH);
    KRATOS_EXPECT_EQ(distance_computed_type, GeometricalProjectionUtilities::DistanceComputed::NO_RADIUS);

    radius = 0.01;
    Point point2(0.0,0.0,0.1);
    distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*line, point2);
    distance_computed_type = GeometricalProjectionUtilities::GeometricalProjectionUtilities::FastMinimalDistanceOnLineWithRadius(distance, *line, point2, radius);
    KRATOS_EXPECT_NEAR(distance, 0.09, TOLERANCE_DISTANCE_PATH);
    KRATOS_EXPECT_EQ(distance_computed_type, GeometricalProjectionUtilities::DistanceComputed::RADIUS_PROJECTED);

    Point point3(-0.1,0.0,0.1);
    distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*line, point3);
    distance_computed_type = GeometricalProjectionUtilities::GeometricalProjectionUtilities::FastMinimalDistanceOnLineWithRadius(distance, *line, point3, radius);
    KRATOS_EXPECT_NEAR(distance, std::sqrt(std::pow(0.09, 2) * 2), TOLERANCE_DISTANCE_PATH);
    KRATOS_EXPECT_EQ(distance_computed_type, GeometricalProjectionUtilities::DistanceComputed::RADIUS_NOT_PROJECTED_OUTSIDE);

    radius = 0.1;
    Point point4(-0.1,0.0,0.09);
    distance = GeometricalProjectionUtilities::FastMinimalDistanceOnLine(*line, point4);
    distance_computed_type = GeometricalProjectionUtilities::GeometricalProjectionUtilities::FastMinimalDistanceOnLineWithRadius(distance, *line, point4, radius);
    KRATOS_EXPECT_NEAR(distance, -(std::sqrt(std::pow(0.01, 2) + std::pow(0.1, 2))), TOLERANCE_DISTANCE_PATH);
    KRATOS_EXPECT_EQ(distance_computed_type, GeometricalProjectionUtilities::DistanceComputed::RADIUS_NOT_PROJECTED_INSIDE);
}

} // namespace Kratos::Testing.
