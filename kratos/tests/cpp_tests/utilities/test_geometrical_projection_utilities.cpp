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

namespace Kratos
{
namespace Testing
{

using NodeType = Node<3>;
using GeometryNodeType = Geometry<NodeType>;
using GeometryPointType = Geometry<Point>;

namespace
{

GeometryNodeType::Pointer CreateLine3D2NForTestNode2D()
{
    GeometryNodeType::PointsArrayType points;
    points.push_back(Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0));
    points.push_back(Kratos::make_intrusive<Node<3>>(2, 2.0, 0.0, 0.0));

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
    points.push_back(Kratos::make_intrusive<Node<3>>(1, 1.0, 3.0, -1.0));
    points.push_back(Kratos::make_intrusive<Node<3>>(2, 3.0, 6.0, 0.0));

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
    points.push_back(Kratos::make_intrusive<Node<3>>(1,0.04, 0.02, 0.0));
    points.push_back(Kratos::make_intrusive<Node<3>>(2,1.1, 0.03, 0.0));
    points.push_back(Kratos::make_intrusive<Node<3>>(3,1.08, 1.0, 0.0));

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

    KRATOS_CHECK_DOUBLE_EQUAL(expected_proj_dist, proj_distance);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.X(), x_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Y(), y_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Z(), 0.0);
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

    KRATOS_CHECK_DOUBLE_EQUAL(expected_proj_dist, proj_distance);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.X(), x_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Y(), y_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Z(), 0.0);
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

    KRATOS_CHECK_DOUBLE_EQUAL(expected_proj_dist, proj_distance);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.X(), x_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Y(), 0.0);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Z(), 0.0);

    proj_distance = GeometricalProjectionUtilities::FastProjectOnLine2D(
        rGeom,
        point_to_proj,
        projected_point);

    KRATOS_CHECK_DOUBLE_EQUAL(expected_proj_dist, proj_distance);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.X(), x_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Y(), 0.0);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Z(), 0.0);
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

    KRATOS_CHECK_DOUBLE_EQUAL(expected_proj_dist, proj_distance);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.X(), expected_coord_x);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Y(), expected_coord_y);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Z(), expected_coord_z);
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

    KRATOS_CHECK_DOUBLE_EQUAL(expected_proj_dist, proj_distance);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.X(), x_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Y(), y_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Z(), 0.0);
}

} // namespace Testing.
} // namespace Kratos.
