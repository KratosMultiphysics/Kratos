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

#include "geometries/triangle_3d_3.h"

namespace Kratos
{
namespace Testing
{

using NodeType = typename GeometricalProjectionUtilities::NodeType;
using GeometryType = typename GeometricalProjectionUtilities::GeometryType;

GeometryType::Pointer CreateTriangle2D3NForTest()
{
    GeometryType::PointsArrayType points;
    points.push_back(Kratos::make_intrusive<NodeType>(1,0.04, 0.02, 0.0));
    points.push_back(Kratos::make_intrusive<NodeType>(2,1.1, 0.03, 0.0));
    points.push_back(Kratos::make_intrusive<NodeType>(3,1.08, 1.0, 0.0));

    return GeometryType::Pointer(new Triangle3D3<NodeType>(points));
}

KRATOS_TEST_CASE_IN_SUITE(GeometricalProjectionUtilitiesFastProjectDirection, KratosCoreFastSuite)
{
    GeometryType::Pointer p_geom = CreateTriangle2D3NForTest();

    double expected_proj_dist = 1.258;

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
        *p_geom,
        point_to_proj,
        projected_point,
        normal_vector,
        dir_vector);

    KRATOS_CHECK_DOUBLE_EQUAL(expected_proj_dist, proj_distance);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.X(), x_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Y(), y_coord);
    KRATOS_CHECK_DOUBLE_EQUAL(projected_point.Z(), 0.0);
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
