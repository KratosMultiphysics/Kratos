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
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"

// // Debugging
// #include "input_output/vtk_output.h"

namespace Kratos::Testing {

// namespace
// {
//     // This method is for debugging
//     void WriteVTK(ModelPart& rModelPart)
//     {
//         auto volume_vtk_output_parameters = Parameters(R"({
//             "model_part_name"                    : "})"+rModelPart.Name()+R"({",
//             "file_format"                        : "ascii",
//             "output_precision"                   : 8,
//             "output_path"                        : ")"+rModelPart.Name()+R"(",
//             "nodal_data_value_variables"         : [],
//             "nodal_flags"                        : [],
//             "condition_data_value_variables"     : [],
//             "condition_flags"                    : []
//         })");
//         VtkOutput volume_vtk_output(rModelPart, volume_vtk_output_parameters);
//         volume_vtk_output.PrintOutput();
//     }
// }

/**
 * @brief Finds the nearest point to the given point on a line segment. 
 * @details It first projects the point into the line. If the projected point is inside the segment boundary 
 * it returns the projected point. If not it returns the nearest end point of the line.
 * @tparam Type of the Point 
 * @tparam TGeometryType The type of the line. Assumes to have [] access and IsInside method
 * @param rPoint The query point which we want to get nearest point to it on the line
 * @param rLine The line in which we want to find the nearest point to rPoint
 * @return The nearest point to rPoint
 */
template<class TPointType, class TGeometryType>
static Point LegacyLineNearestPoint(
    const TPointType& rPoint, 
    const TGeometryType& rLine
    )
{
    KRATOS_DEBUG_ERROR_IF_NOT(rLine.size() == 2) << "This function only accepts Line2D2 as input" << std::endl;
    Point result;
    const Point line_projected_point = GeometricalProjectionUtilities::FastProjectOnLine(rLine, rPoint, result);
    array_1d<double,3> projected_local;
    if(rLine.IsInside(result.Coordinates(), projected_local))
        return result;

    const double distance1 = norm_2(rLine[0] - result);
    const double distance2 = norm_2(rLine[1] - result);

    result = (distance1 < distance2) ? rLine[0] : rLine[1];
    return result;
}


/**
 * @brief Finds the nearest point to the given point on a triangle. 
 * @details It first projects the point into the triangle surface. If the projected point is inside the triangle 
 * it returns the projected point. If not it returns the nearest point on the edges of the triangle.
 * Dividing the plane of the triangle in 7 zones and find the nearest reflecting those zones
 *
 *             \       /
 *              \  6  /
 *               \   /
 *                \ /
 *                 /\ 
 *                /  \ 
 *          2    /    \     3
 *              /   1  \ 
 *             /        \ 
 *    _______ /__________\_____________ 
 *           /            \ 
 *     5    /       4      \    7
 *         /                \ 
 *        /                  \ 
 * @tparam Type of the Point 
 * @tparam TGeometryType The type of the triangle. Assumes to have [] access and IsInside method
 * @param rPoint The query point which we want to get nearest point to it on the line
 * @param rTriangle The triangle in which we want to find the nearest point to rPoint
 * @return The nearest point to rPoint
 */
template<class TPointType, class TGeometryType>
static Point LegacyTriangleNearestPoint(
    const TPointType& rPoint,
    const TGeometryType& rTriangle
    )
{
    constexpr double Tolerance = 1e-12;
    Point result;
    array_1d<double,3> old_coordinates = ZeroVector(3);
    const Point center = rTriangle.Center();
    const array_1d<double, 3> normal = rTriangle.UnitNormal(old_coordinates);
    double distance = 0.0;
    const Point point_projected = GeometricalProjectionUtilities::FastProject( center, rPoint, normal, distance);
    rTriangle.PointLocalCoordinates(old_coordinates, point_projected);
    using line_point_type= typename TGeometryType::PointType;

    if(old_coordinates[0] < -Tolerance) { // case 2,5,6
        if(old_coordinates[1] < -Tolerance) { // case 5
            result = rTriangle[0];
        } else if ((old_coordinates[0] + old_coordinates[1]) > (1.0+Tolerance)) { // case 6
            result = rTriangle[2];
        } else {
            result = LegacyLineNearestPoint(rPoint, Line3D2<line_point_type>(rTriangle.pGetPoint(0), rTriangle.pGetPoint(2)));
        }
    } else if(old_coordinates[1] < -Tolerance) { // case 4,7 (case 5 is already covered in previous if)
        if ((old_coordinates[0] + old_coordinates[1]) > (1.0+Tolerance)) { // case 7
            result = rTriangle[1];
        } else { // case 4
            result = LegacyLineNearestPoint(rPoint, Line3D2<line_point_type>(rTriangle.pGetPoint(0), rTriangle.pGetPoint(1)));
        }
    } else if ((old_coordinates[0] + old_coordinates[1]) > (1.0+Tolerance)) { // case 3
        result = LegacyLineNearestPoint(rPoint, Line3D2<line_point_type>(rTriangle.pGetPoint(1), rTriangle.pGetPoint(2)));
    } else {  // inside
        result = point_projected;
    }

    return result;
}

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
    auto scenario = NearestPointUtilities::TriangleNearestPoint(inside_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, inside_point, 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::INSIDE_TRIANGLE);

    Point border_point(0.2 * length, 0.2 * length, 0.00);
    scenario = NearestPointUtilities::TriangleNearestPoint(border_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, border_point, 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::INSIDE_TRIANGLE);

    Point lower_corner(-2 * length, -0.1 * length, 0.00);
    scenario = NearestPointUtilities::TriangleNearestPoint(lower_corner, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_1), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_0);

    Point below_point(0.2 * length, -0.2 * length, 0.00);
    scenario = NearestPointUtilities::TriangleNearestPoint(below_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.0, 0.00), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_EDGE_20);

    Point below_left_point(-0.1 * length, -0.2 * length, 0.00);
    scenario = NearestPointUtilities::TriangleNearestPoint(below_left_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_1), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_0);

    Point below_point_right(1.2 * length, -0.2 * length, 0.00);
    scenario = NearestPointUtilities::TriangleNearestPoint(below_point_right, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_2), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_2);

    Point right_side(1.2 * length, 0.2 * length, 0.00);
    scenario = NearestPointUtilities::TriangleNearestPoint(right_side, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(length, 0.2 * length, 0.00), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_EDGE_12);

    Point upper_corner(1.2 * length, 1.2 * length, 0.00);
    scenario = NearestPointUtilities::TriangleNearestPoint(upper_corner, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_3), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_1);

    Point left_point(0.2 * length, 0.4 * length, 0.00);
    scenario = NearestPointUtilities::TriangleNearestPoint(left_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.3 * length, 0.3 * length, 0.00), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_EDGE_01);

    Point far_left_point(-0.2 * length, 0.2 * length, 0.00);
    scenario = NearestPointUtilities::TriangleNearestPoint(far_left_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.00, 0.00, 0.00), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_EDGE_01);
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
    auto scenario = NearestPointUtilities::TriangleNearestPoint(inside_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.1 * length, 0.00), 1e-6)
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::INSIDE_TRIANGLE);

    Point border_point(0.2 * length, 0.2 * length, distance);
    scenario = NearestPointUtilities::TriangleNearestPoint(border_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.2 * length, 0.00), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::INSIDE_TRIANGLE);

    Point lower_corner(-2 * length, -0.1 * length, distance);
    scenario = NearestPointUtilities::TriangleNearestPoint(lower_corner, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_1), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_0);

    Point below_point(0.2 * length, -0.2 * length, distance);
    scenario = NearestPointUtilities::TriangleNearestPoint(below_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.0, 0.00), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_EDGE_20);

    Point below_left_point(-0.1 * length, -0.2 * length, distance);
    scenario = NearestPointUtilities::TriangleNearestPoint(below_left_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_1), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_0);

    Point below_point_right(1.2 * length, -0.2 * length, distance);
    scenario = NearestPointUtilities::TriangleNearestPoint(below_point_right, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_2), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_2);

    Point right_side(1.2 * length, 0.2 * length, distance);
    scenario = NearestPointUtilities::TriangleNearestPoint(right_side, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(length, 0.2 * length, 0.00), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_EDGE_12);

    Point upper_corner(1.2 * length, 1.2 * length, distance);
    scenario = NearestPointUtilities::TriangleNearestPoint(upper_corner, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, (*p_point_3), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_1);

    Point left_point(0.2 * length, 0.4 * length, distance);
    scenario = NearestPointUtilities::TriangleNearestPoint(left_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.3 * length, 0.3 * length, 0.00), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_EDGE_01);

    Point far_left_point(-0.2 * length, 0.2 * length, distance);
    scenario = NearestPointUtilities::TriangleNearestPoint(far_left_point, triangle, nearest_point);
    KRATOS_EXPECT_VECTOR_NEAR(nearest_point, Point(0.00, 0.00, 0.00), 1e-6);
    KRATOS_EXPECT_EQ(scenario, TriangleNearestPointLocation::ON_TRIANGLE_EDGE_01);
}

KRATOS_TEST_CASE_IN_SUITE(TriangleCornerCasesNearestPoint1, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    auto p_node_0 = r_model_part.CreateNewNode(0, 0.0150564, 0.0125443, -0.00755447);

    auto p_node_1 = r_model_part.CreateNewNode(1, 0.018141,0.013636,-0.000874817);
    auto p_node_2 = r_model_part.CreateNewNode(2, 0.0158935,0.0173839,-0.00353975);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0155049,0.02,-0.000874817);

    auto p_cond = r_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 2, 3}}, nullptr);

    // Compute
    Point nearest_point = NearestPointUtilities::TriangleNearestPoint(*p_node_0, p_cond->GetGeometry());
    const double distance = norm_2(nearest_point - *p_node_0);
    Point old_nearest_point = LegacyTriangleNearestPoint(*p_node_0, p_cond->GetGeometry());
    const double old_distance = norm_2(old_nearest_point - *p_node_0);

    // Check
    KRATOS_CHECK_LESS_EQUAL(distance, old_distance);
    // KRATOS_EXPECT_VECTOR_NEAR(nearest_point, old_nearest_point, 1e-6);

    // // Debugging
    // WriteVTK(r_model_part);
    // ModelPart& r_model_part_solution = current_model.CreateModelPart("Solution");
    // r_model_part_solution.CreateNewNode(1, nearest_point[0], nearest_point[1], nearest_point[2]);
    // WriteVTK(r_model_part_solution);
    // ModelPart& r_model_part_old_solution = current_model.CreateModelPart("LegacySolution");
    // r_model_part_old_solution.CreateNewNode(1, old_nearest_point[0], old_nearest_point[1], old_nearest_point[2]);
    // WriteVTK(r_model_part_old_solution);
}

KRATOS_TEST_CASE_IN_SUITE(TriangleCornerCasesNearestPoint2, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    auto p_node_0 = r_model_part.CreateNewNode(0, 0.0150564, 0.0125443, -0.00755447);

    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0158935, 0.0173839, -0.00353975);
    auto p_node_2 = r_model_part.CreateNewNode(2, 0.0155049, 0.02, -0.00814909);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0155049, 0.02, -0.000874817);

    auto p_cond = r_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 2, 3}}, nullptr);

    // Compute
    Point nearest_point = NearestPointUtilities::TriangleNearestPoint(*p_node_0, p_cond->GetGeometry());
    const double distance = norm_2(nearest_point - *p_node_0);
    Point old_nearest_point = LegacyTriangleNearestPoint(*p_node_0, p_cond->GetGeometry());
    const double old_distance = norm_2(old_nearest_point - *p_node_0);

    // Check
    KRATOS_CHECK_LESS_EQUAL(distance, old_distance);
    // KRATOS_EXPECT_VECTOR_NEAR(nearest_point, old_nearest_point, 1e-6);

    // // Debugging
    // WriteVTK(r_model_part);
    // ModelPart& r_model_part_solution = current_model.CreateModelPart("Solution");
    // r_model_part_solution.CreateNewNode(1, nearest_point[0], nearest_point[1], nearest_point[2]);
    // WriteVTK(r_model_part_solution);
    // ModelPart& r_model_part_old_solution = current_model.CreateModelPart("LegacySolution");
    // r_model_part_old_solution.CreateNewNode(1, old_nearest_point[0], old_nearest_point[1], old_nearest_point[2]);
    // WriteVTK(r_model_part_old_solution);
}

KRATOS_TEST_CASE_IN_SUITE(TriangleCornerCasesNearestPoint3, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    auto p_node_0 = r_model_part.CreateNewNode(0, 0.0129223, 0.0167283, -0.00755451);

    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0155049, 0.02, -0.00814909);
    auto p_node_2 = r_model_part.CreateNewNode(2, 0.0156193, 0.01857, -0.0129188);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0155632, 0.0189773, -0.015374);

    auto p_cond = r_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, {{1, 2, 3}}, nullptr);

    // Compute
    Point nearest_point = NearestPointUtilities::TriangleNearestPoint(*p_node_0, p_cond->GetGeometry());
    const double distance = norm_2(nearest_point - *p_node_0);
    Point old_nearest_point = LegacyTriangleNearestPoint(*p_node_0, p_cond->GetGeometry());
    const double old_distance = norm_2(old_nearest_point - *p_node_0);

    // Check
    KRATOS_CHECK_LESS_EQUAL(distance, old_distance);
    // KRATOS_EXPECT_VECTOR_NEAR(nearest_point, old_nearest_point, 1e-6);

    // // Debugging
    // WriteVTK(r_model_part);
    // ModelPart& r_model_part_solution = current_model.CreateModelPart("Solution");
    // r_model_part_solution.CreateNewNode(1, nearest_point[0], nearest_point[1], nearest_point[2]);
    // WriteVTK(r_model_part_solution);
    // ModelPart& r_model_part_old_solution = current_model.CreateModelPart("LegacySolution");
    // r_model_part_old_solution.CreateNewNode(1, old_nearest_point[0], old_nearest_point[1], old_nearest_point[2]);
    // WriteVTK(r_model_part_old_solution);
}

} // namespace Kratos::Testing