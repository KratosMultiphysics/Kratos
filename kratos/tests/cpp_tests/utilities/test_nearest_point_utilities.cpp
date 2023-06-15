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
#include "geometries/line_2d_2.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"


namespace Kratos::Testing {


    KRATOS_TEST_CASE_IN_SUITE(LineNearestPoint, KratosCoreFastSuite)
    {
        constexpr double length = 1.2;

        Point::Pointer p_point_1(make_shared<Point>( 0.00, 0.00, 0.00));
        Point::Pointer p_point_2(make_shared<Point>( length, 0.00, 0.00));
        Point::Pointer p_point_3(make_shared<Point>( 0.00, length, 0.00));


        Point nearest_point( 0.00, 0.00, 0.00);

        Line2D2<Point> line_1(p_point_1, p_point_2);
        Line2D2<Point> line_2(p_point_1, p_point_3);

        Point point_1(0.2 * length, 0.1 * length, 0.00);
        nearest_point = NearestPointUtilities::LineNearestPoint(point_1, line_1);
        KRATOS_CHECK_VECTOR_NEAR(nearest_point, Point(0.2 * length, 0.00, 0.00), 1e-6);

        nearest_point = NearestPointUtilities::LineNearestPoint(point_1, line_2);
        KRATOS_CHECK_VECTOR_NEAR(nearest_point, Point(0.00, 0.1 * length, 0.00), 1e-6);

        Point point_2(-0.2 * length, -0.1 * length, 0.00);
        nearest_point = NearestPointUtilities::LineNearestPoint(point_2, line_1);
        KRATOS_CHECK_VECTOR_NEAR(nearest_point,Point(0.00, 0.00, 0.00), 1e-6);

        nearest_point = NearestPointUtilities::LineNearestPoint(point_2, line_2);
        KRATOS_CHECK_VECTOR_NEAR(nearest_point, Point(0.00, 0.00, 0.00), 1e-6);

        Point point_3(1.2 * length, 1.1 * length, 0.00);
        nearest_point = NearestPointUtilities::LineNearestPoint(point_3, line_1);
        KRATOS_CHECK_VECTOR_NEAR(nearest_point,Point(length, 0.00, 0.00), 1e-6);

        nearest_point = NearestPointUtilities::LineNearestPoint(point_3, line_2);
        KRATOS_CHECK_VECTOR_NEAR(nearest_point, Point(0.00, length, 0.00), 1e-6);     
    }

} // namespace Kratos::Testing
