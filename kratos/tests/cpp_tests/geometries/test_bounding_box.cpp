//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/bounding_box.h"
#include "geometries/point.h"

namespace Kratos::Testing
{

/** Checks BoundingBox creation with min max point
*/
KRATOS_TEST_CASE_IN_SUITE(BoundingBoxMinMaxConstruction, KratosCoreGeometriesFastSuite) 
{
    constexpr double tolerance = 1e-12;
    BoundingBox<Point> bounding_box({0, .1, .3}, {2.1, 3.4, 5.6});


    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[0], 0.00, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[1], 0.10, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[2], 0.30, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[0], 2.10, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[1], 3.40, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[2], 5.60, tolerance);
}


/** Checks BoundingBox copy and assingment
*/
KRATOS_TEST_CASE_IN_SUITE(BoundingCopyAndAssignment, KratosCoreGeometriesFastSuite) 
{
    constexpr double tolerance = 1e-12;
    BoundingBox<Point> bounding_box({0, .1, .3}, {2.1, 3.4, 5.6});

    BoundingBox<Point> copied_box(bounding_box);
    BoundingBox<Point> assigned_box({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});

    assigned_box = copied_box;


    KRATOS_EXPECT_NEAR(assigned_box.GetMinPoint()[0], 0.00, tolerance);
    KRATOS_EXPECT_NEAR(assigned_box.GetMinPoint()[1], 0.10, tolerance);
    KRATOS_EXPECT_NEAR(assigned_box.GetMinPoint()[2], 0.30, tolerance);
    KRATOS_EXPECT_NEAR(assigned_box.GetMaxPoint()[0], 2.10, tolerance);
    KRATOS_EXPECT_NEAR(assigned_box.GetMaxPoint()[1], 3.40, tolerance);
    KRATOS_EXPECT_NEAR(assigned_box.GetMaxPoint()[2], 5.60, tolerance);
}

/** Checks BoundingBox creation with iterator of points
*/
KRATOS_TEST_CASE_IN_SUITE(BoundingPointsConstruction, KratosCoreGeometriesFastSuite) 
{
    constexpr double tolerance = 1e-12;

    std::vector<Point> points;
    points.push_back(Point{0.0, 0.4, -0.3});
    points.push_back(Point{0.9, -0.3, 0.1});
    points.push_back(Point{-0.2, 0.4, 0.3});
    points.push_back(Point{0.0, 0.8, -0.5});

    BoundingBox<Point> bounding_box(points.begin(), points.end());

    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[0],-0.20, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[1],-0.30, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[2],-0.50, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[0], 0.90, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[1], 0.80, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[2], 0.30, tolerance);
}

/** Checks BoundingBox set method
*/
KRATOS_TEST_CASE_IN_SUITE(BoundingBoxSet, KratosCoreGeometriesFastSuite) 
{
    constexpr double tolerance = 1e-12;

    std::vector<Point> points;
    points.push_back(Point{0.0, 0.4, -0.3});
    points.push_back(Point{0.9, -0.3, 0.1});
    points.push_back(Point{-0.2, 0.4, 0.3});
    points.push_back(Point{0.0, 0.8, -0.5});

    BoundingBox<Point> bounding_box;

    bounding_box.Set(points.begin(), points.end());

    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[0],-0.20, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[1],-0.30, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[2],-0.50, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[0], 0.90, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[1], 0.80, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[2], 0.30, tolerance);
}

/** Checks BoundingBox extend method
*/
KRATOS_TEST_CASE_IN_SUITE(BoundingBoxExtend, KratosCoreGeometriesFastSuite)
{
    constexpr double tolerance = 1e-12;

    std::vector<Point> points;
    points.push_back(Point{0.0, 0.4, -0.3});
    points.push_back(Point{0.9, -0.3, 0.1});

    BoundingBox<Point> bounding_box(points.begin(), points.end());

    points[0] = Point{-0.2, 0.4, 0.3};
    points[1] = Point{0.0, 0.8, -0.5};

    bounding_box.Extend(points.begin(), points.end());

    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[0],-0.20, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[1],-0.30, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMinPoint()[2],-0.50, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[0], 0.90, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[1], 0.80, tolerance);
    KRATOS_EXPECT_NEAR(bounding_box.GetMaxPoint()[2], 0.30, tolerance);
}

/** Checks BoundingBox IsInside method
*/
KRATOS_TEST_CASE_IN_SUITE(BoundingBoxIsInside, KratosCoreGeometriesFastSuite) 
{
    BoundingBox<Point> bounding_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});

    Point point_inside{0.5, 0.5, 0.5};
    Point point_outside{1.5, 0.5, 0.5};

    KRATOS_CHECK(bounding_box.IsInside(point_inside));
    KRATOS_CHECK_IS_FALSE(bounding_box.IsInside(point_outside));
}

} // namespace Kratos::Testing.
