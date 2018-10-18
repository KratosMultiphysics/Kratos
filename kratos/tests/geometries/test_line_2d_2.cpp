//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_2d_2.h"
#include "tests/geometries/test_geometry.h"


// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos {
namespace Testing {

  /// Factory functions

  /** Generates a point type sample Line2D2N
   * @return  Pointer to a Line2D2N
   */
  Line2D2<Point>::Pointer GeneratePointsUnitXDirectionLine2D2() {
    return Line2D2<Point>::Pointer(new Line2D2<Point>(
      Point::Pointer(new Point(0.0, 0.0, 0.0)),
      Point::Pointer(new Point(1.0, 0.0, 0.0))
    ));
  }

  /** Generates a point type sample Line2D2N
   * @return  Pointer to a Line2D2N
   */
  Line2D2<Point>::Pointer GeneratePointsUnitYDirectionLine2D2() {
    return Line2D2<Point>::Pointer(new Line2D2<Point>(
      Point::Pointer(new Point(0.0, 0.0, 0.0)),
      Point::Pointer(new Point(0.0, 1.0, 0.0))
    ));
  }
  /** Generates a point type sample Line2D2N.
   * @return  Pointer to a Line2D2N
   */
  Line2D2<Point>::Pointer GenerateLine2D2WithPoints(Point::Pointer rPointOne, Point::Pointer rPointTwo ) {
    return Line2D2<Point>::Pointer(new Line2D2<Point>(rPointOne, rPointTwo));
  }


    /**
     * Test an overlaping box and line (line has only one node in the box) HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IntersectionBoxSingleNodeInside, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        Point point_1(0.5, -0.1, 0.0);
        Point point_2(1.5, 0.1, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    /**
     * Test an overlaping box and line (line has both nodes in the box) HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IntersectionBoxTwoNodesInside, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        Point point_1(-0.5, -0.1, 0.0);
        Point point_2(1.5, 0.1, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    /**
     * Test an intersection with another line
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IntersectionWithAnotherLine, KratosCoreGeometriesFastSuite) {
        auto geom_1 = GeneratePointsUnitXDirectionLine2D2();
        Point::Pointer point_1 = Point::Pointer(new Point(0.5, 0.5, 0.0));
        Point::Pointer point_2 = Point::Pointer(new Point(0.5, -0.5, 0.0));
        auto geom_2 = GenerateLine2D2WithPoints(point_1, point_2);
        KRATOS_CHECK(geom_1->HasIntersection(*geom_2));
    }

    /**
     * Test a box inside a line HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IntersectionBoxInsideX, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        Point point_1(0.25, -0.1, 0.0);
        Point point_2(0.75, 0.1, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    /**
     * Test a box inside a line HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IntersectionBoxInsideY, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitYDirectionLine2D2();
        Point point_1(-0.1,0.25, 0.0);
        Point point_2(0.1, 0.75, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }
    /**
     * Test a non overlaping box HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2NoIntersectionBox, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        Point point_1(1, 1, 0.0);
        Point point_2(2, 2, 0.0);
        KRATOS_CHECK_IS_FALSE(geom->HasIntersection(point_1, point_2));
    }

} // namespace Testing.
} // namespace Kratos.
