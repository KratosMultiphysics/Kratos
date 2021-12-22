//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// Project includes
#include "testing/testing.h"
#include "custom_utilities/closest_points.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(PointWithIdBasics, KratosMappingApplicationSerialTestSuite)
{
    const std::size_t id=36;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};
    const double dist=1.236;

    PointWithId point(id, coords, dist);

    KRATOS_CHECK_EQUAL(id, point.GetId());
    KRATOS_CHECK_VECTOR_EQUAL(coords, point.Coordinates());
    KRATOS_CHECK_DOUBLE_EQUAL(dist, point.GetDistance());
}

KRATOS_TEST_CASE_IN_SUITE(PointWithIdEqualComparison, KratosMappingApplicationSerialTestSuite)
{
    const std::size_t id=36;
    const double dist=1.236;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};
    Point::CoordinatesArrayType coords2(coords);
    coords2[2]+=1E-10;

    PointWithId point_1(id, coords, dist);
    PointWithId point_2(id, coords, dist);
    PointWithId point_3(id, coords2, dist);

    // only the position aka the coordinates are used for the equal comparison!
    PointWithId point_4(id+1, coords, dist);
    PointWithId point_5(id, coords, dist+1.0);
    PointWithId point_6(id+1, coords, dist+1.0);

    KRATOS_CHECK_EQUAL(point_1, point_2);
    KRATOS_CHECK_NOT_EQUAL(point_1, point_3);

    KRATOS_CHECK_EQUAL(point_1, point_4);
    KRATOS_CHECK_EQUAL(point_1, point_5);
    KRATOS_CHECK_EQUAL(point_1, point_6);
}

KRATOS_TEST_CASE_IN_SUITE(PointWithIdLessComparison, KratosMappingApplicationSerialTestSuite)
{
    const std::size_t id=36;
    const double dist=1.236;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};
    Point::CoordinatesArrayType coords2(coords);
    coords2[2]+=1E-10;

    PointWithId point_1(id, coords, dist);
    PointWithId point_2(id, coords, dist+0.01);
    PointWithId point_3(id, coords, dist-0.01);

    KRATOS_CHECK_LESS(point_1, point_2);
    KRATOS_CHECK(point_1 < point_2);
    KRATOS_CHECK_IS_FALSE(point_1 < point_3);
}

}  // namespace Testing
}  // namespace Kratos