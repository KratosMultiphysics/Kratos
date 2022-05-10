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
#include "includes/stream_serializer.h"
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

KRATOS_TEST_CASE_IN_SUITE(PointWithIdCopyConstructor, KratosMappingApplicationSerialTestSuite)
{
    const std::size_t id=36;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};
    const double dist=1.236;

    PointWithId point(id, coords, dist);

    KRATOS_CHECK_EQUAL(id, point.GetId());
    KRATOS_CHECK_VECTOR_EQUAL(coords, point.Coordinates());
    KRATOS_CHECK_DOUBLE_EQUAL(dist, point.GetDistance());

    PointWithId copied_point(point);

    KRATOS_CHECK_EQUAL(point.GetId(), copied_point.GetId());
    KRATOS_CHECK_VECTOR_EQUAL(point.Coordinates(), copied_point.Coordinates());
    KRATOS_CHECK_DOUBLE_EQUAL(point.GetDistance(), copied_point.GetDistance());
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

    PointWithId point_1(id, Point(1), dist);
    PointWithId point_2(id, Point(2), dist+0.01);
    PointWithId point_3(id, Point(3), dist-0.01);

    KRATOS_CHECK_LESS(point_1, point_2);
    KRATOS_CHECK(point_1 < point_2);
    KRATOS_CHECK_IS_FALSE(point_1 < point_3);
}

KRATOS_TEST_CASE_IN_SUITE(PointWithIdSerialization, KratosMappingApplicationSerialTestSuite)
{
    const std::size_t id=36;
    const double dist=1.236;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};

    PointWithId point_1(id, coords, dist);
    PointWithId point_new(0, Point(), 0);

    // serializing the object
    StreamSerializer serializer;
    serializer.save("obj", point_1);
    serializer.load("obj", point_new);

    KRATOS_CHECK_EQUAL(point_1.GetId(), point_new.GetId());
    KRATOS_CHECK_VECTOR_EQUAL(point_1.Coordinates(), point_new.Coordinates());
    KRATOS_CHECK_DOUBLE_EQUAL(point_1.GetDistance(), point_new.GetDistance());

    KRATOS_CHECK_EQUAL(point_1, point_new);
}

KRATOS_TEST_CASE_IN_SUITE(ClosestPointsContainerBasics, KratosMappingApplicationSerialTestSuite)
{
    ClosestPointsContainer cont(2);

    KRATOS_CHECK_EQUAL(cont.GetPoints().size(), 0);

    const std::size_t id=36;
    const double dist=1.236;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};

    PointWithId point_1(id, coords, dist);

    cont.Add(point_1);

    KRATOS_CHECK_EQUAL(cont.GetPoints().size(), 1);

    const PointWithId& retrieved_point = *cont.GetPoints().begin();

    KRATOS_CHECK_EQUAL(point_1.GetId(), retrieved_point.GetId());
    KRATOS_CHECK_VECTOR_EQUAL(point_1.Coordinates(), retrieved_point.Coordinates());
    KRATOS_CHECK_DOUBLE_EQUAL(point_1.GetDistance(), retrieved_point.GetDistance());
}

KRATOS_TEST_CASE_IN_SUITE(ClosestPointsContainerAddEqualPoints, KratosMappingApplicationSerialTestSuite)
{
    ClosestPointsContainer cont(2);

    const std::size_t id=36;
    const double dist=1.236;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};

    PointWithId point_1(id, coords, dist);
    PointWithId point_2(id, coords, dist+1);

    cont.Add(point_1);
    cont.Add(point_1); // this must not add the same point again!

    KRATOS_CHECK_EQUAL(cont.GetPoints().size(), 1);

    cont.Add(point_2); // this must also not add the point again!

    KRATOS_CHECK_EQUAL(cont.GetPoints().size(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ClosestPointsContainerAddMultiplePoints, KratosMappingApplicationSerialTestSuite)
{
    constexpr int max_size = 3;
    ClosestPointsContainer cont(max_size);

    const std::size_t id=36;

    PointWithId point_1(id, Point(1), 1);
    PointWithId point_2(id, Point(2), 2);
    PointWithId point_3(id, Point(3), 3);
    PointWithId point_4(id, Point(4), 4);

    cont.Add(point_3);
    cont.Add(point_2);
    cont.Add(point_1);
    cont.Add(point_4);

    KRATOS_CHECK_EQUAL(cont.GetPoints().size(), max_size);

    // check if are points are ordered by distance
    int counter = 1;
    for (const auto& r_point : cont.GetPoints()) {
        KRATOS_CHECK_DOUBLE_EQUAL(counter*1.0, r_point.GetDistance());
        counter++;
    }
}

KRATOS_TEST_CASE_IN_SUITE(ClosestPointsContainerMerge, KratosMappingApplicationSerialTestSuite)
{
    constexpr int max_size = 3;
    ClosestPointsContainer cont(max_size);
    ClosestPointsContainer cont2(max_size);

    const std::size_t id=36;

    PointWithId point_1(id, Point(1), 1);
    PointWithId point_2(id, Point(2), 2);
    PointWithId point_3(id, Point(3), 3);
    PointWithId point_4(id, Point(4), 4);

    PointWithId point_15(id, Point(1.5), 1.5);
    PointWithId point_25(id, Point(2.5), 2.5);
    PointWithId point_35(id, Point(3.5), 3.5);
    PointWithId point_45(id, Point(4.5), 4.5);

    cont.Add(point_3);
    cont.Add(point_2);
    cont.Add(point_1);
    cont.Add(point_4);

    cont2.Add(point_35);
    cont2.Add(point_25);
    cont2.Add(point_15);
    cont2.Add(point_45);

    KRATOS_CHECK_EQUAL(cont.GetPoints().size(), max_size);
    KRATOS_CHECK_EQUAL(cont2.GetPoints().size(), max_size);

    cont.Merge(cont2);
    KRATOS_CHECK_EQUAL(cont.GetPoints().size(), max_size);

    // check if are points are ordered by distance
    int counter = 2;
    for (const auto& r_point : cont.GetPoints()) {
        KRATOS_CHECK_DOUBLE_EQUAL(counter*0.5, r_point.GetDistance());
        counter++;
    }
}

KRATOS_TEST_CASE_IN_SUITE(ClosestPointsContainerSerialization, KratosMappingApplicationSerialTestSuite)
{
    constexpr int max_size = 3;
    ClosestPointsContainer cont(max_size);

    const std::size_t id=36;

    PointWithId point_1(id, Point(1), 1);
    PointWithId point_2(id, Point(2), 2);
    PointWithId point_3(id, Point(3), 3);
    PointWithId point_4(id, Point(4), 4);

    cont.Add(point_3);
    cont.Add(point_2);
    cont.Add(point_1);
    cont.Add(point_4);

    KRATOS_CHECK_EQUAL(cont.GetPoints().size(), max_size);

    // check if are points are ordered by distance
    int counter = 1;
    for (const auto& r_point : cont.GetPoints()) {
        KRATOS_CHECK_DOUBLE_EQUAL(counter*1.0, r_point.GetDistance());
        counter++;
    }

    // serializing the object
    StreamSerializer serializer;
    serializer.save("obj", cont);
    ClosestPointsContainer cont_new(0);
    serializer.load("obj", cont_new);

    KRATOS_CHECK_EQUAL(cont.GetPoints().size(), cont_new.GetPoints().size());
    // check if are points are ordered by distance
    counter = 1;
    for (const auto& r_point : cont.GetPoints()) {
        KRATOS_CHECK_DOUBLE_EQUAL(counter*1.0, r_point.GetDistance());
        counter++;
    }
}

}  // namespace Testing
}  // namespace Kratos