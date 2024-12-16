//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// Project includes
#include "testing/testing.h"
#include "includes/stream_serializer.h"
#include "custom_utilities/closest_points.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(PointWithIdBasics, KratosMappingApplicationSerialTestSuite)
{
    const std::size_t id=36;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};
    const double dist=1.236;

    PointWithId point(id, coords, dist);

    KRATOS_EXPECT_EQ(id, point.GetId());
    KRATOS_EXPECT_VECTOR_EQ(coords, point.Coordinates());
    KRATOS_EXPECT_DOUBLE_EQ(dist, point.GetDistance());
}

KRATOS_TEST_CASE_IN_SUITE(PointWithIdCopyConstructor, KratosMappingApplicationSerialTestSuite)
{
    const std::size_t id=36;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};
    const double dist=1.236;

    PointWithId point(id, coords, dist);

    KRATOS_EXPECT_EQ(id, point.GetId());
    KRATOS_EXPECT_VECTOR_EQ(coords, point.Coordinates());
    KRATOS_EXPECT_DOUBLE_EQ(dist, point.GetDistance());

    PointWithId copied_point(point);

    KRATOS_EXPECT_EQ(point.GetId(), copied_point.GetId());
    KRATOS_EXPECT_VECTOR_EQ(point.Coordinates(), copied_point.Coordinates());
    KRATOS_EXPECT_DOUBLE_EQ(point.GetDistance(), copied_point.GetDistance());
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

    PointWithId point_4(id+1, coords, dist);
    PointWithId point_5(id, coords, dist+1.0);
    PointWithId point_6(id+1, coords, dist+1.0);

    // Only the position (aka the coordinates) are used for the EQ and NE comparison
    KRATOS_EXPECT_EQ(point_1, point_2);
    KRATOS_EXPECT_NE(point_1, point_3);

    KRATOS_EXPECT_EQ(point_1, point_4);
    KRATOS_EXPECT_EQ(point_1, point_5);
    KRATOS_EXPECT_EQ(point_1, point_6);
}

KRATOS_TEST_CASE_IN_SUITE(PointWithIdLessComparison, KratosMappingApplicationSerialTestSuite)
{
    const std::size_t id=36;
    const double dist=1.236;

    PointWithId point_1(id, Point(1), dist);
    PointWithId point_2(id, Point(2), dist+0.01);
    PointWithId point_3(id, Point(3), dist-0.01);

    KRATOS_EXPECT_LT(point_1, point_2);
    KRATOS_EXPECT_TRUE(point_1 < point_2);
    KRATOS_EXPECT_FALSE(point_1 < point_3);
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

    KRATOS_EXPECT_EQ(point_1.GetId(), point_new.GetId());
    KRATOS_EXPECT_VECTOR_EQ(point_1.Coordinates(), point_new.Coordinates());
    KRATOS_EXPECT_DOUBLE_EQ(point_1.GetDistance(), point_new.GetDistance());

    KRATOS_EXPECT_EQ(point_1, point_new);
}

KRATOS_TEST_CASE_IN_SUITE(ClosestPointsContainerBasics, KratosMappingApplicationSerialTestSuite)
{
    ClosestPointsContainer cont(2);

    KRATOS_EXPECT_EQ(cont.GetPoints().size(), 0);

    const std::size_t id=36;
    const double dist=1.236;
    const Point::CoordinatesArrayType coords{1.1,-2.5,31.09};

    PointWithId point_1(id, coords, dist);

    cont.Add(point_1);

    KRATOS_EXPECT_EQ(cont.GetPoints().size(), 1);

    const PointWithId& retrieved_point = *cont.GetPoints().begin();

    KRATOS_EXPECT_EQ(point_1.GetId(), retrieved_point.GetId());
    KRATOS_EXPECT_VECTOR_EQ(point_1.Coordinates(), retrieved_point.Coordinates());
    KRATOS_EXPECT_DOUBLE_EQ(point_1.GetDistance(), retrieved_point.GetDistance());
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

    KRATOS_EXPECT_EQ(cont.GetPoints().size(), 1);

    cont.Add(point_2); // this must also not add the point again!

    KRATOS_EXPECT_EQ(cont.GetPoints().size(), 1);
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

    KRATOS_EXPECT_EQ(cont.GetPoints().size(), max_size);

    // check if are points are ordered by distance
    int counter = 1;
    for (const auto& r_point : cont.GetPoints()) {
        KRATOS_EXPECT_DOUBLE_EQ(counter*1.0, r_point.GetDistance());
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

    KRATOS_EXPECT_EQ(cont.GetPoints().size(), max_size);
    KRATOS_EXPECT_EQ(cont2.GetPoints().size(), max_size);

    cont.Merge(cont2);
    KRATOS_EXPECT_EQ(cont.GetPoints().size(), max_size);

    // check if are points are ordered by distance
    int counter = 2;
    for (const auto& r_point : cont.GetPoints()) {
        KRATOS_EXPECT_DOUBLE_EQ(counter*0.5, r_point.GetDistance());
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

    KRATOS_EXPECT_EQ(cont.GetPoints().size(), max_size);

    // check if are points are ordered by distance
    int counter = 1;
    for (const auto& r_point : cont.GetPoints()) {
        KRATOS_EXPECT_DOUBLE_EQ(counter*1.0, r_point.GetDistance());
        counter++;
    }

    // serializing the object
    StreamSerializer serializer;
    serializer.save("obj", cont);
    ClosestPointsContainer cont_new(0);
    serializer.load("obj", cont_new);

    KRATOS_EXPECT_EQ(cont.GetPoints().size(), cont_new.GetPoints().size());
    // check if are points are ordered by distance
    counter = 1;
    for (const auto& r_point : cont.GetPoints()) {
        KRATOS_EXPECT_DOUBLE_EQ(counter*1.0, r_point.GetDistance());
        counter++;
    }
}

}  // namespace Kratos::Testing