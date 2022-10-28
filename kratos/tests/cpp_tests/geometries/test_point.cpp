//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "geometries/point.h"
#include "testing/testing.h"

namespace Kratos::Testing {

/**
 *  Here the constructors are test
 */
KRATOS_TEST_CASE_IN_SUITE(Constructors, KratosCoreGeometriesFastSuite)
{
    const Point point;

    KRATOS_CHECK_EQUAL(point.Dimension(), 3);
    KRATOS_CHECK_DOUBLE_EQUAL(point.X(), 0.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point.Y(), 0.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point.Z(), 0.0);

    const Point point_3d(1.0, 2.0, 3.0);

    KRATOS_CHECK_DOUBLE_EQUAL(point_3d.X(), 1.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point_3d.Y(), 2.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point_3d.Z(), 3.0);

    const Point copied_point_3d(point_3d);

    KRATOS_CHECK_DOUBLE_EQUAL(copied_point_3d.X(), point_3d.X());
    KRATOS_CHECK_DOUBLE_EQUAL(copied_point_3d.Y(), point_3d.Y());
    KRATOS_CHECK_DOUBLE_EQUAL(copied_point_3d.Z(), point_3d.Z());

    const array_1d<double, 3> coords {1.0, 2.0, 3.0};
    const Point array_point_3d(coords);

    KRATOS_CHECK_DOUBLE_EQUAL(array_point_3d.X(), coords[0]);
    KRATOS_CHECK_DOUBLE_EQUAL(array_point_3d.Y(), coords[1]);
    KRATOS_CHECK_DOUBLE_EQUAL(array_point_3d.Z(), coords[2]);

    const std::vector<double> vector_coords = {0.0,1.0,2.0};
    const Point vector_point_3d(vector_coords);

    KRATOS_CHECK_DOUBLE_EQUAL(vector_point_3d.X(), vector_coords[0]);
    KRATOS_CHECK_DOUBLE_EQUAL(vector_point_3d.Y(), vector_coords[1]);
    KRATOS_CHECK_DOUBLE_EQUAL(vector_point_3d.Z(), vector_coords[2]);
}

/**
 *  Here the accesors are test
 */
KRATOS_TEST_CASE_IN_SUITE(Access, KratosCoreGeometriesFastSuite)
{
    Point point;
    point.X() = 1.0;
    point.Y() = 2.0;
    point.Z() = 3.0;

    KRATOS_CHECK_DOUBLE_EQUAL(point.X(), 1.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point.Y(), 2.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point.Z(), 3.0);

    array_1d<double, 3>& r_point_ref = point.Coordinates();
    KRATOS_CHECK_DOUBLE_EQUAL(r_point_ref[0], 1.0);
    KRATOS_CHECK_DOUBLE_EQUAL(r_point_ref[1], 2.0);
    KRATOS_CHECK_DOUBLE_EQUAL(r_point_ref[2], 3.0);

    r_point_ref[0] = 3.0;
    r_point_ref[2] = 1.0;
    KRATOS_CHECK_DOUBLE_EQUAL(point.X(), 3.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point.Y(), 2.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point.Z(), 1.0);
}

/**
 *  Here the Info is test
 */
KRATOS_TEST_CASE_IN_SUITE(Info, KratosCoreGeometriesFastSuite)
{
    const Point point;
    KRATOS_CHECK_STRING_EQUAL(point.Info(), "Point");
}

/**
 *  Here the Distance is test
 */
KRATOS_TEST_CASE_IN_SUITE(Distance, KratosCoreGeometriesFastSuite)
{
    const Point point;
    const Point point_3d(3.0, 4.0, 0.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point.Distance(point_3d), 5.0);
    KRATOS_CHECK_DOUBLE_EQUAL(point.SquaredDistance(point_3d), 25.0);
}

}  // namespace Kratos::Testing
