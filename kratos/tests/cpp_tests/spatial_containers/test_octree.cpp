//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "testing/testing.h"
#include "spatial_containers/bucket.h"
#include "spatial_containers/octree.h"

namespace Kratos::Testing {

using PointType = Point;
using PointTypePointer = Point::Pointer;
using PointVector = std::vector<PointType::Pointer>;
using PointIterator = std::vector<PointType::Pointer>::iterator;
using DistanceVector = std::vector<double>;
using DistanceIterator = std::vector<double>::iterator;

/// KDtree definitions
using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
using Octree = Tree<OCTreePartition<BucketType>>;

using CoordinateType = Octree::CoordinateType;
using PointerType = Octree::PointerType;
using SearchStructureType = Octree::SearchStructureType;

/**
 * @brief Test that ExistPoint works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(OctreeExistPoint, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    KRATOS_CHECK_EQUAL(testOctree.ExistPoint(PointTypePointer(new PointType(10.0, 10.0, 10.0))), nullptr);
    KRATOS_CHECK_EQUAL(testOctree.ExistPoint(PointTypePointer(new PointType(9.0, 9.0, 9.0))), points[9]);
}

/**
 * @brief Test that SearchNearestPoint works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(OctreeSearchNearestPoint, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    auto point_10 = PointType(10.0, 10.0, 10.0);
    KRATOS_CHECK_EQUAL(testOctree.SearchNearestPoint(point_10), points[9]);
    double distance;
    KRATOS_CHECK_EQUAL(testOctree.SearchNearestPoint(point_10, distance), points[9]);
    KRATOS_CHECK_DOUBLE_EQUAL(distance, 3.0); // NOTE: Should be sqrt of 3, may require to check that
}

/**
 * @brief Test that SearchInRadius works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(OctreeSearchInRadius, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    const std::size_t max_number_results = 10;
    PointVector result_points(max_number_results);
    DistanceVector distances(max_number_results);
    auto point_10 = PointType(10.0, 10.0, 10.0);
    KRATOS_CHECK_EQUAL(testOctree.SearchInRadius(point_10, 1.0, result_points.begin(), distances.begin(), max_number_results), 0);

    KRATOS_CHECK_EQUAL(testOctree.SearchInRadius(point_10, 3.0, result_points.begin(), distances.begin(), max_number_results), 1);
    KRATOS_CHECK_DOUBLE_EQUAL(distances[0], 3.0); // NOTE: Should be sqrt of 3, it is always the quare for performance reasons

    KRATOS_CHECK_EQUAL(testOctree.SearchInRadius(point_10, 4.0, result_points.begin(), max_number_results), 2);
    KRATOS_CHECK_EQUAL(testOctree.SearchInRadius(point_10, 4.0, result_points.begin(), distances.begin(), max_number_results), 2);
    KRATOS_CHECK_DOUBLE_EQUAL(distances[0] + distances[1], 15.0); // NOTE: Should be sqrt of 3 + sqrt of 12, it is always the quare for performance reasons
}

/**
 * @brief Test that SearchInBox works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(OctreeSearchInBox, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    const std::size_t max_number_results = 10;
    PointVector result_points(max_number_results);
    DistanceVector distances(max_number_results);
    auto point_10 = PointType(10.0, 10.0, 10.0);
    auto point_11 = PointType(9.1, 9.1, 9.1);
    KRATOS_CHECK_EQUAL(testOctree.SearchInBox(point_11, point_10, result_points.begin(), max_number_results), 0);

    auto point_12 = PointType(9.0, 9.0, 9.0);
    KRATOS_CHECK_EQUAL(testOctree.SearchInBox(point_12, point_10, result_points.begin(), max_number_results), 1);
}

/**
 * @brief Test that Bounding box points works well works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(OctreeBB, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    KRATOS_CHECK_VECTOR_EQUAL(testOctree.BoundingBoxLowPoint().Coordinates(), points[0]->Coordinates());
    KRATOS_CHECK_VECTOR_EQUAL(testOctree.BoundingBoxHighPoint().Coordinates(), points[9]->Coordinates());
}

} // namespace Kratos::Testing