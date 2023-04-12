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

using OTPointType = Point;
using OTPointPointerType = typename OTPointType::Pointer;
using OTPointVectorType = std::vector<OTPointType::Pointer>;
using OTPointIteratorType = std::vector<OTPointType::Pointer>::iterator;
using OTDistanceVectorType = std::vector<double>;
using OTDistanceIteratorType = std::vector<double>::iterator;

/// Octree definitions
using OTBucketType = Bucket< 3ul, OTPointType, OTPointVectorType, OTPointPointerType, OTPointIteratorType, OTDistanceIteratorType>;
using Octree = Tree<OCTreePartition<OTBucketType>>;

/**
 * @brief Test that ExistPoint works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(OctreeExistPoint, KratosCoreFastSuite)
{
    OTPointVectorType points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(OTPointPointerType(new OTPointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    KRATOS_CHECK_EQUAL(testOctree.ExistPoint(OTPointPointerType(new OTPointType(10.0, 10.0, 10.0))), nullptr);
    KRATOS_CHECK_EQUAL(testOctree.ExistPoint(OTPointPointerType(new OTPointType(9.0, 9.0, 9.0))), points[9]);
}

/**
 * @brief Test that SearchNearestPoint works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(OctreeSearchNearestPoint, KratosCoreFastSuite)
{
    OTPointVectorType points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(OTPointPointerType(new OTPointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    auto point_10 = OTPointType(10.0, 10.0, 10.0);
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
    OTPointVectorType points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(OTPointPointerType(new OTPointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    const std::size_t max_number_results = 10;
    OTPointVectorType result_points(max_number_results);
    OTDistanceVectorType distances(max_number_results);
    auto point_10 = OTPointType(10.0, 10.0, 10.0);
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
    OTPointVectorType points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(OTPointPointerType(new OTPointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    const std::size_t max_number_results = 10;
    OTPointVectorType result_points(max_number_results);
    OTDistanceVectorType distances(max_number_results);
    auto point_10 = OTPointType(10.0, 10.0, 10.0);
    auto point_11 = OTPointType(9.1, 9.1, 9.1);
    KRATOS_CHECK_EQUAL(testOctree.SearchInBox(point_11, point_10, result_points.begin(), max_number_results), 0);

    auto point_12 = OTPointType(9.0, 9.0, 9.0);
    KRATOS_CHECK_EQUAL(testOctree.SearchInBox(point_12, point_10, result_points.begin(), max_number_results), 1);
}

/**
 * @brief Test that Bounding box points works well works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(OctreeBB, KratosCoreFastSuite)
{
    OTPointVectorType points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(OTPointPointerType(new OTPointType(i, i, i)));
    }

    Octree testOctree(points.begin(), points.end(), 100);

    KRATOS_CHECK_VECTOR_EQUAL(testOctree.BoundingBoxLowPoint().Coordinates(), points[0]->Coordinates());
    KRATOS_CHECK_VECTOR_EQUAL(testOctree.BoundingBoxHighPoint().Coordinates(), points[9]->Coordinates());
}

} // namespace Kratos::Testing