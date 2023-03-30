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
#include "spatial_containers/kd_tree.h"

namespace Kratos::Testing {

using PointType = Node<3>;
using PointTypePointer = Node<3>::Pointer;
using PointVector = std::vector<PointType::Pointer>;
using PointIterator = std::vector<PointType::Pointer>::iterator;
using DistanceVector = std::vector<double>;
using DistanceIterator = std::vector<double>::iterator;

/// KDtree definitions
using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
using KDTree = Tree<KDTreePartition<BucketType>>;

using CoordinateType = KDTree::CoordinateType;
using PointerType = KDTree::PointerType;
using SearchStructureType = KDTree::SearchStructureType;

/**
 * @brief Test that ExistPoint works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(KDTreeExistPoint, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    KDTree testKDTree(points.begin(), points.end(), 100);

    KRATOS_CHECK_EQUAL(testKDTree.ExistPoint(PointTypePointer(new PointType(10, 10.0, 10.0, 10.0))), nullptr);
    KRATOS_CHECK_EQUAL(testKDTree.ExistPoint(PointTypePointer(new PointType(9, 9.0, 9.0, 9.0))), points[9]);
}

/**
 * @brief Test that SearchNearestPoint works correctly
 */
KRATOS_TEST_CASE_IN_SUITE(KDTreeSearchNearestPoint, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    KDTree testKDTree(points.begin(), points.end(), 100);

    auto point_10 = PointType(10, 10.0, 10.0, 10.0);
    KRATOS_CHECK_EQUAL(testKDTree.SearchNearestPoint(point_10), points[9]);
    double distance;
    KRATOS_CHECK_EQUAL(testKDTree.SearchNearestPoint(point_10, distance), points[9]);
    KRATOS_CHECK_DOUBLE_EQUAL(distance, 3.0); // NOTE: Should be sqrt of 3, may require to check that
}


} // namespace Kratos::Testing