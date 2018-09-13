//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//
//

// System includes


// External includes


// Project includes
#include "includes/node.h"
#include "testing/testing.h"
#include "spatial_containers/bins_static.h"


namespace Kratos {
namespace Testing {

typedef Node<3>                                     PointType;
typedef Node<3>::Pointer                            PointTypePointer;
typedef std::vector<PointType::Pointer>             PointVector;
typedef std::vector<PointType::Pointer>::iterator   PointIterator;
typedef std::vector<double>                         DistanceVector;
typedef std::vector<double>::iterator               DistanceIterator;

typedef Bins<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator> StaticBins;

typedef StaticBins::CoordinateType                  CoordinateType;
typedef StaticBins::PointerType                     PointerType;
typedef StaticBins::SearchStructureType             SearchStructureType;

/**
 * @brief Test that the bins is constructed correctly
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsDefaultConstructorBoundigBox, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    StaticBins testBins(points.begin(), points.end());

    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[0], 0.0);
    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[1], 0.0);
    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[2], 0.0);

    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[0], 9.0);
    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[1], 9.0);
    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[2], 9.0);
}

/**
 * @brief Test that the number of cells is calculated correctly
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsDefaultConstructorCellNumber, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    StaticBins testBins(points.begin(), points.end());

    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[0], 3);
    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[1], 3);
    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[2], 3);
}

/**
 * @brief Test that the size of the cells is calculated correctly
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsDefaultConstructorCellSize, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    StaticBins testBins(points.begin(), points.end());

    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[0], 3.0);
    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[1], 3.0);
    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[2], 3.0);
}

/**
 * @brief Test that the bins is constructed correctly with the bounding box constructor
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsBBConstructorBoundingBox, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    PointType testMinPoint(10, -10.0, -10.0, -10.0);
    PointType testMaxPoint(11,  10.0,  10.0,  10.0);

    StaticBins testBins(points.begin(), points.end(), testMinPoint, testMaxPoint);

    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[0], testMinPoint[0]);
    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[1], testMinPoint[1]);
    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[2], testMinPoint[2]);

    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[0], testMaxPoint[0]);
    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[1], testMaxPoint[1]);
    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[2], testMaxPoint[2]);
}

/**
 * @brief Test that the number of cells is calculated correctly with the bounding box constructor
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsBBConstructorCellNumber, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    PointType testMinPoint(10, -10.0, -10.0, -10.0);
    PointType testMaxPoint(11,  10.0,  10.0,  10.0);

    StaticBins testBins(points.begin(), points.end(), testMinPoint, testMaxPoint);

    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[0], 3);
    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[1], 3);
    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[2], 3);
}

/**
 * @brief Test that the size of the cells is calculated correctly with the bounding box constructor
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsBBConstructorCellSize, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    PointType testMinPoint(10, -10.0, -10.0, -10.0);
    PointType testMaxPoint(11,  10.0,  10.0,  10.0);

    StaticBins testBins(points.begin(), points.end(), testMinPoint, testMaxPoint);

    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[0], 20.0/3.0);
    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[1], 20.0/3.0);
    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[2], 20.0/3.0);
}

/**
 * @brief Test that the bins is constructed correctly with the cell size constructor
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsCellSizeConstructorBoundingBox, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    CoordinateType testCellSize = 3.2244;

    StaticBins testBins(points.begin(), points.end(), testCellSize);

    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[0], 0.0);
    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[1], 0.0);
    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[2], 0.0);

    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[0], 9.0);
    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[1], 9.0);
    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[2], 9.0);
}

/**
 * @brief Test that the number of cells is calculated correctly with the cell size constructor
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsCellSizeConstructorCellNumber, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    CoordinateType testCellSize = 3.2244;

    StaticBins testBins(points.begin(), points.end(), testCellSize);

    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[0], 3);
    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[1], 3);
    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[2], 3);
}

/**
 * @brief Test that the size of the cells is calculated correctly with the cell size constructor
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsCellSizeConstructorCellSize, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    CoordinateType testCellSize = 3.2244;

    StaticBins testBins(points.begin(), points.end(), testCellSize);

    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[0], testCellSize);
    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[1], testCellSize);
    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[2], testCellSize);
}

/**
 * @brief Searches the nearest point
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsExistPoint, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    StaticBins testBins(points.begin(), points.end());

    PointerType nearestPoint = testBins.ExistPoint(PointerType(new PointType(0, 4.1, 4.1, 4.1)));

    KRATOS_CHECK_EQUAL(nearestPoint->Id(), 4);
}

/**
 * @brief Searches the nearest point (excluding the input)
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsNearestPointInner, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    PointerType pointToSearch = PointerType(new PointType(10, 4.25, 4.25, 4.25));
    points.push_back(pointToSearch);
    
    StaticBins testBins(points.begin(), points.end());

    PointerType nearestPoint = testBins.SearchNearestPointInner(pointToSearch);

    KRATOS_CHECK_EQUAL(nearestPoint->Id(), 4);
}

/**
 * @brief Searches the nearest point with (including the input)
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsNearestPoint, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    PointerType pointToSearch = PointerType(new PointType(10, 4.25, 4.25, 4.25));
    points.push_back(pointToSearch);
    
    StaticBins testBins(points.begin(), points.end());

    PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch);
 
    KRATOS_CHECK_EQUAL(nearestPoint->Id(), 10);
}

/**
 * @brief Searches the nearest point (including the input) with distance
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsNearestPointWithDistance, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    PointerType pointToSearch = PointerType(new PointType(10, 4.25, 4.25, 4.25));
    
    StaticBins testBins(points.begin(), points.end());

    double squaredDistance = 0.0;
    PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch, squaredDistance);
 
    KRATOS_CHECK_EQUAL(nearestPoint->Id(), 4);
    KRATOS_CHECK_EQUAL(squaredDistance, 0.1875);
}

/**
 * @brief Searches the nearest point (including the input) with distance (threadsafe)
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(StaticBinsNearestPointWithDistanceThreadsafe, KratosCoreFastSuite)
{
    PointVector points;

    for(std::size_t i = 0; i < 10; i++) {
        points.push_back(PointTypePointer(new PointType(i, i, i, i)));
    }

    PointerType pointToSearch = PointerType(new PointType(10, 4.25, 4.25, 4.25));
    
    StaticBins testBins(points.begin(), points.end());
    SearchStructureType searchBox;

    double squaredDistance = 0.0;
    PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch, squaredDistance, searchBox);
 
    KRATOS_CHECK_EQUAL(nearestPoint->Id(), 4);
    KRATOS_CHECK_EQUAL(squaredDistance, 0.1875);
}

    
} // namespace Testesing
} // namespace Kratos