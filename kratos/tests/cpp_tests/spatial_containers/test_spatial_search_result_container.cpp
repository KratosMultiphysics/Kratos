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
#include "testing/testing.h"
#include "geometries/line_2d_2.h"
#include "spatial_containers/spatial_search_result.h"
#include "spatial_containers/spatial_search_result_container.h"

namespace Kratos::Testing 
{

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerAddResult, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result(&object);
    result.SetDistance(0.5);

    // Add the result to the container
    container.AddResult(result);

    // Check that the result was added correctly
    auto& r_local_pointers = container.GetLocalPointers();
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), 1);

    // Check distances
    auto& r_distances = container.GetLocalDistances();
    KRATOS_CHECK_EQUAL(r_distances.size(), 1);
    KRATOS_CHECK_EQUAL(r_distances[1], 0.5);

    // Check global pointers
    auto& r_global_pointers = container.GetGlobalPointers();
    KRATOS_CHECK_EQUAL(r_global_pointers.size(), 0); // It should be empty as we have not synchronized
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerClear, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result(&object);
    result.SetDistance(0.5);

    // Add the result to the container
    container.AddResult(result);

    // Clear
    container.Clear();

    // Check that the result was added correctly
    auto& r_local_pointers = container.GetLocalPointers();
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerSynchronizeAll, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    DataCommunicator data_communicator;
    container.SynchronizeAll(data_communicator);

    // Check that the result was added correctly
    auto& r_local_pointers = container.GetLocalPointers();
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), 1);

    // Check global pointers
    auto& r_global_pointers = container.GetGlobalPointers();
    KRATOS_CHECK_EQUAL(r_global_pointers.size(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerGetResultShapeFunctions, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Generate a geometry
    auto p_node1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_node2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Geometry<Node>::Pointer p_geom = Kratos::make_shared<Line2D2<Node>>(p_node1, p_node2);

    // Create a test result
    GeometricalObject object = GeometricalObject(1, p_geom);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    DataCommunicator data_communicator;
    container.SynchronizeAll(data_communicator);

    // Compute shape functions
    Point point = Point(0.5, 0.0, 0.0);
    auto shape_functions = container.GetResultShapeFunctions(point);

    // Check shape functions
    KRATOS_CHECK_EQUAL(shape_functions.size(), 1);
    KRATOS_CHECK_NEAR(shape_functions[0][0], 0.5, 1.0e-12);
    KRATOS_CHECK_NEAR(shape_functions[0][1], 0.5, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerGetResultIndices, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    DataCommunicator data_communicator;
    container.SynchronizeAll(data_communicator);

    // Compute shape functions
    auto indixes = container.GetResultIndices();

    // Check shape functions
    KRATOS_CHECK_EQUAL(indixes.size(), 1);
    KRATOS_CHECK_EQUAL(indixes[0], object.Id());
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerGetResultCoordinates, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Generate a geometry
    auto p_node1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_node2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Geometry<Node>::Pointer p_geom = Kratos::make_shared<Line2D2<Node>>(p_node1, p_node2);

    // Create a test result
    GeometricalObject object = GeometricalObject(1, p_geom);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    DataCommunicator data_communicator;
    container.SynchronizeAll(data_communicator);

    // Compute shape functions
    auto coordinates = container.GetResultCoordinates();

    // Check shape functions
    KRATOS_CHECK_EQUAL(coordinates.size(), 1);
    KRATOS_CHECK_EQUAL(coordinates[0].size(), 2);
    KRATOS_CHECK_VECTOR_NEAR(coordinates[0][2], p_node1->Coordinates(), 1.0e-12);
    KRATOS_CHECK_VECTOR_NEAR(coordinates[0][1], p_node2->Coordinates(), 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerMapInitializeResult, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerMap<GeometricalObject> container_map;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    container_map.InitializeResult(point);

    // Check that the result was added correctly
    KRATOS_CHECK(container_map.HasResult(point));
    Point fake_point = Point(1.5, 0.0, 0.0);
    KRATOS_CHECK_IS_FALSE(container_map.HasResult(fake_point));
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerMapClear, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerMap<GeometricalObject> container_map;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    container_map.InitializeResult(point);

    // Check that the result was added correctly
    KRATOS_CHECK(container_map.HasResult(point));
    container_map.Clear();
    KRATOS_CHECK_IS_FALSE(container_map.HasResult(point));
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerMapOperators, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerMap<GeometricalObject> container_map;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    container_map.InitializeResult(point);

    // Check that the result was added correctly
    auto& r_result = container_map[point];
    auto& r_local_pointers = r_result.GetLocalPointers();
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), 0);
}

}  // namespace Kratos::Testing