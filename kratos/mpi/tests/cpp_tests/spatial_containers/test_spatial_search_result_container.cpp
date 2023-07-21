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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerAddResult, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(r_data_comm.Rank() + 1);
    SpatialSearchResult<GeometricalObject> result(&object);
    result.SetDistance(0.5);

    // Add the result to the container
    container.AddResult(result);

    // Check that the result was added correctly
    auto& r_local_pointers = container.GetLocalResults();
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), 1);
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), container.NumberOfLocalResults());

    // Check distances
    auto distances = container.GetDistances();
    KRATOS_CHECK_EQUAL(static_cast<int>(distances.size()), r_data_comm.Size());
    KRATOS_CHECK_EQUAL(distances[r_data_comm.Rank()], 0.5);

    // Check global pointers
    auto& r_global_pointers = container.GetGlobalResults();
    KRATOS_CHECK_EQUAL(r_global_pointers.size(), 0); // It should be empty as we have not synchronized
    KRATOS_CHECK_EQUAL(r_global_pointers.size(), container.NumberOfGlobalResults()); // It should be empty as we have not synchronized
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerClear, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(r_data_comm.Rank() + 1);
    SpatialSearchResult<GeometricalObject> result(&object);
    result.SetDistance(0.5);

    // Add the result to the container
    container.AddResult(result);

    // Clear
    container.Clear();

    // Check that the result was added correctly
    auto& r_local_pointers = container.GetLocalResults();
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), 0);
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), container.NumberOfLocalResults());
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerSynchronizeAll, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(r_data_comm.Rank() + 1);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    container.SynchronizeAll(r_data_comm);

    // Check that the result was added correctly
    auto& r_local_pointers = container.GetLocalResults();
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), 1);
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), container.NumberOfLocalResults());

    // Check global pointers
    auto& r_global_pointers = container.GetGlobalResults();
    KRATOS_CHECK_EQUAL(static_cast<int>(r_global_pointers.size()), r_data_comm.Size());
    KRATOS_CHECK_EQUAL(r_global_pointers.size(), container.NumberOfGlobalResults()); 
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerGetResultShapeFunctions, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Generate a geometry
    auto p_node1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_node2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Geometry<Node>::Pointer p_geom = Kratos::make_shared<Line2D2<Node>>(p_node1, p_node2);

    // Create a test result
    GeometricalObject object = GeometricalObject(r_data_comm.Rank() + 1, p_geom);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    container.SynchronizeAll(r_data_comm);

    // Compute shape functions
    Point point = Point(0.5, 0.0, 0.0);
    auto shape_functions = container.GetResultShapeFunctions(point);

    // Check shape functions
    KRATOS_CHECK_EQUAL(static_cast<int>(shape_functions.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_CHECK_NEAR(shape_functions[i_rank][0], 0.5, 1.0e-12);
        KRATOS_CHECK_NEAR(shape_functions[i_rank][1], 0.5, 1.0e-12);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerGetResultIndices, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(r_data_comm.Rank() + 1);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    container.SynchronizeAll(r_data_comm);

    // Compute shape functions
    auto indixes = container.GetResultIndices();

    // Check shape functions
    KRATOS_CHECK_EQUAL(static_cast<int>(indixes.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_CHECK_EQUAL(static_cast<int>(indixes[i_rank]), i_rank + 1);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerGetResultCoordinates, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    
    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Generate a geometry
    auto p_node1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_node2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Geometry<Node>::Pointer p_geom = Kratos::make_shared<Line2D2<Node>>(p_node1, p_node2);

    // Create a test result
    GeometricalObject object = GeometricalObject(r_data_comm.Rank() + 1, p_geom);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    container.SynchronizeAll(r_data_comm);

    // Compute shape functions
    auto coordinates = container.GetResultCoordinates();

    // Check shape functions
    KRATOS_CHECK_EQUAL(static_cast<int>(coordinates.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_CHECK_EQUAL(coordinates[i_rank].size(), 2);
        KRATOS_CHECK_VECTOR_NEAR(coordinates[i_rank][0], p_node1->Coordinates(), 1.0e-12);
        KRATOS_CHECK_VECTOR_NEAR(coordinates[i_rank][1], p_node2->Coordinates(), 1.0e-12);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerMapInitializeResult, KratosMPICoreFastSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerMapClear, KratosMPICoreFastSuite)
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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerMapOperators, KratosMPICoreFastSuite)
{    
    // Create a test object
    SpatialSearchResultContainerMap<GeometricalObject> container_map;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    container_map.InitializeResult(point);

    // Check that the result was added correctly
    auto& r_result = container_map[point];
    auto& r_local_pointers = r_result.GetLocalResults();
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), 0);
    KRATOS_CHECK_EQUAL(r_local_pointers.size(), r_result.NumberOfLocalResults());
}

}  // namespace Kratos::Testing