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
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 1);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), container.NumberOfLocalResults());

    // Check distances
    KRATOS_EXPECT_EQ(container[0].GetDistance(), 0.5);

    // Check global pointers
    auto& r_global_pointers = container.GetGlobalResults();
    KRATOS_EXPECT_EQ(r_global_pointers.size(), 0); // It should be empty as we have not synchronized
    KRATOS_EXPECT_EQ(r_global_pointers.size(), container.NumberOfGlobalResults()); // It should be empty as we have not synchronized
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
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 0);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), container.NumberOfLocalResults());
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
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 1);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), container.NumberOfLocalResults());

    // Check global pointers
    auto& r_global_pointers = container.GetGlobalResults();
    KRATOS_EXPECT_EQ(static_cast<int>(r_global_pointers.size()), r_data_comm.Size());
    KRATOS_EXPECT_EQ(r_global_pointers.size(), container.NumberOfGlobalResults());
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
    KRATOS_EXPECT_EQ(static_cast<int>(shape_functions.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_NEAR(shape_functions[i_rank][0], 0.5, 1.0e-12);
        KRATOS_EXPECT_NEAR(shape_functions[i_rank][1], 0.5, 1.0e-12);
    }

    // Check is inside
    auto is_inside_true = container.GetResultIsInside(point, 1.0e-5);
    KRATOS_EXPECT_EQ(static_cast<int>(is_inside_true.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_TRUE(is_inside_true[i_rank]);
    }

    Point point_outside = Point(1.0e6, 1.0e6, 1.0e6);
    auto is_inside_false = container.GetResultIsInside(point_outside, 1.0e-5);
    KRATOS_EXPECT_EQ(static_cast<int>(is_inside_true.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_FALSE(is_inside_false[i_rank]);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerGetResultIsLocal, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(r_data_comm.Rank() + 1);
    SpatialSearchResult<GeometricalObject> result(&object, r_data_comm.Rank());

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    container.SynchronizeAll(r_data_comm);

    // Compute is local
    auto is_local = container.GetResultIsLocal();

    // Check is local
    KRATOS_EXPECT_EQ(static_cast<int>(is_local.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        if (i_rank == r_data_comm.Rank()) {
            KRATOS_EXPECT_TRUE(static_cast<int>(is_local[i_rank]));
        } else {
            KRATOS_EXPECT_FALSE(static_cast<int>(is_local[i_rank]));
        }
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerGetResultIsActive, KratosMPICoreFastSuite)
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

    // Compute is active
    auto is_active = container.GetResultIsActive();

    // Check is active
    KRATOS_EXPECT_EQ(static_cast<int>(is_active.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_TRUE(static_cast<int>(is_active[i_rank]));
    }

    // Deactivate the object
    object.Set(ACTIVE, false);

    // Compute is active
    is_active = container.GetResultIsActive();

    // Check is active
    KRATOS_EXPECT_EQ(static_cast<int>(is_active.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_FALSE(static_cast<int>(is_active[i_rank]));
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

    // Compute indices
    auto indices = container.GetResultIndices();

    // Check indices
    KRATOS_EXPECT_EQ(static_cast<int>(indices.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_EQ(static_cast<int>(indices[i_rank]), i_rank + 1);
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
    KRATOS_EXPECT_EQ(static_cast<int>(coordinates.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_EQ(coordinates[i_rank].size(), 2);
        KRATOS_EXPECT_VECTOR_NEAR(coordinates[i_rank][0], p_node1->Coordinates(), 1.0e-12);
        KRATOS_EXPECT_VECTOR_NEAR(coordinates[i_rank][1], p_node2->Coordinates(), 1.0e-12);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerRemoveResultsFromIndexesList, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    const int rank = r_data_comm.Rank();
    const int world_size = r_data_comm.Size();

    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object_1 = GeometricalObject(3 * rank + 1);
    SpatialSearchResult<GeometricalObject> result_1(&object_1);
    container.AddResult(result_1);
    GeometricalObject object_2 = GeometricalObject(3 * rank + 2);
    SpatialSearchResult<GeometricalObject> result_2(&object_2);
    container.AddResult(result_2);
    GeometricalObject object_3 = GeometricalObject(3 * rank + 3);
    SpatialSearchResult<GeometricalObject> result_3(&object_3);
    container.AddResult(result_3);

    // Check that the result was added correctly
    KRATOS_EXPECT_EQ(container.NumberOfLocalResults(), 3);

    // Synchronize the container between partitions
    container.SynchronizeAll(r_data_comm);

    // Check global pointers
    KRATOS_EXPECT_EQ(container.NumberOfGlobalResults(), static_cast<std::size_t>(3 * world_size));

    // Remove indexes
    std::vector<std::size_t> index_to_remove;
    index_to_remove.reserve(2 * world_size);
    for (int i_rank = 0; i_rank < world_size; ++i_rank) {
        index_to_remove.push_back(3 * i_rank + 2);
        index_to_remove.push_back(3 * i_rank + 3);
    }
    container.RemoveResultsFromIndexesList(index_to_remove);

    // Check that the result was removed correctly
    KRATOS_EXPECT_EQ(container.NumberOfLocalResults(), 1);
    KRATOS_EXPECT_EQ(container.NumberOfGlobalResults(), static_cast<std::size_t>(world_size));

    // Compute indices
    auto indices = container.GetResultIndices();

    // Check indices
    KRATOS_EXPECT_EQ(indices.size(), static_cast<std::size_t>(world_size));
    for (int i_rank = 0; i_rank < world_size; ++i_rank) {
        KRATOS_EXPECT_EQ(indices[i_rank], static_cast<std::size_t>(3 * i_rank + 1));
    }
}

}  // namespace Kratos::Testing