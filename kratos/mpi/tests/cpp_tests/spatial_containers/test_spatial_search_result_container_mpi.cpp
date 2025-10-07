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
#include "mpi/testing/mpi_testing.h"
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
}

}  // namespace Kratos::Testing