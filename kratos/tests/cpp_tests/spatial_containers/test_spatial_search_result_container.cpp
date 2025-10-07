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
    DataCommunicator data_communicator;
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(1);
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
    KRATOS_EXPECT_FALSE(container.IsObjectFound());
    auto& r_global_pointers = container.GetGlobalResults();
    KRATOS_EXPECT_EQ(r_global_pointers.size(), 0); // It should be empty as we have not synchronized
    KRATOS_EXPECT_EQ(r_global_pointers.size(), container.NumberOfGlobalResults()); // It should be empty as we have not synchronized
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerClear, KratosCoreFastSuite)
{
    // Create a test object
    DataCommunicator data_communicator;
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
    auto& r_local_pointers = container.GetLocalResults();
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 0);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), container.NumberOfLocalResults());
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerSynchronizeAll, KratosCoreFastSuite)
{
    // Create a test object
    DataCommunicator data_communicator;
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    container.AddResult(result);

    // Synchronize the container between partitions
    container.SynchronizeAll(data_communicator);

    // Check that the result was added correctly
    auto& r_local_pointers = container.GetLocalResults();
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 1);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), container.NumberOfLocalResults());

    // Check global pointers
    KRATOS_EXPECT_TRUE(container.IsObjectFound());
    auto& r_global_pointers = container.GetGlobalResults();
    KRATOS_EXPECT_EQ(r_global_pointers.size(), 1);
    KRATOS_EXPECT_EQ(r_global_pointers.size(), container.NumberOfGlobalResults());
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerRemoveResultsFromIndexesList, KratosCoreFastSuite)
{
    // Create a test object
    DataCommunicator data_communicator;
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object_1 = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result_1(&object_1);
    container.AddResult(result_1);
    GeometricalObject object_2 = GeometricalObject(2);
    SpatialSearchResult<GeometricalObject> result_2(&object_2);
    container.AddResult(result_2);
    GeometricalObject object_3 = GeometricalObject(3);
    SpatialSearchResult<GeometricalObject> result_3(&object_3);
    container.AddResult(result_3);

    // Check that the result was added correctly
    KRATOS_EXPECT_EQ(container.NumberOfLocalResults(), 3);

    // Synchronize the container between partitions
    container.SynchronizeAll(data_communicator);

    // Check global pointers
    KRATOS_EXPECT_EQ(container.NumberOfGlobalResults(), 3);

    // Remove indexes
    std::vector<std::size_t> index_to_remove({2,3});
    container.RemoveResultsFromIndexesList(index_to_remove);

    // Check that the result was removed correctly
    KRATOS_EXPECT_EQ(container.NumberOfLocalResults(), 1);
    KRATOS_EXPECT_EQ(container.NumberOfGlobalResults(), 1);
}

}  // namespace Kratos::Testing