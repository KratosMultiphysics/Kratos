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
#include "spatial_containers/spatial_search_result_container_vector.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorInitializeResult, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    const std::size_t index = 0;
    DataCommunicator data_communicator;
    container_vector.InitializeResult(data_communicator);

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_vector.HasResult(index));
    const std::size_t fake_index = 1;
    KRATOS_EXPECT_FALSE(container_vector.HasResult(fake_index));
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorInitializeResults, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0,1,2,3,4,5,6,7,8,9};
    const std::vector<const DataCommunicator*> data_communicators(indexes.size(), &data_communicator);
    container_vector.InitializeResults(data_communicators);

    // Check that the result was added correctly
    for (auto index : indexes) {
        KRATOS_EXPECT_TRUE(container_vector.HasResult(index));
    }
    const std::size_t fake_index = 10;
    KRATOS_EXPECT_FALSE(container_vector.HasResult(fake_index));
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorClear, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    const std::size_t index = 0;
    DataCommunicator data_communicator;
    container_vector.InitializeResult(data_communicator);

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_vector.HasResult(index));
    container_vector.Clear();
    KRATOS_EXPECT_FALSE(container_vector.HasResult(index));
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorOperators, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    const std::size_t index = 0;
    DataCommunicator data_communicator;
    container_vector.InitializeResult(data_communicator);

    // Check that the result was added correctly
    auto& r_result = container_vector[index];
    auto& r_local_pointers = r_result.GetLocalResults();
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 0);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), r_result.NumberOfLocalResults());
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorSynchronizeAll, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0,1};
    const std::vector<const DataCommunicator*> data_communicators(indexes.size(), &data_communicator);
    container_vector.InitializeResults(data_communicators);

    // Container 1
    auto& r_container_1 = container_vector[0];

    // Create a test result
    GeometricalObject object_1 = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result_1(&object_1);
    result_1.SetDistance(0.5);
    GeometricalObject object_2 = GeometricalObject(2);
    SpatialSearchResult<GeometricalObject> result_2(&object_2);
    result_2.SetDistance(0.25);

    // Add the result to the container
    r_container_1.AddResult(result_1);
    r_container_1.AddResult(result_2);

    // Container 2
    auto& r_container_2 = container_vector[1];

    // Create a test result
    GeometricalObject object_3 = GeometricalObject(3);
    SpatialSearchResult<GeometricalObject> result_3(&object_3);
    result_3.SetDistance(0.5);

    // Add the result to the container
    r_container_2.AddResult(result_3);

    // SynchronizeAll
    container_vector.SynchronizeAll(data_communicator);

    // Check that the results were added correctly
    KRATOS_EXPECT_EQ(r_container_1.NumberOfLocalResults(), 2);
    KRATOS_EXPECT_EQ(r_container_1.NumberOfGlobalResults(), 2);
    KRATOS_EXPECT_EQ(r_container_2.NumberOfLocalResults(), 1);
    KRATOS_EXPECT_EQ(r_container_2.NumberOfGlobalResults(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetDistances, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0,1};
    const std::vector<const DataCommunicator*> data_communicators(indexes.size(), &data_communicator);
    container_vector.InitializeResults(data_communicators);

    // Container 1
    auto& r_container_1 = container_vector[0];

    // Create a test result
    GeometricalObject object_1 = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result_1(&object_1);
    result_1.SetDistance(0.5);
    GeometricalObject object_2 = GeometricalObject(2);
    SpatialSearchResult<GeometricalObject> result_2(&object_2);
    result_2.SetDistance(0.25);

    // Add the result to the container
    r_container_1.AddResult(result_1);
    r_container_1.AddResult(result_2);

    // Container 2
    auto& r_container_2 = container_vector[1];

    // Create a test result
    GeometricalObject object_3 = GeometricalObject(3);
    SpatialSearchResult<GeometricalObject> result_3(&object_3);
    result_3.SetDistance(0.5);

    // Add the result to the container
    r_container_2.AddResult(result_3);

    // SynchronizeAll
    container_vector.SynchronizeAll(data_communicator);

    // GetDistances
    auto r_distances = container_vector.GetDistances();
    KRATOS_EXPECT_EQ(r_distances.size(), 2);
    KRATOS_EXPECT_EQ(r_distances[0].size(), 2);
    KRATOS_EXPECT_DOUBLE_EQ(r_distances[0][0], 0.5);
    KRATOS_EXPECT_DOUBLE_EQ(r_distances[0][1], 0.25);
    KRATOS_EXPECT_EQ(r_distances[1].size(), 1);
    KRATOS_EXPECT_DOUBLE_EQ(r_distances[1][0], 0.5);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetResultIsLocal, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0,1};
    const std::vector<const DataCommunicator*> data_communicators(indexes.size(), &data_communicator);
    container_vector.InitializeResults(data_communicators);

    // Container 1
    auto& r_container_1 = container_vector[0];

    // Create a test result
    GeometricalObject object_1 = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result_1(&object_1);
    result_1.SetDistance(0.5);
    GeometricalObject object_2 = GeometricalObject(2);
    SpatialSearchResult<GeometricalObject> result_2(&object_2);
    result_2.SetDistance(0.25);

    // Add the result to the container
    r_container_1.AddResult(result_1);
    r_container_1.AddResult(result_2);

    // Container 2
    auto& r_container_2 = container_vector[1];

    // Create a test result
    GeometricalObject object_3 = GeometricalObject(3);
    SpatialSearchResult<GeometricalObject> result_3(&object_3);
    result_3.SetDistance(0.5);

    // Add the result to the container
    r_container_2.AddResult(result_3);

    // SynchronizeAll
    container_vector.SynchronizeAll(data_communicator);

    // GetResultIsLocal
    auto r_is_local = container_vector.GetResultIsLocal();
    KRATOS_EXPECT_EQ(r_is_local.size(), 2);
    KRATOS_EXPECT_EQ(r_is_local[0].size(), 2);
    KRATOS_EXPECT_TRUE(r_is_local[0][0]);
    KRATOS_EXPECT_TRUE(r_is_local[0][1]);
    KRATOS_EXPECT_EQ(r_is_local[1].size(), 1);
    KRATOS_EXPECT_TRUE(r_is_local[1][0]);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetResultRank, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0,1};
    const std::vector<const DataCommunicator*> data_communicators(indexes.size(), &data_communicator);
    container_vector.InitializeResults(data_communicators);

    // Container 1
    auto& r_container_1 = container_vector[0];

    // Create a test result
    GeometricalObject object_1 = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result_1(&object_1);
    result_1.SetDistance(0.5);
    GeometricalObject object_2 = GeometricalObject(2);
    SpatialSearchResult<GeometricalObject> result_2(&object_2);
    result_2.SetDistance(0.25);

    // Add the result to the container
    r_container_1.AddResult(result_1);
    r_container_1.AddResult(result_2);

    // Container 2
    auto& r_container_2 = container_vector[1];

    // Create a test result
    GeometricalObject object_3 = GeometricalObject(3);
    SpatialSearchResult<GeometricalObject> result_3(&object_3);
    result_3.SetDistance(0.5);

    // Add the result to the container
    r_container_2.AddResult(result_3);

    // SynchronizeAll
    container_vector.SynchronizeAll(data_communicator);

    // GetResultRank
    auto r_rank = container_vector.GetResultRank();
    KRATOS_EXPECT_EQ(r_rank.size(), 2);
    KRATOS_EXPECT_EQ(r_rank[0].size(), 2);
    KRATOS_EXPECT_EQ(r_rank[0][0], 0);
    KRATOS_EXPECT_EQ(r_rank[0][1], 0);
    KRATOS_EXPECT_EQ(r_rank[1].size(), 1);
    KRATOS_EXPECT_EQ(r_rank[1][0], 0);
}

}  // namespace Kratos::Testing