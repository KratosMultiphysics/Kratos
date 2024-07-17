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
#include "spatial_containers/spatial_search_result_container_vector.h"

namespace Kratos::Testing
{

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorInitializeResult, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    const std::size_t index = 0;
    container_vector.InitializeResult(r_data_comm);

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_vector.HasResult(index));
    const std::size_t fake_index = 1;
    KRATOS_EXPECT_FALSE(container_vector.HasResult(fake_index));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorInitializeResults, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    const std::vector<std::size_t> indexes{0,1,2,3,4,5,6,7,8,9};
    const std::vector<const DataCommunicator*> data_communicators(indexes.size(), &r_data_comm);
    container_vector.InitializeResults(data_communicators);

    // Check that the result was added correctly
    for (auto index : indexes) {
        KRATOS_EXPECT_TRUE(container_vector.HasResult(index));
    }
    const std::size_t fake_index = 10;
    KRATOS_EXPECT_FALSE(container_vector.HasResult(fake_index));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorClear, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    const std::size_t index = 0;
    container_vector.InitializeResult(r_data_comm);

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_vector.HasResult(index));
    container_vector.Clear();
    KRATOS_EXPECT_FALSE(container_vector.HasResult(index));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorOperators, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    const std::size_t index = 0;
    container_vector.InitializeResult(r_data_comm);

    // Check that the result was added correctly
    auto& r_result = container_vector[index];
    auto& r_local_pointers = r_result.GetLocalResults();
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 0);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), r_result.NumberOfLocalResults());
}

// KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorSynchronizeAllSubDataCommunicators, KratosMPICoreFastSuite)
// {
//     // The data communicator
//     const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

//     // MPI data
//     const int rank = r_data_comm.Rank();
//     const int size = r_data_comm.Size();

//     // Create a test object
//     SpatialSearchResultContainerVector<GeometricalObject, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS> container_vector;

//     // Initialize result
//     std::vector<std::size_t> indexes(rank + 1);
//     std::size_t counter = 0;
//     for (auto& index : indexes) {
//         index = counter;
//         counter++;
//     }
//     std::vector<const DataCommunicator*> data_communicators(indexes.size());
//     for (std::size_t i = 0; i < indexes.size(); i++) {
//         std::vector<int> ranks(size - i);
//         for (std::size_t j = 0; j < ranks.size(); j++) {
//             ranks[j] = size - j - 1;
//         }
//         std::sort(ranks.begin(), ranks.end());
//         const DataCommunicator& r_sub_data_comm = r_data_comm.GetSubDataCommunicator(ranks, std::to_string(i+1));
//         data_communicators[i] = &r_sub_data_comm;
//     }

//     // Initialize results
//     container_vector.InitializeResults(data_communicators);

//     // Create a result
//     GeometricalObject object = GeometricalObject(rank + 1);
//     SpatialSearchResult<GeometricalObject> result(&object);
//     result.SetDistance(0.5*(rank + 1));

//     // Add the result to the containers
//     for (std::size_t i = 0; i < indexes.size(); i++) {
//         // Container
//         auto& r_container = container_vector[i];

//         // Add the result to the container
//         r_container.AddResult(result);
//     }

//     // SynchronizeAll
//     container_vector.SynchronizeAll(r_data_comm);

//     // Check that the results were added correctly
//     KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), static_cast<std::size_t>(rank + 1));
//     for (std::size_t i = 0; i < indexes.size(); i++) {
//         // Container
//         auto& r_container = container_vector[i];
//         KRATOS_EXPECT_EQ(r_container.NumberOfLocalResults(), 1);
//         KRATOS_EXPECT_EQ(r_container.NumberOfGlobalResults(), static_cast<std::size_t>(size - i));
//         const DataCommunicator& r_sub_data_comm = *data_communicators[i];
//         KRATOS_EXPECT_EQ(r_container.NumberOfGlobalResults(), static_cast<std::size_t>(r_sub_data_comm.Size()));
//     }
// }

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorGetDistances, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // MPI data
    const int rank = r_data_comm.Rank();
    const int size = r_data_comm.Size();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    std::vector<const DataCommunicator*> data_communicators(size, &r_data_comm);

    // Initialize results
    container_vector.InitializeResults(data_communicators);

    // Create a result
    GeometricalObject object = GeometricalObject(rank + 1);
    SpatialSearchResult<GeometricalObject> result(&object);
    result.SetDistance(0.5*(rank + 1));

    // Add the result to the containers
    auto& r_container = container_vector[rank];

    // Add the result to the container
    r_container.AddResult(result);

    // SynchronizeAll
    container_vector.SynchronizeAll(r_data_comm);

    // Check that the results were added correctly
    KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), static_cast<std::size_t>(size));

    // GetDistances
    auto r_distances = container_vector.GetDistances();
    KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), r_distances.size());
    for (std::size_t i = 0; i < r_distances.size(); ++i) {
        auto& r_partial_distance = r_distances[i];
        KRATOS_EXPECT_EQ(r_partial_distance.size(), 1);
        KRATOS_EXPECT_DOUBLE_EQ(r_partial_distance[0], 0.5*(i + 1));
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorGetResultIsLocal, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // MPI data
    const int rank = r_data_comm.Rank();
    const int size = r_data_comm.Size();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    std::vector<const DataCommunicator*> data_communicators(size, &r_data_comm);

    // Initialize results
    container_vector.InitializeResults(data_communicators);

    // Create a result
    GeometricalObject object = GeometricalObject(rank + 1);
    SpatialSearchResult<GeometricalObject> result(&object, rank);
    result.SetDistance(0.5*(rank + 1));

    // Add the result to the containers
    auto& r_container = container_vector[rank];

    // Add the result to the container
    r_container.AddResult(result);

    // SynchronizeAll
    container_vector.SynchronizeAll(r_data_comm);

    // Check that the results were added correctly
    KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), static_cast<std::size_t>(size));

    // GetResultIsLocal
    auto r_result_is_local = container_vector.GetResultIsLocal();
    KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), r_result_is_local.size());
    for (int i = 0; i < static_cast<int>(r_result_is_local.size()); ++i) {
        auto& r_partial_result_is_local = r_result_is_local[i];
        KRATOS_EXPECT_EQ(r_partial_result_is_local.size(), 1);
        KRATOS_EXPECT_EQ(r_partial_result_is_local[0], i == rank);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorGetResultRank, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // MPI data
    const int rank = r_data_comm.Rank();
    const int size = r_data_comm.Size();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    std::vector<const DataCommunicator*> data_communicators(size, &r_data_comm);

    // Initialize results
    container_vector.InitializeResults(data_communicators);

    // Create a result
    GeometricalObject object = GeometricalObject(rank + 1);
    SpatialSearchResult<GeometricalObject> result(&object, rank);
    result.SetDistance(0.5*(rank + 1));

    // Add the result to the containers
    auto& r_container = container_vector[rank];

    // Add the result to the container
    r_container.AddResult(result);

    // SynchronizeAll
    container_vector.SynchronizeAll(r_data_comm);

    // Check that the results were added correctly
    KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), static_cast<std::size_t>(size));

    // GetResultRank
    auto r_result_rank = container_vector.GetResultRank();
    KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), r_result_rank.size());
    for (int i = 0; i < static_cast<int>(r_result_rank.size()); ++i) {
        auto& r_partial_result_rank = r_result_rank[i];
        KRATOS_EXPECT_EQ(r_partial_result_rank.size(), 1);
        KRATOS_EXPECT_EQ(r_partial_result_rank[0], i);
    }
}

}  // namespace Kratos::Testing