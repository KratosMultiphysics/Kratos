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
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    const std::size_t index = 0;
    container_vector.InitializeResult();

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_vector.HasResult(index));
    const std::size_t fake_index = 1;
    KRATOS_EXPECT_FALSE(container_vector.HasResult(fake_index));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorInitializeResults, KratosMPICoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    const std::vector<std::size_t> indexes{0,1,2,3,4,5,6,7,8,9};
    container_vector.InitializeResults(indexes.size());

    // Check that the result was added correctly
    for (auto index : indexes) {
        KRATOS_EXPECT_TRUE(container_vector.HasResult(index));
    }
    const std::size_t fake_index = 10;
    KRATOS_EXPECT_FALSE(container_vector.HasResult(fake_index));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorClear, KratosMPICoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    Point point = Point(0.5, 0.0, 0.0);
    const std::size_t index = 0;
    container_vector.InitializeResult();

    // Check that the result was added correctly
    KRATOS_EXPECT_TRUE(container_vector.HasResult(index));
    container_vector.Clear();
    KRATOS_EXPECT_FALSE(container_vector.HasResult(index));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorOperators, KratosMPICoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    const std::size_t index = 0;
    container_vector.InitializeResult();

    // Check that the result was added correctly
    auto& r_result = container_vector[index];
    auto& r_local_pointers = r_result.GetLocalResults();
    KRATOS_EXPECT_EQ(r_local_pointers.size(), 0);
    KRATOS_EXPECT_EQ(r_local_pointers.size(), r_result.NumberOfLocalResults());
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorGetDistances, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // MPI data
    const int rank = r_data_comm.Rank();
    const int size = r_data_comm.Size();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize results
    container_vector.InitializeResults(size);

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

    // Initialize results
    container_vector.InitializeResults(size);

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
    auto r_result_is_local = container_vector.GetResultIsLocal(r_data_comm);
    KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), r_result_is_local.size());
    for (int i = 0; i < static_cast<int>(r_result_is_local.size()); ++i) {
        auto& r_partial_result_is_local = r_result_is_local[i];
        KRATOS_EXPECT_EQ(r_partial_result_is_local.size(), 1);
        KRATOS_EXPECT_EQ(r_partial_result_is_local[0], i == rank);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorGetResultIsActive, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // MPI data
    const int rank = r_data_comm.Rank();
    const int size = r_data_comm.Size();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize results
    container_vector.InitializeResults(size);

    // Add the result to the containers
    auto& r_container = container_vector[rank];

    // Create a result
    GeometricalObject object = GeometricalObject(rank + 1);
    SpatialSearchResult<GeometricalObject> result(&object, rank);

    // Add the result to the container
    r_container.AddResult(result);

    // Synchronize the container between partitions
    container_vector.SynchronizeAll(r_data_comm);

    // Check that the results were added correctly
    KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), static_cast<std::size_t>(size));

    // Compute is active
    auto is_active = container_vector.GetResultIsActive(r_data_comm);

    // Check is active
    KRATOS_EXPECT_EQ(static_cast<int>(is_active.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_TRUE(static_cast<int>(is_active[i_rank][0]));
    }

    // Deactivate the object
    object.Set(ACTIVE, false);

    // Compute is active
    is_active = container_vector.GetResultIsActive(r_data_comm);

    // Check is active
    KRATOS_EXPECT_EQ(static_cast<int>(is_active.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_FALSE(static_cast<int>(is_active[i_rank][0]));
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorGetResultIndices, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // MPI data
    const int rank = r_data_comm.Rank();
    const int size = r_data_comm.Size();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize results
    container_vector.InitializeResults(size);

    // Add the result to the containers
    auto& r_container = container_vector[rank];

    // Create a test object
    SpatialSearchResultContainer<GeometricalObject> container;

    // Create a test result
    GeometricalObject object = GeometricalObject(r_data_comm.Rank() + 1);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    r_container.AddResult(result);

    // Synchronize the container between partitions
    container_vector.SynchronizeAll(r_data_comm);

    // Compute indices
    auto indices = container_vector.GetResultIndices();

    // Check indices
    KRATOS_EXPECT_EQ(static_cast<int>(indices.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_EQ(static_cast<int>(indices[i_rank][0]), i_rank + 1);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISpatialSearchResultContainerVectorGetResultCoordinates, KratosMPICoreFastSuite)
{
    // The data communicator
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // MPI data
    const int rank = r_data_comm.Rank();
    const int size = r_data_comm.Size();

    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize results
    container_vector.InitializeResults(size);

    // Add the result to the containers
    auto& r_container = container_vector[rank];

    // Generate a geometry
    auto p_node1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_node2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Geometry<Node>::Pointer p_geom = Kratos::make_shared<Line2D2<Node>>(p_node1, p_node2);

    // Create a test result
    GeometricalObject object = GeometricalObject(r_data_comm.Rank() + 1, p_geom);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    r_container.AddResult(result);

    // Synchronize the container between partitions
    container_vector.SynchronizeAll(r_data_comm);

    // Compute shape functions
    auto coordinates = container_vector.GetResultCoordinates();

    // Check shape functions
    KRATOS_EXPECT_EQ(static_cast<int>(coordinates.size()), r_data_comm.Size());
    for (int i_rank = 0; i_rank < r_data_comm.Size(); ++i_rank) {
        KRATOS_EXPECT_EQ(coordinates[i_rank][0].size(), 2);
        KRATOS_EXPECT_VECTOR_NEAR(coordinates[i_rank][0][0], p_node1->Coordinates(), 1.0e-12);
        KRATOS_EXPECT_VECTOR_NEAR(coordinates[i_rank][0][1], p_node2->Coordinates(), 1.0e-12);
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

    // Initialize results
    container_vector.InitializeResults(size);

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
    auto r_result_rank = container_vector.GetResultRank(r_data_comm);
    KRATOS_EXPECT_EQ(container_vector.NumberOfSearchResults(), r_result_rank.size());
    for (int i = 0; i < static_cast<int>(r_result_rank.size()); ++i) {
        auto& r_partial_result_rank = r_result_rank[i];
        KRATOS_EXPECT_EQ(r_partial_result_rank.size(), 1);
        KRATOS_EXPECT_EQ(r_partial_result_rank[0], i);
    }
}

}  // namespace Kratos::Testing