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

// System includes
#include <limits>

// External includes

// Project includes
#include "mpi/testing/mpi_testing.h"
#include "containers/model.h"
#include "spatial_containers/parallel_spatial_search.h"
#include "tests/test_utilities/cpp_tests_utilities.h"
#include "mpi/utilities/parallel_fill_communicator.h"

namespace Kratos
{

namespace
{

ModelPart& CreateCubeSkinModelPart(
    Model& rCurrentModel,
    const double HalfX = 0.6,
    const double HalfY = 0.9,
    const double HalfZ = 0.3
    )
{
    // Generate the cube skin
    const DataCommunicator& r_data_communicator = Testing::GetDefaultDataCommunicator();
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(rCurrentModel, HalfX, HalfY, HalfZ, r_data_communicator);

    // Compute communication plan and fill communicator meshes correctly
    ParallelFillCommunicator(r_skin_part, r_data_communicator).Execute();

    // Return the skin model part
    return r_skin_part;
}

ModelPart& CreateCubeModelPart(Model& rCurrentModel)
{
    // Generate the cube
    const DataCommunicator& r_data_communicator = Testing::GetDefaultDataCommunicator();
    ModelPart& r_model_part = CppTestsUtilities::CreateCubeModelPart(rCurrentModel, r_data_communicator);

    // Compute communication plan and fill communicator meshes correctly
    ParallelFillCommunicator(r_model_part, r_data_communicator).Execute();

    // Return the model part
    return r_model_part;
}

// Definition of the geometrical object bins search wrapper
using ParallelSpatialSearchGeometricalObjectsBins = ParallelSpatialSearch<GeometricalObjectsBins>;

/**
 * @brief Test function for searching geometrical objects in bins using a given search wrapper
 * @tparam TParallelSpatialSearch The type of search wrapper to be tested
 */
template<class TParallelSpatialSearch>
void TestParallelSpatialSearchGeometricalObjectsBinsSearchInRadius()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TParallelSpatialSearch parallel_spatial_search_bins(r_skin_part.Elements(), r_data_comm);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    const std::size_t point_id = 1;
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(point_id, 0.0, 0.0, 0.0);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;

    // 0.29 radius
    parallel_spatial_search_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.29, results);
    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }

    // 0.3 radius
    parallel_spatial_search_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }

    // 0.4 radius
    parallel_spatial_search_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.4, results);
    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }

    // 0.6 radius
    parallel_spatial_search_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);
    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }

    // 0.7 radius
    parallel_spatial_search_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.7, results);
    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }

    // 0.9 radius
    parallel_spatial_search_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.9, results);
    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

/**
 * @brief Test function for nearest searching geometrical objects in bins using a given search wrapper
 * @tparam TParallelSpatialSearch The type of search wrapper to be tested
 */
template<class TParallelSpatialSearch>
void TestParallelSpatialSearchGeometricalObjectsBinsSearchNearestInRadius()
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TParallelSpatialSearch parallel_spatial_search_bins(r_skin_part.Elements(), r_data_comm);

    double epsilon = 1.0e-6;
    const std::size_t near_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(near_point_id, epsilon,epsilon,epsilon);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search_bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z - 1.e-4, results);

    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }

    parallel_spatial_search_bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z + 1.e-4, results);

    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

        // Distances are just local
        std::vector<std::vector<double>> distances;
        results.GetDistances(distances);
        KRATOS_EXPECT_NEAR(distances[0][0], (cube_z - epsilon), tolerance);

        // Compute indices
        std::vector<std::vector<IndexType>> indices;
        results.GetResultIndices(indices);
        const std::size_t id = indices[0][0];
        KRATOS_EXPECT_EQ(id, 3);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

/**
 * @brief Test function for nearest searching geometrical objects in bins using a given search wrapper
 * @tparam TParallelSpatialSearch The type of search wrapper to be tested
 */
template<class TParallelSpatialSearch>
void TestParallelSpatialSearchGeometricalObjectsBinsSearchNearest()
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TParallelSpatialSearch parallel_spatial_search_bins(r_skin_part.Elements(), r_data_comm);

    double epsilon = 1.0e-6;
    const std::size_t near_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(near_point_id, epsilon,epsilon,epsilon);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search_bins.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

        // Distances are just local
        std::vector<std::vector<double>> distances;
        results.GetDistances(distances);
        KRATOS_EXPECT_NEAR(distances[0][0], (cube_z - epsilon), tolerance);

        // Compute indices
        std::vector<std::vector<IndexType>> indices;
        results.GetResultIndices(indices);
        const std::size_t id = indices[0][0];
        KRATOS_EXPECT_EQ(id, 3);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

/**
 * @brief Test function for empty nearest searching geometrical objects in bins using a given search wrapper
 * @tparam TParallelSpatialSearch The type of search wrapper to be tested
 */
template<class TParallelSpatialSearch>
void TestParallelSpatialSearchGeometricalObjectsBinsEmptySearchNearest()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TParallelSpatialSearch parallel_spatial_search_bins(r_skin_part.Elements(), r_data_comm);

    const std::size_t point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(point_id, 0.0,0.0,0.0);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search_bins.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

/**
 * @brief Test function for inside searching geometrical objects in bins using a given search wrapper
 * @tparam TParallelSpatialSearch The type of search wrapper to be tested
 */
template<class TParallelSpatialSearch>
void TestParallelSpatialSearchGeometricalObjectsBinsSearchIsInside()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeModelPart(current_model);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TParallelSpatialSearch parallel_spatial_search_bins(r_skin_part.Elements(), r_data_comm);

    const std::size_t inside_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(inside_point_id, 0.5,0.5,0.5);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search_bins.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

/**
 * @brief Test function for non inside searching geometrical objects in bins using a given search wrapper
 * @tparam TParallelSpatialSearch The type of search wrapper to be tested
 */
template<class TParallelSpatialSearch>
void TestParallelSpatialSearchGeometricalObjectsBinsSearchIsNotInside()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeModelPart(current_model);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TParallelSpatialSearch parallel_spatial_search_bins(r_skin_part.Elements(), r_data_comm);

    const std::size_t outside_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(outside_point_id, 100.0,100.0,100.0);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search_bins.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only in first rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

// Definition of the trees search wrapper
using ParallelSpatialSearchKDTreeElement = ParallelSpatialSearch<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>;
using ParallelSpatialSearchOCTreeElement = ParallelSpatialSearch<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>;
using ParallelSpatialSearchStaticBinsTreeElement = ParallelSpatialSearch<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS>;
using ParallelSpatialSearchBinsDynamicElement = ParallelSpatialSearch<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS>;

/**
 * @brief A function to test tree-based search in a specified radius.
 * @details This function tests the tree-based search algorithm within a specified radius.
 * @tparam TParallelSpatialSearch The type of the search wrapper.
 */
template<class TParallelSpatialSearch>
void TestTreeSearchInRadius()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

    // Generate the search wrapper for bins
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    TParallelSpatialSearch parallel_spatial_search(r_skin_part.Elements(), r_data_comm);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    const std::size_t point_id = 1;
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(point_id, 0.0, 0.0, 0.0);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;

    // 0.3 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }

    // 2.0 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 2.0, results);
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

/**
 * @brief A function to test tree-based nearest search in a specified radius.
 * @details This function tests the tree-based nearest search algorithm within a specified radius.
 * @tparam TParallelSpatialSearch The type of the search wrapper.
 */
template<class TParallelSpatialSearch>
void TestTreeSearchNearestInRadius()
{
    constexpr double tolerance = 1e-6;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    // Generate the search wrapper for bins
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    TParallelSpatialSearch parallel_spatial_search(r_skin_part.Elements(), r_data_comm);

    double epsilon = 1.0e-6;
    const std::size_t near_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(near_point_id, epsilon,epsilon,epsilon);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z, results);

    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }

    parallel_spatial_search.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);

    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

        // Distances
        std::vector<std::vector<double>> distances;
        results.GetDistances(distances);
        KRATOS_EXPECT_NEAR(distances[0][0], (cube_z - 0.08 - epsilon), tolerance);

        // Compute indices
        std::vector<std::vector<IndexType>> indices;
        results.GetResultIndices(indices);
        const std::size_t id = indices[0][0];
        KRATOS_EXPECT_EQ(id, 4);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

/**
 * @brief A function to test tree-based nearest search
 * @details This function tests the tree-based nearest search algorithm
 * @tparam TParallelSpatialSearch The type of the search wrapper.
 */
template<class TParallelSpatialSearch>
void TestTreeSearchNearest()
{
    constexpr double tolerance = 1e-6;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    // Generate the search wrapper for bins
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    TParallelSpatialSearch parallel_spatial_search(r_skin_part.Elements(), r_data_comm);

    double epsilon = 1.0e-6;
    const std::size_t near_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(near_point_id, epsilon,epsilon,epsilon);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only one result in first and second rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

        // Distances
        std::vector<std::vector<double>> distances;
        results.GetDistances(distances);
        KRATOS_EXPECT_NEAR(distances[0][0], (cube_z - 0.08 - epsilon), tolerance);

        // Compute indices
        std::vector<std::vector<IndexType>> indices;
        results.GetResultIndices(indices);
        const std::size_t id = indices[0][0];
        KRATOS_EXPECT_EQ(id, 4);
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

/**
 * @brief A function to test tree-based nearest search (empty)
 * @details This function tests the tree-based nearest search algorithm (empty)
 * @tparam TParallelSpatialSearch The type of the search wrapper.
 */
template<class TParallelSpatialSearch>
void TestTreeSearchNearestEmpty()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");

    // Generate the search wrapper for bins
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    TParallelSpatialSearch parallel_spatial_search(r_skin_part.Elements(), r_data_comm);

    const std::size_t point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(point_id, 0.0,0.0,0.0);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    // Only in first rank
    if constexpr (TParallelSpatialSearch::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        KRATOS_ERROR << "Asynchronous spatial search not implemented!" << std::endl;
    }
}

} // namespace

namespace Testing
{

/** Checks parallel_spatial_search_bins search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchGeometricalObjectsBinsSearchInRadius, KratosMPICoreFastSuite)
{
    TestParallelSpatialSearchGeometricalObjectsBinsSearchInRadius<ParallelSpatialSearchGeometricalObjectsBins>();
}

/** Checks parallel_spatial_search_bins search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchGeometricalObjectsBinsSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestParallelSpatialSearchGeometricalObjectsBinsSearchNearestInRadius<ParallelSpatialSearchGeometricalObjectsBins>();
}

/** Checks parallel_spatial_search_bins search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchGeometricalObjectsBinsSearchNearest, KratosMPICoreFastSuite)
{
    TestParallelSpatialSearchGeometricalObjectsBinsSearchNearest<ParallelSpatialSearchGeometricalObjectsBins>();
}

/** Checks parallel_spatial_search_bins empty search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchGeometricalObjectsBinsEmptySearchNearest, KratosMPICoreFastSuite)
{
    TestParallelSpatialSearchGeometricalObjectsBinsEmptySearchNearest<ParallelSpatialSearchGeometricalObjectsBins>();
}

/** Checks parallel_spatial_search_bins search is inside
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchGeometricalObjectsBinsSearchIsInside, KratosMPICoreFastSuite)
{
    TestParallelSpatialSearchGeometricalObjectsBinsSearchIsInside<ParallelSpatialSearchGeometricalObjectsBins>();
}

/** Checks parallel_spatial_search_bins search is inside = not found
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchGeometricalObjectsBinsSearchIsNotInside, KratosMPICoreFastSuite)
{
    TestParallelSpatialSearchGeometricalObjectsBinsSearchIsNotInside<ParallelSpatialSearchGeometricalObjectsBins>();
}

/** Checks ParallelSpatialSearch works for KDTreeElement search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchKDTreeElementSearchInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchInRadius<ParallelSpatialSearchKDTreeElement>();
}

/** Checks ParallelSpatialSearch works for OCTreeElement search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchOCTreeElementSearchInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchInRadius<ParallelSpatialSearchOCTreeElement>();
}

/** Checks ParallelSpatialSearch works for StaticBinsTree search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchStaticBinsTreeElementSearchInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchInRadius<ParallelSpatialSearchStaticBinsTreeElement>();
}

/** Checks ParallelSpatialSearch works for BinsDynamicElement search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchBinsDynamicElementSearchInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchInRadius<ParallelSpatialSearchBinsDynamicElement>();
}

/** Checks ParallelSpatialSearch works for KDTreeElement search nearest in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchKDTreeElementSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestInRadius<ParallelSpatialSearchKDTreeElement>();
}

/** Checks ParallelSpatialSearch works for OCTreeElement search nearest in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchOCTreeElementSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestInRadius<ParallelSpatialSearchOCTreeElement>();
}

/** Checks ParallelSpatialSearch works for StaticBinsTree search nearest in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchStaticBinsTreeElementSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestInRadius<ParallelSpatialSearchStaticBinsTreeElement>();
}

/** Checks ParallelSpatialSearch works for BinsDynamicElement search nearest in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchBinsDynamicElementSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestInRadius<ParallelSpatialSearchBinsDynamicElement>();
}

/** Checks ParallelSpatialSearch works for KDTreeElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchKDTreeElementSearchNearest, KratosMPICoreFastSuite)
{
    TestTreeSearchNearest<ParallelSpatialSearchKDTreeElement>();
}

/** Checks ParallelSpatialSearch works for OCTreeElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchOCTreeElementSearchNearest, KratosMPICoreFastSuite)
{
    TestTreeSearchNearest<ParallelSpatialSearchOCTreeElement>();
}

/** Checks ParallelSpatialSearch works for StaticBinsTree search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchStaticBinsTreeElementSearchNearest, KratosMPICoreFastSuite)
{
    TestTreeSearchNearest<ParallelSpatialSearchStaticBinsTreeElement>();
}

/** Checks ParallelSpatialSearch works for BinsDynamicElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchBinsDynamicElementSearchNearest, KratosMPICoreFastSuite)
{
    TestTreeSearchNearest<ParallelSpatialSearchBinsDynamicElement>();
}

/** Checks ParallelSpatialSearch works for KDTreeElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchKDTreeElementSearchNearestEmpty, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestEmpty<ParallelSpatialSearchKDTreeElement>();
}

/** Checks ParallelSpatialSearch works for OCTreeElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchOCTreeElementSearchNearestEmpty, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestEmpty<ParallelSpatialSearchOCTreeElement>();
}

/** Checks ParallelSpatialSearch works for StaticBinsTree search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchStaticBinsTreeElementSearchNearestEmpty, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestEmpty<ParallelSpatialSearchStaticBinsTreeElement>();
}

/** Checks ParallelSpatialSearch works for BinsDynamicElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIParallelSpatialSearchBinsDynamicElementSearchNearestEmpty, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestEmpty<ParallelSpatialSearchBinsDynamicElement>();
}

} // namespace Testing.

} // namespace Kratos.