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
#include "spatial_containers/search_wrapper.h"
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
// SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS
using SearchWrapperGeometricalObjectsBins = SearchWrapper<GeometricalObjectsBins>;
// SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS
using SearchWrapperGeometricalObjectsBinsHetereogeneous = SearchWrapper<GeometricalObjectsBins, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;

/**
 * @brief Test function for searching geometrical objects in bins using a given search wrapper
 * @tparam TSearchWrapper The type of search wrapper to be tested
 */
template<class TSearchWrapper>
void TestSearchWrapperGeometricalObjectsBinsSearchInRadius()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TSearchWrapper search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;

    // 0.29 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.29, results);
    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }

    // 0.3 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }

    // 0.4 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.4, results);
    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }

    // 0.6 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);
    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }

    // 0.7 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.7, results);
    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }

    // 0.9 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.9, results);
    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

/**
 * @brief Test function for nearest searching geometrical objects in bins using a given search wrapper
 * @tparam TSearchWrapper The type of search wrapper to be tested
 */
template<class TSearchWrapper>
void TestSearchWrapperGeometricalObjectsBinsSearchNearestInRadius()
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TSearchWrapper search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper_bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z - 1.e-4, results);

    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
        }
    }

    search_wrapper_bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z + 1.e-4, results);

    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

        // Distances are just local
        const auto distances = results[0].GetDistances();
        KRATOS_EXPECT_NEAR(distances[0], (cube_z - epsilon), tolerance);

        // Compute indices
        auto indices = results[0].GetResultIndices();
        const std::size_t id = indices[0];
        KRATOS_EXPECT_EQ(id, 3);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

            // Distances are just local
            const auto distances = results[0].GetDistances();
            KRATOS_EXPECT_NEAR(distances[0], (cube_z - epsilon), tolerance);

            // Compute indices
            auto indices = results[0].GetResultIndices();
            const std::size_t id = indices[0];
            KRATOS_EXPECT_EQ(id, 3);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

/**
 * @brief Test function for nearest searching geometrical objects in bins using a given search wrapper
 * @tparam TSearchWrapper The type of search wrapper to be tested
 */
template<class TSearchWrapper>
void TestSearchWrapperGeometricalObjectsBinsSearchNearest()
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TSearchWrapper search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper_bins.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

        // Distances are just local
        const auto distances = results[0].GetDistances();
        KRATOS_EXPECT_NEAR(distances[0], (cube_z - epsilon), tolerance);

        // Compute indices
        auto indices = results[0].GetResultIndices();
        const std::size_t id = indices[0];
        KRATOS_EXPECT_EQ(id, 3);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

            // Distances are just local
            const auto distances = results[0].GetDistances();
            KRATOS_EXPECT_NEAR(distances[0], (cube_z - epsilon), tolerance);

            // Compute indices
            auto indices = results[0].GetResultIndices();
            const std::size_t id = indices[0];
            KRATOS_EXPECT_EQ(id, 3);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

/**
 * @brief Test function for empty nearest searching geometrical objects in bins using a given search wrapper
 * @tparam TSearchWrapper The type of search wrapper to be tested
 */
template<class TSearchWrapper>
void TestSearchWrapperGeometricalObjectsBinsEmptySearchNearest()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TSearchWrapper search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper_bins.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        if (rank == 0) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

/**
 * @brief Test function for inside searching geometrical objects in bins using a given search wrapper
 * @tparam TSearchWrapper The type of search wrapper to be tested
 */
template<class TSearchWrapper>
void TestSearchWrapperGeometricalObjectsBinsSearchIsInside()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeModelPart(current_model);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TSearchWrapper search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper_bins.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

/**
 * @brief Test function for non inside searching geometrical objects in bins using a given search wrapper
 * @tparam TSearchWrapper The type of search wrapper to be tested
 */
template<class TSearchWrapper>
void TestSearchWrapperGeometricalObjectsBinsSearchIsNotInside()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeModelPart(current_model);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    TSearchWrapper search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper_bins.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only in first rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        if (rank == 0) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

// Definition of the trees search wrapper
// SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS
using SearchWrapperKDTreeElement = SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
using SearchWrapperOCTreeElement = SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
using SearchWrapperStaticBinsTreeElement = SearchWrapper<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
using SearchWrapperBinsDynamicElement = SearchWrapper<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;

// SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS
using SearchWrapperKDTreeElementHetereogeneous = SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
using SearchWrapperOCTreeElementHetereogeneous = SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
using SearchWrapperStaticBinsTreeElementHetereogeneous = SearchWrapper<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
using SearchWrapperBinsDynamicElementHetereogeneous = SearchWrapper<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;

/**
 * @brief A function to test tree-based search in a specified radius.
 * @details This function tests the tree-based search algorithm within a specified radius.
 * @tparam TSearchWrapper The type of the search wrapper.
 */
template<class TSearchWrapper>
void TestTreeSearchInRadius()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

    // Generate the search wrapper for bins
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    TSearchWrapper search_wrapper(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;

    // 0.3 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }

    // 2.0 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 2.0, results);
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

/**
 * @brief A function to test tree-based nearest search in a specified radius.
 * @details This function tests the tree-based nearest search algorithm within a specified radius.
 * @tparam TSearchWrapper The type of the search wrapper.
 */
template<class TSearchWrapper>
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
    TSearchWrapper search_wrapper(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z, results);

    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }

    search_wrapper.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);

    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

        // Distances
        auto distances = results[0].GetDistances();
        KRATOS_EXPECT_NEAR(distances[0], (cube_z - 0.08 - epsilon), tolerance);

        // Compute indices
        auto indices = results[0].GetResultIndices();
        const std::size_t id = indices[0];
        KRATOS_EXPECT_EQ(id, 4);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

            // Distances
            auto distances = results[0].GetDistances();
            KRATOS_EXPECT_NEAR(distances[0], (cube_z - 0.08 - epsilon), tolerance);

            // Compute indices
            auto indices = results[0].GetResultIndices();
            const std::size_t id = indices[0];
            KRATOS_EXPECT_EQ(id, 4);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

/**
 * @brief A function to test tree-based nearest search
 * @details This function tests the tree-based nearest search algorithm
 * @tparam TSearchWrapper The type of the search wrapper.
 */
template<class TSearchWrapper>
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
    TSearchWrapper search_wrapper(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    // We expect only one result in first and second rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
        KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

        // Distances
        auto distances = results[0].GetDistances();
        KRATOS_EXPECT_NEAR(distances[0], (cube_z - 0.08 - epsilon), tolerance);

        // Compute indices
        auto indices = results[0].GetResultIndices();
        const std::size_t id = indices[0];
        KRATOS_EXPECT_EQ(id, 4);
    } else {
        if (rank == 0 || rank == 1) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
            KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

            // Distances
            auto distances = results[0].GetDistances();
            KRATOS_EXPECT_NEAR(distances[0], (cube_z - 0.08 - epsilon), tolerance);

            // Compute indices
            auto indices = results[0].GetResultIndices();
            const std::size_t id = indices[0];
            KRATOS_EXPECT_EQ(id, 4);
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

/**
 * @brief A function to test tree-based nearest search (empty)
 * @details This function tests the tree-based nearest search algorithm (empty)
 * @tparam TSearchWrapper The type of the search wrapper.
 */
template<class TSearchWrapper>
void TestTreeSearchNearestEmpty()
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");

    // Generate the search wrapper for bins
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    TSearchWrapper search_wrapper(r_skin_part.Elements(), r_data_comm);

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

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    // Only in first rank
    if constexpr (TSearchWrapper::GetSpatialSearchCommunication() == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS) {
        KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
        KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    } else {
        if (rank == 0) {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
            KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
        } else {
            KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 0);
        }
    }
}

} // namespace


namespace Testing
{


/** Checks search_wrapper_bins search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchInRadius, KratosMPICoreFastSuite)
{
    TestSearchWrapperGeometricalObjectsBinsSearchInRadius<SearchWrapperGeometricalObjectsBins>();
    TestSearchWrapperGeometricalObjectsBinsSearchInRadius<SearchWrapperGeometricalObjectsBinsHetereogeneous>();
}

/** Checks search_wrapper_bins search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestSearchWrapperGeometricalObjectsBinsSearchNearestInRadius<SearchWrapperGeometricalObjectsBins>();
    TestSearchWrapperGeometricalObjectsBinsSearchNearestInRadius<SearchWrapperGeometricalObjectsBinsHetereogeneous>();
}

/** Checks search_wrapper_bins search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchNearest, KratosMPICoreFastSuite)
{
    TestSearchWrapperGeometricalObjectsBinsSearchNearest<SearchWrapperGeometricalObjectsBins>();
    TestSearchWrapperGeometricalObjectsBinsSearchNearest<SearchWrapperGeometricalObjectsBinsHetereogeneous>();
}

/** Checks search_wrapper_bins empty search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsEmptySearchNearest, KratosMPICoreFastSuite)
{
    TestSearchWrapperGeometricalObjectsBinsEmptySearchNearest<SearchWrapperGeometricalObjectsBins>();
    TestSearchWrapperGeometricalObjectsBinsEmptySearchNearest<SearchWrapperGeometricalObjectsBinsHetereogeneous>();
}

/** Checks search_wrapper_bins search is inside
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchIsInside, KratosMPICoreFastSuite)
{
    TestSearchWrapperGeometricalObjectsBinsSearchIsInside<SearchWrapperGeometricalObjectsBins>();
    TestSearchWrapperGeometricalObjectsBinsSearchIsInside<SearchWrapperGeometricalObjectsBinsHetereogeneous>();
}

/** Checks search_wrapper_bins search is inside = not found
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchIsNotInside, KratosMPICoreFastSuite)
{
    TestSearchWrapperGeometricalObjectsBinsSearchIsNotInside<SearchWrapperGeometricalObjectsBins>();
    TestSearchWrapperGeometricalObjectsBinsSearchIsNotInside<SearchWrapperGeometricalObjectsBinsHetereogeneous>();
}

/** Checks SearchWrapper works for KDTreeElement search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperKDTreeElementSearchInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchInRadius<SearchWrapperKDTreeElement>();
    TestTreeSearchInRadius<SearchWrapperKDTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for OCTreeElement search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperOCTreeElementSearchInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchInRadius<SearchWrapperOCTreeElement>();
    TestTreeSearchInRadius<SearchWrapperOCTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for StaticBinsTree search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperStaticBinsTreeElementSearchInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchInRadius<SearchWrapperStaticBinsTreeElement>();
    TestTreeSearchInRadius<SearchWrapperStaticBinsTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for BinsDynamicElement search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperBinsDynamicElementSearchInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchInRadius<SearchWrapperBinsDynamicElement>();
    TestTreeSearchInRadius<SearchWrapperBinsDynamicElementHetereogeneous>();
}

/** Checks SearchWrapper works for KDTreeElement search nearest in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperKDTreeElementSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestInRadius<SearchWrapperKDTreeElement>();
    TestTreeSearchNearestInRadius<SearchWrapperKDTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for OCTreeElement search nearest in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperOCTreeElementSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestInRadius<SearchWrapperOCTreeElement>();
    TestTreeSearchNearestInRadius<SearchWrapperOCTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for StaticBinsTree search nearest in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperStaticBinsTreeElementSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestInRadius<SearchWrapperStaticBinsTreeElement>();
    TestTreeSearchNearestInRadius<SearchWrapperStaticBinsTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for BinsDynamicElement search nearest in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperBinsDynamicElementSearchNearestInRadius, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestInRadius<SearchWrapperBinsDynamicElement>();
    TestTreeSearchNearestInRadius<SearchWrapperBinsDynamicElementHetereogeneous>();
}

/** Checks SearchWrapper works for KDTreeElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperKDTreeElementSearchNearest, KratosMPICoreFastSuite)
{
    TestTreeSearchNearest<SearchWrapperKDTreeElement>();
    TestTreeSearchNearest<SearchWrapperKDTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for OCTreeElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperOCTreeElementSearchNearest, KratosMPICoreFastSuite)
{
    TestTreeSearchNearest<SearchWrapperOCTreeElement>();
    TestTreeSearchNearest<SearchWrapperOCTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for StaticBinsTree search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperStaticBinsTreeElementSearchNearest, KratosMPICoreFastSuite)
{
    TestTreeSearchNearest<SearchWrapperStaticBinsTreeElement>();
    TestTreeSearchNearest<SearchWrapperStaticBinsTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for BinsDynamicElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperBinsDynamicElementSearchNearest, KratosMPICoreFastSuite)
{
    TestTreeSearchNearest<SearchWrapperBinsDynamicElement>();
    TestTreeSearchNearest<SearchWrapperBinsDynamicElementHetereogeneous>();
}

/** Checks SearchWrapper works for KDTreeElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperKDTreeElementSearchNearestEmpty, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestEmpty<SearchWrapperKDTreeElement>();
    TestTreeSearchNearestEmpty<SearchWrapperKDTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for OCTreeElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperOCTreeElementSearchNearestEmpty, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestEmpty<SearchWrapperOCTreeElement>();
    TestTreeSearchNearestEmpty<SearchWrapperOCTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for StaticBinsTree search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperStaticBinsTreeElementSearchNearestEmpty, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestEmpty<SearchWrapperStaticBinsTreeElement>();
    TestTreeSearchNearestEmpty<SearchWrapperStaticBinsTreeElementHetereogeneous>();
}

/** Checks SearchWrapper works for BinsDynamicElement search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperBinsDynamicElementSearchNearestEmpty, KratosMPICoreFastSuite)
{
    TestTreeSearchNearestEmpty<SearchWrapperBinsDynamicElement>();
    TestTreeSearchNearestEmpty<SearchWrapperBinsDynamicElementHetereogeneous>();
}

} // namespace Testing.

} // namespace Kratos.