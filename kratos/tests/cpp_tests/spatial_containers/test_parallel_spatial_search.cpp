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
#include "testing/testing.h"
#include "spatial_containers/parallel_spatial_search.h"
#include "includes/data_communicator.h"
#include "containers/model.h"
#include "geometries/triangle_3d_3.h"
#include "tests/test_utilities/cpp_tests_utilities.h"

namespace Kratos::Testing
{

// Definition of the geometrical object bins search wrapper
using ParallelSpatialSearchGeometricalObjectsBins = ParallelSpatialSearch<GeometricalObjectsBins>;

/** Checks ParallelSpatialSearch works for GeometricalObjectBins search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchGeometricalObjectsBinsSearchInRadius, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    ParallelSpatialSearchGeometricalObjectsBins parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_point = r_point_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    ParallelSpatialSearchGeometricalObjectsBins::ResultContainerVectorType results;

    // 0.29 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.29, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);

    // 0.3 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);

    // 0.4 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.4, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);

    // 0.6 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);

    // 0.7 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.7, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);

    // 0.9 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.9, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
}

/** Checks ParallelSpatialSearch works for GeometricalObjectBins search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchGeometricalObjectsBinsSearchNearestInRadius, KratosCoreFastSuite)
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    ParallelSpatialSearchGeometricalObjectsBins parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    double epsilon = 1.0e-6;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_near_point = r_point_model_part.CreateNewNode(1, epsilon,epsilon,epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    ParallelSpatialSearchGeometricalObjectsBins::ResultContainerVectorType results;
    parallel_spatial_search.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z - 1.e-4, results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());

    parallel_spatial_search.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z + 1.e-4, results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

    // Distances
    std::vector<std::vector<double>> distances;
    results.GetDistances(distances);
    KRATOS_EXPECT_NEAR(distances[0][0], (cube_z - epsilon), tolerance);

    // Compute indices
    std::vector<std::vector<std::size_t>> indices;
    results.GetResultIndices(indices);
    const std::size_t id = indices[0][0];
    KRATOS_EXPECT_EQ(id, 3);
}

/** Checks ParallelSpatialSearch works for GeometricalObjectBins search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchGeometricalObjectsBinsSearchNearest, KratosCoreFastSuite)
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    ParallelSpatialSearchGeometricalObjectsBins parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    double epsilon = 1.0e-6;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_near_point = r_point_model_part.CreateNewNode(1, epsilon,epsilon,epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    ParallelSpatialSearchGeometricalObjectsBins::ResultContainerVectorType results;
    parallel_spatial_search.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

    // Distances
    std::vector<std::vector<double>> distances;
    results.GetDistances(distances);
    KRATOS_EXPECT_NEAR(distances[0][0], (cube_z - epsilon), tolerance);

    // Compute indices
    std::vector<std::vector<std::size_t>> indices;
    results.GetResultIndices(indices);
    const std::size_t id = indices[0][0];
    KRATOS_EXPECT_EQ(id, 3);
}

/** Checks ParallelSpatialSearch works for GeometricalObjectBins empty search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchGeometricalObjectsBinsEmptySearchNearest, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    ParallelSpatialSearchGeometricalObjectsBins parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_point = r_point_model_part.CreateNewNode(1, 0.0,0.0,0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    ParallelSpatialSearchGeometricalObjectsBins::ResultContainerVectorType results;
    parallel_spatial_search.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
}

/** Checks ParallelSpatialSearch works for GeometricalObjectBins search is inside 
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchGeometricalObjectsBinsSearchIsInside, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeModelPart(current_model);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    ParallelSpatialSearchGeometricalObjectsBins parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_inside_point = r_point_model_part.CreateNewNode(1, 0.5,0.5,0.5);
    auto& r_array_nodes = r_point_model_part.Nodes();

    ParallelSpatialSearchGeometricalObjectsBins::ResultContainerVectorType results;
    parallel_spatial_search.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);
}

/** Checks ParallelSpatialSearch works for GeometricalObjectBins search is inside = not found
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchGeometricalObjectsBinsSearchIsNotInside, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeModelPart(current_model);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    ParallelSpatialSearchGeometricalObjectsBins parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_outside_point = r_point_model_part.CreateNewNode(1, 100.0,100.0,100.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    ParallelSpatialSearchGeometricalObjectsBins::ResultContainerVectorType results;
    parallel_spatial_search.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
}

// Definition of the trees search wrapper
using ParallelSpatialSearchKDTreeElement = ParallelSpatialSearch<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>>;
using ParallelSpatialSearchOCTreeElement = ParallelSpatialSearch<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>>;
using ParallelSpatialSearchStaticBinsTreeElement = ParallelSpatialSearch<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>;
using ParallelSpatialSearchBinsDynamicElement = ParallelSpatialSearch<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>;

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
    DataCommunicator serial_communicator;
    TParallelSpatialSearch parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_point = r_point_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;

    // 0.3 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);

    // 2.0 radius
    parallel_spatial_search.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 2.0, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
}

/** Checks ParallelSpatialSearch works for KDTreeElement search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchKDTreeElementSearchInRadius, KratosCoreFastSuite)
{
    TestTreeSearchInRadius<ParallelSpatialSearchKDTreeElement>();
}

/** Checks ParallelSpatialSearch works for OCTreeElement search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchOCTreeElementSearchInRadius, KratosCoreFastSuite)
{
    TestTreeSearchInRadius<ParallelSpatialSearchOCTreeElement>();
}

/** Checks ParallelSpatialSearch works for StaticBinsTree search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchStaticBinsTreeElementSearchInRadius, KratosCoreFastSuite)
{
    TestTreeSearchInRadius<ParallelSpatialSearchStaticBinsTreeElement>();
}

/** Checks ParallelSpatialSearch works for BinsDynamicElement search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchBinsDynamicElementSearchInRadius, KratosCoreFastSuite)
{
    TestTreeSearchInRadius<ParallelSpatialSearchBinsDynamicElement>();
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
    DataCommunicator serial_communicator;
    TParallelSpatialSearch parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    double epsilon = 1.0e-6;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_near_point = r_point_model_part.CreateNewNode(1, epsilon,epsilon,epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z, results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());

    parallel_spatial_search.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

    // Distances
    std::vector<std::vector<double>> distances;
    results.GetDistances(distances);
    KRATOS_EXPECT_NEAR(distances[0][0], (cube_z - 0.08 - epsilon), tolerance);

    // Compute indices
    std::vector<std::vector<std::size_t>> indices;
    results.GetResultIndices(indices);
    const std::size_t id = indices[0][0];
    KRATOS_EXPECT_EQ(id, 4);
}

/** Checks ParallelSpatialSearch works for KDTreeElement search nearest in radius
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchKDTreeElementSearchNearestInRadius, KratosCoreFastSuite)
{
    TestTreeSearchNearestInRadius<ParallelSpatialSearchKDTreeElement>();
}

/** Checks ParallelSpatialSearch works for OCTreeElement search nearest in radius
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchOCTreeElementSearchNearestInRadius, KratosCoreFastSuite)
{
    TestTreeSearchNearestInRadius<ParallelSpatialSearchOCTreeElement>();
}

/** Checks ParallelSpatialSearch works for StaticBinsTree search nearest in radius
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchStaticBinsTreeElementSearchNearestInRadius, KratosCoreFastSuite)
{
    TestTreeSearchNearestInRadius<ParallelSpatialSearchStaticBinsTreeElement>();
}

/** Checks ParallelSpatialSearch works for BinsDynamicElement search nearest in radius
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchBinsDynamicElementSearchNearestInRadius, KratosCoreFastSuite)
{
    TestTreeSearchNearestInRadius<ParallelSpatialSearchBinsDynamicElement>();
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
    DataCommunicator serial_communicator;
    TParallelSpatialSearch parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    double epsilon = 1.0e-6;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_near_point = r_point_model_part.CreateNewNode(1, epsilon,epsilon,epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

    // Distances
    std::vector<std::vector<double>> distances;
    results.GetDistances(distances);
    KRATOS_EXPECT_NEAR(distances[0][0], (cube_z - 0.08 - epsilon), tolerance);

    // Compute indices
    std::vector<std::vector<std::size_t>> indices;
    results.GetResultIndices(indices);
    const std::size_t id = indices[0][0];
    KRATOS_EXPECT_EQ(id, 4);
}

/** Checks ParallelSpatialSearch works for KDTreeElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchKDTreeElementSearchNearest, KratosCoreFastSuite)
{
    TestTreeSearchNearest<ParallelSpatialSearchKDTreeElement>();
}

/** Checks ParallelSpatialSearch works for OCTreeElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchOCTreeElementSearchNearest, KratosCoreFastSuite)
{
    TestTreeSearchNearest<ParallelSpatialSearchOCTreeElement>();
}

/** Checks ParallelSpatialSearch works for StaticBinsTree search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchStaticBinsTreeElementSearchNearest, KratosCoreFastSuite)
{
    TestTreeSearchNearest<ParallelSpatialSearchStaticBinsTreeElement>();
}

/** Checks ParallelSpatialSearch works for BinsDynamicElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchBinsDynamicElementSearchNearest, KratosCoreFastSuite)
{
    TestTreeSearchNearest<ParallelSpatialSearchBinsDynamicElement>();
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
    DataCommunicator serial_communicator;
    TParallelSpatialSearch parallel_spatial_search(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_point = r_point_model_part.CreateNewNode(1, 0.0,0.0,0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TParallelSpatialSearch::ResultContainerVectorType results;
    parallel_spatial_search.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
}

/** Checks ParallelSpatialSearch works for KDTreeElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchKDTreeElementSearchNearestEmpty, KratosCoreFastSuite)
{
    TestTreeSearchNearestEmpty<ParallelSpatialSearchKDTreeElement>();
}

/** Checks ParallelSpatialSearch works for OCTreeElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchOCTreeElementSearchNearestEmpty, KratosCoreFastSuite)
{
    TestTreeSearchNearestEmpty<ParallelSpatialSearchOCTreeElement>();
}

/** Checks ParallelSpatialSearch works for StaticBinsTree search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchStaticBinsTreeElementSearchNearestEmpty, KratosCoreFastSuite)
{
    TestTreeSearchNearestEmpty<ParallelSpatialSearchStaticBinsTreeElement>();
}

/** Checks ParallelSpatialSearch works for BinsDynamicElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(ParallelSpatialSearchBinsDynamicElementSearchNearestEmpty, KratosCoreFastSuite)
{
    TestTreeSearchNearestEmpty<ParallelSpatialSearchBinsDynamicElement>();
}

} // namespace Kratos::Testing.
