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
#include "spatial_containers/search_wrapper.h"
#include "includes/data_communicator.h"
#include "containers/model.h"
#include "geometries/triangle_3d_3.h"
#include "tests/test_utilities/cpp_tests_utilities.h"

namespace Kratos::Testing
{

// Definition of the geometrical object bins search wrapper
using SearchWrapperGeometricalObjectsBins = SearchWrapper<GeometricalObjectsBins>;

/** Checks SearchWrapper works for GeometricalObjectBins search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperGeometricalObjectsBinsSearchInRadius, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    SearchWrapperGeometricalObjectsBins search_wrapper(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_point = r_point_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;

    // 0.29 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.29, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);

    // 0.3 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);

    // 0.4 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.4, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 4);

    // 0.6 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);

    // 0.7 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.7, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 8);

    // 0.9 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.9, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
}

/** Checks SearchWrapper works for GeometricalObjectBins search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperGeometricalObjectsBinsSearchNearestInRadius, KratosCoreFastSuite)
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    SearchWrapperGeometricalObjectsBins search_wrapper(r_skin_part.Elements(), serial_communicator);

    double epsilon = 1.0e-6;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_near_point = r_point_model_part.CreateNewNode(1, epsilon,epsilon,epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z - 1.e-4, results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());

    search_wrapper.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z + 1.e-4, results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

    // Distances
    auto distances = results[0].GetDistances();
    KRATOS_EXPECT_NEAR(distances[0], (cube_z - epsilon), tolerance);

    // Compute indices
    auto indices = results[0].GetResultIndices();
    const std::size_t id = indices[0];
    KRATOS_EXPECT_EQ(id, 3);
}

/** Checks SearchWrapper works for GeometricalObjectBins search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperGeometricalObjectsBinsSearchNearest, KratosCoreFastSuite)
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    SearchWrapperGeometricalObjectsBins search_wrapper(r_skin_part.Elements(), serial_communicator);

    double epsilon = 1.0e-6;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_near_point = r_point_model_part.CreateNewNode(1, epsilon,epsilon,epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);

    // Distances
    auto distances = results[0].GetDistances();
    KRATOS_EXPECT_NEAR(distances[0], (cube_z - epsilon), tolerance);

    // Compute indices
    auto indices = results[0].GetResultIndices();
    const std::size_t id = indices[0];
    KRATOS_EXPECT_EQ(id, 3);
}

/** Checks SearchWrapper works for GeometricalObjectBins empty search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperGeometricalObjectsBinsEmptySearchNearest, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    SearchWrapperGeometricalObjectsBins search_wrapper(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_point = r_point_model_part.CreateNewNode(1, 0.0,0.0,0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
}

/** Checks SearchWrapper works for GeometricalObjectBins search is inside 
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperGeometricalObjectsBinsSearchIsInside, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeModelPart(current_model);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    SearchWrapperGeometricalObjectsBins search_wrapper(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_inside_point = r_point_model_part.CreateNewNode(1, 0.5,0.5,0.5);
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 1);
}

/** Checks SearchWrapper works for GeometricalObjectBins search is inside = not found
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperGeometricalObjectsBinsSearchIsNotInside, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeModelPart(current_model);

    // Generate the search wrapper for bins
    DataCommunicator serial_communicator;
    SearchWrapperGeometricalObjectsBins search_wrapper(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_outside_point = r_point_model_part.CreateNewNode(1, 100.0,100.0,100.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
}

// Definition of the trees search wrapper
using SearchWrapperKDTreeElement = SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>>;
using SearchWrapperOCTreeElement = SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>>;
using SearchWrapperStaticBinsTreeElement = SearchWrapper<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>;
using SearchWrapperBinsDynamicElement = SearchWrapper<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>;

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
    DataCommunicator serial_communicator;
    TSearchWrapper search_wrapper(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_point = r_point_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TSearchWrapper::ResultContainerVectorType results;

    // 0.3 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 0);

    // 2.0 radius
    search_wrapper.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 2.0, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[0].IsObjectFound());
    KRATOS_EXPECT_EQ(results[0].NumberOfGlobalResults(), 12);
}

/** Checks SearchWrapper works for KDTreeElement search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperKDTreeElementSearchInRadius, KratosCoreFastSuite)
{
    TestTreeSearchInRadius<SearchWrapperKDTreeElement>();
}

/** Checks SearchWrapper works for OCTreeElement search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperOCTreeElementSearchInRadius, KratosCoreFastSuite)
{
    TestTreeSearchInRadius<SearchWrapperOCTreeElement>();
}

/** Checks SearchWrapper works for StaticBinsTree search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperStaticBinsTreeElementSearchInRadius, KratosCoreFastSuite)
{
    TestTreeSearchInRadius<SearchWrapperStaticBinsTreeElement>();
}

/** Checks SearchWrapper works for BinsDynamicElement search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperBinsDynamicElementSearchInRadius, KratosCoreFastSuite)
{
    TestTreeSearchInRadius<SearchWrapperBinsDynamicElement>();
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
    DataCommunicator serial_communicator;
    TSearchWrapper search_wrapper(r_skin_part.Elements(), serial_communicator);

    double epsilon = 1.0e-6;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_near_point = r_point_model_part.CreateNewNode(1, epsilon,epsilon,epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z, results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());

    search_wrapper.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);

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
}

/** Checks SearchWrapper works for KDTreeElement search nearest in radius
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperKDTreeElementSearchNearestInRadius, KratosCoreFastSuite)
{
    TestTreeSearchNearestInRadius<SearchWrapperKDTreeElement>();
}

/** Checks SearchWrapper works for OCTreeElement search nearest in radius
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperOCTreeElementSearchNearestInRadius, KratosCoreFastSuite)
{
    TestTreeSearchNearestInRadius<SearchWrapperOCTreeElement>();
}

/** Checks SearchWrapper works for StaticBinsTree search nearest in radius
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperStaticBinsTreeElementSearchNearestInRadius, KratosCoreFastSuite)
{
    TestTreeSearchNearestInRadius<SearchWrapperStaticBinsTreeElement>();
}

/** Checks SearchWrapper works for BinsDynamicElement search nearest in radius
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperBinsDynamicElementSearchNearestInRadius, KratosCoreFastSuite)
{
    TestTreeSearchNearestInRadius<SearchWrapperBinsDynamicElement>();
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
    DataCommunicator serial_communicator;
    TSearchWrapper search_wrapper(r_skin_part.Elements(), serial_communicator);

    double epsilon = 1.0e-6;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_near_point = r_point_model_part.CreateNewNode(1, epsilon,epsilon,epsilon);
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

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
}

/** Checks SearchWrapper works for KDTreeElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperKDTreeElementSearchNearest, KratosCoreFastSuite)
{
    TestTreeSearchNearest<SearchWrapperKDTreeElement>();
}

/** Checks SearchWrapper works for OCTreeElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperOCTreeElementSearchNearest, KratosCoreFastSuite)
{
    TestTreeSearchNearest<SearchWrapperOCTreeElement>();
}

/** Checks SearchWrapper works for StaticBinsTree search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperStaticBinsTreeElementSearchNearest, KratosCoreFastSuite)
{
    TestTreeSearchNearest<SearchWrapperStaticBinsTreeElement>();
}

/** Checks SearchWrapper works for BinsDynamicElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperBinsDynamicElementSearchNearest, KratosCoreFastSuite)
{
    TestTreeSearchNearest<SearchWrapperBinsDynamicElement>();
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
    DataCommunicator serial_communicator;
    TSearchWrapper search_wrapper(r_skin_part.Elements(), serial_communicator);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    auto p_point = r_point_model_part.CreateNewNode(1, 0.0,0.0,0.0);
    auto& r_array_nodes = r_point_model_part.Nodes();

    typename TSearchWrapper::ResultContainerVectorType results;
    search_wrapper.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[0].IsObjectFound());
}

/** Checks SearchWrapper works for KDTreeElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperKDTreeElementSearchNearestEmpty, KratosCoreFastSuite)
{
    TestTreeSearchNearestEmpty<SearchWrapperKDTreeElement>();
}

/** Checks SearchWrapper works for OCTreeElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperOCTreeElementSearchNearestEmpty, KratosCoreFastSuite)
{
    TestTreeSearchNearestEmpty<SearchWrapperOCTreeElement>();
}

/** Checks SearchWrapper works for StaticBinsTree search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperStaticBinsTreeElementSearchNearestEmpty, KratosCoreFastSuite)
{
    TestTreeSearchNearestEmpty<SearchWrapperStaticBinsTreeElement>();
}

/** Checks SearchWrapper works for BinsDynamicElement search nearest
*/
KRATOS_TEST_CASE_IN_SUITE(SearchWrapperBinsDynamicElementSearchNearestEmpty, KratosCoreFastSuite)
{
    TestTreeSearchNearestEmpty<SearchWrapperBinsDynamicElement>();
}

} // namespace Kratos::Testing.
