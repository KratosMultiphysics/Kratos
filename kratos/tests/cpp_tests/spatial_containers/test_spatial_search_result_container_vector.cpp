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
#include "containers/model.h"
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
    container_vector.InitializeResult();

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
    container_vector.InitializeResults(indexes.size());

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
    container_vector.InitializeResult();

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
    container_vector.InitializeResult();

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
    container_vector.InitializeResults(indexes.size());

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
    container_vector.InitializeResults(indexes.size());

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
    std::vector<std::vector<double>> distances;
    container_vector.GetDistances(distances);
    KRATOS_EXPECT_EQ(distances.size(), 2);
    KRATOS_EXPECT_EQ(distances[0].size(), 2);
    KRATOS_EXPECT_DOUBLE_EQ(distances[0][0], 0.5);
    KRATOS_EXPECT_DOUBLE_EQ(distances[0][1], 0.25);
    KRATOS_EXPECT_EQ(distances[1].size(), 1);
    KRATOS_EXPECT_DOUBLE_EQ(distances[1][0], 0.5);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetResultIsLocal, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0,1};
    container_vector.InitializeResults(indexes.size());

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
    std::vector<std::vector<bool>> is_local;
    container_vector.GetResultIsLocal(is_local, data_communicator);
    KRATOS_EXPECT_EQ(is_local.size(), 2);
    KRATOS_EXPECT_EQ(is_local[0].size(), 2);
    KRATOS_EXPECT_TRUE(is_local[0][0]);
    KRATOS_EXPECT_TRUE(is_local[0][1]);
    KRATOS_EXPECT_EQ(is_local[1].size(), 1);
    KRATOS_EXPECT_TRUE(is_local[1][0]);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetResultIsActive, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0};
    container_vector.InitializeResults(indexes.size());

    // Container 1
    auto& r_container_1 = container_vector[0];

    // Create a test result
    GeometricalObject object = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    r_container_1.AddResult(result);

    // SynchronizeAll
    container_vector.SynchronizeAll(data_communicator);

    // Compute is active
    std::vector<std::vector<bool>> is_active;
    container_vector.GetResultIsActive(is_active, data_communicator);

    // Check is active
    KRATOS_EXPECT_EQ(is_active.size(), 1);
    KRATOS_EXPECT_TRUE(is_active[0][0]);

    // Deactivate the object
    object.Set(ACTIVE, false);

    // Compute is active
    container_vector.GetResultIsActive(is_active, data_communicator);

    // Check is active
    KRATOS_EXPECT_EQ(is_active.size(), 1);
    KRATOS_EXPECT_FALSE(is_active[0][0]);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetResultShapeFunctions, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0};
    container_vector.InitializeResults(indexes.size());

    // Container 1
    auto& r_container_1 = container_vector[0];

    // Generate a geometry
    auto p_node1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_node2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Geometry<Node>::Pointer p_geom = Kratos::make_shared<Line2D2<Node>>(p_node1, p_node2);

    // Create a test result
    GeometricalObject object = GeometricalObject(1, p_geom);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    r_container_1.AddResult(result);
    r_container_1.SetLocalIndex(0);

    // Synchronize the container between partitions
    container_vector.SynchronizeAll(data_communicator);

    // Define the model
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    auto p_test_node = r_model_part.CreateNewNode(1, 0.5, 0.0, 0.0);

    // Compute shape functions
    std::vector<std::vector<Vector>> shape_functions;
    container_vector.GetResultShapeFunctions(shape_functions, r_model_part.Nodes(), data_communicator);

    // Check shape functions
    KRATOS_EXPECT_EQ(shape_functions[0].size(), 1);
    KRATOS_EXPECT_NEAR(shape_functions[0][0][0], 0.5, 1.0e-12);
    KRATOS_EXPECT_NEAR(shape_functions[0][0][1], 0.5, 1.0e-12);

    // Check is inside
    std::vector<std::vector<bool>> is_inside_true;
    container_vector.GetResultIsInside(is_inside_true, r_model_part.Nodes(), data_communicator, 1.0e-5);
    KRATOS_EXPECT_EQ(is_inside_true[0].size(), 1);
    KRATOS_EXPECT_TRUE(is_inside_true[0][0]);

    // Check is outside
    p_test_node->X() = 1.0e6;
    p_test_node->Y() = 1.0e6;
    p_test_node->Z() = 1.0e6;
    std::vector<std::vector<bool>> is_inside_false;
    container_vector.GetResultIsInside(is_inside_false, r_model_part.Nodes(), data_communicator, 1.0e-5);
    KRATOS_EXPECT_EQ(is_inside_false[0].size(), 1);
    KRATOS_EXPECT_FALSE(is_inside_false[0][0]);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetResultIndices, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0};
    container_vector.InitializeResults(indexes.size());

    // Container
    auto& r_container = container_vector[0];

    // Create a test result
    GeometricalObject object = GeometricalObject(1);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    r_container.AddResult(result);

    // SynchronizeAll
    container_vector.SynchronizeAll(data_communicator);

    // Compute indices
    std::vector<std::vector<std::size_t>> indices;
    container_vector.GetResultIndices(indices);

    // Check indices
    KRATOS_EXPECT_EQ(indices[0].size(), 1);
    KRATOS_EXPECT_EQ(indices[0][0], object.Id());
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetResultNodeIndices, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0};
    container_vector.InitializeResults(indexes.size());

    // Container
    auto& r_container = container_vector[0];

    // Generate a geometry
    auto p_node1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_node2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Geometry<Node>::Pointer p_geom = Kratos::make_shared<Line2D2<Node>>(p_node1, p_node2);

    // Create a test result
    GeometricalObject object = GeometricalObject(1, p_geom);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    r_container.AddResult(result);

    // SynchronizeAll
    container_vector.SynchronizeAll(data_communicator);

    // Compute indices
    std::vector<std::vector<std::vector<IndexType>>> indices;
    container_vector.GetResultNodeIndices(indices);

    // Check indices
    KRATOS_EXPECT_EQ(indices[0].size(), 1);
    KRATOS_EXPECT_EQ(indices[0][0][0], 1);
    KRATOS_EXPECT_EQ(indices[0][0][1], 2);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetResultCoordinates, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0};
    container_vector.InitializeResults(indexes.size());

    // Container
    auto& r_container = container_vector[0];

    // Generate a geometry
    auto p_node1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_node2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Geometry<Node>::Pointer p_geom = Kratos::make_shared<Line2D2<Node>>(p_node1, p_node2);

    // Create a test result
    GeometricalObject object = GeometricalObject(1, p_geom);
    SpatialSearchResult<GeometricalObject> result(&object);

    // Add the result to the container
    r_container.AddResult(result);

    // SynchronizeAll
    container_vector.SynchronizeAll(data_communicator);

    // Compute result coordinates
    std::vector<std::vector<std::vector<array_1d<double, 3>>>> coordinates;
    container_vector.GetResultCoordinates(coordinates);

    // Check result coordinates
    KRATOS_EXPECT_EQ(coordinates[0].size(), 1);
    KRATOS_EXPECT_EQ(coordinates[0][0].size(), 2);
    KRATOS_EXPECT_VECTOR_NEAR(coordinates[0][0][0], p_node1->Coordinates(), 1.0e-12);
    KRATOS_EXPECT_VECTOR_NEAR(coordinates[0][0][1], p_node2->Coordinates(), 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSearchResultContainerVectorGetResultRank, KratosCoreFastSuite)
{
    // Create a test object
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Initialize result
    DataCommunicator data_communicator;
    const std::vector<std::size_t> indexes{0,1};
    container_vector.InitializeResults(indexes.size());

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
    std::vector<std::vector<int>> rank;
    container_vector.GetResultRank(rank);
    KRATOS_EXPECT_EQ(rank.size(), 2);
    KRATOS_EXPECT_EQ(rank[0].size(), 2);
    KRATOS_EXPECT_EQ(rank[0][0], 0);
    KRATOS_EXPECT_EQ(rank[0][1], 0);
    KRATOS_EXPECT_EQ(rank[1].size(), 1);
    KRATOS_EXPECT_EQ(rank[1][0], 0);
}

}  // namespace Kratos::Testing