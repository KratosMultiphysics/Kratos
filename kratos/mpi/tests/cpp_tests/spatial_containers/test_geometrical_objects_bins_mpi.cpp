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
#include "mpi/spatial_containers/geometrical_objects_bins_mpi.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "containers/model.h"
#include "geometries/triangle_3d_3.h"

// TODO: Update serial test in order to test with the new containers (we can test both interfaces)
namespace Kratos::Testing {

ModelPart& CreateCubeSkinModelPart(
    Model& rCurrentModel,
    const double HalfX,
    const double HalfY,
    const double HalfZ
    )
{
    // Generate the cube skin
    ModelPart& r_skin_part = rCurrentModel.CreateModelPart("Skin");
    r_skin_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Create properties
    auto p_properties = r_skin_part.CreateNewProperties(1, 0);

    // Set partitions
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
    const int rank =  r_data_comm.Rank();
    const int world_size = r_data_comm.Size();
    if (world_size == 1) {
        r_skin_part.CreateNewNode(1, -HalfX, -HalfY, -HalfZ);
        r_skin_part.CreateNewNode(2,  HalfX, -HalfY, -HalfZ);
        r_skin_part.CreateNewNode(3,  HalfX,  HalfY, -HalfZ);
        r_skin_part.CreateNewNode(4, -HalfX,  HalfY, -HalfZ);
        r_skin_part.CreateNewNode(5, -HalfX, -HalfY,  HalfZ);
        r_skin_part.CreateNewNode(6,  HalfX, -HalfY,  HalfZ);
        r_skin_part.CreateNewNode(7,  HalfX,  HalfY,  HalfZ);
        r_skin_part.CreateNewNode(8, -HalfX,  HalfY,  HalfZ);
        for (auto& r_node : r_skin_part.Nodes()) {
            r_node.FastGetSolutionStepValue(PARTITION_INDEX) = 0;
        }
        // Create elements
        r_skin_part.CreateNewElement("Element3D3N",  1, { 1,2,3 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  2, { 1,3,4 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  3, { 5,6,7 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  4, { 5,7,8 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  5, { 3,6,2 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  6, { 3,7,6 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  7, { 4,5,1 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  8, { 4,8,5 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  9, { 3,4,8 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N", 10, { 3,8,7 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N", 11, { 2,1,5 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N", 12, { 2,5,6 }, p_properties);
    } else { // Assuming always two partitions
        if (rank == 0) {
            // Create nodes
            auto p_node1 = r_skin_part.CreateNewNode(1, -HalfX, -HalfY, -HalfZ);
            auto p_node2 = r_skin_part.CreateNewNode(2,  HalfX, -HalfY, -HalfZ);
            auto p_node3 = r_skin_part.CreateNewNode(3,  HalfX,  HalfY, -HalfZ);
            auto p_node4 = r_skin_part.CreateNewNode(4, -HalfX,  HalfY, -HalfZ);

            // Set partitions
            p_node1->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node2->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node3->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node4->FastGetSolutionStepValue(PARTITION_INDEX) = 0;

            // Create elements
            r_skin_part.CreateNewElement("Element3D3N",  1, { 1,2,3 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  2, { 1,3,4 }, p_properties);
        } else if (rank == 1) {
            // Create nodes
            auto p_node1 = r_skin_part.CreateNewNode(1, -HalfX, -HalfY, -HalfZ);
            auto p_node2 = r_skin_part.CreateNewNode(2,  HalfX, -HalfY, -HalfZ);
            auto p_node3 = r_skin_part.CreateNewNode(3,  HalfX,  HalfY, -HalfZ);
            auto p_node4 = r_skin_part.CreateNewNode(4, -HalfX,  HalfY, -HalfZ);
            auto p_node5 = r_skin_part.CreateNewNode(5, -HalfX, -HalfY,  HalfZ);
            auto p_node6 = r_skin_part.CreateNewNode(6,  HalfX, -HalfY,  HalfZ);
            auto p_node7 = r_skin_part.CreateNewNode(7,  HalfX,  HalfY,  HalfZ);
            auto p_node8 = r_skin_part.CreateNewNode(8, -HalfX,  HalfY,  HalfZ);

            // Set partitions
            p_node1->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node2->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node3->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node4->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node5->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
            p_node6->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
            p_node7->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
            p_node8->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

            // Create elements
            r_skin_part.CreateNewElement("Element3D3N",  3, { 5,6,7 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  4, { 5,7,8 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  5, { 3,6,2 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  6, { 3,7,6 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  7, { 4,5,1 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  8, { 4,8,5 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  9, { 3,4,8 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N", 10, { 3,8,7 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N", 11, { 2,1,5 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N", 12, { 2,5,6 }, p_properties);
        }
    }

    // Compute communicaton plan and fill communicator meshes correctly
    ParallelFillCommunicator(r_skin_part, r_data_comm).Execute();

    return r_skin_part;
}

/** Checks bins bounding box
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(GeometricalObjectsBinsMPIBoundingBox, KratosMPICoreFastSuite)
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    const double cube_x = 0.6;
    const double cube_y = 0.9;
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    GeometricalObjectsBinsMPI bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd(), r_data_comm);

    auto bounding_box = bins.GetBoundingBox();

    KRATOS_CHECK_NEAR(bounding_box.GetMinPoint()[0],-cube_x, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMinPoint()[1],-cube_y, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMinPoint()[2],-cube_z, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMaxPoint()[0], cube_x, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMaxPoint()[1], cube_y, tolerance);
    KRATOS_CHECK_NEAR(bounding_box.GetMaxPoint()[2], cube_z, tolerance);
}

/** Checks bins search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(GeometricalObjectsBinsMPISearchInRadius, KratosMPICoreFastSuite)
{
    Model current_model;

    const double cube_x = 0.6;
    const double cube_y = 0.9;
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    GeometricalObjectsBinsMPI bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd(), r_data_comm);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    Point point(0.0, 0.0, 0.0);
    if (rank == 0) {
        r_point_model_part.CreateNewNode(1, point.X(), point.Y(), point.Z());
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    GeometricalObjectsBinsMPI::ResultTypeContainerMap results;

    // 0.29 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.29, results);
    // Like the point is on the first rank the results should be empty in all partitions except the first one
    KRATOS_CHECK_EQUAL(results.NumberOfPointsResults(), 1);
    KRATOS_CHECK_EQUAL(results[point].NumberOfGlobalResults(), 0);

    // 0.3 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    // Like the point is on the first rank the results should be empty in all partitions except the first one
    KRATOS_CHECK_EQUAL(results.NumberOfPointsResults(), 1);
    KRATOS_CHECK_EQUAL(results[point].NumberOfGlobalResults(), 4);

    // 0.4 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.4, results);
    // Like the point is on the first rank the results should be empty in all partitions except the first one
    KRATOS_CHECK_EQUAL(results.NumberOfPointsResults(), 1);
    KRATOS_CHECK_EQUAL(results[point].NumberOfGlobalResults(), 4);

    // 0.6 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);
    // Like the point is on the first rank the results should be empty in all partitions except the first one
    KRATOS_CHECK_EQUAL(results.NumberOfPointsResults(), 1);
    KRATOS_CHECK_EQUAL(results[point].NumberOfGlobalResults(), 8);

    // 0.7 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.7, results);
    // Like the point is on the first rank the results should be empty in all partitions except the first one
    KRATOS_CHECK_EQUAL(results.NumberOfPointsResults(), 1);
    KRATOS_CHECK_EQUAL(results[point].NumberOfGlobalResults(), 8);

    // 0.9 radius
    bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.9, results);
    // Like the point is on the first rank the results should be empty in all partitions except the first one
    KRATOS_CHECK_EQUAL(results.NumberOfPointsResults(), 1);
    KRATOS_CHECK_EQUAL(results[point].NumberOfGlobalResults(), 12);
}

/** Checks bins search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(GeometricalObjectsBinsMPISearchNearestInRadius, KratosMPICoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    const double cube_x = 0.6;
    const double cube_y = 0.9;
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    GeometricalObjectsBinsMPI bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd(), r_data_comm);

    double epsilon = 1.0e-6;
    Point near_point{epsilon,epsilon,epsilon};

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        r_point_model_part.CreateNewNode(1, epsilon,epsilon,epsilon);
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    GeometricalObjectsBinsMPI::ResultTypeContainerMap results;
    bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z - 1.e-4, results);

    KRATOS_CHECK_EQUAL(results.NumberOfPointsResults(), 1);
    KRATOS_CHECK_EQUAL(results[near_point].NumberOfGlobalResults(), 0);

    bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z + 1.e-4, results);

    KRATOS_CHECK_EQUAL(results.NumberOfPointsResults(), 1);
    KRATOS_CHECK_EQUAL(results[near_point].NumberOfGlobalResults(), 1);

    // Distances are just local
    if (rank == 1) {
        auto& r_distances = results[near_point].GetLocalDistances();
        KRATOS_CHECK_NEAR(r_distances.begin()->second, (cube_z - epsilon), tolerance);
    }

    // Compute indices
    auto indices = results[near_point].GetResultIndices();
    const std::size_t id = indices[0];
    KRATOS_CHECK_EQUAL(id, 3);
}

// /** Checks bins search nearest
// */
// KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(GeometricalObjectsBinsMPISearchNearest, KratosMPICoreFastSuite) 
// {
//     constexpr double tolerance = 1e-12;

//     Model current_model;

//     const double cube_x = 0.6;
//     const double cube_y = 0.9;
//     const double cube_z = 0.3;

//     // Generate the cube skin
//     ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, cube_x, cube_y, cube_z);
//     const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

//     GeometricalObjectsBinsMPI bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd(), r_data_comm);

//     double epsilon = 1.0e-6;
//     Point near_point{epsilon,epsilon,epsilon};
//     auto result = bins.SearchNearest(near_point);

//     KRATOS_CHECK_NEAR(result.GetDistance(), (cube_z - epsilon), tolerance);

//     std::size_t id = result.Get()->Id();
//     KRATOS_CHECK(id == 3); 
// }

// /** Checks bins empty search nearest 
// */
// KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(GeometricalObjectsBinsMPIEmptySearchNearest, KratosMPICoreFastSuite) 
// {
//     Model current_model;

//     // Generate the cube skin
//     ModelPart& r_skin_part = current_model.CreateModelPart("Skin");
//     const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

//     GeometricalObjectsBinsMPI bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd(), r_data_comm);

//     Point center_point{0.0,0.0,0.0};
//     auto result = bins.SearchNearest(center_point);

//     KRATOS_CHECK_IS_FALSE(result.IsObjectFound());
// }

// /** Checks bins search is inside 
// */
// KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(GeometricalObjectsBinsMPISearchIsInside, KratosMPICoreFastSuite) 
// {
//     constexpr double tolerance = 1e-12;

//     Model current_model;

//     // Generate the cube skin
//     ModelPart& r_skin_part = current_model.CreateModelPart("Skin");
//     const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
//     r_skin_part.CreateNewNode(1, 0.0, 0.0, 0.0);
//     r_skin_part.CreateNewNode(2, 1.0, 0.0, 0.0);
//     r_skin_part.CreateNewNode(3, 0.0, 1.0, 0.0);
//     r_skin_part.CreateNewNode(4, 0.0, 0.0, 1.0);
//     Properties::Pointer p_properties(new Properties(0));
//     r_skin_part.CreateNewElement("Element3D4N",  1, { 1,2,3,4 }, p_properties);

//     GeometricalObjectsBinsMPI bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd(), r_data_comm);

//     Point near_point{0.1,0.1,0.1};
//     auto result = bins.SearchIsInside(near_point);

//     KRATOS_CHECK(result.IsObjectFound());
//     KRATOS_CHECK_NEAR(result.GetDistance(), 0.0, tolerance);
// }

// /** Checks bins search is inside = not found
// */
// KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(GeometricalObjectsBinsMPISearchIsNotInside, KratosMPICoreFastSuite) 
// {
//     Model current_model;

//     // Generate the cube skin
//     ModelPart& r_skin_part = current_model.CreateModelPart("Skin");
//     const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();
//     r_skin_part.CreateNewNode(1, 0.0, 0.0, 0.0);
//     r_skin_part.CreateNewNode(2, 1.0, 0.0, 0.0);
//     r_skin_part.CreateNewNode(3, 0.0, 1.0, 0.0);
//     r_skin_part.CreateNewNode(4, 0.0, 0.0, 1.0);
//     Properties::Pointer p_properties(new Properties(0));
//     r_skin_part.CreateNewElement("Element3D4N",  1, { 1,2,3,4 }, p_properties);

//     GeometricalObjectsBinsMPI bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd(), r_data_comm);

//     Point near_point{0.5,0.5,0.5};
//     auto result = bins.SearchIsInside(near_point);

//     KRATOS_CHECK_IS_FALSE(result.IsObjectFound());
// }

} // namespace Kratos::Testing.
