//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "mapping_application_variables.h"
#include "custom_utilities/mapper_utilities.h"
#include "custom_mappers/nearest_neighbor_mapper.h"

namespace Kratos {
namespace Testing {

void CreateNodesForMapping(ModelPart& rModelPart, const int NumNodes)
{
    const int rank = rModelPart.GetCommunicator().MyPID();
    const int size = rModelPart.GetCommunicator().TotalProcesses();

    const int start_id = NumNodes * rank + 1;

    // creating nodes with random coordinates
    for (int i=0; i< NumNodes; ++i)
        rModelPart.CreateNewNode(i+start_id, i*0.1*rank*size+0.134,
                                             i*0.2+rank*3.48*size,
                                             i*0.3*rank*6.13*size);
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_AssignInterfaceEquationIds, KratosMappingApplicationSerialTestSuite)
{
    const int num_nodes = 11;
    ModelPart model_part("ForTest");

    CreateNodesForMapping(model_part, num_nodes);

    MapperUtilities::AssignInterfaceEquationIds(model_part.GetCommunicator());

    int idx = 0;

    for (const auto& r_node : model_part/*.GetCommunicator().LocalMesh()*/.Nodes())
    {
        KRATOS_CHECK_EQUAL(idx, r_node.GetValue(INTERFACE_EQUATION_ID));
        idx += 1;
    }
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_ComputeBoundingBox, KratosMappingApplicationSerialTestSuite)
{
    ModelPart model_part("ForTest");
    model_part.CreateNewNode(1, 0.2, 5.3, -8.3);
    model_part.CreateNewNode(2, 8.2, 25.3, 16.4);
    model_part.CreateNewNode(3, -9.2, -17.13, 1.5);
    model_part.CreateNewNode(4, 12.6, 5.3, -8.3);

    const auto bbox = MapperUtilities::ComputeLocalBoundingBox(model_part);

    // std::cout << MapperUtilities::BoundingBoxStringStream(bbox) << std::endl;

    KRATOS_CHECK_EQUAL(bbox.size(), 6);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[0], 12.6);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[1], -9.2);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[2], 25.3);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[3], -17.13);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[4], 16.4);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[5], -8.3);
}

double GetBBoxValue(const int Index, const double Factor, const double Offset)
{
    return static_cast<double>(Index)*Factor - Offset;
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_ComputeBoundingBoxWithTol, KratosMappingApplicationSerialTestSuite)
{
    std::vector<double> bboxes_wrong_size(5);
    std::vector<double> bboxes_with_tol;

    KRATOS_CHECK_EXCEPTION_IS_THROWN(MapperUtilities::ComputeBoundingBoxesWithTolerance(bboxes_wrong_size, 1.235, bboxes_with_tol),
        "Error: Bounding Boxes size has to be a multiple of 6!");

    // Cretae a vector containing the fake bboxes
    const int num_entries = 24;
    std::vector<double> bboxes(num_entries);

    const double factor = 1.2589;
    const double offset = 8.4;

    for (int i=0; i<num_entries; ++i)
        bboxes[i] = GetBBoxValue(i, factor, offset);

    const double tolerance = 5.478;

    MapperUtilities::ComputeBoundingBoxesWithTolerance(bboxes,
                                                       tolerance,
                                                       bboxes_with_tol);

    for (int i=0; i<num_entries; i+=2)
        KRATOS_CHECK_NEAR(bboxes_with_tol[i], (GetBBoxValue(i, factor, offset) + tolerance), 1e-12);

    for (int i=1; i<num_entries; i+=2)
        KRATOS_CHECK_NEAR(bboxes_with_tol[i], (GetBBoxValue(i, factor, offset) - tolerance), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_PointIsInsideBoundingBox, KratosMappingApplicationSerialTestSuite)
{
    const std::vector<double> bounding_box {10.5, -2.8, 3.89, -77.6, 4.64, 2.3};
    // xmax, xmin,  ymax, ymin,  zmax, zmin

    const Point p_out_x(10.6, 1.0, 3.8);
    const Point p_out_y(10.1, -80.0, 3.8);
    const Point p_out_z(10.1, 1.0, -3.8);
    const Point p_in(10.0, -30.78, 3.7);

    KRATOS_CHECK_IS_FALSE(MapperUtilities::PointIsInsideBoundingBox(bounding_box, p_out_x));
    KRATOS_CHECK_IS_FALSE(MapperUtilities::PointIsInsideBoundingBox(bounding_box, p_out_y));
    KRATOS_CHECK_IS_FALSE(MapperUtilities::PointIsInsideBoundingBox(bounding_box, p_out_z));

    KRATOS_CHECK(MapperUtilities::PointIsInsideBoundingBox(bounding_box, p_in));
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_MapperInterfaceInfoSerializer, KratosMappingApplicationSerialTestSuite)
{
    // this test checks if the serialization/deserialization of the Helper-Class that
    // serializes/deserializes the MapperInterfaceInfos works correctly
    // => this is needed to transfer the data btw the ranks in MPI
    // The test is rather large, this is needed to cover a complete example
    // Note that the same checks are performed before and after ther serialization
    // to make sure that the objects are properly initialized

    using MapperInterfaceInfoUniquePointerType = Kratos::unique_ptr<MapperInterfaceInfo>;

    using MapperInterfaceInfoPointerType = Kratos::shared_ptr<MapperInterfaceInfo>;
    using MapperInterfaceInfoPointerVectorType = std::vector<std::vector<MapperInterfaceInfoPointerType>>;
    using MapperInterfaceInfoPointerVectorPointerType = Kratos::unique_ptr<MapperInterfaceInfoPointerVectorType>;


    // A "NearestNeighborInterfaceInfo" is being used since "MapperInterfaceInfo" is a pure virtual class
    Point coords_1(1.0, 2.45, 33.8);
    Point coords_2(10.0, 20.45, 100.0);
    Point coords_3(2.0, 2.45, -2.38);

    std::size_t source_local_sys_idx_1 = 123;
    std::size_t source_local_sys_idx_2 = 1235214;
    std::size_t source_local_sys_idx_3 = 8;

    MapperInterfaceInfoPointerType p_nearest_neighbor_info_1(
        Kratos::make_shared<NearestNeighborInterfaceInfo>(coords_1, source_local_sys_idx_1, 0));
    MapperInterfaceInfoPointerType p_nearest_neighbor_info_2(
        Kratos::make_shared<NearestNeighborInterfaceInfo>(coords_2, source_local_sys_idx_2, 0));
    MapperInterfaceInfoPointerType p_nearest_neighbor_info_3(
        Kratos::make_shared<NearestNeighborInterfaceInfo>(coords_3, source_local_sys_idx_3, 0));

    // Auxiliary objects to fill the NearestNeighborInterfaceInfos with values that can be checked afterwards
    auto node_1(Kratos::make_shared<Node<3>>(1, 1.0, 2.5, 30.0));
    auto node_2(Kratos::make_shared<Node<3>>(3, 10.5, 20.0, 96.8));
    auto node_3(Kratos::make_shared<Node<3>>(15, 2.3, 1.9, -2.5));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3));

    const int expected_id_found_1 = 108;
    const int expected_id_found_2 = 18;
    const int expected_id_found_3 = 896;

    node_1->SetValue(INTERFACE_EQUATION_ID, expected_id_found_1);
    node_2->SetValue(INTERFACE_EQUATION_ID, expected_id_found_2);
    node_3->SetValue(INTERFACE_EQUATION_ID, expected_id_found_3);

    // We compute the real distance bcs this would also be computed by the search
    const double dist_1_1 = MapperUtilities::ComputeDistance(coords_1, *interface_node_1);
    const double dist_2_1 = MapperUtilities::ComputeDistance(coords_1, *interface_node_2);
    const double dist_3_1 = MapperUtilities::ComputeDistance(coords_1, *interface_node_3);

    p_nearest_neighbor_info_1->ProcessSearchResult(interface_node_1, dist_1_1);
    p_nearest_neighbor_info_1->ProcessSearchResult(interface_node_2, dist_2_1);
    p_nearest_neighbor_info_1->ProcessSearchResult(interface_node_3, dist_3_1);

    // Now some the checks are performed to make sure the objects are correctly initialized
    int found_id;
    p_nearest_neighbor_info_1->GetValue(found_id);
    KRATOS_CHECK_EQUAL(found_id, expected_id_found_1);
    double neighbor_dist;
    p_nearest_neighbor_info_1->GetValue(neighbor_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_1_1);

    const double dist_1_2 = MapperUtilities::ComputeDistance(coords_2, *interface_node_1);
    const double dist_2_2 = MapperUtilities::ComputeDistance(coords_2, *interface_node_2);
    const double dist_3_2 = MapperUtilities::ComputeDistance(coords_2, *interface_node_3);

    p_nearest_neighbor_info_2->ProcessSearchResult(interface_node_1, dist_1_2);
    p_nearest_neighbor_info_2->ProcessSearchResult(interface_node_2, dist_2_2);
    p_nearest_neighbor_info_2->ProcessSearchResult(interface_node_3, dist_3_2);

    // Now some the checks are performed to make sure the objects are correctly initialized
    p_nearest_neighbor_info_2->GetValue(found_id);
    KRATOS_CHECK_EQUAL(found_id, expected_id_found_2);
    p_nearest_neighbor_info_2->GetValue(neighbor_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_2_2);

    const double dist_1_3 = MapperUtilities::ComputeDistance(coords_3, *interface_node_1);
    const double dist_2_3 = MapperUtilities::ComputeDistance(coords_3, *interface_node_2);
    const double dist_3_3 = MapperUtilities::ComputeDistance(coords_3, *interface_node_3);

    p_nearest_neighbor_info_3->ProcessSearchResult(interface_node_1, dist_1_3);
    p_nearest_neighbor_info_3->ProcessSearchResult(interface_node_2, dist_2_3);
    p_nearest_neighbor_info_3->ProcessSearchResult(interface_node_3, dist_3_3);

    // Now some the checks are performed to make sure the objects are correctly initialized
    p_nearest_neighbor_info_3->GetValue(found_id);
    KRATOS_CHECK_EQUAL(found_id, expected_id_found_3);
    p_nearest_neighbor_info_3->GetValue(neighbor_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_3_3);

    // Now finally we can construct the container
    MapperInterfaceInfoPointerVectorPointerType p_interface_info_container
        = Kratos::make_unique<MapperInterfaceInfoPointerVectorType>();

    p_interface_info_container->resize(2);

    // Note: the order is choosen intentionally
    (*p_interface_info_container)[0].push_back(p_nearest_neighbor_info_3);
    (*p_interface_info_container)[1].push_back(p_nearest_neighbor_info_1);
    (*p_interface_info_container)[1].push_back(p_nearest_neighbor_info_2);

    // Construct a reference obj, needed to create the correct objects while loading/deserializing
    auto p_ref_nearest_neighbor_info = p_nearest_neighbor_info_1->Create();

    MapperUtilities::MapperInterfaceInfoSerializer serializer_helper(
        p_interface_info_container, p_ref_nearest_neighbor_info );

    // serializing the object (=> happens on the partition that sends the objects)
    Serializer serializer;
    serializer.save("obj", serializer_helper);

    MapperInterfaceInfoPointerVectorPointerType p_interface_info_container_new
        = Kratos::make_unique<MapperInterfaceInfoPointerVectorType>();

    MapperUtilities::MapperInterfaceInfoSerializer serializer_helper_new(
        p_interface_info_container_new, p_ref_nearest_neighbor_info );

    // deserializing the object (=> happens on the partition that receives the objects)
    serializer.load("obj", serializer_helper_new);

    // Checking for the sizes of the container
    KRATOS_CHECK_EQUAL(p_interface_info_container_new->size(), 2);
    KRATOS_CHECK_EQUAL((*p_interface_info_container_new)[0].size(), 1);
    KRATOS_CHECK_EQUAL((*p_interface_info_container_new)[1].size(), 2);

    // Checking the objects inside the container
    const auto& r_info_1 = (*p_interface_info_container_new)[1][0];
    const auto& r_info_2 = (*p_interface_info_container_new)[1][1];
    const auto& r_info_3 = (*p_interface_info_container_new)[0][0];

    r_info_1->GetValue(found_id);
    KRATOS_CHECK_EQUAL(found_id, expected_id_found_1);
    r_info_2->GetValue(found_id);
    KRATOS_CHECK_EQUAL(found_id, expected_id_found_2);
    r_info_3->GetValue(found_id);
    KRATOS_CHECK_EQUAL(found_id, expected_id_found_3);

    r_info_1->GetValue(neighbor_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_1_1);
    r_info_2->GetValue(neighbor_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_2_2);
    r_info_3->GetValue(neighbor_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_3_3);
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_FillBufferBeforeLocalSearch, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_CreateMapperInterfaceInfosFromBuffer, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_SelectInterfaceInfosSuccessfulSearch, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_AssignInterfaceInfosAfterRemoteSearch, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

}  // namespace Testing
}  // namespace Kratos