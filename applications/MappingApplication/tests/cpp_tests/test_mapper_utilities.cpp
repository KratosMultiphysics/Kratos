//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/stream_serializer.h"
#include "tests/test_utilities/cpp_tests_utilities.h"
#include "utilities/variable_utils.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "mapping_application_variables.h"
#include "custom_utilities/mapper_utilities.h"
#include "custom_mappers/nearest_neighbor_mapper.h"

namespace Kratos::Testing {

typedef std::size_t IndexType;
typedef std::size_t SizeType;

typedef Kratos::unique_ptr<MapperInterfaceInfo> MapperInterfaceInfoUniquePointerType;

typedef Kratos::shared_ptr<MapperInterfaceInfo> MapperInterfaceInfoPointerType;
typedef std::vector<std::vector<MapperInterfaceInfoPointerType>> MapperInterfaceInfoPointerVectorType;

namespace {

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

} //empty namespace

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_AssignInterfaceEquationIds, KratosMappingApplicationSerialTestSuite)
{
    const int num_nodes = 11;

    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");

    CreateNodesForMapping(model_part, num_nodes);

    MapperUtilities::AssignInterfaceEquationIds(model_part.GetCommunicator());

    int idx = 0;

    for (const auto& r_node : model_part/*.GetCommunicator().LocalMesh()*/.Nodes())
    {
        KRATOS_EXPECT_EQ(idx, r_node.GetValue(INTERFACE_EQUATION_ID));
        idx += 1;
    }
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_ComputeLocalBoundingBox, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");

    model_part.CreateNewNode(1, 0.2, 5.3, -8.3);
    model_part.CreateNewNode(2, 8.2, 25.3, 16.4);
    model_part.CreateNewNode(3, -9.2, -17.13, 1.5);
    model_part.CreateNewNode(4, 12.6, 5.3, -8.3);

    const auto bbox = MapperUtilities::ComputeLocalBoundingBox(model_part);

    // std::cout << MapperUtilities::BoundingBoxStringStream(bbox) << std::endl;

    KRATOS_EXPECT_EQ(bbox.size(), 6);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[0], 12.6);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[1], -9.2);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[2], 25.3);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[3], -17.13);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[4], 16.4);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[5], -8.3);
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_ComputeGlobalBoundingBox, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");

    model_part.CreateNewNode(1, 0.2, 5.3, -8.3);
    model_part.CreateNewNode(2, 8.2, 25.3, 16.4);
    model_part.CreateNewNode(3, -9.2, -17.13, 1.5);
    model_part.CreateNewNode(4, 12.6, 5.3, -8.3);

    const auto bbox = MapperUtilities::ComputeGlobalBoundingBox(model_part);

    // std::cout << MapperUtilities::BoundingBoxStringStream(bbox) << std::endl;

    KRATOS_EXPECT_EQ(bbox.size(), 6);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[0], 12.6);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[1], -9.2);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[2], 25.3);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[3], -17.13);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[4], 16.4);
    KRATOS_EXPECT_DOUBLE_EQ(bbox[5], -8.3);
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_FillBufferBeforeLocalSearch, KratosMappingApplicationSerialTestSuite)
{
    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;

    MapperLocalSystemPointerVector local_systems;

    const SizeType buffer_size_estimate = 3;

    const int comm_size = 3;

    std::vector<std::vector<double>> send_buffer(comm_size);
    std::vector<int> send_sizes(comm_size);

    // each row is one bbox
    const std::vector<double> bounding_boxes {10.5, -2.8, 3.89, -77.6, 4.64, 2.3,
                                              -1.5, -10.3, 20.6, 3.4, 3.77, -20.8,
                                              25.998, 6.4, 50.6, 15.2, 10.88, 4.12};

    KRATOS_EXPECT_EQ(bounding_boxes.size(), (6*comm_size)); // ensure the test is set up correctly

    // xmax, xmin,  ymax, ymin,  zmax, zmin
    const std::vector<double> missized_bounding_boxes {10.5, -2.8, 3.89};

    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(MapperUtilities::FillBufferBeforeLocalSearch(
        local_systems,
        missized_bounding_boxes,
        buffer_size_estimate,
        send_buffer,
        send_sizes),
        "Error: Bounding Boxes size has to be a multiple of 6!");

    // Node-ids do not matter here
    auto node_local_sys_1(Kratos::make_shared<Node>(87, -2.0, 3.5, 3.0)); // in bbox 1&2
    auto node_local_sys_2(Kratos::make_shared<Node>(26, 10.0, -25.0, 3.0)); // in bbox 1
    auto node_local_sys_3(Kratos::make_shared<Node>(36, -10.0, 15.5, -5.0)); // in bbox 2
    auto node_local_sys_4(Kratos::make_shared<Node>(46, 12.6, 50.1, 5.0)); // in bbox 3

    local_systems.reserve(7);

    local_systems.push_back(Kratos::make_unique<NearestNeighborLocalSystem>(node_local_sys_1.get()));
    local_systems.push_back(Kratos::make_unique<NearestNeighborLocalSystem>(node_local_sys_2.get()));
    local_systems.push_back(Kratos::make_unique<NearestNeighborLocalSystem>(node_local_sys_3.get()));
    local_systems.push_back(Kratos::make_unique<NearestNeighborLocalSystem>(node_local_sys_4.get()));

    MapperUtilities::FillBufferBeforeLocalSearch(
        local_systems,
        bounding_boxes,
        buffer_size_estimate,
        send_buffer,
        send_sizes);

    // preforming checks on the buffer
    KRATOS_EXPECT_EQ(send_buffer.size(), comm_size); // equal to comm_size, should not have changed
    KRATOS_EXPECT_EQ(send_sizes.size(), comm_size); // equal to comm_size, should not have changed

    KRATOS_EXPECT_EQ(send_sizes[0], 8); // two objs with each 4 doubles are sent to this bbox
    KRATOS_EXPECT_EQ(send_sizes[1], 8); // two objs with each 4 doubles are sent to this bbox
    KRATOS_EXPECT_EQ(send_sizes[2], 4); // one obj with 4 doubles are sent to this bbox

    KRATOS_EXPECT_EQ(send_buffer[0].size(), 8); // two objs with each 4 doubles are sent to this bbox
    KRATOS_EXPECT_EQ(send_buffer[1].size(), 8); // two objs with each 4 doubles are sent to this bbox
    KRATOS_EXPECT_EQ(send_buffer[2].size(), 4); // one obj with 4 doubles are sent to this bbox

    // set up expected send buffer
    std::vector<double> exp_send_buffer_rank_1 {0.0, -2.0, 3.5, 3.0,
                                                1.0, 10.0, -25.0, 3.0};
    std::vector<double> exp_send_buffer_rank_2 {0.0, -2.0, 3.5, 3.0,
                                                2.0, -10.0, 15.5, -5.0};
    std::vector<double> exp_send_buffer_rank_3 {3.0, 12.6, 50.1, 5.0};

    std::vector<std::vector<double>> exp_send_buffer {exp_send_buffer_rank_1,
                                                      exp_send_buffer_rank_2,
                                                      exp_send_buffer_rank_3};

    // compare the send buffers
    for (IndexType i=0; i<exp_send_buffer.size(); ++i)
        for (IndexType j=0; j<exp_send_buffer[i].size(); ++j)
            KRATOS_EXPECT_DOUBLE_EQ(send_buffer[i][j], exp_send_buffer[i][j]);

    /////
    // now we "update" the Interface and then check the buffers again
    /////

    Point new_coords_node_1(1008.3, -89.123, 125.7);
    noalias(node_local_sys_1->Coordinates()) = new_coords_node_1;
    // => now "node_local_sys_1" will not fall into any bbox any more

    auto node_local_sys_5(Kratos::make_shared<Node>(50, -8.301, 17.75, -15.18)); // in bbox 2
    auto node_local_sys_6(Kratos::make_shared<Node>(416, 13.5, 44.58, 7.5)); // in bbox 3
    auto node_local_sys_7(Kratos::make_shared<Node>(417, 13.5125, 44.68, 8.5)); // in bbox 3

    local_systems.push_back(Kratos::make_unique<NearestNeighborLocalSystem>(node_local_sys_5.get()));
    local_systems.push_back(Kratos::make_unique<NearestNeighborLocalSystem>(node_local_sys_6.get()));
    local_systems.push_back(Kratos::make_unique<NearestNeighborLocalSystem>(node_local_sys_7.get()));

    MapperUtilities::FillBufferBeforeLocalSearch(
        local_systems,
        bounding_boxes,
        buffer_size_estimate,
        send_buffer,
        send_sizes);

    // preforming checks on the buffer
    KRATOS_EXPECT_EQ(send_buffer.size(), comm_size); // equal to comm_size, should not have changed
    KRATOS_EXPECT_EQ(send_sizes.size(), comm_size); // equal to comm_size, should not have changed

    KRATOS_EXPECT_EQ(send_sizes[0], 4); // one obj with 4 doubles are sent to this bbox
    KRATOS_EXPECT_EQ(send_sizes[1], 8); // two objs with each 4 doubles are sent to this bbox
    KRATOS_EXPECT_EQ(send_sizes[2], 12); // three objs with each 4 doubles are sent to this bbox

    KRATOS_EXPECT_EQ(send_buffer[0].size(), 4); // one obj with 4 doubles are sent to this bbox
    KRATOS_EXPECT_EQ(send_buffer[1].size(), 8); // two objs with each 4 doubles are sent to this bbox
    KRATOS_EXPECT_EQ(send_buffer[2].size(), 12); // three objs with each 4 doubles are sent to this bbox

    // set up expected send buffer
    std::vector<double> exp_send_buffer_rank_1_2 {1.0, 10.0, -25.0, 3.0};
    std::vector<double> exp_send_buffer_rank_2_2 {2.0, -10.0, 15.5, -5.0,
                                                  4.0, -8.301, 17.75, -15.18};
    std::vector<double> exp_send_buffer_rank_3_2 {3.0, 12.6, 50.1, 5.0,
                                                  5.0, 13.5, 44.58, 7.5,
                                                  6.0, 13.5125, 44.68, 8.5};

    std::vector<std::vector<double>> exp_send_buffer_2 {exp_send_buffer_rank_1_2,
                                                        exp_send_buffer_rank_2_2,
                                                        exp_send_buffer_rank_3_2};

    // compare the send buffers
    for (IndexType i=0; i<exp_send_buffer_2.size(); ++i)
        for (IndexType j=0; j<exp_send_buffer_2[i].size(); ++j)
            KRATOS_EXPECT_DOUBLE_EQ(send_buffer[i][j], exp_send_buffer_2[i][j]);
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_CreateMapperInterfaceInfosFromBuffer, KratosMappingApplicationSerialTestSuite)
{
    // set up receive buffer
    std::vector<double> recv_buffer_rank_1 {0.0, -2.0, 3.5, 3.0,
                                            15.0, 10.0, -25.0, 3.0};
    std::vector<double> recv_buffer_rank_2 {0.0, -21.0, 73.5, 35.89,
                                            2.0, -10.0, 15.5, -5.0,
                                            44.0, -89.2, 25.4, 78.7};
    std::vector<double> recv_buffer_rank_3 {}; // nothing received from this rank
    std::vector<double> recv_buffer_rank_4 {2.9999999999999, 12.6, 50.1, 5.0}; // => check in case of truncation

    std::vector<std::vector<double>> recv_buffer {recv_buffer_rank_1,
                                                  recv_buffer_rank_2,
                                                  recv_buffer_rank_3,
                                                  recv_buffer_rank_4};

    std::vector<double> recv_buffer_rank_3_wrong {3.0, 50.1, 5.0}; // has wrong size
    std::vector<double> recv_buffer_rank_3_wrong_2 {3.11, 12.6, 50.1, 5.0}; // 3.11 cannot stem from casting an int to double

    std::vector<std::vector<double>> recv_buffer_wrong {recv_buffer_rank_1,
                                                        recv_buffer_rank_2,
                                                        recv_buffer_rank_3_wrong,
                                                        recv_buffer_rank_4};

    std::vector<std::vector<double>> recv_buffer_wrong_2 {recv_buffer_rank_1,
                                                          recv_buffer_rank_2,
                                                          recv_buffer_rank_3_wrong_2,
                                                          recv_buffer_rank_4};

    const SizeType comm_size = recv_buffer.size();
    const int comm_rank = 16; // only needed for error printing

    MapperInterfaceInfoPointerVectorType interface_info_container;

    MapperInterfaceInfoUniquePointerType p_ref_interface_info(Kratos::make_unique<NearestNeighborInterfaceInfo>());

    // throws bcs "interface_info_container" has the wrong size
    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(MapperUtilities::CreateMapperInterfaceInfosFromBuffer(
        recv_buffer,
        p_ref_interface_info,
        comm_rank,
        interface_info_container),
        "Error: Buffer-size mismatch!");

    interface_info_container.resize(comm_size);

    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(MapperUtilities::CreateMapperInterfaceInfosFromBuffer(
        recv_buffer_wrong,
        p_ref_interface_info,
        comm_rank,
        interface_info_container),
        "Error: Rank 16 received a wrong buffer-size from rank 3!");

    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(MapperUtilities::CreateMapperInterfaceInfosFromBuffer(
        recv_buffer_wrong_2,
        p_ref_interface_info,
        comm_rank,
        interface_info_container),
        "Error: Buffer contains a double (3.11) that was not casted from an int, i.e. it contains a fractional part of 0.11!");

    MapperUtilities::CreateMapperInterfaceInfosFromBuffer(
        recv_buffer,
        p_ref_interface_info,
        comm_rank,
        interface_info_container);

    // Check the created InterfaceInfos
    KRATOS_EXPECT_EQ(interface_info_container[0].size(), 2);
    KRATOS_EXPECT_EQ(interface_info_container[1].size(), 3);
    KRATOS_EXPECT_EQ(interface_info_container[2].size(), 0);
    KRATOS_EXPECT_EQ(interface_info_container[3].size(), 1);

    KRATOS_EXPECT_EQ(interface_info_container[0][0]->GetLocalSystemIndex(), 0);
    KRATOS_EXPECT_EQ(interface_info_container[0][1]->GetLocalSystemIndex(), 15);
    KRATOS_EXPECT_EQ(interface_info_container[1][0]->GetLocalSystemIndex(), 0);
    KRATOS_EXPECT_EQ(interface_info_container[1][1]->GetLocalSystemIndex(), 2);
    KRATOS_EXPECT_EQ(interface_info_container[1][2]->GetLocalSystemIndex(), 44);
    KRATOS_EXPECT_EQ(interface_info_container[3][0]->GetLocalSystemIndex(), 3);

    const auto coords_to_check = interface_info_container[1][0]->Coordinates();
    Point coords_exp(-21.0, 73.5, 35.89);

    for (IndexType i=0; i<3; ++i)
        KRATOS_EXPECT_DOUBLE_EQ(coords_to_check[i], coords_exp[i]);

    // Test if the "Create" function returns the correct object
    const auto& r_arg_1 = *p_ref_interface_info;
    const auto& r_arg_2 = *interface_info_container[0][0];
    KRATOS_EXPECT_EQ(typeid(r_arg_1), typeid(r_arg_2));

    /////
    // now we "update" the Interface and then check again
    /////

    std::vector<double> recv_buffer_rank_1_2 {0.0, -2.0, 3.5, 3.0,
                                              15.0, 10.0, -25.0, 3.0,
                                              18, 12.6, 50.1, 5.0};
    std::vector<double> recv_buffer_rank_2_2 {2.0, -10.0, 15.5, -5.0,
                                              44.0, -89.2, 25.4, 78.7};
    std::vector<double> recv_buffer_rank_3_2 {0.0, -21.0, 73.5, 35.89};
    std::vector<double> recv_buffer_rank_4_2 {};

    std::vector<std::vector<double>> recv_buffer_2 {recv_buffer_rank_1_2,
                                                  recv_buffer_rank_2_2,
                                                  recv_buffer_rank_3_2,
                                                  recv_buffer_rank_4_2};

    MapperUtilities::CreateMapperInterfaceInfosFromBuffer(
        recv_buffer_2,
        p_ref_interface_info,
        comm_rank,
        interface_info_container);

    // Check the created InterfaceInfos
    KRATOS_EXPECT_EQ(interface_info_container[0].size(), 3);
    KRATOS_EXPECT_EQ(interface_info_container[1].size(), 2);
    KRATOS_EXPECT_EQ(interface_info_container[2].size(), 1);
    KRATOS_EXPECT_EQ(interface_info_container[3].size(), 0);

    KRATOS_EXPECT_EQ(interface_info_container[0][0]->GetLocalSystemIndex(), 0);
    KRATOS_EXPECT_EQ(interface_info_container[0][1]->GetLocalSystemIndex(), 15);
    KRATOS_EXPECT_EQ(interface_info_container[0][2]->GetLocalSystemIndex(), 18);
    KRATOS_EXPECT_EQ(interface_info_container[1][0]->GetLocalSystemIndex(), 2);
    KRATOS_EXPECT_EQ(interface_info_container[1][1]->GetLocalSystemIndex(), 44);
    KRATOS_EXPECT_EQ(interface_info_container[2][0]->GetLocalSystemIndex(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_MapperInterfaceInfoSerializer, KratosMappingApplicationSerialTestSuite)
{
    // this test checks if the serialization/deserialization of the Helper-Class that
    // serializes/deserializes the MapperInterfaceInfos works correctly
    // => this is needed to transfer the data btw the ranks in MPI
    // The test is rather large, this is needed to cover a complete example
    // Note that the same checks are performed before and after the serialization
    // to make sure that the objects are properly initialized

    typedef Kratos::shared_ptr<MapperInterfaceInfo> MapperInterfaceInfoPointerType;
    typedef std::vector<std::vector<MapperInterfaceInfoPointerType>> MapperInterfaceInfoPointerVectorType;

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
    auto node_1(Kratos::make_shared<Node>(1, 1.0, 2.5, 30.0));
    auto node_2(Kratos::make_shared<Node>(3, 10.5, 20.0, 96.8));
    auto node_3(Kratos::make_shared<Node>(15, 2.3, 1.9, -2.5));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3.get()));

    const int expected_id_found_1 = 108;
    const int expected_id_found_2 = 18;
    const int expected_id_found_3 = 896;

    node_1->SetValue(INTERFACE_EQUATION_ID, expected_id_found_1);
    node_2->SetValue(INTERFACE_EQUATION_ID, expected_id_found_2);
    node_3->SetValue(INTERFACE_EQUATION_ID, expected_id_found_3);

    // We compute the real distance bcs this would also be computed by the search
    const double dist_1_1 = coords_1.Distance(*interface_node_1);

    p_nearest_neighbor_info_1->ProcessSearchResult(*interface_node_1);
    p_nearest_neighbor_info_1->ProcessSearchResult(*interface_node_2);
    p_nearest_neighbor_info_1->ProcessSearchResult(*interface_node_3);

    // Now some the checks are performed to make sure the objects are correctly initialized
    std::vector<int> found_id(1);
    p_nearest_neighbor_info_1->GetValue(found_id, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_id[0], expected_id_found_1);
    double neighbor_dist;
    p_nearest_neighbor_info_1->GetValue(neighbor_dist, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_1_1);

    const double dist_2_2 = coords_2.Distance(*interface_node_2);

    p_nearest_neighbor_info_2->ProcessSearchResult(*interface_node_1);
    p_nearest_neighbor_info_2->ProcessSearchResult(*interface_node_2);
    p_nearest_neighbor_info_2->ProcessSearchResult(*interface_node_3);

    // Now some the checks are performed to make sure the objects are correctly initialized
    p_nearest_neighbor_info_2->GetValue(found_id, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_id[0], expected_id_found_2);
    p_nearest_neighbor_info_2->GetValue(neighbor_dist, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_2_2);

    const double dist_3_3 = coords_3.Distance(*interface_node_3);

    p_nearest_neighbor_info_3->ProcessSearchResult(*interface_node_1);
    p_nearest_neighbor_info_3->ProcessSearchResult(*interface_node_2);
    p_nearest_neighbor_info_3->ProcessSearchResult(*interface_node_3);

    // Now some the checks are performed to make sure the objects are correctly initialized
    p_nearest_neighbor_info_3->GetValue(found_id, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_id[0], expected_id_found_3);
    p_nearest_neighbor_info_3->GetValue(neighbor_dist, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_3_3);

    // Now finally we can construct the container
    MapperInterfaceInfoPointerVectorType interface_info_container(2);

    // Note: the order is chosen intentionally
    interface_info_container[0].push_back(p_nearest_neighbor_info_3);
    interface_info_container[1].push_back(p_nearest_neighbor_info_1);
    interface_info_container[1].push_back(p_nearest_neighbor_info_2);

    // Construct a reference obj, needed to create the correct objects while loading/deserializing
    const Kratos::unique_ptr<MapperInterfaceInfo> p_ref_nearest_neighbor_info(Kratos::make_unique<NearestNeighborInterfaceInfo>());

    MapperUtilities::MapperInterfaceInfoSerializer serializer_helper_0(
        interface_info_container[0], p_ref_nearest_neighbor_info );
    MapperUtilities::MapperInterfaceInfoSerializer serializer_helper_1(
        interface_info_container[1], p_ref_nearest_neighbor_info );

    // serializing the object (=> happens on the partition that sends the objects)
    StreamSerializer serializer_0;
    serializer_0.save("obj", serializer_helper_0);
    StreamSerializer serializer_1;
    serializer_1.save("obj", serializer_helper_1);

    MapperInterfaceInfoPointerVectorType interface_info_container_new(2);

    MapperUtilities::MapperInterfaceInfoSerializer serializer_helper_new_0(
        interface_info_container_new[0], p_ref_nearest_neighbor_info );

    MapperUtilities::MapperInterfaceInfoSerializer serializer_helper_new_1(
        interface_info_container_new[1], p_ref_nearest_neighbor_info );

    // deserializing the object (=> happens on the partition that receives the objects)
    serializer_0.load("obj", serializer_helper_new_0);
    serializer_1.load("obj", serializer_helper_new_1);

    // Checking for the sizes of the container
    KRATOS_EXPECT_EQ(interface_info_container_new[0].size(), 1);
    KRATOS_EXPECT_EQ(interface_info_container_new[1].size(), 2);

    // Checking the objects inside the container
    const auto& r_info_1 = interface_info_container_new[1][0];
    const auto& r_info_2 = interface_info_container_new[1][1];
    const auto& r_info_3 = interface_info_container_new[0][0];

    r_info_1->GetValue(found_id, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_id[0], expected_id_found_1);
    r_info_2->GetValue(found_id, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_id[0], expected_id_found_2);
    r_info_3->GetValue(found_id, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_id[0], expected_id_found_3);

    r_info_1->GetValue(neighbor_dist, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_1_1);
    r_info_2->GetValue(neighbor_dist, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_2_2);
    r_info_3->GetValue(neighbor_dist, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_3_3);

    // Test if the correct object type was created
    const auto& r_arg_1 = *r_info_1;
    const auto& r_arg_2 = *p_ref_nearest_neighbor_info;
    KRATOS_EXPECT_EQ(typeid(r_arg_1), typeid(r_arg_2));
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_CreateMapperLocalSystemsFromNodes, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");
    CppTestsUtilities::Create2DGeometry(model_part, "Element2D3N", false);

    KRATOS_EXPECT_GE(model_part.NumberOfNodes(), 0);

    std::vector<Kratos::unique_ptr<MapperLocalSystem>> mapper_local_systems;

    MapperUtilities::CreateMapperLocalSystemsFromNodes(
        NearestNeighborLocalSystem(nullptr),
        model_part.GetCommunicator(),
        mapper_local_systems);

    KRATOS_EXPECT_EQ(model_part.NumberOfNodes(), mapper_local_systems.size());
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_EraseNodalVariable, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");
    CppTestsUtilities::Create2DGeometry(model_part, "Element2D3N", false);

    KRATOS_EXPECT_GE(model_part.NumberOfNodes(), 0);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.Has(DISPLACEMENT_X));
        r_node[DISPLACEMENT_X] = 15.3;
        KRATOS_EXPECT_TRUE(r_node.Has(DISPLACEMENT_X));
    }

    MapperUtilities::EraseNodalVariable(model_part, DISPLACEMENT_X);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.Has(DISPLACEMENT_X));
    }
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_SaveCurrentConfiguration, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");
    CppTestsUtilities::Create2DGeometry(model_part, "Element2D3N", false);

    KRATOS_EXPECT_GE(model_part.NumberOfNodes(), 0);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.Has(CURRENT_COORDINATES));
    }

    MapperUtilities::SaveCurrentConfiguration(model_part);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.Has(CURRENT_COORDINATES));
        KRATOS_EXPECT_DOUBLE_EQ(r_node.X(), r_node[CURRENT_COORDINATES][0]);
        KRATOS_EXPECT_DOUBLE_EQ(r_node.Y(), r_node[CURRENT_COORDINATES][1]);
        KRATOS_EXPECT_DOUBLE_EQ(r_node.Z(), r_node[CURRENT_COORDINATES][2]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_RestoreCurrentConfiguration, KratosMappingApplicationSerialTestSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");
    CppTestsUtilities::Create2DGeometry(model_part, "Element2D3N", false);

    KRATOS_EXPECT_GE(model_part.NumberOfNodes(), 0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MapperUtilities::RestoreCurrentConfiguration(model_part), "Nodes do not have CURRENT_COORDINATES for restoring the current configuration!");

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.Has(CURRENT_COORDINATES));
        r_node.X() += 0.1;
        r_node.Y() -= 0.125;
        r_node.Z() += 0.33;
    }

    MapperUtilities::SaveCurrentConfiguration(model_part);

    // X = X0
    VariableUtils().UpdateCurrentToInitialConfiguration(model_part.Nodes());

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.Has(CURRENT_COORDINATES));
        KRATOS_EXPECT_DOUBLE_EQ(r_node.X(), r_node.X0());
        KRATOS_EXPECT_DOUBLE_EQ(r_node.Y(), r_node.Y0());
        KRATOS_EXPECT_DOUBLE_EQ(r_node.Z(), r_node.Z0());
    }

    MapperUtilities::RestoreCurrentConfiguration(model_part);

    for (auto& r_node : model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.Has(CURRENT_COORDINATES));
        KRATOS_EXPECT_DOUBLE_EQ(r_node.X(), (r_node.X0()+0.1));
        KRATOS_EXPECT_DOUBLE_EQ(r_node.Y(), (r_node.Y0()-0.125));
        KRATOS_EXPECT_DOUBLE_EQ(r_node.Z(), (r_node.Z0()+0.33));
    }
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_PointsAreCollinear, KratosMappingApplicationSerialTestSuite)
{
    Point p1(0,0,0);
    Point p2(1,0,0);
    Point p3(2,0,0);
    Point p4(2,1,0);

    KRATOS_EXPECT_TRUE(MapperUtilities::PointsAreCollinear(p1,p2,p3));
    KRATOS_EXPECT_TRUE(MapperUtilities::PointsAreCollinear(p2,p3,p1));
    KRATOS_EXPECT_FALSE(MapperUtilities::PointsAreCollinear(p1,p2,p4));
    KRATOS_EXPECT_FALSE(MapperUtilities::PointsAreCollinear(p1,p3,p4));
    KRATOS_EXPECT_FALSE(MapperUtilities::PointsAreCollinear(p2,p3,p4));
    KRATOS_EXPECT_FALSE(MapperUtilities::PointsAreCollinear(p2,p3,p4));
}

}  // namespace Kratos::Testing