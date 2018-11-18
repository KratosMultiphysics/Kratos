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
#include "includes/serializer.h"
#include "testing/testing.h"
#include "custom_mappers/nearest_neighbor_mapper.h"
#include "includes/stream_serializer.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos {
namespace Testing {

typedef typename MapperLocalSystem::MatrixType MatrixType;
typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;

typedef Node<3> NodeType;

KRATOS_TEST_CASE_IN_SUITE(MapperInterfaceInfo_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    // This test covers the basic functionalities provided by the "MapperInterfaceInfo"
    // A "NearestNeighborInterfaceInfo" is being used since "MapperInterfaceInfo" is a pure virtual class

    Point coords_1(1.0, 2.45, -23.8);

    std::size_t source_local_sys_idx = 55;
    std::size_t source_rank = 23;

    NearestNeighborInterfaceInfo nearest_neighbor_info(coords_1, source_local_sys_idx, source_rank);
    const auto nearest_neighbor_info_2(nearest_neighbor_info.Create(coords_1, source_local_sys_idx));

    KRATOS_CHECK_EQUAL(nearest_neighbor_info.GetLocalSystemIndex(), source_local_sys_idx);

    KRATOS_CHECK_EQUAL(nearest_neighbor_info.GetSourceRank(), source_rank);
    KRATOS_CHECK_EQUAL(nearest_neighbor_info_2->GetSourceRank(), 0); // Testing against default

    KRATOS_CHECK_IS_FALSE(nearest_neighbor_info.GetLocalSearchWasSuccessful());
    KRATOS_CHECK_IS_FALSE(nearest_neighbor_info.GetIsApproximation());

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(nearest_neighbor_info.Coordinates()[i], coords_1[i]);
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfo_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    const Point coords(1.0, 2.45, -23.8);

    const std::size_t source_local_sys_idx = 123;
    const std::size_t dummy_rank = 78;

    NearestNeighborInterfaceInfo nearest_neighbor_info(coords, source_local_sys_idx, 0);

    const auto nearest_neighbor_info_1(nearest_neighbor_info.Create());
    const auto nearest_neighbor_info_2(nearest_neighbor_info.Create(coords, source_local_sys_idx));
    const auto nearest_neighbor_info_3(nearest_neighbor_info.Create(coords, source_local_sys_idx, dummy_rank));

    // Test if the "Create" function returns the correct object
    KRATOS_CHECK_EQUAL(typeid(nearest_neighbor_info), typeid(*nearest_neighbor_info_1));
    KRATOS_CHECK_EQUAL(typeid(nearest_neighbor_info), typeid(*nearest_neighbor_info_2));
    KRATOS_CHECK_EQUAL(typeid(nearest_neighbor_info), typeid(*nearest_neighbor_info_3));
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfo_NeighborsFound, KratosMappingApplicationSerialTestSuite)
{
    // NOTE: In this tests "node_3" is the closest node

    Point coords(1.0, 2.5, -3.0);

    std::size_t source_local_sys_idx = 123;

    NearestNeighborInterfaceInfo nearest_neighbor_info(coords, source_local_sys_idx, 0);

    auto node_1(Kratos::make_shared<NodeType>(1, 1.0, 2.5, 30.0));
    auto node_2(Kratos::make_shared<NodeType>(3, 10.5, 20.0, 96.8));
    auto node_3(Kratos::make_shared<NodeType>(15, 2.3, 1.9, -2.5));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3));

    const int expected_id_found = 108;

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, expected_id_found);

    // We compute the real distance bcs this would also be computed by the search
    const double dist_1 = MapperUtilities::ComputeDistance(coords, *interface_node_1);
    const double dist_2 = MapperUtilities::ComputeDistance(coords, *interface_node_2);
    const double dist_3 = MapperUtilities::ComputeDistance(coords, *interface_node_3);

    nearest_neighbor_info.ProcessSearchResult(interface_node_1, dist_1);
    nearest_neighbor_info.ProcessSearchResult(interface_node_2, dist_2);
    nearest_neighbor_info.ProcessSearchResult(interface_node_3, dist_3);

    KRATOS_CHECK(nearest_neighbor_info.GetLocalSearchWasSuccessful());
    // this function should never return true for this class!
    // It already works with nodes, which could be approximations for other mappers
    KRATOS_CHECK_IS_FALSE(nearest_neighbor_info.GetIsApproximation());

    int found_id;
    nearest_neighbor_info.GetValue(found_id);
    KRATOS_CHECK_EQUAL(found_id, expected_id_found);

    double neighbor_dist;
    nearest_neighbor_info.GetValue(neighbor_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_3);
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfo_MatchingNeighborFound, KratosMappingApplicationSerialTestSuite)
{
    // NOTE: In this tests "node_2" is a matching Node i.e. it has the same coordinates

    Point coords(1.0, 2.5, -3.0);

    std::size_t source_local_sys_idx = 123;

    NearestNeighborInterfaceInfo nearest_neighbor_info(coords, source_local_sys_idx, 0);

    auto node_1(Kratos::make_shared<NodeType>(1, 18.0, 2.7, 30.0));
    auto node_2(Kratos::make_shared<NodeType>(3, 1.0, 2.5, -3.0));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2));

    const int expected_id_found = 67;

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, expected_id_found);

    // We compute the real distance bcs this would also be computed by the search
    const double dist_1 = MapperUtilities::ComputeDistance(coords, *interface_node_1);
    const double dist_2 = MapperUtilities::ComputeDistance(coords, *interface_node_2);

    KRATOS_CHECK_IS_FALSE(nearest_neighbor_info.GetLocalSearchWasSuccessful()); // this is the default

    nearest_neighbor_info.ProcessSearchResult(interface_node_1, dist_1);
    nearest_neighbor_info.ProcessSearchResult(interface_node_2, dist_2);

    KRATOS_CHECK(nearest_neighbor_info.GetLocalSearchWasSuccessful());

    int found_id;
    nearest_neighbor_info.GetValue(found_id);
    KRATOS_CHECK_EQUAL(found_id, expected_id_found);

    double neighbor_dist;
    nearest_neighbor_info.GetValue(neighbor_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_2);
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfo_Serialization, KratosMappingApplicationSerialTestSuite)
{
    // NOTE: In this tests "node_3" is the closest node

    Point coords(1.0, 2.5, -3.0);

    std::size_t source_local_sys_idx = 123;

    NearestNeighborInterfaceInfo nearest_neighbor_info(coords, source_local_sys_idx, 0);

    auto node_2(Kratos::make_shared<NodeType>(3, 10.5, 20.0, 96.8));
    auto node_3(Kratos::make_shared<NodeType>(15, 2.3, 1.9, -2.5));

    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3));

    const int expected_id_found = 108;

    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, expected_id_found);

    // We compute the real distance bcs this would also be computed by the search
    const double dist_2 = MapperUtilities::ComputeDistance(coords, *interface_node_2);
    const double dist_3 = MapperUtilities::ComputeDistance(coords, *interface_node_3);

    nearest_neighbor_info.ProcessSearchResult(interface_node_2, dist_2);
    nearest_neighbor_info.ProcessSearchResult(interface_node_3, dist_3);

    // serializing the object
    StreamSerializer serializer;
    serializer.save("nearest_neighbor_interface_info", nearest_neighbor_info);
    // deserializing the object => this happens if the remote search was successful and
    // sending back of the object to the partition where it came from is required
    NearestNeighborInterfaceInfo nearest_neighbor_info_new;
    serializer.load("nearest_neighbor_interface_info", nearest_neighbor_info_new);

    KRATOS_CHECK_EQUAL(nearest_neighbor_info_new.GetLocalSystemIndex(), source_local_sys_idx);

    int found_id;
    nearest_neighbor_info_new.GetValue(found_id);
    KRATOS_CHECK_EQUAL(found_id, expected_id_found);

    double neighbor_dist;
    nearest_neighbor_info_new.GetValue(neighbor_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_3);
}

KRATOS_TEST_CASE_IN_SUITE(MapperLocalSystem_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    // This test covers the basic functionalities provided by the "MapperLocalSystem"
    // A "NearestNeighborLocalSystem" is being used since "MapperLocalSystem" is a pure virtual class

    Point coords_1(1.0, 2.45, -23.8);
    auto node_local_sys(Kratos::make_shared<NodeType>(5, coords_1));

    NearestNeighborLocalSystem local_sys(node_local_sys);

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(local_sys.Coordinates()[i], coords_1[i]);

    KRATOS_CHECK_IS_FALSE(local_sys.HasInterfaceInfo());
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborLocalSystem_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    NearestNeighborLocalSystem local_sys_dummy;

    auto node_local_sys(Kratos::make_shared<NodeType>(8, 1.0, 2.5, -5.0));

    auto local_sys(local_sys_dummy.Create(node_local_sys));

    // Test if the "Create" function returns the correct object
    KRATOS_CHECK_EQUAL(typeid(local_sys_dummy), typeid(*local_sys));

    KRATOS_CHECK(local_sys->UseNodesAsBasis());

    // Computing the local system
    // this should return nothing since no InterfaceInfos are available
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType origin_ids2;
    EquationIdVectorType destination_ids;
    EquationIdVectorType destination_ids2;

    local_sys->EquationIdVectors(origin_ids, destination_ids);

    KRATOS_CHECK_EQUAL(origin_ids.size(), 0);
    KRATOS_CHECK_EQUAL(destination_ids.size(), 0);

    local_sys->CalculateLocalSystem(local_mapping_matrix, origin_ids2, destination_ids2);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size1(), 0);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size2(), 0);
    KRATOS_CHECK_EQUAL(origin_ids2.size(), 0);
    KRATOS_CHECK_EQUAL(destination_ids2.size(), 0);

    KRATOS_CHECK_C_STRING_EQUAL((local_sys->PairingInfo(2,23)).c_str(),
        "NearestNeighborLocalSystem based on Node #8 at Coodinates 1 | 2.5 | -5 in rank 23");
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborLocalSystem_ComputeLocalSystem, KratosMappingApplicationSerialTestSuite)
{
    NearestNeighborLocalSystem local_sys_dummy;

    const int dest_id = 13;

    auto node_local_sys(Kratos::make_shared<NodeType>(5, 1.0, 2.5, -5.0));
    node_local_sys->SetValue(INTERFACE_EQUATION_ID, dest_id);

    auto local_sys(local_sys_dummy.Create(node_local_sys));

    // Test if the "Create" function returns the correct object
    KRATOS_CHECK_EQUAL(typeid(local_sys_dummy), typeid(*local_sys));

    // Create the NearestNeighborInfos to be used by the NearestNeighborLocalSystem
    auto node_1(Kratos::make_shared<NodeType>(1, 18.0, 2.7, 30.0));
    auto node_2(Kratos::make_shared<NodeType>(3, 1.0, 2.5, -3.0)); // this is the nearest neighbor

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2));

    const int expected_id_found = 67;

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, expected_id_found);

    // We compute the real distance bcs this would also be computed by the search
    const double dist_1 = MapperUtilities::ComputeDistance(local_sys->Coordinates(), *interface_node_1);
    const double dist_2 = MapperUtilities::ComputeDistance(local_sys->Coordinates(), *interface_node_2);

    MapperInterfaceInfo::Pointer p_nearest_neighbor_info_1(Kratos::make_shared<NearestNeighborInterfaceInfo>(local_sys->Coordinates(), 0, 0));
    MapperInterfaceInfo::Pointer p_nearest_neighbor_info_2(Kratos::make_shared<NearestNeighborInterfaceInfo>(local_sys->Coordinates(), 0, 0));

    p_nearest_neighbor_info_1->ProcessSearchResult(interface_node_1, dist_1);
    p_nearest_neighbor_info_2->ProcessSearchResult(interface_node_2, dist_2);

    local_sys->AddInterfaceInfo(p_nearest_neighbor_info_1);
    local_sys->AddInterfaceInfo(p_nearest_neighbor_info_2);

    // Computing the local system
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType origin_ids2;
    EquationIdVectorType destination_ids;
    EquationIdVectorType destination_ids2;

    local_sys->EquationIdVectors(origin_ids, destination_ids);

    KRATOS_CHECK_EQUAL(origin_ids.size(), 1);
    KRATOS_CHECK_EQUAL(destination_ids.size(), 1);

    local_sys->CalculateLocalSystem(local_mapping_matrix, origin_ids2, destination_ids2);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size1(), 1);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size2(), 1);
    KRATOS_CHECK_EQUAL(origin_ids2.size(), 1);
    KRATOS_CHECK_EQUAL(destination_ids2.size(), 1);

    KRATOS_CHECK_DOUBLE_EQUAL(local_mapping_matrix(0,0), 1.0);
    KRATOS_CHECK_EQUAL(origin_ids2[0], expected_id_found);
    KRATOS_CHECK_EQUAL(destination_ids2[0], dest_id);
}

}  // namespace Testing
}  // namespace Kratos