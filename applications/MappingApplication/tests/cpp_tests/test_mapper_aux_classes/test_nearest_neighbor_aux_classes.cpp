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
#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(InterfaceInfoBasicTests, KratosMappingApplicationSerialTestSuite)
{
    // This test covers the basic functionalities provided by the "MapperInterfaceInfo"
    // A "NearestNeigborInterfaceInfo" is being used since "MapperInterfaceInfo" is a pure virtual class

    Point coords_1(1.0, 2.45, -23.8);

    std::size_t source_local_sys_idx = 55;
    std::size_t source_rank = 23;

    NearestNeigborInterfaceInfo nearest_neighbor_info(coords_1, source_local_sys_idx, source_rank);
    const auto nearest_neighbor_info_2(nearest_neighbor_info.Create(coords_1, source_local_sys_idx));

    KRATOS_CHECK_EQUAL(nearest_neighbor_info.GetLocalSystemIndex(), source_local_sys_idx);

    KRATOS_CHECK_EQUAL(nearest_neighbor_info.GetSourceRank(), source_rank);
    KRATOS_CHECK_EQUAL(nearest_neighbor_info_2->GetSourceRank(), 0); // Testing against default

    KRATOS_CHECK_IS_FALSE(nearest_neighbor_info.GetLocalSearchWasSuccessful());
    KRATOS_CHECK_IS_FALSE(nearest_neighbor_info.GetIsApproximation());

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(nearest_neighbor_info.GetCoordinates()[i], coords_1[i]);
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfo_NoNeighbor, KratosMappingApplicationSerialTestSuite)
{
    Point coords(1.0, 2.45, -23.8);

    std::size_t source_local_sys_idx = 123;

    NearestNeigborInterfaceInfo nearest_neighbor_info(coords, source_local_sys_idx, 0);
    const auto nearest_neighbor_info_2(nearest_neighbor_info.Create(coords, source_local_sys_idx));

    // Test if the "Create" function returns the correct object
    KRATOS_CHECK_EQUAL(typeid(nearest_neighbor_info), typeid(*nearest_neighbor_info_2));

    KRATOS_CHECK_EQUAL(nearest_neighbor_info_2->GetLocalSystemIndex(), source_local_sys_idx);
    KRATOS_CHECK_EQUAL(nearest_neighbor_info_2->GetSourceRank(), 0);

    KRATOS_CHECK_IS_FALSE(nearest_neighbor_info_2->GetLocalSearchWasSuccessful());
    KRATOS_CHECK_IS_FALSE(nearest_neighbor_info_2->GetIsApproximation());

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(nearest_neighbor_info_2->GetCoordinates()[i], coords[i]);
}

KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfo_NeighborsFound, KratosMappingApplicationSerialTestSuite)
{
    // NOTE: In this tests "node_3" is the closest node

    Point coords(1.0, 2.5, -3.0);

    std::size_t source_local_sys_idx = 123;

    NearestNeigborInterfaceInfo nearest_neighbor_info(coords, source_local_sys_idx, 0);

    auto node_1(Kratos::make_shared<Node<3>>(1, 1.0, 2.5, 30.0));
    auto node_2(Kratos::make_shared<Node<3>>(3, 10.5, 20.0, 96.8));
    auto node_3(Kratos::make_shared<Node<3>>(15, 2.3, 1.9, -2.5));

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

    NearestNeigborInterfaceInfo nearest_neighbor_info(coords, source_local_sys_idx, 0);

    auto node_1(Kratos::make_shared<Node<3>>(1, 18.0, 2.7, 30.0));
    auto node_2(Kratos::make_shared<Node<3>>(3, 1.0, 2.5, -3.0));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2));

    const int expected_id_found = 67;

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, expected_id_found);

    // We compute the real distance bcs this would also be computed by the search
    const double dist_1 = MapperUtilities::ComputeDistance(coords, *interface_node_1);
    const double dist_2 = MapperUtilities::ComputeDistance(coords, *interface_node_2);

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

    NearestNeigborInterfaceInfo nearest_neighbor_info(coords, source_local_sys_idx, 0);

    auto node_2(Kratos::make_shared<Node<3>>(3, 10.5, 20.0, 96.8));
    auto node_3(Kratos::make_shared<Node<3>>(15, 2.3, 1.9, -2.5));

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

    // SERIALIZE

    // DESERIALIZE

    // int found_id;
    // nearest_neighbor_info_new.GetValue(found_id);
    // KRATOS_CHECK_EQUAL(found_id, expected_id_found);

    // double neighbor_dist;
    // nearest_neighbor_info_new.GetValue(neighbor_dist);
    // KRATOS_CHECK_DOUBLE_EQUAL(neighbor_dist, dist_3);
}

}  // namespace Testing
}  // namespace Kratos