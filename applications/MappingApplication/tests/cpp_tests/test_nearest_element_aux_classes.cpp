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
#include "geometries/triangle_3d_3.h"
#include "testing/testing.h"
#include "custom_mappers/nearest_element_mapper.h"
#include "includes/stream_serializer.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos {
namespace Testing {

typedef typename MapperLocalSystem::MatrixType MatrixType;
typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;

typedef Node<3> NodeType;

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    const Point coords(1.0, 2.45, -23.8);

    const std::size_t source_local_sys_idx = 123;
    const std::size_t dummy_rank = 78;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);

    const auto nearest_element_info_1(nearest_element_info.Create());
    const auto nearest_element_info_2(nearest_element_info.Create(coords, source_local_sys_idx));
    const auto nearest_element_info_3(nearest_element_info.Create(coords, source_local_sys_idx, dummy_rank));

    // Test if the "Create" function returns the correct object
    KRATOS_CHECK_EQUAL(typeid(nearest_element_info), typeid(*nearest_element_info_1));
    KRATOS_CHECK_EQUAL(typeid(nearest_element_info), typeid(*nearest_element_info_2));
    KRATOS_CHECK_EQUAL(typeid(nearest_element_info), typeid(*nearest_element_info_3));
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_ValidProjectionExists, KratosMappingApplicationSerialTestSuite)
{
    // NOTE: In this tests "tria_2" has a valid projection

    const double dist = 1.1;
    Point coords(0.3, -0.3, dist);

    std::size_t source_local_sys_idx = 123;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);

    auto node_1(Kratos::make_shared<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_shared<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_shared<NodeType>(3, 1.0, 1.0, 0.0));
    auto node_4(Kratos::make_shared<NodeType>(4, 0.0, -1.0, 0.0));
    auto node_5(Kratos::make_shared<NodeType>(5, 2.0, -1.0, 0.0));

    const Geometry<NodeType>::Pointer tria_1(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));
    const Geometry<NodeType>::Pointer tria_2(Kratos::make_shared<Triangle3D3<NodeType>>(node_4, node_2, node_1));
    const Geometry<NodeType>::Pointer tria_3(Kratos::make_shared<Triangle3D3<NodeType>>(node_4, node_5, node_2));

    InterfaceObject::Pointer interface_geom_obj_1(Kratos::make_shared<InterfaceGeometryObject>(tria_1));
    InterfaceObject::Pointer interface_geom_obj_2(Kratos::make_shared<InterfaceGeometryObject>(tria_2));
    InterfaceObject::Pointer interface_geom_obj_3(Kratos::make_shared<InterfaceGeometryObject>(tria_3));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);
    node_4->SetValue(INTERFACE_EQUATION_ID, 61);
    node_5->SetValue(INTERFACE_EQUATION_ID, 899);

    // Distances do not matter bcs only one projection is valid!
    nearest_element_info.ProcessSearchResult(interface_geom_obj_1, 10.0);
    nearest_element_info.ProcessSearchResult(interface_geom_obj_2, 11.0);
    nearest_element_info.ProcessSearchResult(interface_geom_obj_3, 33.5);

    KRATOS_CHECK(nearest_element_info.GetLocalSearchWasSuccessful());
    KRATOS_CHECK_IS_FALSE(nearest_element_info.GetIsApproximation());

    double proj_dist;
    nearest_element_info.GetValue(proj_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(proj_dist, dist);

    // the nodes found are 4,2,1, aka (61,18,35)
    std::vector<int> found_ids;
    nearest_element_info.GetValue(found_ids);
    KRATOS_CHECK_EQUAL(found_ids.size(), 3);
    KRATOS_CHECK_EQUAL(found_ids[0], 61);
    KRATOS_CHECK_EQUAL(found_ids[1], 18);
    KRATOS_CHECK_EQUAL(found_ids[2], 35);

    std::vector<double> sf_values;
    nearest_element_info.GetValue(sf_values);
    KRATOS_CHECK_EQUAL(sf_values.size(), 3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[0], 0.3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[1], 0.3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[2], 0.4);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_Approximation, KratosMappingApplicationSerialTestSuite)
{
    const double dist = 1.33;
    Point coords(0.0, -dist, 0.0);

    std::size_t source_local_sys_idx = 123;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);

    auto node_1(Kratos::make_shared<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_shared<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_shared<NodeType>(3, 1.0, 1.0, 0.0));

    const Geometry<NodeType>::Pointer tria_1(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));

    InterfaceObject::Pointer interface_geom_obj_1(Kratos::make_shared<InterfaceGeometryObject>(tria_1));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);

    // Distances do not matter bcs only one projection is valid!
    nearest_element_info.ProcessSearchResult(interface_geom_obj_1, 10.0);
    KRATOS_CHECK_IS_FALSE(nearest_element_info.GetLocalSearchWasSuccessful());
    // since the no valid projection could be found we try to get an approximation
    nearest_element_info.ProcessSearchResultForApproximation(interface_geom_obj_1, 10.0);

    KRATOS_CHECK(nearest_element_info.GetLocalSearchWasSuccessful());
    KRATOS_CHECK(nearest_element_info.GetIsApproximation());

    double approximation_dist;
    nearest_element_info.GetValue(approximation_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(approximation_dist, dist);

    std::vector<int> found_ids;
    nearest_element_info.GetValue(found_ids);
    KRATOS_CHECK_EQUAL(found_ids.size(), 1);
    KRATOS_CHECK_EQUAL(found_ids[0], 35);

    std::vector<double> sf_values;
    nearest_element_info.GetValue(sf_values);
    KRATOS_CHECK_EQUAL(sf_values.size(), 1);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[0], 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_NothingFound, KratosMappingApplicationSerialTestSuite)
{
    // here no search-results are being processed, therfore nothing is found
    const Point coords(1.0, 2.45, -23.8);

    const std::size_t source_local_sys_idx = 123;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);

    KRATOS_CHECK_IS_FALSE(nearest_element_info.GetLocalSearchWasSuccessful());
    KRATOS_CHECK_IS_FALSE(nearest_element_info.GetIsApproximation());

    // Testing the uninitialized Values
    std::vector<int> found_ids;
    nearest_element_info.GetValue(found_ids);
    KRATOS_CHECK_EQUAL(found_ids.size(), 0);

    std::vector<double> sf_values;
    nearest_element_info.GetValue(sf_values);
    KRATOS_CHECK_EQUAL(sf_values.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_Serialization, KratosMappingApplicationSerialTestSuite)
{
    // NOTE: In this tests "tria_2" has a valid projection, the rest is being removed to test only the serialization

    const double dist = 1.1;
    Point coords(0.3, -0.3, dist);

    std::size_t source_local_sys_idx = 123;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);

    auto node_1(Kratos::make_shared<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_shared<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_4(Kratos::make_shared<NodeType>(4, 0.0, -1.0, 0.0));

    const Geometry<NodeType>::Pointer tria_2(Kratos::make_shared<Triangle3D3<NodeType>>(node_4, node_2, node_1));

    InterfaceObject::Pointer interface_geom_obj_2(Kratos::make_shared<InterfaceGeometryObject>(tria_2));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_4->SetValue(INTERFACE_EQUATION_ID, 61);

    // Distances do not matter bcs only one projection is valid!
    nearest_element_info.ProcessSearchResult(interface_geom_obj_2, 11.0);

    KRATOS_CHECK(nearest_element_info.GetLocalSearchWasSuccessful());

    // serializing the object
    StreamSerializer serializer;
    serializer.save("nearest_element_interface_info", nearest_element_info);
    // deserializing the object => this happens if the remote search was successful and
    // sending back of the object to the partition where it came from is required
    NearestElementInterfaceInfo nearest_element_info_new;
    serializer.load("nearest_element_interface_info", nearest_element_info_new);

    KRATOS_CHECK_EQUAL(nearest_element_info_new.GetLocalSystemIndex(), source_local_sys_idx);

    KRATOS_CHECK_IS_FALSE(nearest_element_info_new.GetIsApproximation());

    double proj_dist;
    nearest_element_info_new.GetValue(proj_dist);
    KRATOS_CHECK_DOUBLE_EQUAL(proj_dist, dist);

    // the nodes found are 4,2,1, aka (61,18,35)
    std::vector<int> found_ids;
    nearest_element_info_new.GetValue(found_ids);
    KRATOS_CHECK_EQUAL(found_ids.size(), 3);
    KRATOS_CHECK_EQUAL(found_ids[0], 61);
    KRATOS_CHECK_EQUAL(found_ids[1], 18);
    KRATOS_CHECK_EQUAL(found_ids[2], 35);

    std::vector<double> sf_values;
    nearest_element_info_new.GetValue(sf_values);
    KRATOS_CHECK_EQUAL(sf_values.size(), 3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[0], 0.3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[1], 0.3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[2], 0.4);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    NearestElementLocalSystem local_sys_dummy;

    auto node_local_sys(Kratos::make_shared<Node<3>>(8, 1.0, 2.5, -5.0));

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
        "NearestElementLocalSystem based on Node #8 at Coodinates 1 | 2.5 | -5 in rank 23");
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Line, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_ERROR <<  "This test is not yet implemented!" << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Triangle, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_ERROR <<  "This test is not yet implemented!" << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Quad, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_ERROR <<  "This test is not yet implemented!" << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Tetra, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_ERROR <<  "This test is not yet implemented!" << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Hexa, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_ERROR <<  "This test is not yet implemented!" << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystemWithApproximation, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_ERROR <<  "This test is not yet implemented!" << std::endl;
}

}  // namespace Testing
}  // namespace Kratos