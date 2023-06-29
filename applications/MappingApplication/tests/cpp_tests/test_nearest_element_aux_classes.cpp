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
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "testing/testing.h"
#include "custom_mappers/nearest_element_mapper.h"
#include "includes/stream_serializer.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos::Testing {

typedef typename MapperLocalSystem::MatrixType MatrixType;
typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;

typedef std::vector<std::vector<double>> MatrixResultsType;
typedef std::vector<int> EqIDVectorResultsType;

typedef Node NodeType;

namespace {

void TestNearestElementLocalSystem(const MatrixResultsType& rExpMatrix,
                                   const EqIDVectorResultsType& rExpOriginIds,
                                   const Geometry<NodeType>::Pointer pGeom)
{
    // validate the input
    KRATOS_CHECK_EQUAL(rExpMatrix[0].size(), rExpOriginIds.size());
    double sum_entries = 0.0;
    for (auto& n : rExpMatrix[0]){
        sum_entries += n;
    }
    KRATOS_CHECK_DOUBLE_EQUAL(sum_entries, 1.0);

    auto node_local_sys(Kratos::make_intrusive<NodeType>(5, 0.5, 0.3, 0.2));
    const int dest_id = 13;
    node_local_sys->SetValue(INTERFACE_EQUATION_ID, dest_id);

    NearestElementLocalSystem local_sys(node_local_sys.get());


    InterfaceObject::Pointer interface_geom_obj(Kratos::make_shared<InterfaceGeometryObject>(pGeom.get()));
    MapperInterfaceInfo::Pointer p_nearest_elem_info(Kratos::make_shared<NearestElementInterfaceInfo>(local_sys.Coordinates(), 0, 0));
    p_nearest_elem_info->ProcessSearchResult(*interface_geom_obj);

    if (!local_sys.HasInterfaceInfo()) { // in case we test an approximation
        p_nearest_elem_info->ProcessSearchResultForApproximation(*interface_geom_obj);
    }

    local_sys.AddInterfaceInfo(p_nearest_elem_info);

    // Computing the local system
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType origin_ids2;
    EquationIdVectorType destination_ids;
    EquationIdVectorType destination_ids2;

    local_sys.EquationIdVectors(origin_ids, destination_ids);

    KRATOS_CHECK_EQUAL(origin_ids.size(), rExpOriginIds.size());
    for (std::size_t i=0; i<rExpOriginIds.size(); ++i) {
        KRATOS_CHECK_EQUAL(origin_ids[i], rExpOriginIds[i]);
    }
    KRATOS_CHECK_EQUAL(destination_ids.size(), 1);
    KRATOS_CHECK_EQUAL(destination_ids[0], dest_id);


    local_sys.CalculateLocalSystem(local_mapping_matrix, origin_ids2, destination_ids2);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size1(), 1);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size2(), rExpOriginIds.size());
    KRATOS_CHECK_EQUAL(origin_ids2.size(), rExpOriginIds.size());
    KRATOS_CHECK_EQUAL(destination_ids2.size(), 1);

    for (std::size_t i=0; i<local_mapping_matrix.size1(); ++i) {
        for (std::size_t j=0; j<local_mapping_matrix.size2(); ++j) {
            KRATOS_CHECK_NEAR(local_mapping_matrix(i,j), rExpMatrix[i][j], 1e-14);
        }
    }
    for (std::size_t i=0; i<rExpOriginIds.size(); ++i) {
        KRATOS_CHECK_EQUAL(origin_ids[i], rExpOriginIds[i]);
    }
    KRATOS_CHECK_EQUAL(destination_ids2[0], dest_id);
}

}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    const Point coords(1.0, 2.45, -23.8);

    const std::size_t source_local_sys_idx = 123;
    const std::size_t dummy_rank = 78;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);

    const auto nearest_element_info_1(nearest_element_info.Create());
    const auto nearest_element_info_2(nearest_element_info.Create(coords, source_local_sys_idx, dummy_rank));

    // Test if the "Create" function returns the correct object
    const auto& r_arg_1 = *nearest_element_info_1;
    const auto& r_arg_2 = *nearest_element_info_2;
    KRATOS_CHECK_EQUAL(typeid(nearest_element_info), typeid(r_arg_1));
    KRATOS_CHECK_EQUAL(typeid(nearest_element_info), typeid(r_arg_2));
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_ValidProjectionExists, KratosMappingApplicationSerialTestSuite)
{
    // NOTE: In this tests "tria_2" has a valid projection

    const double dist = 1.1;
    Point coords(0.3, -0.3, dist);

    std::size_t source_local_sys_idx = 123;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);

    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.0, -1.0, 0.0));
    auto node_5(Kratos::make_intrusive<NodeType>(5, 2.0, -1.0, 0.0));

    const Geometry<NodeType>::Pointer tria_1(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));
    const Geometry<NodeType>::Pointer tria_2(Kratos::make_shared<Triangle3D3<NodeType>>(node_4, node_2, node_1));
    const Geometry<NodeType>::Pointer tria_3(Kratos::make_shared<Triangle3D3<NodeType>>(node_4, node_5, node_2));

    InterfaceObject::Pointer interface_geom_obj_1(Kratos::make_shared<InterfaceGeometryObject>(tria_1.get()));
    InterfaceObject::Pointer interface_geom_obj_2(Kratos::make_shared<InterfaceGeometryObject>(tria_2.get()));
    InterfaceObject::Pointer interface_geom_obj_3(Kratos::make_shared<InterfaceGeometryObject>(tria_3.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);
    node_4->SetValue(INTERFACE_EQUATION_ID, 61);
    node_5->SetValue(INTERFACE_EQUATION_ID, 899);

    // Distances do not matter bcs only one projection is valid!
    nearest_element_info.ProcessSearchResult(*interface_geom_obj_1);
    nearest_element_info.ProcessSearchResult(*interface_geom_obj_2);
    nearest_element_info.ProcessSearchResult(*interface_geom_obj_3);

    KRATOS_CHECK(nearest_element_info.GetLocalSearchWasSuccessful());
    KRATOS_CHECK_IS_FALSE(nearest_element_info.GetIsApproximation());

    double proj_dist;
    nearest_element_info.GetValue(proj_dist, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_CHECK_DOUBLE_EQUAL(proj_dist, dist);

    // the nodes found are 4,2,1, aka (61,18,35)
    std::vector<int> found_ids;
    nearest_element_info.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_CHECK_EQUAL(found_ids.size(), 3);
    KRATOS_CHECK_EQUAL(found_ids[0], 61);
    KRATOS_CHECK_EQUAL(found_ids[1], 18);
    KRATOS_CHECK_EQUAL(found_ids[2], 35);

    std::vector<double> sf_values;
    nearest_element_info.GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_CHECK_EQUAL(sf_values.size(), 3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[0], 0.3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[1], 0.3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[2], 0.4);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_WithoutApproximation, KratosMappingApplicationSerialTestSuite)
{
    const double dist = 1.1;
    Point coords(-0.01, -0.3, dist);

    std::size_t source_local_sys_idx = 123;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);
    NearestElementInterfaceInfo nearest_element_info_without_approx(coords, source_local_sys_idx, 0, {false});

    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 0.0, -1.0, 0.0));

    const Geometry<NodeType>::Pointer tria(Kratos::make_shared<Triangle3D3<NodeType>>(node_3, node_2, node_1));

    InterfaceObject::Pointer interface_geom_obj(Kratos::make_shared<InterfaceGeometryObject>(tria.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, 61);

    nearest_element_info.ProcessSearchResultForApproximation(*interface_geom_obj);
    nearest_element_info_without_approx.ProcessSearchResultForApproximation(*interface_geom_obj);

    KRATOS_CHECK(nearest_element_info.GetLocalSearchWasSuccessful());
    KRATOS_CHECK(nearest_element_info.GetIsApproximation());

    KRATOS_CHECK_IS_FALSE(nearest_element_info_without_approx.GetLocalSearchWasSuccessful());
}

// KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_Approximation, KratosMappingApplicationSerialTestSuite)
// {
//     const double z_coord = 0.2;
//     Point coords_1(1.1, 0.5, z_coord);
//     Point coords_2(-0.1, 1.33, 0.0);
//     Point coords_3(10.0, 10.0, 0.0);

//     std::size_t source_local_sys_idx = 123;

//     NearestElementInterfaceInfo nearest_element_info_1(coords_1, source_local_sys_idx, 0, 0.2);
//     NearestElementInterfaceInfo nearest_element_info_2(coords_2, source_local_sys_idx, 0);
//     NearestElementInterfaceInfo nearest_element_info_3(coords_3, source_local_sys_idx, 0, 0.2);

//     auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
//     auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
//     auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));

//     const Geometry<NodeType>::Pointer tria_1(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));

//     InterfaceObject::Pointer interface_geom_obj(Kratos::make_shared<InterfaceGeometryObject>(tria_1.get()));

//     node_1->SetValue(INTERFACE_EQUATION_ID, 35);
//     node_2->SetValue(INTERFACE_EQUATION_ID, 18);
//     node_3->SetValue(INTERFACE_EQUATION_ID, 108);

//     // Testing different projections
//     double approximation_dist;
//     int pairing_index;
//     std::vector<int> found_ids;
//     std::vector<double> sf_values;

//     // Surface Outside //
//     nearest_element_info_1.ProcessSearchResult(*interface_geom_obj, 10.0);
//     KRATOS_CHECK_IS_FALSE(nearest_element_info_1.GetLocalSearchWasSuccessful());
//     // since the no valid projection could be found we try to get an approximation
//     nearest_element_info_1.ProcessSearchResultForApproximation(*interface_geom_obj, 10.0);

//     KRATOS_CHECK(nearest_element_info_1.GetLocalSearchWasSuccessful());
//     KRATOS_CHECK(nearest_element_info_1.GetIsApproximation());

//     nearest_element_info_1.GetValue(approximation_dist, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_DOUBLE_EQUAL(approximation_dist, z_coord);

//     nearest_element_info_1.GetValue(pairing_index, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_EQUAL((ProjectionUtilities::PairingIndex)pairing_index, ProjectionUtilities::PairingIndex::Surface_Outside);

//     nearest_element_info_1.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_EQUAL(found_ids.size(), 3);
//     KRATOS_CHECK_EQUAL(found_ids[0], 35);
//     KRATOS_CHECK_EQUAL(found_ids[1], 18);
//     KRATOS_CHECK_EQUAL(found_ids[2], 108);

//     nearest_element_info_1.GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_EQUAL(sf_values.size(), 3);
//     KRATOS_CHECK_DOUBLE_EQUAL(sf_values[0], -0.1);
//     KRATOS_CHECK_DOUBLE_EQUAL(sf_values[1], 0.6);
//     KRATOS_CHECK_DOUBLE_EQUAL(sf_values[2], 0.5);

//     // Line Inside //
//     nearest_element_info_2.ProcessSearchResult(*interface_geom_obj, 10.0);
//     KRATOS_CHECK_IS_FALSE(nearest_element_info_2.GetLocalSearchWasSuccessful());
//     // since the no valid projection could be found we try to get an approximation
//     nearest_element_info_2.ProcessSearchResultForApproximation(*interface_geom_obj, 10.0);

//     KRATOS_CHECK(nearest_element_info_2.GetLocalSearchWasSuccessful());
//     KRATOS_CHECK(nearest_element_info_2.GetIsApproximation());

//     nearest_element_info_2.GetValue(approximation_dist, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_NEAR(approximation_dist, 1.01116, 1E-4);

//     nearest_element_info_2.GetValue(pairing_index, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_EQUAL((ProjectionUtilities::PairingIndex)pairing_index, ProjectionUtilities::PairingIndex::Line_Inside);

//     nearest_element_info_2.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_EQUAL(found_ids.size(), 2);
//     KRATOS_CHECK_EQUAL(found_ids[0], 108);
//     KRATOS_CHECK_EQUAL(found_ids[1], 35);

//     nearest_element_info_2.GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_EQUAL(sf_values.size(), 2);
//     KRATOS_CHECK_NEAR(sf_values[0], 0.615, 1E-4);
//     KRATOS_CHECK_NEAR(sf_values[1], 0.385, 1E-4);

//     // Closest Point //
//     nearest_element_info_3.ProcessSearchResult(*interface_geom_obj, 10.0);
//     KRATOS_CHECK_IS_FALSE(nearest_element_info_3.GetLocalSearchWasSuccessful());
//     // since the no valid projection could be found we try to get an approximation
//     nearest_element_info_3.ProcessSearchResultForApproximation(*interface_geom_obj, 10.0);

//     KRATOS_CHECK(nearest_element_info_3.GetLocalSearchWasSuccessful());
//     KRATOS_CHECK(nearest_element_info_3.GetIsApproximation());

//     nearest_element_info_3.GetValue(approximation_dist, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_NEAR(approximation_dist, 12.727922061, 1E-9);

//     nearest_element_info_3.GetValue(pairing_index, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_EQUAL((ProjectionUtilities::PairingIndex)pairing_index, ProjectionUtilities::PairingIndex::Closest_Point);

//     nearest_element_info_3.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_EQUAL(found_ids.size(), 1);
//     KRATOS_CHECK_EQUAL(found_ids[0], 108);

//     nearest_element_info_3.GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);
//     KRATOS_CHECK_EQUAL(sf_values.size(), 1);
//     KRATOS_CHECK_DOUBLE_EQUAL(sf_values[0], 1.0);
// }

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
    nearest_element_info.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_CHECK_EQUAL(found_ids.size(), 0);

    std::vector<double> sf_values;
    nearest_element_info.GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_CHECK_EQUAL(sf_values.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_Serialization, KratosMappingApplicationSerialTestSuite)
{
    // NOTE: In this tests "tria_2" has a valid projection, the rest is being removed to test only the serialization

    const double dist = 1.1;
    Point coords(0.3, -0.3, dist);

    std::size_t source_local_sys_idx = 123;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);

    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.0, -1.0, 0.0));

    const Geometry<NodeType>::Pointer tria_2(Kratos::make_shared<Triangle3D3<NodeType>>(node_4, node_2, node_1));

    InterfaceObject::Pointer interface_geom_obj_2(Kratos::make_shared<InterfaceGeometryObject>(tria_2.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_4->SetValue(INTERFACE_EQUATION_ID, 61);

    // Distances do not matter bcs only one projection is valid!
    nearest_element_info.ProcessSearchResult(*interface_geom_obj_2);

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
    nearest_element_info_new.GetValue(proj_dist, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_CHECK_DOUBLE_EQUAL(proj_dist, dist);

    // int pairing_index;
    // nearest_element_info_new.GetValue(pairing_index, MapperInterfaceInfo::InfoType::Dummy);
    // KRATOS_CHECK_EQUAL((ProjectionUtilities::PairingIndex)pairing_index, ProjectionUtilities::PairingIndex::Surface_Inside);

    // the nodes found are 4,2,1, aka (61,18,35)
    std::vector<int> found_ids;
    nearest_element_info_new.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_CHECK_EQUAL(found_ids.size(), 3);
    KRATOS_CHECK_EQUAL(found_ids[0], 61);
    KRATOS_CHECK_EQUAL(found_ids[1], 18);
    KRATOS_CHECK_EQUAL(found_ids[2], 35);

    std::vector<double> sf_values;
    nearest_element_info_new.GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_CHECK_EQUAL(sf_values.size(), 3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[0], 0.3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[1], 0.3);
    KRATOS_CHECK_DOUBLE_EQUAL(sf_values[2], 0.4);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    auto node_local_sys(Kratos::make_shared<Node>(8, 1.0, 2.5, -5.0));

    NearestElementLocalSystem local_sys(node_local_sys.get());

    // Computing the local system
    // this should return nothing since no InterfaceInfos are available
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType origin_ids2;
    EquationIdVectorType destination_ids;
    EquationIdVectorType destination_ids2;

    local_sys.EquationIdVectors(origin_ids, destination_ids);

    KRATOS_CHECK_EQUAL(origin_ids.size(), 0);
    KRATOS_CHECK_EQUAL(destination_ids.size(), 0);

    local_sys.CalculateLocalSystem(local_mapping_matrix, origin_ids2, destination_ids2);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size1(), 0);
    KRATOS_CHECK_EQUAL(local_mapping_matrix.size2(), 0);
    KRATOS_CHECK_EQUAL(origin_ids2.size(), 0);
    KRATOS_CHECK_EQUAL(destination_ids2.size(), 0);

    std::stringstream str_steam;
    local_sys.PairingInfo(str_steam, 4);
    KRATOS_CHECK_STRING_EQUAL(str_steam.str(),
        "NearestElementLocalSystem based on Node #8 at Coordinates 1 | 2.5 | -5");
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Line, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);

    const Geometry<NodeType>::Pointer p_geom(Kratos::make_shared<Line2D2<NodeType>>(node_1, node_2));

    MatrixResultsType exp_loc_matrix {{0.5, 0.5}};
    EqIDVectorResultsType exp_origin_ids {35, 18};

    TestNearestElementLocalSystem(exp_loc_matrix, exp_origin_ids, p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Triangle, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);

    const Geometry<NodeType>::Pointer p_geom(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));

    MatrixResultsType exp_loc_matrix {{0.5, 0.2, 0.3}};
    EqIDVectorResultsType exp_origin_ids {35, 18, 108};

    TestNearestElementLocalSystem(exp_loc_matrix, exp_origin_ids, p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Quad, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.0, 1.0, 0.0));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);
    node_4->SetValue(INTERFACE_EQUATION_ID, 95);

    const Geometry<NodeType>::Pointer p_geom(Kratos::make_shared<Quadrilateral3D4<NodeType>>(node_1, node_2, node_3, node_4));

    MatrixResultsType exp_loc_matrix {{0.35, 0.35, 0.15, 0.15}};
    EqIDVectorResultsType exp_origin_ids {35, 18, 108, 95};

    TestNearestElementLocalSystem(exp_loc_matrix, exp_origin_ids, p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Tetra, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 0.5, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.5, 0.3, 1.0));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);
    node_4->SetValue(INTERFACE_EQUATION_ID, 95);

    const Geometry<NodeType>::Pointer p_geom(Kratos::make_shared<Tetrahedra3D4<NodeType>>(node_1, node_2, node_3, node_4));

    MatrixResultsType exp_loc_matrix {{0.28, 0.28, 0.24, 0.2}};
    EqIDVectorResultsType exp_origin_ids {35, 18, 108, 95};

    TestNearestElementLocalSystem(exp_loc_matrix, exp_origin_ids, p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem_Hexa, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.0, 1.0, 0.0));
    auto node_5(Kratos::make_intrusive<NodeType>(5, 0.0, 0.0, 1.0));
    auto node_6(Kratos::make_intrusive<NodeType>(6, 1.0, 0.0, 1.0));
    auto node_7(Kratos::make_intrusive<NodeType>(7, 1.0, 1.0, 1.0));
    auto node_8(Kratos::make_intrusive<NodeType>(8, 0.0, 1.0, 1.0));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);
    node_4->SetValue(INTERFACE_EQUATION_ID, 95);
    node_5->SetValue(INTERFACE_EQUATION_ID, 12);
    node_6->SetValue(INTERFACE_EQUATION_ID, 14);
    node_7->SetValue(INTERFACE_EQUATION_ID, 19);
    node_8->SetValue(INTERFACE_EQUATION_ID, 22);

    const Geometry<NodeType>::Pointer p_geom(Kratos::make_shared<Hexahedra3D8<NodeType>>(node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8));

    MatrixResultsType exp_loc_matrix {{0.28, 0.28, 0.12, 0.12, 0.07, 0.07, 0.03, 0.03}};
    EqIDVectorResultsType exp_origin_ids {35, 18, 108, 95, 12, 14, 19, 22};

    TestNearestElementLocalSystem(exp_loc_matrix, exp_origin_ids, p_geom);
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystemWithApproximation, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, -3.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, -2.0, 0.0, 0.0));

    node_1->SetValue(INTERFACE_EQUATION_ID, 35);
    node_2->SetValue(INTERFACE_EQUATION_ID, 18);

    const Geometry<NodeType>::Pointer p_geom(Kratos::make_shared<Line2D2<NodeType>>(node_1, node_2));

    MatrixResultsType exp_loc_matrix {{1.0}};
    EqIDVectorResultsType exp_origin_ids {18};

    TestNearestElementLocalSystem(exp_loc_matrix, exp_origin_ids, p_geom);
}

}  // namespace Kratos::Testing