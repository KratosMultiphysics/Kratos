//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <limits>

// Project includes
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "testing/testing.h"
#include "custom_mappers/barycentric_mapper.h"
#include "includes/stream_serializer.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos::Testing {

typedef typename MapperLocalSystem::MatrixType MatrixType;
typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;
typedef Kratos::shared_ptr<MapperInterfaceInfo> MapperInterfaceInfoPointerType;

typedef Node NodeType;

KRATOS_TEST_CASE_IN_SUITE(BarycentricInterfaceInfo_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    const Point coords(1.0, 2.45, -23.8);

    const std::size_t source_local_sys_idx = 123;
    const std::size_t dummy_rank = 78;

    BarycentricInterfaceInfo barycentric_info(coords, source_local_sys_idx, 0, BarycentricInterpolationType::LINE);

    const auto barycentric_info_1(barycentric_info.Create());
    const auto barycentric_info_2(barycentric_info.Create(coords, source_local_sys_idx, dummy_rank));

    // Test if the "Create" function returns the correct object
    const auto& r_arg_1 = *barycentric_info_1;
    const auto& r_arg_2 = *barycentric_info_2;
    KRATOS_EXPECT_EQ(typeid(barycentric_info), typeid(r_arg_1));
    KRATOS_EXPECT_EQ(typeid(barycentric_info), typeid(r_arg_2));
}

KRATOS_TEST_CASE_IN_SUITE(BarycentricInterfaceInfo_simple_line_interpolation, KratosMappingApplicationSerialTestSuite)
{
    Point coords(0.4, 0.0, 0.0);

    std::size_t source_local_sys_idx = 123;

    BarycentricInterfaceInfo barycentric_info(coords, source_local_sys_idx, 0, BarycentricInterpolationType::LINE);

    auto node_1(Kratos::make_intrusive<NodeType>(1,  3.3, 0.0, 0.0)); // third closest (not used)
    auto node_2(Kratos::make_intrusive<NodeType>(3,  1.0, 0.1, -0.2)); // second closest
    auto node_3(Kratos::make_intrusive<NodeType>(15, 0.3, 0.0, 0.0)); // closest

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 13);
    node_2->SetValue(INTERFACE_EQUATION_ID, 5);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);

    barycentric_info.ProcessSearchResult(*interface_node_1);
    barycentric_info.ProcessSearchResult(*interface_node_2);
    barycentric_info.ProcessSearchResult(*interface_node_3);

    KRATOS_EXPECT_TRUE(barycentric_info.GetLocalSearchWasSuccessful());
    KRATOS_EXPECT_FALSE(barycentric_info.GetIsApproximation());

    ClosestPointsContainer exp_closest_points(2);
    exp_closest_points.Add(PointWithId(108, Point(0.3, 0.0, 0.0), 0.1));
    exp_closest_points.Add(PointWithId(5, Point(1.0, 0.1, -0.2), coords.Distance(*node_2)));

    KRATOS_EXPECT_EQ(barycentric_info.GetClosestPoints(), exp_closest_points);
}

KRATOS_TEST_CASE_IN_SUITE(BarycentricInterfaceInfo_line_only_one_point, KratosMappingApplicationSerialTestSuite)
{
    Point coords(0.4, 0.0, 0.0);

    std::size_t source_local_sys_idx = 123;

    BarycentricInterfaceInfo barycentric_info(coords, source_local_sys_idx, 0, BarycentricInterpolationType::LINE);

    auto node_1(Kratos::make_intrusive<NodeType>(1,  3.3, 0.0, 0.0));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 13);

    barycentric_info.ProcessSearchResult(*interface_node_1);

    KRATOS_EXPECT_TRUE(barycentric_info.GetLocalSearchWasSuccessful());
    KRATOS_EXPECT_TRUE(barycentric_info.GetIsApproximation());

    ClosestPointsContainer exp_closest_points(2);
    exp_closest_points.Add(PointWithId(13, Point(3.3, 0.0, 0.0), 2.9));

    KRATOS_EXPECT_EQ(barycentric_info.GetClosestPoints(), exp_closest_points);
}

KRATOS_TEST_CASE_IN_SUITE(BarycentricInterfaceInfo_line_duplicated_point, KratosMappingApplicationSerialTestSuite)
{
    // test to make sure that if several points share the same coords only one of them is used!

    Point coords(0.4, 0.0, 0.0);

    std::size_t source_local_sys_idx = 123;

    BarycentricInterfaceInfo barycentric_info(coords, source_local_sys_idx, 0, BarycentricInterpolationType::LINE);

    auto node_1(Kratos::make_intrusive<NodeType>(1,  3.3, 0.0, 0.0)); // third closest (not used)
    auto node_2(Kratos::make_intrusive<NodeType>(3,  1.0, 0.1, -0.2)); // second closest
    auto node_3(Kratos::make_intrusive<NodeType>(15, 0.3, 0.0, 0.0)); // closest
    auto node_4(Kratos::make_intrusive<NodeType>(16, 0.3, 0.0, 0.0)); // closest (but not used bcs duplicate)

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3.get()));
    InterfaceObject::Pointer interface_node_4(Kratos::make_shared<InterfaceNode>(node_4.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 13);
    node_2->SetValue(INTERFACE_EQUATION_ID, 5);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);
    node_4->SetValue(INTERFACE_EQUATION_ID, 32);

    barycentric_info.ProcessSearchResult(*interface_node_1);
    barycentric_info.ProcessSearchResult(*interface_node_2);
    barycentric_info.ProcessSearchResult(*interface_node_3);
    barycentric_info.ProcessSearchResult(*interface_node_4);

    KRATOS_EXPECT_TRUE(barycentric_info.GetLocalSearchWasSuccessful());
    KRATOS_EXPECT_FALSE(barycentric_info.GetIsApproximation());

    ClosestPointsContainer exp_closest_points(2);
    exp_closest_points.Add(PointWithId(108, Point(0.3, 0.0, 0.0), 0.1));
    exp_closest_points.Add(PointWithId(5, Point(1.0, 0.1, -0.2), coords.Distance(*node_2)));

    KRATOS_EXPECT_EQ(barycentric_info.GetClosestPoints(), exp_closest_points);
}

KRATOS_TEST_CASE_IN_SUITE(BarycentricInterfaceInfo_Serialization, KratosMappingApplicationSerialTestSuite)
{
    Point coords(0.4, 0.0, 0.0);

    std::size_t source_local_sys_idx = 123;

    BarycentricInterfaceInfo barycentric_info(coords, source_local_sys_idx, 0, BarycentricInterpolationType::LINE);

    auto node_1(Kratos::make_intrusive<NodeType>(1,  3.3, 0.0, 0.0)); // third closest (not used)
    auto node_2(Kratos::make_intrusive<NodeType>(3,  1.0, 0.1, -0.2)); // second closest
    auto node_3(Kratos::make_intrusive<NodeType>(15, 0.3, 0.0, 0.0)); // closest

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 13);
    node_2->SetValue(INTERFACE_EQUATION_ID, 5);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);

    barycentric_info.ProcessSearchResult(*interface_node_1);
    barycentric_info.ProcessSearchResult(*interface_node_2);
    barycentric_info.ProcessSearchResult(*interface_node_3);

    KRATOS_EXPECT_TRUE(barycentric_info.GetLocalSearchWasSuccessful());
    KRATOS_EXPECT_FALSE(barycentric_info.GetIsApproximation());

    ClosestPointsContainer exp_closest_points(2);
    exp_closest_points.Add(PointWithId(108, Point(0.3, 0.0, 0.0), 0.1));
    exp_closest_points.Add(PointWithId(5, Point(1.0, 0.1, -0.2), coords.Distance(*node_2)));

    KRATOS_EXPECT_EQ(barycentric_info.GetClosestPoints(), exp_closest_points);

    // serializing the object
    StreamSerializer serializer;
    serializer.save("barycentric_interface_info", barycentric_info);
    // deserializing the object => this happens if the remote search was successful and
    // sending back of the object to the partition where it came from is required
    BarycentricInterfaceInfo barycentric_info_new(BarycentricInterpolationType::LINE);
    serializer.load("barycentric_interface_info", barycentric_info_new);

    KRATOS_EXPECT_EQ(barycentric_info_new.GetLocalSystemIndex(), source_local_sys_idx);

    KRATOS_EXPECT_EQ(barycentric_info_new.GetClosestPoints(), exp_closest_points);
}

KRATOS_TEST_CASE_IN_SUITE(BarycentricLocalSystem_simple_line_interpolation, KratosMappingApplicationSerialTestSuite)
{
    Point coords(0.4, 0.0, 0.0);
    auto node_orig(Kratos::make_intrusive<NodeType>(1, coords[0], coords[1], coords[2]));
    node_orig->SetValue(INTERFACE_EQUATION_ID, 8);

    std::size_t source_local_sys_idx = 123;

    MapperInterfaceInfoPointerType p_interface_info(Kratos::make_shared<BarycentricInterfaceInfo>(coords, source_local_sys_idx, 0, BarycentricInterpolationType::LINE));

    auto node_1(Kratos::make_intrusive<NodeType>(3,  1.0, 0.0, 0.0)); // second closest
    auto node_2(Kratos::make_intrusive<NodeType>(15, 0.0, 0.0, 0.0)); // closest

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 5);
    node_2->SetValue(INTERFACE_EQUATION_ID, 13);

    p_interface_info->ProcessSearchResult(*interface_node_1);
    p_interface_info->ProcessSearchResult(*interface_node_2);

    KRATOS_EXPECT_TRUE(p_interface_info->GetLocalSearchWasSuccessful());
    KRATOS_EXPECT_FALSE(p_interface_info->GetIsApproximation());

    BarycentricLocalSystem local_sys(node_orig.get());
    local_sys.AddInterfaceInfo(p_interface_info);

    MatrixType exp_matrix;
    exp_matrix.resize(1,2);
    exp_matrix(0,0) = 0.6;
    exp_matrix(0,1) = 0.4;

    const EquationIdVectorType exp_origin_ids{13, 5};
    const int exp_destination_id = 8;

    // Computing the local system
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType origin_ids2;
    EquationIdVectorType destination_ids;
    EquationIdVectorType destination_ids2;

    local_sys.EquationIdVectors(origin_ids, destination_ids);

    KRATOS_EXPECT_EQ(origin_ids.size(), exp_origin_ids.size());
    for (std::size_t i=0; i<exp_origin_ids.size(); ++i) {
        KRATOS_EXPECT_EQ(origin_ids[i], exp_origin_ids[i]);
    }
    KRATOS_EXPECT_EQ(destination_ids.size(), 1);
    KRATOS_EXPECT_EQ(destination_ids[0], exp_destination_id);

    local_sys.CalculateLocalSystem(local_mapping_matrix, origin_ids2, destination_ids2);

    KRATOS_EXPECT_EQ(local_mapping_matrix.size1(), 1);
    KRATOS_EXPECT_EQ(local_mapping_matrix.size2(), exp_origin_ids.size());
    KRATOS_EXPECT_EQ(origin_ids2.size(), exp_origin_ids.size());
    KRATOS_EXPECT_EQ(destination_ids2.size(), 1);

    KRATOS_EXPECT_MATRIX_NEAR(local_mapping_matrix, exp_matrix, 1e-14)

    for (std::size_t i=0; i<exp_origin_ids.size(); ++i) {
        KRATOS_EXPECT_EQ(origin_ids[i], exp_origin_ids[i]);
    }
    KRATOS_EXPECT_EQ(destination_ids2[0], exp_destination_id);
}

KRATOS_TEST_CASE_IN_SUITE(BarycentricInterfaceInfo_simple_triangle_interpolation, KratosMappingApplicationSerialTestSuite)
{
    Point coords(0.2, 0.2, 0.0);

    std::size_t source_local_sys_idx = 123;

    BarycentricInterfaceInfo barycentric_info(coords, source_local_sys_idx, 0, BarycentricInterpolationType::TRIANGLE);

    auto node_1(Kratos::make_intrusive<NodeType>(1,  0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(3,  1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(15, 0.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(18, 1.0, 1.0, 0.0));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3.get()));
    InterfaceObject::Pointer interface_node_4(Kratos::make_shared<InterfaceNode>(node_4.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 13);
    node_2->SetValue(INTERFACE_EQUATION_ID, 5);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);
    node_4->SetValue(INTERFACE_EQUATION_ID, 41);

    barycentric_info.ProcessSearchResult(*interface_node_1);
    barycentric_info.ProcessSearchResult(*interface_node_2);
    barycentric_info.ProcessSearchResult(*interface_node_3);
    barycentric_info.ProcessSearchResult(*interface_node_4);

    KRATOS_EXPECT_TRUE(barycentric_info.GetLocalSearchWasSuccessful());
    KRATOS_EXPECT_TRUE(barycentric_info.GetIsApproximation());

    KRATOS_SKIP_TEST << "this test needs some updates" << std::endl;

    std::vector<int> found_ids;
    barycentric_info.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_ids.size(), 3);
    KRATOS_EXPECT_EQ(found_ids[0], 13);
    KRATOS_EXPECT_EQ(found_ids[1], 5);
    KRATOS_EXPECT_EQ(found_ids[2], 108);

    std::vector<double> neighbor_coords;
    barycentric_info.GetValue(neighbor_coords, MapperInterfaceInfo::InfoType::Dummy);
    const std::vector<double> exp_results {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0
    };
    KRATOS_EXPECT_VECTOR_EQ(exp_results, neighbor_coords)
}

KRATOS_TEST_CASE_IN_SUITE(BarycentricInterfaceInfo_triangle_collinear_nodes_interpolation, KratosMappingApplicationSerialTestSuite)
{
    Point coords(0.2, 0.2, 0.0);

    std::size_t source_local_sys_idx = 123;

    BarycentricInterfaceInfo barycentric_info(coords, source_local_sys_idx, 0, BarycentricInterpolationType::TRIANGLE);

    auto node_1(Kratos::make_intrusive<NodeType>(1,  0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(3,  1.0, 0.0, 0.0)); // this is closer than the next one but collinear and hence forbidden
    auto node_3(Kratos::make_intrusive<NodeType>(15, 0.0, 1.1, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(18, 1.0, 1.0, 0.0));
    auto node_5(Kratos::make_intrusive<NodeType>(21, 0.7, 0.0, 0.0));
    auto node_6(Kratos::make_intrusive<NodeType>(22, 0.5, 0.0, 0.0));
    auto node_7(Kratos::make_intrusive<NodeType>(23, 0.75, 0.0, 0.0));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3.get()));
    InterfaceObject::Pointer interface_node_4(Kratos::make_shared<InterfaceNode>(node_4.get()));
    InterfaceObject::Pointer interface_node_5(Kratos::make_shared<InterfaceNode>(node_5.get()));
    InterfaceObject::Pointer interface_node_6(Kratos::make_shared<InterfaceNode>(node_6.get()));
    InterfaceObject::Pointer interface_node_7(Kratos::make_shared<InterfaceNode>(node_7.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 13);
    node_2->SetValue(INTERFACE_EQUATION_ID, 5);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);
    node_4->SetValue(INTERFACE_EQUATION_ID, 41);
    node_5->SetValue(INTERFACE_EQUATION_ID, 63);
    node_6->SetValue(INTERFACE_EQUATION_ID, 64);
    node_7->SetValue(INTERFACE_EQUATION_ID, 65);

    barycentric_info.ProcessSearchResult(*interface_node_1);
    barycentric_info.ProcessSearchResult(*interface_node_2);
    barycentric_info.ProcessSearchResult(*interface_node_3);
    barycentric_info.ProcessSearchResult(*interface_node_4);
    barycentric_info.ProcessSearchResult(*interface_node_5);
    barycentric_info.ProcessSearchResult(*interface_node_6);
    barycentric_info.ProcessSearchResult(*interface_node_7);

    KRATOS_EXPECT_TRUE(barycentric_info.GetLocalSearchWasSuccessful());
    KRATOS_EXPECT_TRUE(barycentric_info.GetIsApproximation());

    KRATOS_SKIP_TEST << "this test needs some updates" << std::endl;

    std::vector<int> found_ids;
    barycentric_info.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_ids.size(), 3);
    KRATOS_EXPECT_EQ(found_ids[0], 13);
    KRATOS_EXPECT_EQ(found_ids[1], 64);
    KRATOS_EXPECT_EQ(found_ids[2], 108);

    std::vector<double> neighbor_coords;
    barycentric_info.GetValue(neighbor_coords, MapperInterfaceInfo::InfoType::Dummy);
    const std::vector<double> exp_results {
        0.0, 0.0, 0.0,
        0.5, 0.0, 0.0,
        0.0, 1.1, 0.0
    };
    KRATOS_EXPECT_VECTOR_EQ(exp_results, neighbor_coords)
}

KRATOS_TEST_CASE_IN_SUITE(BarycentricInterfaceInfo_triangle_only_collinear_nodes_approx, KratosMappingApplicationSerialTestSuite)
{
    Point coords(0.2, 0.2, 0.0);

    std::size_t source_local_sys_idx = 123;

    BarycentricInterfaceInfo barycentric_info(coords, source_local_sys_idx, 0, BarycentricInterpolationType::TRIANGLE);

    auto node_1(Kratos::make_intrusive<NodeType>(1,  0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(3,  1.0, 0.0, 0.0)); // this should NOT be saved since it is collinear!
    auto node_3(Kratos::make_intrusive<NodeType>(22, 0.5, 0.0, 0.0));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 13);
    node_2->SetValue(INTERFACE_EQUATION_ID, 5);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);

    barycentric_info.ProcessSearchResult(*interface_node_1);
    barycentric_info.ProcessSearchResult(*interface_node_2);
    barycentric_info.ProcessSearchResult(*interface_node_3);

    KRATOS_EXPECT_TRUE(barycentric_info.GetLocalSearchWasSuccessful());
    KRATOS_EXPECT_TRUE(barycentric_info.GetIsApproximation());

    KRATOS_SKIP_TEST << "this test needs some updates" << std::endl;

    std::vector<int> found_ids;
    barycentric_info.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_ids.size(), 3);
    KRATOS_EXPECT_EQ(found_ids[0], 13);
    KRATOS_EXPECT_EQ(found_ids[1], 108);
    KRATOS_EXPECT_EQ(found_ids[2], -1);

    std::vector<double> neighbor_coords;
    barycentric_info.GetValue(neighbor_coords, MapperInterfaceInfo::InfoType::Dummy);
    const std::vector<double> exp_results {
        0.0, 0.0, 0.0,
        0.5, 0.0, 0.0,
        std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()
    };
    KRATOS_EXPECT_VECTOR_EQ(exp_results, neighbor_coords)
}

KRATOS_TEST_CASE_IN_SUITE(BarycentricInterfaceInfo_simple_tetra_interpolation, KratosMappingApplicationSerialTestSuite)
{
    Point coords(0.3, 0.2, 0.0);

    std::size_t source_local_sys_idx = 123;

    BarycentricInterfaceInfo barycentric_info(coords, source_local_sys_idx, 0, BarycentricInterpolationType::TETRAHEDRA);

    auto node_1(Kratos::make_intrusive<NodeType>(1,  0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(3,  1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(15, 0.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(18, 0.0, 0.0, 1.1));

    InterfaceObject::Pointer interface_node_1(Kratos::make_shared<InterfaceNode>(node_1.get()));
    InterfaceObject::Pointer interface_node_2(Kratos::make_shared<InterfaceNode>(node_2.get()));
    InterfaceObject::Pointer interface_node_3(Kratos::make_shared<InterfaceNode>(node_3.get()));
    InterfaceObject::Pointer interface_node_4(Kratos::make_shared<InterfaceNode>(node_4.get()));

    node_1->SetValue(INTERFACE_EQUATION_ID, 13);
    node_2->SetValue(INTERFACE_EQUATION_ID, 5);
    node_3->SetValue(INTERFACE_EQUATION_ID, 108);
    node_4->SetValue(INTERFACE_EQUATION_ID, 41);

    barycentric_info.ProcessSearchResult(*interface_node_1);
    barycentric_info.ProcessSearchResult(*interface_node_2);
    barycentric_info.ProcessSearchResult(*interface_node_3);
    barycentric_info.ProcessSearchResult(*interface_node_4);

    KRATOS_EXPECT_TRUE(barycentric_info.GetLocalSearchWasSuccessful());
    KRATOS_EXPECT_TRUE(barycentric_info.GetIsApproximation());

    KRATOS_SKIP_TEST << "this test needs some updates" << std::endl;

    std::vector<int> found_ids;
    barycentric_info.GetValue(found_ids, MapperInterfaceInfo::InfoType::Dummy);
    KRATOS_EXPECT_EQ(found_ids.size(), 4);
    KRATOS_EXPECT_EQ(found_ids[0], 13);
    KRATOS_EXPECT_EQ(found_ids[1], 5);
    KRATOS_EXPECT_EQ(found_ids[2], 108);
    KRATOS_EXPECT_EQ(found_ids[3], 41);

    std::vector<double> neighbor_coords;
    barycentric_info.GetValue(neighbor_coords, MapperInterfaceInfo::InfoType::Dummy);
    const std::vector<double> exp_results {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.1
    };
    KRATOS_EXPECT_VECTOR_EQ(exp_results, neighbor_coords)
}

}  // namespace Kratos::Testing