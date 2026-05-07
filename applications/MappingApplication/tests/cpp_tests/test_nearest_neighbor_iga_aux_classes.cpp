//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti
//

// Project includes
#include "includes/serializer.h"
#include "testing/testing.h"
#include "custom_mappers/nearest_neighbor_mapper_iga.h"
#include "includes/stream_serializer.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos::Testing {

    typedef typename MapperLocalSystem::MatrixType MatrixType;
    typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;

    typedef Node NodeType;
    typedef Geometry<Node> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    /** Generates a point type sample Triangle2D3 to serve as the quadrature point geometry parent.
         * @return  Pointer to a Triangle2D3
         */
    GeometryPointerType GeneratePointsTriangle2D3_Element(
        const NodeType::Pointer& p_node_1,
        const NodeType::Pointer& p_node_2,
        const NodeType::Pointer& p_node_3)
    {
        return Kratos::make_shared<Triangle2D3<NodeType>>(p_node_1, p_node_2, p_node_3);
    }

    KRATOS_TEST_CASE_IN_SUITE(MapperInterfaceInfoIGA_BasicTests, KratosMappingApplicationSerialTestSuite)
    {
        // This test covers the basic functionalities provided by the "MapperInterfaceInfo"
        // A "NearestNeighborInterfaceInfo" is being used since "MapperInterfaceInfo" is a pure virtual class

        Point coords_1(1.0, 2.45, -23.8);

        std::size_t source_local_sys_idx = 55;
        std::size_t source_rank = 23;

        NearestNeighborInterfaceInfoIGA nearest_neighbor_iga_info(coords_1, source_local_sys_idx, source_rank);
        const auto nearest_neighbor_iga_info_2(nearest_neighbor_iga_info.Create(coords_1, source_local_sys_idx, 0));

        KRATOS_EXPECT_EQ(nearest_neighbor_iga_info.GetLocalSystemIndex(), source_local_sys_idx);

        KRATOS_EXPECT_EQ(nearest_neighbor_iga_info.GetSourceRank(), source_rank);
        KRATOS_EXPECT_EQ(nearest_neighbor_iga_info_2->GetSourceRank(), 0); // Testing against default

        KRATOS_EXPECT_FALSE(nearest_neighbor_iga_info.GetLocalSearchWasSuccessful());
        KRATOS_EXPECT_FALSE(nearest_neighbor_iga_info.GetIsApproximation());

        for (std::size_t i=0; i<3; ++i)
            KRATOS_EXPECT_DOUBLE_EQ(nearest_neighbor_iga_info.Coordinates()[i], coords_1[i]);
    }

    KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfoIGA_BasicTests, KratosMappingApplicationSerialTestSuite)
    {
        const Point coords(1.0, 2.45, -23.8);

        const std::size_t source_local_sys_idx = 123;
        const std::size_t dummy_rank = 78;

        NearestNeighborInterfaceInfoIGA nearest_neighbor_iga_info(coords, source_local_sys_idx, 0);

        const auto nearest_neighbor_iga_info_1(nearest_neighbor_iga_info.Create());
        const auto nearest_neighbor_iga_info_2(nearest_neighbor_iga_info.Create(coords, source_local_sys_idx, dummy_rank));

        // Test if the "Create" function returns the correct object
        const auto& r_arg_1 = *nearest_neighbor_iga_info_1;
        const auto& r_arg_2 = *nearest_neighbor_iga_info_2;
        KRATOS_EXPECT_EQ(typeid(nearest_neighbor_iga_info), typeid(r_arg_1));
        KRATOS_EXPECT_EQ(typeid(nearest_neighbor_iga_info), typeid(r_arg_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfoIGA_NeighborsFound, KratosMappingApplicationSerialTestSuite)
    {
        // NOTE: In this tests "qp_geom_2" is the closest quadrature point

        Point coords(1.5, 0.5, 0.1);

        std::size_t source_local_sys_idx = 123;

        NearestNeighborInterfaceInfoIGA nearest_neighbor_iga_info(coords, source_local_sys_idx, 0);

        // Create 2 triangular elements to serve as the quadrature point geometry parent
        auto p_node_1 = Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0);
        auto p_node_2 = Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0);
        auto p_node_3 = Kratos::make_intrusive<NodeType>(3, 0.0, 1.0, 0.0);
        auto p_node_4 = Kratos::make_intrusive<NodeType>(4, 1.0, 1.0, 0.0);
        auto element_1 = GeneratePointsTriangle2D3_Element(p_node_1, p_node_2, p_node_3);
        auto element_2 = GeneratePointsTriangle2D3_Element(p_node_3, p_node_2, p_node_4);

        // Obtain the element shape functions
        auto N_1 = element_1->ShapeFunctionsValues();
        auto N_2 = element_2->ShapeFunctionsValues();

        Matrix N_1_mat(1, 3), N_2_mat(1, 3);
        for (std::size_t i = 0; i < 3; ++i){
            N_1_mat(0, i) = N_1(0, i);
            N_2_mat(0, i) = N_2(0, i);
        }

        // Obtain the element derivatives
        Matrix DN_De_1 = element_1->ShapeFunctionLocalGradient(0);
        Matrix DN_De_2 = element_2->ShapeFunctionLocalGradient(0);

        // Create the shape functions container
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> container_1(
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            element_1->IntegrationPoints()[0],
            N_1_mat,
            DN_De_1);
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> container_2(
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            element_2->IntegrationPoints()[0],
            N_2_mat,
            DN_De_2);

        // Create the quadrature point geometries for the mapper
        auto qp_geom_1 = Kratos::make_shared<QuadraturePointGeometry<Node, 2, 2>>(
            element_1->Points(), container_1, element_1.get());
        auto qp_geom_2 = Kratos::make_shared<QuadraturePointGeometry<Node, 2, 2>>(
            element_2->Points(), container_2, element_2.get());

        InterfaceObject::Pointer interface_quadrature_point_1(Kratos::make_shared<InterfaceGeometryObject>((qp_geom_1.get())));
        InterfaceObject::Pointer interface_quadrature_point_2(Kratos::make_shared<InterfaceGeometryObject>((qp_geom_2.get())));

        const std::vector<int> expected_ids = {2, 1, 3};
        
        // Set the INTERFACE_EQUATION_ID for the nodes
        p_node_1->SetValue(INTERFACE_EQUATION_ID, 0);
        p_node_2->SetValue(INTERFACE_EQUATION_ID, 1);
        p_node_3->SetValue(INTERFACE_EQUATION_ID, 2);
        p_node_4->SetValue(INTERFACE_EQUATION_ID, 3);

        // We compute the real distance bcs this would also be computed by the search
        const double dist_2 = coords.Distance(*interface_quadrature_point_2);

        nearest_neighbor_iga_info.ProcessSearchResult(*interface_quadrature_point_1);
        nearest_neighbor_iga_info.ProcessSearchResult(*interface_quadrature_point_2);

        KRATOS_EXPECT_TRUE(nearest_neighbor_iga_info.GetLocalSearchWasSuccessful());
        // this function should never return true for this class!
        // It already works with nodes, which could be approximations for other mappers
        KRATOS_EXPECT_FALSE(nearest_neighbor_iga_info.GetIsApproximation());

        std::vector<int> found_id(1);
        nearest_neighbor_iga_info.GetValue(found_id, MapperInterfaceInfo::InfoType::Dummy);
        KRATOS_EXPECT_VECTOR_EQ(found_id, expected_ids);

        double neighbor_dist;
        nearest_neighbor_iga_info.GetValue(neighbor_dist, MapperInterfaceInfo::InfoType::Dummy);
        KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_2);
    }


    KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfoIGA_SameDistanceFound, KratosMappingApplicationSerialTestSuite)
    {
        Point coords(0.5, 0.5, 0.1);

        std::size_t source_local_sys_idx = 123;

        NearestNeighborInterfaceInfoIGA nearest_neighbor_iga_info(coords, source_local_sys_idx, 0);

       // Create 2 triangular elements to serve as the quadrature point geometry parent
        auto p_node_1 = Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0);
        auto p_node_2 = Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0);
        auto p_node_3 = Kratos::make_intrusive<NodeType>(3, 0.0, 1.0, 0.0);
        auto p_node_4 = Kratos::make_intrusive<NodeType>(4, 1.0, 1.0, 0.0);
        auto element_1 = GeneratePointsTriangle2D3_Element(p_node_1, p_node_2, p_node_3);
        auto element_2 = GeneratePointsTriangle2D3_Element(p_node_3, p_node_2, p_node_4);

        // Obtain the element shape functions
        auto N_1 = element_1->ShapeFunctionsValues();
        auto N_2 = element_2->ShapeFunctionsValues();

        Matrix N_1_mat(1, 3), N_2_mat(1, 3);
        for (std::size_t i = 0; i < 3; ++i){
            N_1_mat(0, i) = N_1(0, i);
            N_2_mat(0, i) = N_2(0, i);
        }

        // Obtain the element derivatives
        Matrix DN_De_1 = element_1->ShapeFunctionLocalGradient(0);
        Matrix DN_De_2 = element_2->ShapeFunctionLocalGradient(0);

        // Create the shape functions container
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> container_1(
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            element_1->IntegrationPoints()[0],
            N_1_mat,
            DN_De_1);
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> container_2(
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            element_2->IntegrationPoints()[0],
            N_2_mat,
            DN_De_2);

        // Create the quadrature point geometries for the mapper
        auto qp_geom_1 = Kratos::make_shared<QuadraturePointGeometry<Node, 2, 2>>(
            element_1->Points(), container_1, element_1.get());
        auto qp_geom_2 = Kratos::make_shared<QuadraturePointGeometry<Node, 2, 2>>(
            element_2->Points(), container_2, element_2.get());
        
        InterfaceObject::Pointer interface_quadrature_point_1(Kratos::make_shared<InterfaceGeometryObject>((qp_geom_1.get())));
        InterfaceObject::Pointer interface_quadrature_point_2(Kratos::make_shared<InterfaceGeometryObject>((qp_geom_2.get())));

        const std::vector<int> expected_ids = {0, 1, 2, 2, 1, 3};

        // Set the INTERFACE_EQUATION_ID for the nodes
        p_node_1->SetValue(INTERFACE_EQUATION_ID, 0);
        p_node_2->SetValue(INTERFACE_EQUATION_ID, 1);
        p_node_3->SetValue(INTERFACE_EQUATION_ID, 2);
        p_node_4->SetValue(INTERFACE_EQUATION_ID, 3);

        // We compute the real distance bcs this would also be computed by the search
        const double dist_1 = coords.Distance(*interface_quadrature_point_1);
        const double dist_2 = coords.Distance(*interface_quadrature_point_2);
        KRATOS_EXPECT_DOUBLE_EQ(dist_1, dist_2);

        KRATOS_EXPECT_FALSE(nearest_neighbor_iga_info.GetLocalSearchWasSuccessful()); // this is the default

        nearest_neighbor_iga_info.ProcessSearchResult(*interface_quadrature_point_1);
        nearest_neighbor_iga_info.ProcessSearchResult(*interface_quadrature_point_2);

        KRATOS_EXPECT_TRUE(nearest_neighbor_iga_info.GetLocalSearchWasSuccessful());

        std::vector<int> found_id(2);
        nearest_neighbor_iga_info.GetValue(found_id, MapperInterfaceInfo::InfoType::Dummy);
        KRATOS_EXPECT_VECTOR_EQ(found_id, expected_ids);

        double neighbor_dist;
        nearest_neighbor_iga_info.GetValue(neighbor_dist, MapperInterfaceInfo::InfoType::Dummy);
        KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_1);
        KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_2);
    }

    KRATOS_TEST_CASE_IN_SUITE(NearestNeighborInterfaceInfoIGA_Serialization, KratosMappingApplicationSerialTestSuite)
    {
        // NOTE: In this tests "qp_geom_2" is the closest quadrature point

        Point coords(1.5, 0.5, 0.1);

        std::size_t source_local_sys_idx = 123;

        NearestNeighborInterfaceInfoIGA nearest_neighbor_iga_info(coords, source_local_sys_idx, 0);

        // Create 2 triangular elements to serve as the quadrature point geometry parent
        auto p_node_1 = Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0);
        auto p_node_2 = Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0);
        auto p_node_3 = Kratos::make_intrusive<NodeType>(3, 0.0, 1.0, 0.0);
        auto p_node_4 = Kratos::make_intrusive<NodeType>(4, 1.0, 1.0, 0.0);
        auto element_1 = GeneratePointsTriangle2D3_Element(p_node_1, p_node_2, p_node_3);
        auto element_2 = GeneratePointsTriangle2D3_Element(p_node_3, p_node_2, p_node_4);

        // Obtain the element shape functions
        auto N_1 = element_1->ShapeFunctionsValues();
        auto N_2 = element_2->ShapeFunctionsValues();

        Matrix N_1_mat(1, 3), N_2_mat(1, 3);
        for (std::size_t i = 0; i < 3; ++i){
            N_1_mat(0, i) = N_1(0, i);
            N_2_mat(0, i) = N_2(0, i);
        }

        // Obtain the element derivatives
        Matrix DN_De_1 = element_1->ShapeFunctionLocalGradient(0);
        Matrix DN_De_2 = element_2->ShapeFunctionLocalGradient(0);

        // Create the shape functions container
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> container_1(
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            element_1->IntegrationPoints()[0],
            N_1_mat,
            DN_De_1);
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> container_2(
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            element_2->IntegrationPoints()[0],
            N_2_mat,
            DN_De_2);

        // Create the quadrature point geometries for the mapper
        auto qp_geom_1 = Kratos::make_shared<QuadraturePointGeometry<Node, 2, 2>>(
            element_1->Points(), container_1, element_1.get());
        auto qp_geom_2 = Kratos::make_shared<QuadraturePointGeometry<Node, 2, 2>>(
            element_2->Points(), container_2, element_2.get());
        
        InterfaceObject::Pointer interface_quadrature_point_1(Kratos::make_shared<InterfaceGeometryObject>((qp_geom_1.get())));
        InterfaceObject::Pointer interface_quadrature_point_2(Kratos::make_shared<InterfaceGeometryObject>((qp_geom_2.get())));

        const std::vector<int> expected_ids = {2, 1, 3};

        // Set the INTERFACE_EQUATION_ID for the nodes
        p_node_1->SetValue(INTERFACE_EQUATION_ID, 0);
        p_node_2->SetValue(INTERFACE_EQUATION_ID, 1);
        p_node_3->SetValue(INTERFACE_EQUATION_ID, 2);
        p_node_4->SetValue(INTERFACE_EQUATION_ID, 3);


        // We compute the real distance bcs this would also be computed by the search
        const double dist_2 = coords.Distance(*interface_quadrature_point_2);

        nearest_neighbor_iga_info.ProcessSearchResult(*interface_quadrature_point_1);
        nearest_neighbor_iga_info.ProcessSearchResult(*interface_quadrature_point_2);

        // serializing the object
        StreamSerializer serializer;
        serializer.save("nearest_neighbor_interface_info", nearest_neighbor_iga_info);
        // deserializing the object => this happens if the remote search was successful and
        // sending back of the object to the partition where it came from is required
        NearestNeighborInterfaceInfoIGA nearest_neighbor_iga_info_new;
        serializer.load("nearest_neighbor_interface_info", nearest_neighbor_iga_info_new);

        KRATOS_EXPECT_EQ(nearest_neighbor_iga_info_new.GetLocalSystemIndex(), source_local_sys_idx);

        std::vector<int> found_id(1);
        nearest_neighbor_iga_info_new.GetValue(found_id, MapperInterfaceInfo::InfoType::Dummy);
        KRATOS_EXPECT_VECTOR_EQ(found_id, expected_ids);

        double neighbor_dist;
        nearest_neighbor_iga_info_new.GetValue(neighbor_dist, MapperInterfaceInfo::InfoType::Dummy);
        KRATOS_EXPECT_DOUBLE_EQ(neighbor_dist, dist_2);
    }

    KRATOS_TEST_CASE_IN_SUITE(MapperLocalSystemIGA_BasicTests, KratosMappingApplicationSerialTestSuite)
    {
        // This test covers the basic functionalities provided by the "MapperLocalSystem"
        // A "NearestNeighborLocalSystem" is being used since "MapperLocalSystem" is a pure virtual class

        Point coords_1(1.5, 0.5, 0.1);
        auto node_local_sys(Kratos::make_intrusive<NodeType>(5, coords_1));

        NearestNeighborLocalSystemIGA local_sys(node_local_sys.get());

        for (std::size_t i=0; i<3; ++i)
            KRATOS_EXPECT_DOUBLE_EQ(local_sys.Coordinates()[i], coords_1[i]);

        KRATOS_EXPECT_FALSE(local_sys.HasInterfaceInfo());
    }

    KRATOS_TEST_CASE_IN_SUITE(NearestNeighborLocalSystemIGA_BasicTests, KratosMappingApplicationSerialTestSuite)
    {
        auto node_local_sys(Kratos::make_intrusive<NodeType>(8, 1.5, 0.5, 0.1));

        NearestNeighborLocalSystemIGA local_sys(node_local_sys.get());

        // Computing the local system
        // this should return nothing since no InterfaceInfos are available
        MatrixType local_mapping_matrix;
        EquationIdVectorType origin_ids;
        EquationIdVectorType origin_ids2;
        EquationIdVectorType destination_ids;
        EquationIdVectorType destination_ids2;

        local_sys.EquationIdVectors(origin_ids, destination_ids);

        KRATOS_EXPECT_EQ(origin_ids.size(), 0);
        KRATOS_EXPECT_EQ(destination_ids.size(), 0);

        local_sys.CalculateLocalSystem(local_mapping_matrix, origin_ids2, destination_ids2);
        KRATOS_EXPECT_EQ(local_mapping_matrix.size1(), 0);
        KRATOS_EXPECT_EQ(local_mapping_matrix.size2(), 0);
        KRATOS_EXPECT_EQ(origin_ids2.size(), 0);
        KRATOS_EXPECT_EQ(destination_ids2.size(), 0);

        std::stringstream str_steam;
        local_sys.PairingInfo(str_steam, 4);
        KRATOS_EXPECT_EQ(str_steam.str(),
            "NearestNeighborLocalSystem based on Node #8 at Coordinates 1.5 | 0.5 | 0.1");
    }

    KRATOS_TEST_CASE_IN_SUITE(NearestNeighborLocalSystemIGA_ComputeLocalSystem, KratosMappingApplicationSerialTestSuite)
    {
        const int dest_id = 13;

        auto node_local_sys(Kratos::make_intrusive<NodeType>(5, 1.5, 0.5, 0.1));
        node_local_sys->SetValue(INTERFACE_EQUATION_ID, dest_id);

        NearestNeighborLocalSystemIGA local_sys(node_local_sys.get());

        // Create 2 triangular elements to serve as the quadrature point geometry parent
        auto p_node_1 = Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0);
        auto p_node_2 = Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0);
        auto p_node_3 = Kratos::make_intrusive<NodeType>(3, 0.0, 1.0, 0.0);
        auto p_node_4 = Kratos::make_intrusive<NodeType>(4, 1.0, 1.0, 0.0);
        auto element_1 = GeneratePointsTriangle2D3_Element(p_node_1, p_node_2, p_node_3);
        auto element_2 = GeneratePointsTriangle2D3_Element(p_node_3, p_node_2, p_node_4);

        // Obtain the element shape functions
        auto N_1 = element_1->ShapeFunctionsValues();
        auto N_2 = element_2->ShapeFunctionsValues();

        Matrix N_1_mat(1, 3), N_2_mat(1, 3);
        for (std::size_t i = 0; i < 3; ++i){
            N_1_mat(0, i) = N_1(0, i);
            N_2_mat(0, i) = N_2(0, i);
        }

        // Obtain the element derivatives
        Matrix DN_De_1 = element_1->ShapeFunctionLocalGradient(0);
        Matrix DN_De_2 = element_2->ShapeFunctionLocalGradient(0);

        // Create the shape functions container
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> container_1(
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            element_1->IntegrationPoints()[0],
            N_1_mat,
            DN_De_1);
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> container_2(
            GeometryData::IntegrationMethod::GI_GAUSS_1,
            element_2->IntegrationPoints()[0],
            N_2_mat,
            DN_De_2);

        // Create the quadrature point geometries for the mapper
        auto qp_geom_1 = Kratos::make_shared<QuadraturePointGeometry<Node, 2, 2>>(
            element_1->Points(), container_1, element_1.get());
        auto qp_geom_2 = Kratos::make_shared<QuadraturePointGeometry<Node, 2, 2>>(
            element_2->Points(), container_2, element_2.get());
        
        InterfaceObject::Pointer interface_quadrature_point_1(Kratos::make_shared<InterfaceGeometryObject>((qp_geom_1.get())));
        InterfaceObject::Pointer interface_quadrature_point_2(Kratos::make_shared<InterfaceGeometryObject>((qp_geom_2.get())));

        const std::vector<int> expected_ids = {2, 1, 3};

        // Set the INTERFACE_EQUATION_ID for the nodes
        p_node_1->SetValue(INTERFACE_EQUATION_ID, 0);
        p_node_2->SetValue(INTERFACE_EQUATION_ID, 1);
        p_node_3->SetValue(INTERFACE_EQUATION_ID, 2);
        p_node_4->SetValue(INTERFACE_EQUATION_ID, 3);

        MapperInterfaceInfo::Pointer p_nearest_neighbor_iga_info(Kratos::make_shared<NearestNeighborInterfaceInfoIGA>(local_sys.Coordinates(), 0, 0));

        p_nearest_neighbor_iga_info->ProcessSearchResult(*interface_quadrature_point_1);
        p_nearest_neighbor_iga_info->ProcessSearchResult(*interface_quadrature_point_2);

        local_sys.AddInterfaceInfo(p_nearest_neighbor_iga_info);

        // Computing the local system
        MatrixType local_mapping_matrix;
        EquationIdVectorType origin_ids;
        EquationIdVectorType origin_ids2;
        EquationIdVectorType destination_ids;
        EquationIdVectorType destination_ids2;

        local_sys.EquationIdVectors(origin_ids, destination_ids);

        KRATOS_EXPECT_EQ(origin_ids.size(), 3);
        KRATOS_EXPECT_EQ(destination_ids.size(), 1);

        local_sys.CalculateLocalSystem(local_mapping_matrix, origin_ids2, destination_ids2);
        KRATOS_EXPECT_EQ(local_mapping_matrix.size1(), 1);
        KRATOS_EXPECT_EQ(local_mapping_matrix.size2(), 3);
        KRATOS_EXPECT_EQ(origin_ids2.size(), 3);
        KRATOS_EXPECT_EQ(destination_ids2.size(), 1);

        KRATOS_EXPECT_DOUBLE_EQ(local_mapping_matrix(0,0), 1.0/3.0);
        KRATOS_EXPECT_DOUBLE_EQ(local_mapping_matrix(0,1), 1.0/3.0);
        KRATOS_EXPECT_DOUBLE_EQ(local_mapping_matrix(0,2), 1.0/3.0);
        KRATOS_EXPECT_VECTOR_EQ(origin_ids2, expected_ids);
        KRATOS_EXPECT_EQ(destination_ids2[0], dest_id);
    }
}

// namespace Kratos::Testing