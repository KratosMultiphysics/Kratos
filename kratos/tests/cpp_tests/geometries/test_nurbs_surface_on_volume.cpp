//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

// System includes

// External includes

// Project includes
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "testing/testing.h"
#include "containers/pointer_vector.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_surface_on_volume_geometry.h"
#include "geometries/quadrature_point_surface_in_volume_geometry.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"

namespace Kratos {
namespace Testing {

    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    // Helper Functions
    Kratos::shared_ptr<NurbsVolumeGeometry<PointerVector<NodeType>>> GenerateCubeNurbsVolume() {

        // Define control points
        PointerVector<NodeType> points;
        std::vector<double> z_direction = {0.0, 1.0, 2.0};
        std::vector<double> y_direction = {0.0, 2.0/9.0, 2.0/3.0, 4.0/3.0, 16.0/9.0, 2.0};
        std::vector<double> x_direction = {0.0, 1.0, 2.0};
        std::size_t id = 1;
        for( auto i : z_direction){
            for( auto j : y_direction) {
                for( auto k : x_direction) {
                    points.push_back( Kratos::make_intrusive<NodeType>(id, k, j, i) );
                    id++;
                }
            }
        }

        // Polynomial orders.
        SizeType polynomial_degree_u = 2;
        SizeType polynomial_degree_v = 3;
        SizeType polynomial_degree_w = 1;

        // Assign knots of the basis along u.
        Vector knot_vector_u(4);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 1.0;
        knot_vector_u[3] = 1.0;

        // Assign knots of the basis along v.
        Vector knot_vector_v(8);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 0.0;
        knot_vector_v[2] = 0.0;
        knot_vector_v[3] = 1.0/3.0;
        knot_vector_v[4] = 2.0/3.0;
        knot_vector_v[5] = 1.0;
        knot_vector_v[6] = 1.0;
        knot_vector_v[7] = 1.0;

        // Assign knots of the basis along w.
        Vector knot_vector_w(3);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 0.5;
        knot_vector_w[2] = 1.0;

        auto volume = NurbsVolumeGeometry<PointerVector<NodeType>>(points, polynomial_degree_u,
            polynomial_degree_v, polynomial_degree_w, knot_vector_u, knot_vector_v, knot_vector_w);

        return Kratos::make_shared<NurbsVolumeGeometry<PointerVector<NodeType>>>(volume);
    }

    Kratos::shared_ptr<NurbsVolumeGeometry<PointerVector<NodeType>>> GenerateCuboidNurbsVolume() {

        // Helper Functions
        PointerVector<NodeType> points;
        std::vector<double> z_direction = {0.0, 2.0/3.0, 4.0/3.0, 2.0};
        std::vector<double> y_direction = {-1.0, 1.0};
        std::vector<double> x_direction = {-0.5, 0.375, 2.125, 3.0};
        std::size_t id = 1;
        for( auto i : z_direction){
            for( auto j : y_direction) {
                for( auto k : x_direction) {
                    points.push_back( Kratos::make_intrusive<NodeType>(id, k, j, i) );
                    id++;
                }
            }
        }

        // Polynomial orders
        SizeType polynomial_degree_u = 2;
        SizeType polynomial_degree_v = 1;
        SizeType polynomial_degree_w = 3;

        // Assign knots of the basis along u.
        Vector knot_vector_u(5);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 0.5;
        knot_vector_u[3] = 1.0;
        knot_vector_u[4] = 1.0;

        // Assign knots of the basis along v.
        Vector knot_vector_v(2);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 1.0;

        // Assign knots of the basis along w.
        Vector knot_vector_w(6);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 0.0;
        knot_vector_w[2] = 0.0;
        knot_vector_w[3] = 1.0;
        knot_vector_w[4] = 1.0;
        knot_vector_w[5] = 1.0;

        auto volume = NurbsVolumeGeometry<PointerVector<NodeType>>(points, polynomial_degree_u,
            polynomial_degree_v, polynomial_degree_w, knot_vector_u, knot_vector_v, knot_vector_w);

        return Kratos::make_shared<NurbsVolumeGeometry<PointerVector<NodeType>>>(volume);
    }

    // Tests
    KRATOS_TEST_CASE_IN_SUITE(SurfaceInVolumeGeometryTriangleInCubeTest, KratosCoreNurbsGeometriesFastSuite)
    {
        // Get nurbs cube
        auto p_nurbs_cube = GenerateCubeNurbsVolume();

        // Create surface
        auto p_node_1 = Kratos::make_intrusive<NodeType>(0, 0.0, 0.0, 0.0); // nodes are in local space.
        auto p_node_2 = Kratos::make_intrusive<NodeType>(1, 0.5, 0.5, 0.0);
        auto p_node_3 = Kratos::make_intrusive<NodeType>(2, 0.0, 0.5, 0.0);
        auto p_triangle = Kratos::make_shared<Triangle3D3<NodeType>>( p_node_1, p_node_2, p_node_3);

        // Create surface in nurbs volume
        NurbsSurfaceOnVolumeGeometry<3, PointerVector<NodeType>> surface_in_volume(p_nurbs_cube, p_triangle);

        // Create integration points and quadrature point geometries.
        std::vector<IntegrationPoint<3>> integration_points_created;
        IntegrationInfo integration_info = surface_in_volume.GetDefaultIntegrationInfo();
        surface_in_volume.CreateIntegrationPoints(integration_points_created, integration_info);

        PointerVector<Geometry<NodeType>> quad_geometries;
        surface_in_volume.CreateQuadraturePointGeometries(quad_geometries, 2, integration_points_created, integration_info);

        // Check get functions
        KRATOS_EXPECT_EQ(surface_in_volume.HasGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX ), 1);
        auto p_nurbs_volume = surface_in_volume.pGetGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX );
        auto p_surface = surface_in_volume.pGetSurface();

        // Check kratos geometry families
        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Surface_On_Volume;
        KRATOS_EXPECT_EQ(surface_in_volume.GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(surface_in_volume.GetGeometryType(), geometry_type);

        // Check dimension
        KRATOS_EXPECT_EQ( p_nurbs_volume->LocalSpaceDimension(), 3);
        KRATOS_EXPECT_EQ( p_surface->LocalSpaceDimension(), 2);
        KRATOS_EXPECT_EQ( surface_in_volume.LocalSpaceDimension(), 3 );
        KRATOS_EXPECT_EQ( quad_geometries[0].LocalSpaceDimension(), 3 );
        KRATOS_EXPECT_EQ( surface_in_volume.WorkingSpaceDimension(), 3 );
        KRATOS_EXPECT_EQ( quad_geometries[0].WorkingSpaceDimension(), 3 );

        // Check geometrical information
        double global_area_triangle =  surface_in_volume.Area();
        KRATOS_EXPECT_NEAR( global_area_triangle, 0.5, 1e-10);

        CoordinatesArrayType test_point; // Center in Dimension Space
        test_point[0] = 1.0/3.0;
        test_point[1] = 1.0/3.0;
        test_point[2] = 0.0;

        CoordinatesArrayType global_coord;
        global_coord = surface_in_volume.GlobalCoordinates(global_coord, test_point);
        std::vector<double> global_center_ref = {1.0/3.0, 2.0/3.0, 0.0};
        KRATOS_EXPECT_VECTOR_NEAR( global_coord, global_center_ref, 1e-10);
        KRATOS_EXPECT_VECTOR_NEAR(surface_in_volume.Center(), global_center_ref, 1e-10);

        std::vector<double> normal_ref = {0, 0, 1};
        auto integration_method = quad_geometries[0].GetDefaultIntegrationMethod();
        KRATOS_EXPECT_VECTOR_NEAR(quad_geometries[0].Normal( 0, integration_method), normal_ref, 1e-10 );

        // Check integration
        double global_area_triangle_ref = 0.5;
        KRATOS_EXPECT_NEAR( surface_in_volume.DeterminantOfJacobian(test_point), 2.0*global_area_triangle_ref, 1e-10);
        auto integration_points = quad_geometries[0].IntegrationPoints();
        KRATOS_EXPECT_EQ(integration_points.size(), 1);
        double weight = integration_points[0].Weight();
        KRATOS_EXPECT_VECTOR_NEAR(integration_points[0].Coordinates(), test_point, 1e-10);
        KRATOS_EXPECT_NEAR( weight, 0.5, 1e-10);
        KRATOS_EXPECT_NEAR( quad_geometries[0].DeterminantOfJacobian(0, integration_method )*weight, 0.5, 1e-10);

        // Check kratos geometry families
        const auto geometry_family_qp = GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry;
        const auto geometry_type_qp = GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry;
        KRATOS_EXPECT_EQ(quad_geometries[0].GetGeometryFamily(), geometry_family_qp);
        KRATOS_EXPECT_EQ(quad_geometries[0].GetGeometryType(), geometry_type_qp);
    }

    KRATOS_TEST_CASE_IN_SUITE(SurfaceInVolumeGeometryQuadInCubeTest, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_nurbs_cube = GenerateCubeNurbsVolume();

        auto p_node_1 = Kratos::make_intrusive<NodeType>(0, 0.0, 0.0, 0.0);
        auto p_node_2 = Kratos::make_intrusive<NodeType>(1, 0.5, 0.5, 0.0);
        auto p_node_3 = Kratos::make_intrusive<NodeType>(2, 0.5, 0.5, 1.0);
        auto p_node_4 = Kratos::make_intrusive<NodeType>(3, 0.0, 0.0, 1.0);

        auto p_quad = Kratos::make_shared<Quadrilateral3D4<NodeType>>( p_node_1, p_node_2, p_node_3, p_node_4 );

        NurbsSurfaceOnVolumeGeometry<3, PointerVector<NodeType>> surface_in_volume(p_nurbs_cube, p_quad);

        // Create integration points and quadrature point geometries.
        std::vector<IntegrationPoint<3>> integration_points_created;
        IntegrationInfo integration_info = surface_in_volume.GetDefaultIntegrationInfo();
        surface_in_volume.CreateIntegrationPoints(integration_points_created, integration_info);

        PointerVector<Geometry<NodeType>> quad_geometries;
        surface_in_volume.CreateQuadraturePointGeometries(quad_geometries, 2, integration_points_created, integration_info);

        // Check kratos geometry families
        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Surface_On_Volume;
        KRATOS_EXPECT_EQ(surface_in_volume.GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(surface_in_volume.GetGeometryType(), geometry_type);

        // Check get functions
        KRATOS_EXPECT_EQ(surface_in_volume.HasGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX ), 1);
        auto p_nurbs_volume = surface_in_volume.pGetGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX );
        auto p_surface = surface_in_volume.pGetSurface();

        // Check dimensions
        KRATOS_EXPECT_EQ( p_nurbs_volume->LocalSpaceDimension(), 3);
        KRATOS_EXPECT_EQ( p_surface->LocalSpaceDimension(), 2);
        KRATOS_EXPECT_EQ( surface_in_volume.LocalSpaceDimension(), 3 );
        KRATOS_EXPECT_EQ( quad_geometries[0].LocalSpaceDimension(), 3 );
        KRATOS_EXPECT_EQ( surface_in_volume.WorkingSpaceDimension(), 3 );
        KRATOS_EXPECT_EQ( quad_geometries[0].WorkingSpaceDimension(), 3 );

        // Check geometrical information
        CoordinatesArrayType test_point; // Center in Dimension Space
        test_point[0] = 0.0;
        test_point[1] = 0.0;
        test_point[2] = 0.0;

        double global_area_triangle =  surface_in_volume.Area();
        KRATOS_EXPECT_NEAR( global_area_triangle, std::sqrt(2.0)*2.0, 1e-10);

        CoordinatesArrayType global_coord;
        global_coord = surface_in_volume.GlobalCoordinates(global_coord, test_point);
        std::vector<double> global_center_ref = {0.5, 0.5, 1.0};
        KRATOS_EXPECT_VECTOR_NEAR( global_coord, global_center_ref, 1e-10);
        KRATOS_EXPECT_VECTOR_NEAR(surface_in_volume.Center(), global_center_ref, 1e-10);

        std::vector<double> normal_ref = {0.5, -0.5, 0};
        auto integration_method = quad_geometries[0].GetDefaultIntegrationMethod();
        KRATOS_EXPECT_VECTOR_NEAR(quad_geometries[0].Normal( 0, integration_method), normal_ref, 1e-10 );

        // Check integration
        KRATOS_EXPECT_EQ( quad_geometries.size(), 4);
        global_area_triangle = 0.0;
        for( SizeType i = 0; i < quad_geometries.size(); ++i ){
            global_area_triangle += quad_geometries[i].DeterminantOfJacobian(0, integration_method ) *
                quad_geometries[i].IntegrationPoints()[0].Weight();
        }

        KRATOS_EXPECT_NEAR( global_area_triangle, std::sqrt(2.0)*2.0, 1e-10);

        // Check kratos geometry families
        const auto geometry_family_qp = GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry;
        const auto geometry_type_qp = GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry;
        for( SizeType i = 0; i < quad_geometries.size(); ++i ){
            KRATOS_EXPECT_EQ(quad_geometries[i].GetGeometryFamily(), geometry_family_qp);
            KRATOS_EXPECT_EQ(quad_geometries[i].GetGeometryType(), geometry_type_qp);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(SurfaceInVolumeGeometryQuadInCuboidTest, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_nurbs_cuboid = GenerateCuboidNurbsVolume();

        PointerVector<NodeType> points_global_space;
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(0, 0.0, 0.0, 0.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(1, 1.0, 0.0, 0.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(2, 2.0, 0.0, 0.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(3, 2.5, 0.5, 1.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(4, 3.0, 1.0, 2.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(5, 2.0, 1.0, 2.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(6, 1.0, 1.0, 2.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(7, 0.5, 0.5, 1.0) );

        auto p_quad_global_space = Kratos::make_shared<Quadrilateral3D8<NodeType>>(points_global_space);

        double area_ref =  p_quad_global_space->Area();
        auto center_ref = p_quad_global_space->Center();
        CoordinatesArrayType test_point; // Center in Dimension Space
        test_point[0] = 0.0;
        test_point[1] = 0.0;
        test_point[2] = 0.0;
        auto normal_ref =  p_quad_global_space->Normal(test_point);
        CoordinatesArrayType coord_ref;
        p_quad_global_space->GlobalCoordinates(coord_ref, test_point);

        std::vector<double> z_direction = {0.0, 0.5, 1.1, 2.0};
        std::vector<double> y_direction = {-1.0, 1.0};
        std::vector<double> x_direction = {-0.5, 0.5,  3.0};

        // Transform nodes into parameter space of nurbs volume
        PointerVector<NodeType> points_local_space;
        for( IndexType i = 0; i < points_global_space.size(); ++i){
            double x_coord_local = (points_global_space[i].Coordinates()[0] + 0.5) / 3.5;
            double y_coord_local = (points_global_space[i].Coordinates()[1] + 1.0) / 2.0;
            double z_coord_local = (points_global_space[i].Coordinates()[2] + 0.0) / 2.0;
            points_local_space.push_back(Kratos::make_intrusive<NodeType>(i, x_coord_local, y_coord_local, z_coord_local));
        }
        auto p_quad_local_space = Kratos::make_shared<Quadrilateral3D8<NodeType>>(points_local_space);

        NurbsSurfaceOnVolumeGeometry<3, PointerVector<NodeType>> surface_in_volume(p_nurbs_cuboid, p_quad_local_space);

        // Create integration points and quadrature point geometries.
        std::vector<IntegrationPoint<3>> integration_points_created;
        IntegrationInfo integration_info = surface_in_volume.GetDefaultIntegrationInfo();
        surface_in_volume.CreateIntegrationPoints(integration_points_created, integration_info);

        PointerVector<Geometry<NodeType>> quad_geometries;
        surface_in_volume.CreateQuadraturePointGeometries(quad_geometries, 2, integration_points_created, integration_info);

        // Check kratos geometry families
        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Surface_On_Volume;
        KRATOS_EXPECT_EQ(surface_in_volume.GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(surface_in_volume.GetGeometryType(), geometry_type);

        // Check Get functions
        KRATOS_EXPECT_EQ(surface_in_volume.HasGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX ), 1);
        auto p_nurbs_volume = surface_in_volume.pGetGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX );
        auto p_surface = surface_in_volume.pGetSurface();

        // Check local tangents
        Matrix local_tangents;
        BoundedMatrix<double,3,2> local_tangents_ref;
        local_tangents_ref(0,0) = 1.10963416786212;
        local_tangents_ref(0,1) = -0.15090726181280;
        local_tangents_ref(1,0) = 0.81368437741389;
        local_tangents_ref(1,1) = 0.34503521010352;
        local_tangents_ref(2,0) = 1.62736875482778;
        local_tangents_ref(2,1) = 0.69007042020704;

        quad_geometries[0].Calculate( LOCAL_TANGENT_MATRIX, local_tangents);
        KRATOS_EXPECT_EQ(local_tangents.size1(), 3);
        KRATOS_EXPECT_EQ(local_tangents.size2(), 2);
        KRATOS_EXPECT_MATRIX_NEAR( local_tangents, local_tangents_ref, 1e-10);

        // Check dimensions
        KRATOS_EXPECT_EQ( p_nurbs_volume->LocalSpaceDimension(), 3);
        KRATOS_EXPECT_EQ( p_surface->LocalSpaceDimension(), 2);
        KRATOS_EXPECT_EQ( surface_in_volume.LocalSpaceDimension(), 3 );
        KRATOS_EXPECT_EQ( quad_geometries[0].LocalSpaceDimension(), 3 );
        KRATOS_EXPECT_EQ( surface_in_volume.WorkingSpaceDimension(), 3 );
        KRATOS_EXPECT_EQ( quad_geometries[0].WorkingSpaceDimension(), 3 );

        // Check geometrical information
        auto center = surface_in_volume.Center();
        KRATOS_EXPECT_VECTOR_NEAR( center, center_ref, 1e-10);
        double global_area_triangle = surface_in_volume.Area();
        KRATOS_EXPECT_NEAR(global_area_triangle, area_ref, 1e-10);

        // Check integration
        KRATOS_EXPECT_EQ( quad_geometries.size(), 9);
        auto integration_method = quad_geometries[0].GetDefaultIntegrationMethod();
        global_area_triangle = 0.0;
        for( SizeType i = 0; i < quad_geometries.size(); ++i ){
            global_area_triangle += quad_geometries[i].DeterminantOfJacobian(0, integration_method ) *
                quad_geometries[i].IntegrationPoints()[0].Weight();
            auto normal = quad_geometries[i].Normal( 0, integration_method);
            auto normal_ref = p_quad_global_space->Normal( quad_geometries[i].IntegrationPoints()[0].Coordinates() );
            KRATOS_EXPECT_VECTOR_NEAR( normal, normal_ref, 1e-10);

            double det_J_quad = quad_geometries[i].DeterminantOfJacobian(0, integration_method );
            double det_J = surface_in_volume.DeterminantOfJacobian( quad_geometries[i].IntegrationPoints()[0].Coordinates() );
            KRATOS_EXPECT_NEAR(det_J_quad, det_J, 1e-10);
        }
        KRATOS_EXPECT_NEAR( global_area_triangle, area_ref, 1e-10);

        // Check kratos geometry families
        const auto geometry_family_qp = GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry;
        const auto geometry_type_qp = GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry;
        for( SizeType i = 0; i < quad_geometries.size(); ++i ){
            KRATOS_EXPECT_EQ(quad_geometries[i].GetGeometryFamily(), geometry_family_qp);
            KRATOS_EXPECT_EQ(quad_geometries[i].GetGeometryType(), geometry_type_qp);
        }
    }
} // End namespace Testsing
} // End namespace Kratos