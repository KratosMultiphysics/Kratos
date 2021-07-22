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
#include "testing/testing.h"
#include "containers/pointer_vector.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/surface_in_nurbs_volume_geometry.h"
#include "geometries/quadrature_point_surface_in_volume_geometry.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    Kratos::shared_ptr<NurbsVolumeGeometry<PointerVector<NodeType>>> GenerateCubeNurbsVolume() {
        PointerVector<NodeType> points(18);

        std::vector<double> z_direction = {0.0, 2.0};
        std::vector<double> y_direction = {0.0, 1.0, 2.0};
        std::vector<double> x_direction = {0.0, 1.0, 2.0};
        std::size_t id = 1;
        for( auto i : z_direction){
            for( auto j : y_direction) {
                for( auto k : x_direction) {
                    double x = k;
                    double y = j;
                    double z = i;
                    points(id-1) = Kratos::make_intrusive<NodeType>(id, x, y, z);
                    id++;
                }
            }
        }
        // Polynomial orders
        SizeType polynomial_degree_u = 2;
        SizeType polynomial_degree_v = 2;
        SizeType polynomial_degree_w = 1;

        // Assign knots of the basis along u.
        Vector knot_vector_u(4);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 1.0;
        knot_vector_u[3] = 1.0;

        // Assign knots of the basis along v.
        Vector knot_vector_v(4);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 0.0;
        knot_vector_v[2] = 1.0;
        knot_vector_v[3] = 1.0;

        // Assign knots of the basis along v
        Vector knot_vector_w(2);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 1.0;

        auto volume = NurbsVolumeGeometry<PointerVector<NodeType>>(points, polynomial_degree_u,
            polynomial_degree_v, polynomial_degree_w, knot_vector_u, knot_vector_v, knot_vector_w);

        return Kratos::make_shared<NurbsVolumeGeometry<PointerVector<NodeType>>>(volume);
    }

    Kratos::shared_ptr<NurbsVolumeGeometry<PointerVector<NodeType>>> GenerateCuboidNurbsVolume() {
        PointerVector<NodeType> points(8);
        double t = 1.0;
        std::vector<double> z_direction = {0.0, 0.5, 1.1, 2.0};
        std::vector<double> y_direction = {-1.0, 1.0};
        std::vector<double> x_direction = {-0.5, 0.5, 1.0, 3.0};
        std::size_t id = 1;
        for( auto i : z_direction){
            for( auto j : y_direction) {
                for( auto k : x_direction) {
                    double x = k;
                    double y = j;
                    double z = i;
                    points(id-1) = Kratos::make_intrusive<NodeType>(id, x, y, z);
                    id++;
                }
            }
            t += 0.8/6.0;
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

        // Assign knots of the basis along v
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

    KRATOS_TEST_CASE_IN_SUITE(SurfaceInVolumeGeometryTriangleInCubeTest, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_cube = GenerateCubeNurbsVolume();

        PointerVector<Geometry<NodeType>> quad_geometries(1);
        auto p_node_1 = Kratos::make_intrusive<NodeType>(0, 0.0, 0.0, 0.0);
        auto p_node_2 = Kratos::make_intrusive<NodeType>(1, 0.5, 0.5, 0.0);
        auto p_node_3 = Kratos::make_intrusive<NodeType>(2, 0.0, 0.5, 0.0);

        auto p_triangle = Kratos::make_shared<Triangle3D3<NodeType>>( p_node_1, p_node_2, p_node_3);
        //auto quad = Kratos::make_shared<Quadrilateral3D4<NodeType>>( p_node_1, p_node_2, p_node_3, p_node_4);
        // Create quadrature point geometries in the background.

        SurfaceInNurbsVolumeGeometry<3, PointerVector<NodeType>> surface_in_volume(p_cube, p_triangle);

        std::vector<IntegrationPoint<3>> integration_points_created;
        surface_in_volume.CreateIntegrationPoints(integration_points_created);

        surface_in_volume.CreateQuadraturePointGeometries(quad_geometries, 2, integration_points_created);


        double local_area_triangle =  surface_in_volume.Area();
        KRATOS_CHECK_NEAR( local_area_triangle, 0.5*0.5*0.5, 1e-10);

        std::vector<double> local_center_ref = {1.0/6.0, 1.0/3.0, 0.0};
        KRATOS_CHECK_VECTOR_NEAR(surface_in_volume.Center(), local_center_ref, 1e-10);
        // TODO: Rename this
        CoordinatesArrayType test_point;
        test_point[0] = 1.0/3.0;
        test_point[1] = 1.0/3.0;
        test_point[2] = 0.0;

        // Check get functions
        KRATOS_CHECK_EQUAL(surface_in_volume.HasGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX ), 1);
        auto p_nurbs_volume = surface_in_volume.pGetGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX );
        auto p_surface = surface_in_volume.pGetSurface();

        // Check dimension
        KRATOS_CHECK_EQUAL( p_nurbs_volume->LocalSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL( p_surface->LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL( surface_in_volume.LocalSpaceDimension(), 3 );
        KRATOS_CHECK_EQUAL( quad_geometries[0].LocalSpaceDimension(), 3 );

        // Check geometrical information
        CoordinatesArrayType global_coord;
        global_coord = surface_in_volume.GlobalCoordinates(global_coord, test_point);
        std::vector<double> global_center_ref = {1.0/3.0, 2.0/3.0, 0.0};
        KRATOS_CHECK_VECTOR_NEAR( global_coord, global_center_ref, 1e-10);
        std::vector<double> normal_ref = {0, 0, 1};
        auto integration_method = quad_geometries[0].GetDefaultIntegrationMethod();
        KRATOS_CHECK_VECTOR_NEAR(quad_geometries[0].Normal( 0, integration_method), normal_ref, 1e-10 );

        // Check integration
        double global_area_triangle_ref = 0.5;
        KRATOS_CHECK_NEAR( surface_in_volume.DeterminantOfJacobian(test_point), 2.0*global_area_triangle_ref, 1e-10);
        auto integration_points = quad_geometries[0].IntegrationPoints();
        KRATOS_CHECK_EQUAL(integration_points.size(), 1);
        double weight = integration_points[0].Weight();
        KRATOS_CHECK_VECTOR_NEAR(integration_points[0].Coordinates(), test_point, 1e-10);
        KRATOS_CHECK_NEAR( weight, 0.5, 1e-10);
        KRATOS_CHECK_NEAR( quad_geometries[0].DeterminantOfJacobian(0, integration_method )*weight, 0.5, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(SurfaceInVolumeGeometrySquareInCubeTest, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_cube = GenerateCubeNurbsVolume();

        PointerVector<Geometry<NodeType>> quad_geometries(1);
        auto p_node_1 = Kratos::make_intrusive<NodeType>(0, 0.0, 0.0, 0.0);
        auto p_node_2 = Kratos::make_intrusive<NodeType>(1, 0.5, 0.5, 0.0);
        auto p_node_3 = Kratos::make_intrusive<NodeType>(2, 0.5, 0.5, 1.0);
        auto p_node_4 = Kratos::make_intrusive<NodeType>(3, 0.0, 0.0, 1.0);

        auto p_quad = Kratos::make_shared<Quadrilateral3D4<NodeType>>( p_node_1, p_node_2, p_node_3, p_node_4 );

        SurfaceInNurbsVolumeGeometry<3, PointerVector<NodeType>> surface_in_volume(p_cube, p_quad);

        std::vector<IntegrationPoint<3>> integration_points_created;
        surface_in_volume.CreateIntegrationPoints(integration_points_created);

        surface_in_volume.CreateQuadraturePointGeometries(quad_geometries, 2, integration_points_created);

        double local_area_triangle =  surface_in_volume.Area();
        KRATOS_CHECK_NEAR( local_area_triangle, sqrt(0.5*0.5 + 0.5*0.5), 1e-10);

        std::vector<double> local_center_ref = {0.25, 0.25, 0.5};
        KRATOS_CHECK_VECTOR_NEAR(surface_in_volume.Center(), local_center_ref, 1e-10);
        CoordinatesArrayType test_point;
        test_point[0] = 0.0;
        test_point[1] = 0.0;
        test_point[2] = 0.0;

        // Check Get functions
        KRATOS_CHECK_EQUAL(surface_in_volume.HasGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX ), 1);
        auto p_nurbs_volume = surface_in_volume.pGetGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX );
        auto p_surface = surface_in_volume.pGetSurface();

        // Check Dimesions
        KRATOS_CHECK_EQUAL( p_nurbs_volume->LocalSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL( p_surface->LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL( surface_in_volume.LocalSpaceDimension(), 3 );
        KRATOS_CHECK_EQUAL( quad_geometries[0].LocalSpaceDimension(), 3 );

        // Check geometrical information
        CoordinatesArrayType global_coord;
        global_coord = surface_in_volume.GlobalCoordinates(global_coord, test_point);
        std::vector<double> global_center_ref = {0.5, 0.5, 1.0};
        KRATOS_CHECK_VECTOR_NEAR( global_coord, global_center_ref, 1e-10);
        std::vector<double> normal_ref = {0.5, -0.5, 0};
        auto integration_method = quad_geometries[0].GetDefaultIntegrationMethod();
        KRATOS_CHECK_VECTOR_NEAR(quad_geometries[0].Normal( 0, integration_method), normal_ref, 1e-10 );

        // Check integration
        KRATOS_CHECK_EQUAL( quad_geometries.size(), 4);
        double global_area_triangle = 0.0;
        for( SizeType i = 0; i < quad_geometries.size(); ++i ){
            global_area_triangle += quad_geometries[i].DeterminantOfJacobian(0, integration_method ) *
                quad_geometries[i].IntegrationPoints()[0].Weight();
        }

        KRATOS_CHECK_NEAR( global_area_triangle, std::sqrt(2.0)*2.0, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(SurfaceInVolumeGeometryTriangleInCubeTest_2, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_cube = GenerateCubeNurbsVolume();

        PointerVector<NodeType> points_global_space;
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(0, 0.0, 0.0, 0.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(1, 1.0, 0.0, 0.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(2, 2.0, 0.0, 1.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(3, 1.0, 0.5, 1.5) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(4, 3.0, 1.0, 2.0) );
        points_global_space.push_back( Kratos::make_intrusive<NodeType>(5, 1.5, 0.5, 1.0) );

        auto p_triangle_global_space = Kratos::make_shared<Triangle3D6<NodeType>>(points_global_space);

        // KRATOS_WATCH( p_triangle->Area() );
        // CoordinatesArrayType test_point;
        // test_point[0] = 1.0/3.0;
        // test_point[1] = 1.0/3.0;
        // test_point[2] = 0.0;
        // KRATOS_WATCH( p_triangle->Normal(test_point) );

        PointerVector<NodeType> points_local_space;
        for( IndexType i = 0; i < points_global_space.size(); ++i){
            double x_coord_local = (points_global_space[i].Coordinates()[0] + 0.5) / 3.5;
            double y_coord_local = (points_global_space[i].Coordinates()[1] + 1.0) / 2.0;
            double z_coord_local = (points_global_space[i].Coordinates()[2] + 0.0) / 2.0;
            points_local_space.push_back(Kratos::make_intrusive<NodeType>(i, x_coord_local, y_coord_local, z_coord_local));
        }
        auto p_triangle_local_space = Kratos::make_shared<Triangle3D6<NodeType>>(points_local_space);
        // std::vector<Kratos_intrusive_ptr<NodeType> nodes_global = {p_node_global_1, p_node_global_2 , p_node_global_3
        //      p_node_global_1, p_node_global_1, p_node_global_1,p_node_global_1}
        // std::vector<double> z_direction = {0.0, 0.5, 1.1, 2.0};
        // std::vector<double> y_direction = {-1.0, 1.0};
        // std::vector<double> x_direction = {-0.5, 0.5, 1.0, 3.0};


        // Transform nodes into parameter space of nurbs volume


        // SurfaceInNurbsVolumeGeometry<3, PointerVector<NodeType>> surface_in_volume(p_cube, p_quad);

        // std::vector<IntegrationPoint<3>> integration_points_created;
        // surface_in_volume.CreateIntegrationPoints(integration_points_created);

        // surface_in_volume.CreateQuadraturePointGeometries(quad_geometries, 2, integration_points_created);

        // double local_area_triangle =  surface_in_volume.Area();
        // KRATOS_CHECK_NEAR( local_area_triangle, sqrt(0.5*0.5 + 0.5*0.5), 1e-10);

        // std::vector<double> local_center_ref = {0.25, 0.25, 0.5};
        // KRATOS_CHECK_VECTOR_NEAR(surface_in_volume.Center(), local_center_ref, 1e-10);
        // CoordinatesArrayType test_point;
        // test_point[0] = 0.0;
        // test_point[1] = 0.0;
        // test_point[2] = 0.0;

        // // Check Get functions
        // KRATOS_CHECK_EQUAL(surface_in_volume.HasGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX ), 1);
        // auto p_nurbs_volume = surface_in_volume.pGetGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX );
        // auto p_surface = surface_in_volume.pGetSurface();

        // // Check Dimesions
        // KRATOS_CHECK_EQUAL( p_nurbs_volume->LocalSpaceDimension(), 3);
        // KRATOS_CHECK_EQUAL( p_surface->LocalSpaceDimension(), 2);
        // KRATOS_CHECK_EQUAL( surface_in_volume.LocalSpaceDimension(), 3 );
        // KRATOS_CHECK_EQUAL( quad_geometries[0].LocalSpaceDimension(), 3 );

        // // Check geometrical information
        // CoordinatesArrayType global_coord;
        // global_coord = surface_in_volume.GlobalCoordinates(global_coord, test_point);
        // std::vector<double> global_center_ref = {0.5, 0.5, 1.0};
        // KRATOS_CHECK_VECTOR_NEAR( global_coord, global_center_ref, 1e-10);
        // std::vector<double> normal_ref = {0.5, -0.5, 0};
        // auto integration_method = quad_geometries[0].GetDefaultIntegrationMethod();
        // KRATOS_CHECK_VECTOR_NEAR(quad_geometries[0].Normal( 0, integration_method), normal_ref, 1e-10 );

        // // Check integration
        // KRATOS_CHECK_EQUAL( quad_geometries.size(), 4);
        // double global_area_triangle = 0.0;
        // for( SizeType i = 0; i < quad_geometries.size(); ++i ){
        //     global_area_triangle += quad_geometries[i].DeterminantOfJacobian(0, integration_method ) *
        //         quad_geometries[i].IntegrationPoints()[0].Weight();
        // }

        // KRATOS_CHECK_NEAR( global_area_triangle, std::sqrt(2.0)*2.0, 1e-10);
    }
} // End namespace Testsing
} // End namespace Kratos