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
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    Kratos::shared_ptr<NurbsVolumeGeometry<PointerVector<NodeType>>> GenerateCubeNurbsVolume() {
        PointerVector<NodeType> points(8);
        double t = 1.0;
        std::vector<double> z_direction = {0.0, 1.0};
        std::vector<double> y_direction = {0.0, 1.0};
        std::vector<double> x_direction = {0.0, 1.0};
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
        SizeType polynomial_degree_u = 1;
        SizeType polynomial_degree_v = 1;
        SizeType polynomial_degree_w = 1;

        // Assign knots of the basis along u.
        Vector knot_vector_u(2);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 1.0;

        // Assign knots of the basis along v.
        Vector knot_vector_v(2);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 1.0;

        // Assign knots of the basis along v
        Vector knot_vector_w(2);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 1.0;

        auto volume = NurbsVolumeGeometry<PointerVector<NodeType>>(points, polynomial_degree_u,
            polynomial_degree_v, polynomial_degree_w, knot_vector_u, knot_vector_v, knot_vector_w);

        return Kratos::make_shared<NurbsVolumeGeometry<PointerVector<NodeType>>>(volume);
    }

    KRATOS_TEST_CASE_IN_SUITE(QuadraturePointSurfaceInVolumeGeometryTest, KratosCoreNurbsGeometriesFastSuite)
    {
        auto cube = GenerateCubeNurbsVolume();

        IntegrationPoint<3> integration_point;
        integration_point[0] = 0.5;
        integration_point[1] = 0.5;
        integration_point[2] = 0.5;
        integration_point.Weight() = 0.8;

        std::vector<IntegrationPoint<3>> integration_points(1);
        integration_points[0] = integration_point;

        PointerVector<Geometry<NodeType>> quad_geometries(1);

        auto p_node_1 = Kratos::make_intrusive<NodeType>(0, 0.0, 0.0, 0.0);
        auto p_node_2 = Kratos::make_intrusive<NodeType>(1, 1.0, 1.0, 0.0);
        auto p_node_3 = Kratos::make_intrusive<NodeType>(2, 0.0, 1.0, 0.0);
        NodeType node1(0, 0.2, 0.2, 0.2);
        auto triangle = Kratos::make_shared<Triangle3D3<NodeType>>( p_node_1, p_node_2, p_node_3);
        // Create quadrature point geometries in the background.
        SurfaceInNurbsVolumeGeometry<3, PointerVector<NodeType>> surface_in_volume(cube, triangle);
        //cube.CreateQuadraturePointGeometries(quad_geometries, 2, integration_points);

        // KRATOS_CHECK_EQUAL(surface_in_volume.HasGeometryPart( GeometryType::BACKGROUND_GEOMETRY_INDEX ), 1);
        KRATOS_CHECK_NEAR(surface_in_volume.Area(), 0.5, 1e-10);
        std::vector<double> center_ref = {1.0/3.0, 2.0/3.0, 0.0};
        KRATOS_CHECK_VECTOR_NEAR(surface_in_volume.Center(), center_ref, 1e-10);
        CoordinatesArrayType test_point;
        test_point[0] = 1.0/3.0;
        test_point[1] = 1.0/3.0;
        test_point[2] = 0.0;
        // Check what should determinatn be??
        KRATOS_WATCH( surface_in_volume.DeterminantOfJacobian(test_point) );


        // Check everthing after double DeterminantOfJacobian(


    //     auto quad_surface_in_volume = QuadraturePointSurfaceInVolumeGeometry<NodeType>(
    //         casted, 1.0, 1.0, 1.0);

    //     KRATOS_CHECK_EQUAL(quad_surface_in_volume.Dimension(), 2);
    //     KRATOS_CHECK_EQUAL(quad_surface_in_volume.WorkingSpaceDimension(), 3);
    //     KRATOS_CHECK_EQUAL(quad_surface_in_volume.LocalSpaceDimension(), 3);
    //     KRATOS_CHECK_EQUAL(quad_surface_in_volume.PointsNumber(), 8);

    //     KRATOS_WATCH(quad_surface_in_volume.DeterminantOfJacobian(0));
        }

} // End namespace Testsing
} // End namespace Kratos