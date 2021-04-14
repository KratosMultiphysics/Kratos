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
#include "geometries/quadrature_point_surface_in_volume_geometry.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;

    NurbsVolumeGeometry<PointerVector<NodeType>> GenerateCubeForQPSurfaceInVolume() {
        // Construct Truncated Pyramid with: lower_base = 2x2; uper_base = 1.8x1.8; heigth = 4
        PointerVector<NodeType> points(196);
        double t = 1.0;
        std::vector<double> z_direction = {0.0, 1.0};
        std::vector<double> y_direction = {0.0, 1.0};
        std::vector<double> x_direction = {0.0, 1.0};
        std::size_t id = 1;
        for( auto i : z_direction){
            for( auto j : y_direction) {
                for( auto k : x_direction) {
                    double x = k*t;
                    double y = j*t;
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

        return NurbsVolumeGeometry<PointerVector<NodeType>>(
            points, polynomial_degree_u, polynomial_degree_v, polynomial_degree_w,
                knot_vector_u, knot_vector_v, knot_vector_w);
    }

    KRATOS_TEST_CASE_IN_SUITE(QuadraturePointSurfaceInVolumeGeometryTest, KratosCoreNurbsGeometriesFastSuite)
    {
        NurbsVolumeGeometry<PointerVector<NodeType>> cube = GenerateCubeForQPSurfaceInVolume();

        IntegrationPoint<3> integration_point;
        integration_point[0] = 0.5;
        integration_point[1] = 0.5;
        integration_point[2] = 0.5;
        integration_point.Weight() = 0.8;

        std::vector<IntegrationPoint<3>> integration_points(1);
        integration_points[0] = integration_point;

        PointerVector<Geometry<NodeType>> quad_geometries(1);

        // Create quadrature point geometries in the background.
        cube.CreateQuadraturePointGeometries(quad_geometries, 2, integration_points);

        auto quad_surface_in_volume = QuadraturePointSurfaceInVolumeGeometry<NodeType>(
            static_cast<QuadraturePointGeometry<NodeType,3,3,3>>(quad_geometries[0]), 1.0, 1.0, 1.0);

        KRATOS_CHECK_EQUAL(quad_surface_in_volume.Dimension(), 2);
        KRATOS_CHECK_EQUAL(quad_surface_in_volume.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(quad_surface_in_volume.LocalSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(quad_surface_in_volume.PointsNumber(), 8);

        KRATOS_WATCH(quad_surface_in_volume.DeterminantOfJacobian(0));
    }

} // End namespace Testsing
} // End namespace Kratos