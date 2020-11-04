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
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions.h"

#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/point.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;

    NurbsVolumeGeometry<3, PointerVector<Point>> GenerateUniformTruncatedPyramid() {

        // Construct Truncated Pyramid with: lower_base = 2x2; uper_base = 1.8x1.8; heigth = 4
        PointerVector<Point> points;
        double t = 0.8;
        for( int i = 0; i <=4; ++i){
            t += 0.2;
            for( int j = -1; j <= 1; ++j) {
                for( int k = -1; k <= 1; ++k) {
                    double x = k*t;
                    double y = j*t;
                    double z = i;

                    points.push_back(Point::Pointer(new Point(x, y, z)));
                }
            }
        }
        // Polynomial orders
        SizeType polynomial_degree_u = 2;
        SizeType polynomial_degree_v = 2;
        SizeType polynomial_degree_w = 4;

        // TODO: Add inner knots
        // Assign knots of the basis along u
        SizeType number_of_knots_u = 6;
        Vector knot_vector_u(number_of_knots_u);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 0.0;
        knot_vector_u[3] = 1.0;
        knot_vector_u[4] = 1.0;
        knot_vector_u[5] = 1.0;

        // Assign knots of the basis along v
        SizeType number_of_knots_v = 6;
        Vector knot_vector_v(number_of_knots_v);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 0.0;
        knot_vector_v[2] = 0.0;
        knot_vector_v[3] = 1.0;
        knot_vector_v[4] = 1.0;
        knot_vector_v[5] = 1.0;

        // Assign knots of the basis along v
        SizeType number_of_knots_w = 10;
        Vector knot_vector_w(number_of_knots_w);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 0.0;
        knot_vector_w[2] = 0.0;
        knot_vector_w[3] = 0.0;
        knot_vector_w[4] = 0.0;
        knot_vector_w[5] = 1.0;
        knot_vector_w[6] = 1.0;
        knot_vector_w[7] = 1.0;
        knot_vector_w[8] = 1.0;
        knot_vector_w[9] = 1.0;

        Vector weights = ZeroVector(360);

        return NurbsVolumeGeometry<3, PointerVector<Point>>(
            points, polynomial_degree_u, polynomial_degree_v, polynomial_degree_w,
                knot_vector_u, knot_vector_v, knot_vector_w);
    }

    NurbsVolumeGeometry<3, PointerVector<Point>> GenerateDistortedCube() {

        // Construct Truncated Pyramid with: lower_base = 2x2; uper_base = 1.8x1.8; heigth = 4
        PointerVector<Point> points;
        double t = 0.8;

        for( int i = 0; i <=4; ++i){
            t += 0.2;
            for( double j = -1; j < 1.001; j = j + 2.0/3.0) {
                for( int k = -2; k <=2; ++k ) {
                    double x = k;
                    double y = j;
                    double z = i;
                    if( k == 0){
                        y = j*std::max(i,2);
                    }
                    if( j == 0){
                        x = k*std::max(i,2);
                    }

                    points.push_back(Point::Pointer(new Point(x, y, z)));
                }
            }
        }

        // Polynomial orders
        SizeType polynomial_degree_u = 2;
        SizeType polynomial_degree_v = 2;
        SizeType polynomial_degree_w = 4;

        // Assign knots of the basis along u
        SizeType number_of_knots_u = 8;
        Vector knot_vector_u(number_of_knots_u);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 0.0;
        knot_vector_u[3] = 0.25;
        knot_vector_u[4] = 0.75;
        knot_vector_u[5] = 1.0;
        knot_vector_u[6] = 1.0;
        knot_vector_u[7] = 1.0;
        //knot_vector_u[7] = 1.0;

        // Assign knots of the basis along v
        SizeType number_of_knots_v = 7;
        Vector knot_vector_v(number_of_knots_v);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 0.0;
        knot_vector_v[2] = 0.0;
        knot_vector_v[3] = 0.5;
        knot_vector_v[4] = 1.0;
        knot_vector_v[5] = 1.0;
        knot_vector_v[6] = 1.0;

        // Assign knots of the basis along v
        SizeType number_of_knots_w = 10;
        Vector knot_vector_w(number_of_knots_w);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 0.0;
        knot_vector_w[2] = 0.0;
        knot_vector_w[3] = 0.0;
        knot_vector_w[4] = 0.0;
        knot_vector_w[5] = 1.0;
        knot_vector_w[6] = 1.0;
        knot_vector_w[7] = 1.0;
        knot_vector_w[8] = 1.0;
        knot_vector_w[9] = 1.0;

        Vector weights = ZeroVector(360);

        return NurbsVolumeGeometry<3, PointerVector<Point>>(
            points, polynomial_degree_u, polynomial_degree_v, polynomial_degree_w,
                knot_vector_u, knot_vector_v, knot_vector_w);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeGeometryPyramid, KratosCoreNurbsGeometriesFastSuite) {
        NurbsVolumeGeometry<3, PointerVector<Point>> TruncatedPyramid = GenerateUniformTruncatedPyramid();

        typename Geometry<Node<3>>::IntegrationPointsArrayType integration_points;
        TruncatedPyramid.CreateIntegrationPoints(integration_points);
        //KRATOS_CHECK_EQUAL( integration_points.size(), 45);
        // Compute and check volume
        double volume = 0;
        for (IndexType i = 0; i < integration_points.size(); ++i) {
            volume += integration_points[i].Weight() * TruncatedPyramid.DeterminantOfJacobian(integration_points[i].Coordinates());
        }
        KRATOS_CHECK_NEAR(volume, 32.21333333333333, TOLERANCE);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 1.0;
        parameter[1] = 1.0;
        parameter[2] = 1.0;
        array_1d<double, 3> global_coordinates(0.0);

        TruncatedPyramid.GlobalCoordinates(global_coordinates, parameter);
        // Check Coordinates
        KRATOS_CHECK_NEAR(global_coordinates[0],1.8,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[1],1.8,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[2],4,TOLERANCE);

        std::vector<Point::CoordinatesArrayType> derivatives;
        TruncatedPyramid.GlobalSpaceDerivatives(derivatives, parameter, 2);
        // Check dN/dx
        KRATOS_CHECK_NEAR(derivatives[1][0], 3.6, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[1][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[1][2], 0.0, TOLERANCE);
        // Check dN/dy
        KRATOS_CHECK_NEAR(derivatives[2][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[2][1], 3.6, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[2][2], 0.0, TOLERANCE);
        // Check dN/dz
        KRATOS_CHECK_NEAR(derivatives[3][0], 0.8, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[3][1], 0.8, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[3][2], 4.0, TOLERANCE);
        // Check dN2/dx2
        KRATOS_CHECK_NEAR(derivatives[4][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[4][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[4][2], 0.0, TOLERANCE);

        parameter[0] = 0.5;
        parameter[1] = 1.0;
        parameter[2] = 0.5;

        TruncatedPyramid.GlobalCoordinates(global_coordinates, parameter);
        // Check Coordinates
        KRATOS_CHECK_NEAR(global_coordinates[0],0.0,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[1],1.4,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[2],2.0,TOLERANCE);

        TruncatedPyramid.GlobalSpaceDerivatives(derivatives, parameter, 2);
        // Check dN/dx
        KRATOS_CHECK_NEAR(derivatives[1][0], 2.8, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[1][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[1][2], 0.0, TOLERANCE);
        // Check dN/dy
        KRATOS_CHECK_NEAR(derivatives[2][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[2][1], 2.8, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[2][2], 0.0, TOLERANCE);
        // Check dN/dz
        KRATOS_CHECK_NEAR(derivatives[3][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[3][1], 0.8, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[3][2], 4.0, TOLERANCE);
        // Check dN2/dx2
        KRATOS_CHECK_NEAR(derivatives[4][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[4][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[4][2], 0.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeGeometryDistortedCube, KratosCoreNurbsGeometriesFastSuite) {
        NurbsVolumeGeometry<3, PointerVector<Point>> DistortedCube = GenerateDistortedCube();

        typename Geometry<Node<3>>::IntegrationPointsArrayType integration_points;
        DistortedCube.CreateIntegrationPoints(integration_points);
        KRATOS_CHECK_EQUAL(DistortedCube.Dimension(), 3);
        KRATOS_CHECK_EQUAL(DistortedCube.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(DistortedCube.LocalSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(DistortedCube.IsRational(), false);
        KRATOS_CHECK_EQUAL(DistortedCube.PointsNumber(), 100 );
        // Compute and check volume
        double volume = 0;
        for (IndexType i = 0; i < integration_points.size(); ++i) {
            volume += integration_points[i].Weight() * DistortedCube.DeterminantOfJacobian(integration_points[i].Coordinates());
        }
        KRATOS_CHECK_NEAR(volume, 44.3259259259, TOLERANCE);

        // Check the local to global space mapping.
        array_1d<double, 3> parameter(0.0);
        parameter[0] = 0.89;
        parameter[1] = 0.213;
        parameter[2] = 0.33;
        array_1d<double, 3> global_coordinates(0.0);
        DistortedCube.GlobalCoordinates(global_coordinates, parameter);
        // Check Coordinates
        KRATOS_CHECK_NEAR(global_coordinates[0],1.2490666666666668,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[1],-0.5280889485640089,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[2],1.32,TOLERANCE);

        parameter[0] = 0.5;
        parameter[1] = 1.0;
        parameter[2] = 0.5;
        DistortedCube.GlobalCoordinates(global_coordinates, parameter);
        // Check Coordinates
        KRATOS_CHECK_NEAR(global_coordinates[0],0.0,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[1],1.9166666666667,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[2],1.9999999999999,TOLERANCE);

        // Check derivatives
        std::vector<Point::CoordinatesArrayType> derivatives;
        DistortedCube.GlobalSpaceDerivatives(derivatives, parameter, 2);

        // First Order
        // Check dN/dx
        KRATOS_CHECK_NEAR(derivatives[1][0], 2.666666666667, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[1][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[1][2], 0.0, TOLERANCE);
        // Check dN/dy
        KRATOS_CHECK_NEAR(derivatives[2][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[2][1], 5.111111111111, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[2][2], 0.0, TOLERANCE);
        // Check dN/dz
        KRATOS_CHECK_NEAR(derivatives[3][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[3][1], 1.333333333333, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[3][2], 4.0, TOLERANCE);
        // Second Order
        // Check dN2/dx2
        KRATOS_CHECK_NEAR(derivatives[4][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[4][1], -14.66666666667, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[4][2], 0.0, TOLERANCE);
        // Check dN2/ dx dy
        KRATOS_CHECK_NEAR(derivatives[5][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[5][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[5][2], 0.0, TOLERANCE);
        // Check dN2/ dx dz
        KRATOS_CHECK_NEAR(derivatives[6][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[6][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[6][2], 0.0, TOLERANCE);
        // Check dN2/ dy2
        KRATOS_CHECK_NEAR(derivatives[7][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[7][1], 5.1111111111, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[7][2], 0.0, TOLERANCE);
        // Check dN2/ dy dz
        KRATOS_CHECK_NEAR(derivatives[8][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[8][1], 3.55555555556, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[8][2], 0.0, TOLERANCE);
        // Check dN2/ dz2
        KRATOS_CHECK_NEAR(derivatives[9][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[9][1], 4.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[9][2], 0.0, TOLERANCE);
    }


} // End namespace Testsing
} // End namespace Kratos