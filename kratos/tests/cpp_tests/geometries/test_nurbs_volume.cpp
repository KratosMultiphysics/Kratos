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
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/point.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;

    NurbsVolumeGeometry<PointerVector<NodeType>> GenerateTruncatedPyramid() {
        // Construct Truncated Pyramid with: lower_base = 2x2; uper_base = 1.8x1.8; heigth = 4
        PointerVector<NodeType> points(196);
        double t = 1.0;
        std::vector<double> z_direction = {0, 2.0/3.0, 4.0/3.0, 2.0, 8.0/3.0, 10/3.0, 4.0};
        std::vector<double> y_direction = {-1.0, -1.0/3.0, 1.0/3.0, 1.0};
        std::vector<double> x_direction = {-1.0, -2.0/3.0, -1.0/3.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0};
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
        SizeType polynomial_degree_u = 3;
        SizeType polynomial_degree_v = 2;
        SizeType polynomial_degree_w = 4;

        // Assign knots of the basis along u.
        SizeType number_of_knots_u = 11;
        Vector knot_vector_u(number_of_knots_u);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 0.0;
        knot_vector_u[3] = 0.0;
        knot_vector_u[4] = 0.25;
        knot_vector_u[5] = 0.6;
        knot_vector_u[6] = 0.85;
        knot_vector_u[7] = 1.0;
        knot_vector_u[8] = 1.0;
        knot_vector_u[9] = 1.0;
        knot_vector_u[10] = 1.0;

        // Assign knots of the basis along v.
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
        SizeType number_of_knots_w = 12;
        Vector knot_vector_w(number_of_knots_w);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 0.0;
        knot_vector_w[2] = 0.0;
        knot_vector_w[3] = 0.0;
        knot_vector_w[4] = 0.0;
        knot_vector_w[5] = 1.0/3.0;
        knot_vector_w[6] = 2.0/3.0;
        knot_vector_w[7] = 1.0;
        knot_vector_w[8] = 1.0;
        knot_vector_w[9] = 1.0;
        knot_vector_w[10] = 1.0;
        knot_vector_w[11] = 1.0;

        return NurbsVolumeGeometry<PointerVector<NodeType>>(
            points, polynomial_degree_u, polynomial_degree_v, polynomial_degree_w,
                knot_vector_u, knot_vector_v, knot_vector_w);
    }

    NurbsVolumeGeometry<PointerVector<Point>> GenerateDistortedCube() {
        // Construct a distroted cube.
        PointerVector<Point> points(100);
        std::vector<double> y_direction = {-1.0, -1.0/3.0, 1.0/3.0, 1.0};
        double t = 0.8;
        int index = 0;
        for( int i = 0; i <=4; ++i){
            t += 0.2;
            for( auto j : y_direction) {
                for( int k = -2; k <=2; ++k ) {
                    double x = k;
                    double y = j;
                    double z = i;
                    if( k == 0)
                        y = j*std::max(i,2);
                    if( j == 0)
                        x = k*std::max(i,2);
                    points(index) = Kratos::make_shared<Point>(x, y, z);
                    index++;
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

        return NurbsVolumeGeometry<PointerVector<Point>>(
            points, polynomial_degree_u, polynomial_degree_v, polynomial_degree_w,
                knot_vector_u, knot_vector_v, knot_vector_w);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeGeometryIntegrationPoints1, KratosCoreNurbsGeometriesFastSuite) {
            NurbsVolumeGeometry<PointerVector<NodeType>> TruncatedPyramid = GenerateTruncatedPyramid();

            KRATOS_CHECK_EQUAL(TruncatedPyramid.Dimension(), 3);
            KRATOS_CHECK_EQUAL(TruncatedPyramid.WorkingSpaceDimension(), 3);
            KRATOS_CHECK_EQUAL(TruncatedPyramid.LocalSpaceDimension(), 3);
            KRATOS_CHECK_EQUAL(TruncatedPyramid.IsRational(), false);
            KRATOS_CHECK_EQUAL(TruncatedPyramid.PointsNumber(), 196 );

            typename Geometry<NodeType>::IntegrationPointsArrayType integration_points;
            IntegrationInfo integration_info = TruncatedPyramid.GetDefaultIntegrationInfo();
            TruncatedPyramid.CreateIntegrationPoints(integration_points, integration_info);
            KRATOS_CHECK_EQUAL( integration_points.size(), 1440);
            // Compute and check volume
            double volume = 0;
            for (IndexType i = 0; i < integration_points.size(); ++i) {
                volume += integration_points[i].Weight() * TruncatedPyramid.DeterminantOfJacobian(integration_points[i].Coordinates());
            }
            KRATOS_CHECK_NEAR(volume, 32.21333333333333, TOLERANCE);

            const int geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
            const int geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Volume;
            KRATOS_CHECK_EQUAL(TruncatedPyramid.GetGeometryFamily(), geometry_family);
            KRATOS_CHECK_EQUAL(TruncatedPyramid.GetGeometryType(), geometry_type);
        }

    KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeGeometryEvaluation1, KratosCoreNurbsGeometriesFastSuite) {
        NurbsVolumeGeometry<PointerVector<NodeType>> TruncatedPyramid = GenerateTruncatedPyramid();

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 1.0;
        parameter[1] = 1.0;
        parameter[2] = 1.0;
        array_1d<double, 3> global_coordinates(0.0);

        TruncatedPyramid.GlobalCoordinates(global_coordinates, parameter);
        // Check results. All values are compared to the Python NURBS library geomdl.
        // https://nurbs-python.readthedocs.io/en/latest/module_nurbs.html
        // Check Coordinates
        KRATOS_CHECK_NEAR(global_coordinates[0],1.8,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[1],1.8,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[2],4,TOLERANCE);

        std::vector<Point::CoordinatesArrayType> derivatives;
        TruncatedPyramid.GlobalSpaceDerivatives(derivatives, parameter, 4);
        // First order
        // Check dN/dx
        KRATOS_CHECK_NEAR(derivatives[1][0], 12.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[1][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[1][2], 0.0, TOLERANCE);
        // Check dN/dy
        KRATOS_CHECK_NEAR(derivatives[2][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[2][1], 4.8, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[2][2], 0.0, TOLERANCE);
        // Check dN/dz
        KRATOS_CHECK_NEAR(derivatives[3][0], 1.6, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[3][1], 1.6, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[3][2], 8.0, TOLERANCE);
        // Second Order
        // Check dN2/dx2
        KRATOS_CHECK_NEAR(derivatives[4][0], 100, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[4][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[4][2], 0.0, TOLERANCE);
        // Check dN2/dx dy
        KRATOS_CHECK_NEAR(derivatives[5][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[5][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[5][2], 0.0, TOLERANCE);
        // Check dN2/dx dz
        KRATOS_CHECK_NEAR(derivatives[6][0], 10.6666666667, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[6][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[6][2], 0.0, TOLERANCE);
        // Check dN2/dy2
        KRATOS_CHECK_NEAR(derivatives[7][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[7][1], 4.8, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[7][2], 0.0, TOLERANCE);
        // Check dN2/dy dz
        KRATOS_CHECK_NEAR(derivatives[8][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[8][1], 4.266666666667, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[8][2], 0.0, TOLERANCE);
        // Check dN2/dz2
        KRATOS_CHECK_NEAR(derivatives[9][0], 7.2, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[9][1], 7.2, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[9][2], 36.0, TOLERANCE);
        // Third order
        // Check dN3/dx3
        KRATOS_CHECK_NEAR(derivatives[10][0], 596.666666666667, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[10][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[10][2], 0.0, TOLERANCE);
        // Check dN3/dx2 dy
        KRATOS_CHECK_NEAR(derivatives[11][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[11][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[11][2], 0.0, TOLERANCE);
        // Check dN3/dx2 dz
        KRATOS_CHECK_NEAR(derivatives[12][0], 88.8888888888887, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[12][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[12][2], 0.0, TOLERANCE);
        // Check dN3/dx dz2
        KRATOS_CHECK_NEAR(derivatives[15][0], 48.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[15][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[15][2], 0.0, TOLERANCE);
        // Check dN3/dy2 dz
        KRATOS_CHECK_NEAR(derivatives[17][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[17][1], 4.2666666666667, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[17][2], 0.0, TOLERANCE);
        // Check dN3/dy dz2
        KRATOS_CHECK_NEAR(derivatives[18][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[18][1], 19.2, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[18][2], 0.0, TOLERANCE);
        // Check dN3/dz3
        KRATOS_CHECK_NEAR(derivatives[19][0], 36.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[19][1], 36.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[19][2], 180.0, TOLERANCE);
        // Remaining derivative on this level must be zero:
        std::vector<int> third_order_zero_derivatives = {13, 14, 16};
        for( auto index : third_order_zero_derivatives){
            for( std::size_t i = 0; i < 3; ++i){
                KRATOS_CHECK_NEAR(derivatives[index][i], 0.0, TOLERANCE);
            }
        }
        // Fourth Order
        // Check dN4/dx4
        KRATOS_CHECK_NEAR(derivatives[20][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[20][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[20][2], 0.0, TOLERANCE);
        // Check dN4/dx3 dz
        KRATOS_CHECK_NEAR(derivatives[22][0], 530.37037037037, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[22][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[22][2], 0.0, TOLERANCE);
        // Check dN4/dx2 dz2
        KRATOS_CHECK_NEAR(derivatives[25][0], 400.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[25][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[25][2], 0.0, TOLERANCE);
        // Check dN4/dx dz3
        KRATOS_CHECK_NEAR(derivatives[29][0], 240.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[29][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[29][2], 0.0, TOLERANCE);
        // Check dN4/dy2 dz2
        KRATOS_CHECK_NEAR(derivatives[32][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[32][1], 19.2, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[32][2], 0.0, TOLERANCE);
        // Check dN4/dy dz3
        KRATOS_CHECK_NEAR(derivatives[33][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[33][1], 96.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[33][2], 0.0, TOLERANCE);
        // Check dN4/dz4
        KRATOS_CHECK_NEAR(derivatives[34][0], 97.2, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[34][1], 97.2, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives[34][2], 486.0, TOLERANCE);
        // Remaining derivative on this level must be zero:
        std::vector<int> fourth_order_zero_derivatives = {23, 24, 26, 27, 28, 30, 31};
        for( auto index : fourth_order_zero_derivatives){
            for( std::size_t i = 0; i < 3; ++i){
                KRATOS_CHECK_NEAR(derivatives[index][i], 0.0, TOLERANCE);
            }
        }
        // Check another point inside the volume
        parameter[0] = 0.88672;
        parameter[1] = 0.231;
        parameter[2] = 0.664;

        TruncatedPyramid.GlobalCoordinates(global_coordinates, parameter);
        // Check Coordinates
        KRATOS_CHECK_NEAR(global_coordinates[0], 0.7776905789441982,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[1], -0.6794661290038271,TOLERANCE);
        KRATOS_CHECK_NEAR(global_coordinates[2], 2.4642328320000004,TOLERANCE);

        // Check kratos geometry family
        const int geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const int geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Volume;
        KRATOS_CHECK_EQUAL(TruncatedPyramid.GetGeometryFamily(), geometry_family);
        KRATOS_CHECK_EQUAL(TruncatedPyramid.GetGeometryType(), geometry_type);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeGeometryIntegrationPoints2, KratosCoreNurbsGeometriesFastSuite) {
        NurbsVolumeGeometry<PointerVector<Point>> DistortedCube = GenerateDistortedCube();

        typename Geometry<Point>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = DistortedCube.GetDefaultIntegrationInfo();
        DistortedCube.CreateIntegrationPoints(integration_points, integration_info);
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

        // Check kratos geometry family
        const int geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const int geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Volume;
        KRATOS_CHECK_EQUAL(DistortedCube.GetGeometryFamily(), geometry_family);
        KRATOS_CHECK_EQUAL(DistortedCube.GetGeometryType(), geometry_type);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeGeometryEvaluation2, KratosCoreNurbsGeometriesFastSuite) {
        NurbsVolumeGeometry<PointerVector<Point>> DistortedCube = GenerateDistortedCube();

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

        // Check kratos geometry family
        const int geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const int geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Volume;
        KRATOS_CHECK_EQUAL(DistortedCube.GetGeometryFamily(), geometry_family);
        KRATOS_CHECK_EQUAL(DistortedCube.GetGeometryType(), geometry_type);
    }

    /// Check quadrature point geometries of nurbs volume.
    KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeQuadraturePointGeometries, KratosCoreNurbsGeometriesFastSuite) {
        NurbsVolumeGeometry<PointerVector<NodeType>> pyramid = GenerateTruncatedPyramid();

        // Check general information, input to ouput
        typename Geometry<NodeType>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = pyramid.GetDefaultIntegrationInfo();
        pyramid.CreateIntegrationPoints(integration_points, integration_info);

        typename Geometry<NodeType>::GeometriesArrayType quadrature_points;
        pyramid.CreateQuadraturePointGeometries(quadrature_points, 3, integration_points, integration_info);

        KRATOS_CHECK_EQUAL(quadrature_points.size(), 1440);
        double sum = 0;
        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                sum += quadrature_points[i].IntegrationPoints()[j].Weight();
            }
        }
        KRATOS_CHECK_NEAR(sum, 1.0, TOLERANCE);

        auto element = Element(0, quadrature_points(2));

        // Check shape functions
        KRATOS_CHECK_MATRIX_NEAR(
            element.pGetGeometry()->ShapeFunctionsValues(),
            quadrature_points(2)->ShapeFunctionsValues(),
            TOLERANCE);

        // Check first derivatives
        KRATOS_CHECK_MATRIX_NEAR(
            element.GetGeometry().ShapeFunctionDerivatives(1, 0),
            quadrature_points(2)->ShapeFunctionLocalGradient(0),
            TOLERANCE);

        // Check second derivatives
        KRATOS_CHECK_MATRIX_NEAR(
            element.GetGeometry().ShapeFunctionDerivatives(2, 0),
            quadrature_points(2)->ShapeFunctionDerivatives(2, 0),
            TOLERANCE);

        // Check kratos geometry family
        const int geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const int geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Volume;
        KRATOS_CHECK_EQUAL(pyramid.GetGeometryFamily(), geometry_family);
        KRATOS_CHECK_EQUAL(pyramid.GetGeometryType(), geometry_type);
    }

} // End namespace Testsing
} // End namespace Kratos