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

    NurbsVolumeGeometry<3, PointerVector<Point>> GenerateVolume() {
        
        PointerVector<Point> points;
        int discretization = 3;
        int length = 1;
        double scale = (double)length / ((double)discretization - 1.0);
        for( int i = 0; i < discretization; ++i){
            for( int j = 0; j < discretization; ++j) {
                for( int k = 0; k < discretization; ++k) {
                    if( k== 1){
                        points.push_back(Point::Pointer(new Point(i*scale, (j+0.2)*scale, k*scale)));
                    }
                    else {
                        points.push_back(Point::Pointer(new Point(i*scale, j*scale, k*scale)));
                    }
                }
            }
        }
        // Polynomial orders
        SizeType polynomial_degree_u = 2;
        SizeType polynomial_degree_v = 2;
        SizeType polynomial_degree_w = 2;

        // Assign the derivative order of the basis
        SizeType derivative_order = 0;

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
        SizeType number_of_knots_w = 6;
        Vector knot_vector_w(number_of_knots_w);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 0.0;
        knot_vector_w[2] = 0.0;
        knot_vector_w[3] = 1.0;
        knot_vector_w[4] = 1.0;
        knot_vector_w[5] = 1.0;
        Vector weights = ZeroVector(36);

        return NurbsVolumeGeometry<3, PointerVector<Point>>(
            points, polynomial_degree_u, polynomial_degree_v, polynomial_degree_w,
                knot_vector_u, knot_vector_v, knot_vector_w);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeGeometry, KratosCoreNurbsGeometriesFastSuite) {
        NurbsVolumeGeometry<3, PointerVector<Point>> volume = GenerateVolume();

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 0.223;
        parameter[1] = 0.88;
        parameter[2] = 0.12;
        array_1d<double, 3> result(0.0);

        volume.GlobalCoordinates(result, parameter);
        std::cout << "volume: " << result[0] << " " << result[1] << " " << result[2] << std::endl; 
        
    }

} // End namespace Testsing
} // End namespace Kratos