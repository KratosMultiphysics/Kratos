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
#include "geometries/nurbs_shape_function_utilities/knot_insertion.h"

#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/point.h"

namespace Kratos {
namespace Testing {

typedef Node<3> NodeType;

    NurbsVolumeGeometry<PointerVector<NodeType>> GenerateVolume() {
        // Construct Truncated Pyramid with: lower_base = 2x2; uper_base = 1.8x1.8; heigth = 4
        PointerVector<NodeType> points;
        double t = 1.0;
        std::vector<double> z_direction = {0, 1.0};
        std::vector<double> y_direction = {0, 1.0};
        std::vector<double> x_direction = {0, 0.333, 0.444, 0.555, 0.666, 1.0};
        std::size_t id = 1;
        for( auto i : z_direction){
            for( auto j : y_direction) {
                for( auto k : x_direction) {
                    points.push_back(NodeType::Pointer(new NodeType(id, k, j, i)));
                    id++;
                }
            }
        }
        // Polynomial orders
        SizeType polynomial_degree_u = 3;
        SizeType polynomial_degree_v = 1;
        SizeType polynomial_degree_w = 1;

        // TODO: Add inner knots
        // Assign knots of the basis along u
        SizeType number_of_knots_u = 10;
        Vector knot_vector_u(number_of_knots_u);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 0.0;
        knot_vector_u[3] = 0.0;
        knot_vector_u[4] = 0.2;
        knot_vector_u[5] = 0.5;
        knot_vector_u[6] = 1.0;
        knot_vector_u[7] = 1.0;
        knot_vector_u[8] = 1.0;
        knot_vector_u[9] = 1.0;


        // Assign knots of the basis along v
        SizeType number_of_knots_v = 4;
        Vector knot_vector_v(number_of_knots_v);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 0.0;
        knot_vector_v[2] = 1.0;
        knot_vector_v[3] = 1.0;

        // Assign knots of the basis along v
        SizeType number_of_knots_w = 4;
        Vector knot_vector_w(number_of_knots_w);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 0.0;
        knot_vector_w[2] = 1.0;
        knot_vector_w[3] = 1.0;

        return NurbsVolumeGeometry<PointerVector<NodeType>>(
            points, polynomial_degree_u, polynomial_degree_v, polynomial_degree_w,
                knot_vector_u, knot_vector_v, knot_vector_w);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsKnotRefinement, KratosCoreNurbsGeometriesFastSuite) {
        auto volume = GenerateVolume();

        KnotInsertion<PointerVector<NodeType>>::Insert_knot_u(volume, 0.1, 2);
        KnotInsertion<PointerVector<NodeType>>::Insert_knot_u(volume, 0.22, 2);
        KnotInsertion<PointerVector<NodeType>>::Insert_knot_u(volume, 0.7, 1);
        // (0.1, 0.1, 0.22, 0.22, 0.7)
        for( int i = 0; i < volume.size(); ++i){
            std::cout << "cps: " << i << " - " << volume[i] << std::endl;
        }
        std::cout << "Knotvector: " << volume.KnotsU() << std::endl;
        std::cout << "\n" << std::endl;
        // Check positions of CP's
        std::size_t k = 0;
        for( std::size_t i = 0; i < 4; ++i){
            KRATOS_CHECK_NEAR(volume[k+0][0], 0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+1][0], 0.1665, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+2][0], 0.26085, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+3][0], 0.37518, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+4][0], 0.39915599999999996, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+5][0], 0.41993519999999995, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+6][0], 0.470653875, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+7][0], 0.524266875, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+8][0], 0.624375, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+9][0], 0.7996, TOLERANCE);
            KRATOS_CHECK_NEAR(volume[k+10][0], 1.0, TOLERANCE);
            k += 11;
        }

        auto volume_2 = GenerateVolume();
        std::vector<double> u = {0.1, 0.1, 0.22, 0.22, 0.7};
        KnotInsertion<PointerVector<NodeType>>::KnotRefinement(volume_2, u);
        // Check positions of CP's
        k = 0;
        for( std::size_t i = 0; i < 4; ++i){
            KRATOS_CHECK_NEAR(volume_2[k+0][0], 0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+1][0], 0.1665, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+2][0], 0.26085, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+3][0], 0.37518, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+4][0], 0.39915599999999996, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+5][0], 0.41993519999999995, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+6][0], 0.470653875, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+7][0], 0.524266875, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+8][0], 0.624375, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+9][0], 0.7996, TOLERANCE);
            KRATOS_CHECK_NEAR(volume_2[k+10][0], 1.0, TOLERANCE);
            k += 11;
        }

        for( int i = 0; i < volume_2.size(); ++i){
            std::cout << "cps: " << i << " - " << volume_2[i] << " ref: " << volume[i] << std::endl;
        }
        std::cout << "Knotvector: " << volume_2.KnotsU() << std::endl;
    }
}
}