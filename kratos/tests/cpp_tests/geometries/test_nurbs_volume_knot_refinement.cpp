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
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_utilities.h"

#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/point.h"

namespace Kratos {
namespace Testing {

typedef Node<3> NodeType;
typedef NurbsVolumeGeometry<PointerVector<NodeType>> NurbsVolumeGeometryType;
typedef typename NurbsVolumeGeometry<PointerVector<NodeType>>::Pointer NurbsVolumeGeometryPointerType;

NurbsVolumeGeometryPointerType GenerateVolume(IndexType Direction) {
    PointerVector<NodeType> points;
    std::vector<double> x_direction = {0, 1.0};
    std::vector<double> y_direction = {0, 1.0};
    std::vector<double> z_direction = {0, 1.0};
    if( Direction == 0){
        x_direction = {0, 0.333, 0.444, 0.555, 0.666, 1.0};
    } else if( Direction == 1){
        y_direction = {0, 0.333, 0.444, 0.555, 0.666, 1.0};
    } else {
        z_direction = {0, 0.333, 0.444, 0.555, 0.666, 1.0};
    }

    std::size_t id = 1;
    for( auto i : z_direction){
        for( auto j : y_direction) {
            for( auto k : x_direction) {
                points.push_back(Kratos::make_intrusive<NodeType>(id, k, j, i));
                id++;
            }
        }
    }
    // Polynomial orders
    SizeType polynomial_degree_1 = 3;
    SizeType polynomial_degree_2 = 1;

    // Knots
    SizeType number_of_knots_1 = 10;
    Vector knot_vector_1(number_of_knots_1);
    knot_vector_1[0] = 0.0;
    knot_vector_1[1] = 0.0;
    knot_vector_1[2] = 0.0;
    knot_vector_1[3] = 0.0;
    knot_vector_1[4] = 0.2;
    knot_vector_1[5] = 0.5;
    knot_vector_1[6] = 1.0;
    knot_vector_1[7] = 1.0;
    knot_vector_1[8] = 1.0;
    knot_vector_1[9] = 1.0;


    // Assign knots of the basis along v
    SizeType number_of_knots_2 = 4;
    Vector knot_vector_2(number_of_knots_2);
    knot_vector_2[0] = 0.0;
    knot_vector_2[1] = 0.0;
    knot_vector_2[2] = 1.0;
    knot_vector_2[3] = 1.0;

    if( Direction == 0){
        return Kratos::make_shared<NurbsVolumeGeometryType>(
            points, polynomial_degree_1, polynomial_degree_2, polynomial_degree_2,
            knot_vector_1, knot_vector_2, knot_vector_2);
    } else if( Direction == 1){
        return Kratos::make_shared<NurbsVolumeGeometryType>(
            points, polynomial_degree_2, polynomial_degree_1, polynomial_degree_2,
            knot_vector_2, knot_vector_1, knot_vector_2);
    } else {
        return Kratos::make_shared<NurbsVolumeGeometryType>(
            points, polynomial_degree_2, polynomial_degree_2, polynomial_degree_1,
            knot_vector_2, knot_vector_2, knot_vector_1);
    }

}

KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeKnotRefinementU, KratosCoreNurbsGeometriesFastSuite) {

    auto volume_u = GenerateVolume(0);
    std::vector<double> insert_knots = {0.1, 0.1, 0.22, 0.22, 0.7};
    auto refined_volume_u = NurbsVolumeUtilities::KnotRefinementU(*volume_u, insert_knots);

    // Check polynomial degree
    KRATOS_CHECK_EQUAL(refined_volume_u->PolynomialDegreeU(),3);
    KRATOS_CHECK_EQUAL(refined_volume_u->PolynomialDegreeV(),1);
    KRATOS_CHECK_EQUAL(refined_volume_u->PolynomialDegreeW(),1);
    // Check knot vectors
    std::vector<double> knots_u_ref = {0, 0, 0, 0.1, 0.1, 0.2, 0.22, 0.22, 0.5, 0.7, 1.0, 1.0, 1.0};
    std::vector<double> knots_v_ref = {0, 1.0};
    std::vector<double> knots_w_ref = {0, 1.0};
    KRATOS_CHECK_VECTOR_NEAR(refined_volume_u->KnotsU(), knots_u_ref, TOLERANCE);
    KRATOS_CHECK_VECTOR_NEAR(refined_volume_u->KnotsV(), knots_v_ref, TOLERANCE);
    KRATOS_CHECK_VECTOR_NEAR(refined_volume_u->KnotsW(), knots_w_ref, TOLERANCE);
    // Check number of CP's
    KRATOS_CHECK_EQUAL(refined_volume_u->NumberOfControlPointsU(),11);
    KRATOS_CHECK_EQUAL(refined_volume_u->NumberOfControlPointsV(),2);
    KRATOS_CHECK_EQUAL(refined_volume_u->NumberOfControlPointsW(),2);

    // Check positions of CP's
    IndexType k = 0;
    for( std::size_t i = 0; i < 4; ++i){
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+0][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+1][0], 0.1665, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+2][0], 0.26085, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+3][0], 0.37518, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+4][0], 0.39915599999999996, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+5][0], 0.41993519999999995, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+6][0], 0.470653875, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+7][0], 0.524266875, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+8][0], 0.624375, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+9][0], 0.7996, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_u)[k+10][0], 1.0, TOLERANCE);
        k += 11;
    }
    // Check ID's
    for( long unsigned int i = 0; i < refined_volume_u->size(); ++i){
        KRATOS_CHECK_EQUAL((*refined_volume_u)[i].Id(), i+1);
    }
}

KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeKnotRefinementV, KratosCoreNurbsGeometriesFastSuite) {
    auto volume_v = GenerateVolume(1);
    std::vector<double> insert_knots = {0.1, 0.1, 0.22, 0.22, 0.7};
    auto refined_volume_v = NurbsVolumeUtilities::KnotRefinementV(*volume_v, insert_knots);

    // Check polynomial degree
    KRATOS_CHECK_EQUAL(refined_volume_v->PolynomialDegreeU(),1);
    KRATOS_CHECK_EQUAL(refined_volume_v->PolynomialDegreeV(),3);
    KRATOS_CHECK_EQUAL(refined_volume_v->PolynomialDegreeW(),1);
    // Check knot vectors
    std::vector<double> knots_u_ref = {0, 1.0};
    std::vector<double> knots_v_ref = {0, 0, 0, 0.1, 0.1, 0.2, 0.22, 0.22, 0.5, 0.7, 1.0, 1.0, 1.0};
    std::vector<double> knots_w_ref = {0, 1.0};
    KRATOS_CHECK_VECTOR_NEAR(refined_volume_v->KnotsU(), knots_u_ref, TOLERANCE);
    KRATOS_CHECK_VECTOR_NEAR(refined_volume_v->KnotsV(), knots_v_ref, TOLERANCE);
    KRATOS_CHECK_VECTOR_NEAR(refined_volume_v->KnotsW(), knots_w_ref, TOLERANCE);
    // Check number of CP's
    KRATOS_CHECK_EQUAL(refined_volume_v->NumberOfControlPointsU(),2);
    KRATOS_CHECK_EQUAL(refined_volume_v->NumberOfControlPointsV(),11);
    KRATOS_CHECK_EQUAL(refined_volume_v->NumberOfControlPointsW(),2);
    // Check position of CP's
    IndexType k = 0;
    for( std::size_t i = 0; i < 2; ++i){
        for( std::size_t j = 0; j < 2; ++j){
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+0][1], 0.0, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+2][1], 0.1665, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+4][1], 0.26085, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+6][1], 0.37518, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+8][1], 0.39915599999999996, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+10][1], 0.41993519999999995, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+12][1], 0.470653875, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+14][1], 0.524266875, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+16][1], 0.624375, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+18][1], 0.7996, TOLERANCE);
            KRATOS_CHECK_NEAR((*refined_volume_v)[k+j+20][1], 1.0, TOLERANCE);
        }
        k += 22;
    }
    // Check ID's
    for( long unsigned int i = 0; i < refined_volume_v->size(); ++i){
        KRATOS_CHECK_EQUAL((*refined_volume_v)[i].Id(), i+1);
    }
}

KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeKnotRefinementW, KratosCoreNurbsGeometriesFastSuite) {
    auto volume_w = GenerateVolume(2);
    std::vector<double> insert_knots = {0.1, 0.1, 0.22, 0.22, 0.7};
    auto refined_volume_w = NurbsVolumeUtilities::KnotRefinementW(*volume_w, insert_knots);

    // Check polynomial degree
    KRATOS_CHECK_EQUAL(refined_volume_w->PolynomialDegreeU(),1);
    KRATOS_CHECK_EQUAL(refined_volume_w->PolynomialDegreeV(),1);
    KRATOS_CHECK_EQUAL(refined_volume_w->PolynomialDegreeW(),3);
    // Check knot vectors
    std::vector<double> knots_u_ref = {0, 1.0};
    std::vector<double> knots_v_ref = {0, 1.0};
    std::vector<double> knots_w_ref = {0, 0, 0, 0.1, 0.1, 0.2, 0.22, 0.22, 0.5, 0.7, 1.0, 1.0, 1.0};
    KRATOS_CHECK_VECTOR_NEAR(refined_volume_w->KnotsU(), knots_u_ref, TOLERANCE);
    KRATOS_CHECK_VECTOR_NEAR(refined_volume_w->KnotsV(), knots_v_ref, TOLERANCE);
    KRATOS_CHECK_VECTOR_NEAR(refined_volume_w->KnotsW(), knots_w_ref, TOLERANCE);
    // Check number of CP's
    KRATOS_CHECK_EQUAL(refined_volume_w->NumberOfControlPointsU(),2);
    KRATOS_CHECK_EQUAL(refined_volume_w->NumberOfControlPointsV(),2);
    KRATOS_CHECK_EQUAL(refined_volume_w->NumberOfControlPointsW(),11);

    // Check position of CP's
    for( std::size_t i = 0; i < 4; ++i){
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+0][2], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+4][2], 0.1665, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+8][2], 0.26085, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+12][2], 0.37518, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+16][2], 0.39915599999999996, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+20][2], 0.41993519999999995, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+24][2], 0.470653875, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+28][2], 0.524266875, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+32][2], 0.624375, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+36][2], 0.7996, TOLERANCE);
        KRATOS_CHECK_NEAR((*refined_volume_w)[i+40][2], 1.0, TOLERANCE);
    }

    // Check ID's
    for( long unsigned int i = 0; i < refined_volume_w->size(); ++i){
        KRATOS_CHECK_EQUAL((*refined_volume_w)[i].Id(), i+1);
    }
}

} // End namespace testing
} // End namespace kratos