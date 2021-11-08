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
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    /// Factory functions
    NurbsVolumeShapeFunction GenerateVolumeShapeFunctions()
    {
        // Polynomial orders
        SizeType polynomial_degree_u = 3;
        SizeType polynomial_degree_v = 2;
        SizeType polynomial_degree_w = 2;

        // Assign the derivative order of the basis
        SizeType derivative_order = 2;

        // Assign knots of the basis along u
        SizeType number_of_knots_u = 11;
        Vector knot_vector_u(number_of_knots_u);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 0.0;
        knot_vector_u[3] = 0.0;
        knot_vector_u[4] = 0.5;
        knot_vector_u[5] = 0.5;
        knot_vector_u[6] = 0.75;
        knot_vector_u[7] = 1.0;
        knot_vector_u[8] = 1.0;
        knot_vector_u[9] = 1.0;
        knot_vector_u[10] = 1.0;

        // Assign knots of the basis along v
        SizeType number_of_knots_v = 7;
        Vector knot_vector_v(number_of_knots_v);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 0.0;
        knot_vector_v[2] = 0.0;
        knot_vector_v[3] = 0.33;
        knot_vector_v[4] = 1.0;
        knot_vector_v[5] = 1.0;
        knot_vector_v[6] = 1.0;

        // Assign knots of the basis along v
        SizeType number_of_knots_w = 7;
        Vector knot_vector_w(number_of_knots_w);
        knot_vector_w[0] = 0.0;
        knot_vector_w[1] = 0.0;
        knot_vector_w[2] = 0.0;
        knot_vector_w[3] = 0.4;
        knot_vector_w[4] = 1.0;
        knot_vector_w[5] = 1.0;
        knot_vector_w[6] = 1.0;

        // Parameter on the surface where to compute the basis functions and their derivatives
        array_1d<double, 3> parameter(0.0);
        parameter[0] = 0.5;
        parameter[1] = 0.5;
        parameter[2] = 0.5;

        // Create the shape function container
        NurbsVolumeShapeFunction shape_functions(polynomial_degree_u, polynomial_degree_v, polynomial_degree_w, derivative_order);

        // Get the knot span where parameter parameterT lies in along u
        const int knot_span_u = NurbsUtilities::GetLowerSpan(shape_functions.PolynomialDegreeU(), knot_vector_u, parameter[0]);

        // Get the knot span where parameter parameterT lies in along v
        const int knot_span_v = NurbsUtilities::GetLowerSpan(shape_functions.PolynomialDegreeV(), knot_vector_v, parameter[1]);

        // Get the knot span where parameter parameterT lies in along w
        const int knot_span_w = NurbsUtilities::GetLowerSpan(shape_functions.PolynomialDegreeW(), knot_vector_w, parameter[2]);

        // Compute the shape functions at the given parametric location
        shape_functions.ComputeBSplineShapeFunctionValuesAtSpan(knot_vector_u, knot_vector_v, knot_vector_w,
            knot_span_u, knot_span_v, knot_span_w, parameter[0], parameter[1], parameter[2]);

        // Return the object
        return shape_functions;
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeShapeFunctions, KratosCoreNurbsGeometriesFastSuite) {
        NurbsVolumeShapeFunction shape_functions = GenerateVolumeShapeFunctions();

        KRATOS_CHECK_EQUAL( shape_functions.NumberOfNonzeroControlPoints(), 36);
        // Check Values
        double sum = 0.0;
        for( std::size_t i = 0; i < shape_functions.NumberOfNonzeroControlPoints(); ++i){
            sum += shape_functions(i,0);
        }
        KRATOS_CHECK_NEAR( sum, 1.0, TOLERANCE);

        // Derivatives
        //dN/du
        KRATOS_CHECK_NEAR( shape_functions(1,1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(4,1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(7,1), 0.9374767951, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(10,1), -0.1072993243, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(13,1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(16,1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(19,1), 1.2499690602, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(22,1), -0.1430657657, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(25,1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(28,1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(31,1), 0.0624984530, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(34,1), -0.0071532883, TOLERANCE);
        //dN/dv
        KRATOS_CHECK_NEAR( shape_functions(1,2), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(4,2), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(7,2), 0.2042028663, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(10,2), 0.1051954160, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(13,2), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(16,2), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(19,2), 0.2722704884, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(22,2), 0.1402605546, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(25,2), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(28,2), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(31,2), 0.0136135244, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(34,2), 0.0070130277, TOLERANCE);
        //dN/dw
        KRATOS_CHECK_NEAR( shape_functions(1,3), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(4,3), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(7,3), -0.6249845301, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(10,3), -0.0357664414, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(13,3), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(16,3), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(19,3), 0.4166563534, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(22,3), 0.0238442943, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(25,3), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(28,3), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(31,3), 0.2083281767, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(34,3), 0.0119221471, TOLERANCE);
        //dN2/du2
        KRATOS_CHECK_NEAR( shape_functions(1,4), 3.7313432836, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(4,4), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(7,4), 3.7499071805, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(10,4), -1.0729932427, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(13,4), 4.9751243781, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(16,4), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(19,4), 4.9998762407, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(22,4), -1.4306576570, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(25,4), 0.2487562189, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(28,4), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(31,4), 0.2499938120, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(34,4), -0.0715328828, TOLERANCE);
        //dN2/du dv
        KRATOS_CHECK_NEAR( shape_functions(1,5), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(4,5), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(7,5), 1.2252171976, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(10,5), -1.2623449915, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(13,5), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(16,5), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(19,5), 1.6336229301, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(22,5), -1.6831266553, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(25,5), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(28,5), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(31,5), 0.0816811465, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(34,5), -0.0841563328, TOLERANCE);
        //dN2/du dw
        KRATOS_CHECK_NEAR( shape_functions(1,6), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(4,6), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(7,6), -3.7499071805, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(10,6), 0.4291972971, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(13,6), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(16,6), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(19,6), 2.4999381203, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(22,6), -0.2861315314, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(25,6), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(28,6), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(31,6), 1.2499690602, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(34,6), -0.1430657657, TOLERANCE);
        //dN2/dv2
        KRATOS_CHECK_NEAR( shape_functions(1,7), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(4,7), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(7,7), -2.0667805252, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(10,7), 0.6187965644, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(13,7), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(16,7), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(19,7), -2.7557073670, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(22,7), 0.8250620859, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(25,7), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(28,7), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(31,7), -0.1377853683, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(34,7), 0.0412531043, TOLERANCE);
        //dN2/dv dw
        KRATOS_CHECK_NEAR( shape_functions(1,8), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(4,8), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(7,8), -0.8168114651, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(10,8), -0.4207816638, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(13,8), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(16,8), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(19,8), 0.5445409767, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(22,8), 0.2805211092, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(25,8), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(28,8), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(31,8), 0.2722704884, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(34,8), 0.1402605546, TOLERANCE);
        //dN2/dw2
        KRATOS_CHECK_NEAR( shape_functions(1,9), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(4,9), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(7,9), 1.2499690602, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(10,9), 0.0715328828, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(13,9), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(16,9), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(19,9), -3.3332508271, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(22,9), -0.1907543543, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(25,9), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(28,9), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(31,9), 2.0832817670, TOLERANCE);
        KRATOS_CHECK_NEAR( shape_functions(34,9), 0.1192214714, TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.