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
#include <limits>

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "testing/testing.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    /// Factory functions
    NurbsVolumeShapeFunction GenerateVolumeShapeFunctions()
    {
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

        // Parameter on the surface where to compute the basis functions and their derivatives
        array_1d<double, 3> parameter(0.0);
        parameter[0] = 0.5;
        parameter[1] = 0.5;
        parameter[2] = 0.5;

        SizeType number_of_control_points = 36;
        Vector weights(number_of_control_points);

        // Create the shape function container
        NurbsVolumeShapeFunction shape_functions(polynomial_degree_u, polynomial_degree_v, polynomial_degree_w, derivative_order);

        // Get the knot span where parameter parameterT lies in along u
        const int knot_span_u = NurbsUtilities::GetLowerSpan(shape_functions.PolynomialDegreeU(), knot_vector_u, parameter[0]);

        // Get the knot span where parameter parameterT lies in along v
        const int knot_span_v = NurbsUtilities::GetLowerSpan(shape_functions.PolynomialDegreeV(), knot_vector_v, parameter[1]);

        // Get the knot span where parameter parameterT lies in along v
        const int knot_span_w = NurbsUtilities::GetLowerSpan(shape_functions.PolynomialDegreeW(), knot_vector_w, parameter[2]);

        // Computer the shape functions at the given parametric location
        shape_functions.ComputeNurbsShapeFunctionValuesAtSpan(knot_vector_u, knot_vector_v, knot_vector_w,
            knot_span_u, knot_span_v, knot_span_w, weights, parameter[0], parameter[1], parameter[2]);

        // Return the object
        return shape_functions;
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsBasisFunctions3dCube, KratosCoreNurbsGeometriesFastSuite) {
        NurbsVolumeShapeFunction shape_functions = GenerateVolumeShapeFunctions();
        // TODO: Add test here...
    }

} // namespace Testing.
} // namespace Kratos.