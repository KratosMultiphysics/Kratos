//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Apostolatos
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/pointer_vector.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_utilities.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

        ///// Tests
        KRATOS_TEST_CASE_IN_SUITE(nurbsUtilities, KratosCoreNurbsGeometriesFastSuite) {
            
            // Dummy polynomial order
            int polynomial_degree = 1;

            // Dummy number of knots
            int number_of_knots = 2;

            // Dummy number of control points
            int number_of_control_points = 2;
            
            // Dummy knot vector
            Vector knot_vector = ZeroVector(number_of_knots);
            knot_vector[0] = 0.0;
            knot_vector[1] = 1.0;

            // Dummy curve parameter
            double parameter_t = 0.5;

            // Call the NURBS utility functions
            int number_of_spans = NurbsUtilities::GetNumberOfSpans(polynomial_degree, number_of_knots);
            polynomial_degree = NurbsUtilities::GetPolynomialDegree(number_of_knots, number_of_control_points);
            number_of_control_points = NurbsUtilities::GetNumberOfControlPoints(polynomial_degree, number_of_knots);
            int lower_span = NurbsUtilities::GetLowerSpan(polynomial_degree, knot_vector, parameter_t);
            int upper_span = NurbsUtilities::GetUpperSpan(polynomial_degree, knot_vector, parameter_t);
            number_of_knots = NurbsUtilities::GetNumberOfKnots(polynomial_degree, number_of_control_points);

            // Verify the results
            KRATOS_CHECK_EQUAL(number_of_spans, 1);
            KRATOS_CHECK_EQUAL(polynomial_degree, 1);
            KRATOS_CHECK_EQUAL(number_of_control_points, 2);
            KRATOS_CHECK_EQUAL(lower_span, 0);
            KRATOS_CHECK_EQUAL(upper_span, 0);
            KRATOS_CHECK_EQUAL(number_of_knots, 2);

        }

} // namespace Testing.
} // namespace Kratos.
