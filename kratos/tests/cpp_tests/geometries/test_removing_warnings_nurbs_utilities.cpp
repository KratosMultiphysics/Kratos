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
        KRATOS_TEST_CASE_IN_SUITE(removingWarningsNurbsUtilities, KratosCoreNurbsGeometriesFastSuite) {
            
            // Dummy polynomial order
            const int polynomial_degree = 1;

            // Dummy number of knots
            const int number_of_knots = 2;

            // Dummy number of control points
            const int number_of_control_points = 2;
            
            // Dummy knot vector
            Vector knot_vector = ZeroVector(number_of_knots);
            knot_vector[0] = 0.0;
            knot_vector[1] = 1.0;

            // Dummy curve parameter
            double parameter_t = .5;

            // Call the functions from the geometries/nurbs_shape_function_utilities/nurbs_utilities.h to remove the warning
            NurbsUtilities::GetNumberOfSpans(polynomial_degree, number_of_knots);
            NurbsUtilities::GetNumberOfControlPoints(polynomial_degree, number_of_knots);
            NurbsUtilities::GetPolynomialDegree(number_of_knots, number_of_control_points);
            NurbsUtilities::GetLowerSpan(polynomial_degree, knot_vector, parameter_t);
            NurbsUtilities::GetUpperSpan(polynomial_degree, knot_vector, parameter_t);
            NurbsUtilities::GetNumberOfKnots(polynomial_degree, number_of_control_points);
            
        }

} // namespace Testing.
} // namespace Kratos.
