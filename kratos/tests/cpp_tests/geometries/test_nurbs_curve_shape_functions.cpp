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
//                   Tobias Teschemacher
//

// System includes
#include <limits>

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "testing/testing.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Factory functions
    NurbsCurveShapeFunction GenerateReferenceNurbsCurveShapeFunction()
    {
        // Assign the polynomial order of the basis
        SizeType polynomial_degree = 4;

        // Assign the derivative order of the basis
        SizeType derivative_order = 6;

        // Assign knots of the basis
        SizeType number_of_knots = 11;
        Vector knots(number_of_knots);
        knots[0] = -1.000000000000000;
        knots[1] = -1.000000000000000;
        knots[2] = -1.000000000000000;
        knots[3] = -1.000000000000000;
        knots[4] = 4.500000000000000;
        knots[5] = 4.500000000000000;
        knots[6] = 4.500000000000000;
        knots[7] = 9.000000000000000;
        knots[8] = 9.000000000000000;
        knots[9] = 9.000000000000000;
        knots[10] = 9.000000000000000;

        // Assign the weights
        SizeType number_of_control_points = 8;
        Vector weights(number_of_control_points);
        weights[0] = 1.000000000000000;
        weights[1] = 2.450000000000000;
        weights[2] = 3.618333333333333;
        weights[3] = 4.505000000000000;
        weights[4] = 5.605000000000000;
        weights[5] = 5.085000000000000;
        weights[6] = 3.550000000000000;
        weights[7] = 1.000000000000000;

        // Parameter on the curve where to compute the basis functions and their derivatives
        double parameter_t = 6.5;

        // Create the shape function container
        NurbsCurveShapeFunction shape_functions(polynomial_degree, derivative_order);

        // Get the knot span where parameter parameterT lies in
        const int knot_span = NurbsUtilities::GetLowerSpan(shape_functions.PolynomialDegree(), knots, parameter_t);

        // Computer the shape functions at the given parametric location
        shape_functions.ComputeNurbsShapeFunctionValuesAtSpan(knots, knot_span, weights, parameter_t);

        // Return the object
        return shape_functions;
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsBasisFunctions1d, KratosCoreNurbsGeometriesFastSuite) {
        
        // Get reference Nurbs curve shape functions
        NurbsCurveShapeFunction shape_functions = GenerateReferenceNurbsCurveShapeFunction();

        // Check the Nurbs basis functions
        KRATOS_CHECK_NEAR(shape_functions(0, 0), 0.040341356830719, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1, 0), 0.418263605320343, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2, 0), 0.388566516226091, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3, 0), 0.144677667136296, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4, 0), 0.008150854486552, TOLERANCE);

        // Check the first derivatives of the Nurbs basis functions
        KRATOS_CHECK_NEAR(shape_functions(0, 1), -0.058116523341430, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1, 1), -0.281332029339492, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2, 1), 0.139643440454173, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3, 1), 0.182204311552113, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4, 1), 0.017600800674635, TOLERANCE);

        // Check the second derivatives of the Nurbs basis functions
        KRATOS_CHECK_NEAR(shape_functions(0, 2), 0.063998870285300, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1, 2), -0.005001367942227, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2, 2), -0.209746254260036, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3, 2), 0.119661582154355, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4, 2), 0.031087169762608, TOLERANCE);

        // Check the third derivatives of the Nurbs basis functions
        KRATOS_CHECK_NEAR(shape_functions(0, 3), -0.053270373094405, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1, 3), 0.174241510932662, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2, 3), -0.140906951084536, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3, 3), -0.026015381549528, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4, 3), 0.045951194795806, TOLERANCE);

        // Check the fourth derivatives of the Nurbs basis functions
        KRATOS_CHECK_NEAR(shape_functions(0, 4), 0.039072547258617, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1, 4), -0.129109960704594, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2, 4), 0.125067122234968, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3, 4), -0.099987226672343, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4, 4), 0.064957517883351, TOLERANCE);

        // Check the fifth derivatives of the Nurbs basis functions
        KRATOS_CHECK_NEAR(shape_functions(0, 5), -0.035796024935348, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1, 5), 0.116042517778655, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2, 5), -0.077380079426900, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3, 5), -0.112368141670365, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4, 5), 0.109501728253958, TOLERANCE);

        // Check the sixth derivatives of the Nurbs basis functions
        KRATOS_CHECK_NEAR(shape_functions(0, 6), 0.039409546074973, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1, 6), -0.132366199073095, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2, 6), 0.161718941270974, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3, 6), -0.295903990157367, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4, 6), 0.227141701884515, TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.