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
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    /// Factory functions
    NurbsSurfaceShapeFunction GenerateReferenceHemisphereSurfaceShapeFunctions()
    {
        // Polynomial orders
        SizeType polynomial_degree_u = 2;
        SizeType polynomial_degree_v = 3;

        // Assign the derivative order of the basis
        SizeType derivative_order = 4;

        // Assign knots of the basis along u
        SizeType number_of_knots_u = 10;
        Vector knot_vector_u(number_of_knots_u);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 0.25;
        knot_vector_u[3] = 0.25;
        knot_vector_u[4] = 0.5;
        knot_vector_u[5] = 0.5;
        knot_vector_u[6] = 0.75;
        knot_vector_u[7] = 0.75;
        knot_vector_u[8] = 1.0;
        knot_vector_u[9] = 1.0;

        // Assign knots of the basis along v
        SizeType number_of_knots_v = 6;
        Vector knot_vector_v(number_of_knots_v);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 0.0;
        knot_vector_v[2] = 0.0;
        knot_vector_v[3] = 1.0;
        knot_vector_v[4] = 1.0;
        knot_vector_v[5] = 1.0;

        // Assign the weights
        SizeType number_of_control_points = 36;
        Vector weights(number_of_control_points);
        weights[0] = 1.0;
        weights[1] = 7.071067811865476e-01;
        weights[2] = 1.0;
        weights[3] = 7.071067811865476e-01;
        weights[4] = 1.0;
        weights[5] = 7.071067811865476e-01;
        weights[6] = 1.0;
        weights[7] = 7.071067811865476e-01;
        weights[8] = 1.0;

        weights[9] = 8.047378541243650e-01;
        weights[10] = 5.690355937288492e-01;
        weights[11] = 8.047378541243650e-01;
        weights[12] = 5.690355937288492e-01;
        weights[13] = 8.047378541243650e-01;
        weights[14] = 5.690355937288492e-01;
        weights[15] = 8.047378541243650e-01;
        weights[16] = 5.690355937288492e-01;
        weights[17] = 8.047378541243650e-01;

        weights[18] = 8.047378541243650e-01;
        weights[19] = 5.690355937288492e-01;
        weights[20] = 8.047378541243650e-01;
        weights[21] = 5.690355937288492e-01;
        weights[22] = 8.047378541243650e-01;
        weights[23] = 5.690355937288492e-01;
        weights[24] = 8.047378541243650e-01;
        weights[25] = 5.690355937288492e-01;
        weights[26] = 8.047378541243650e-01;

        weights[27] = 1.0;
        weights[28] = 7.071067811865476e-01;
        weights[29] = 1.0;
        weights[30] = 7.071067811865476e-01;
        weights[31] = 1.0;
        weights[32] = 7.071067811865476e-01;
        weights[33] = 1.0;
        weights[34] = 7.071067811865476e-01;
        weights[35] = 1.0;

        // Parameter on the surface where to compute the basis functions and their derivatives
        array_1d<double, 3> parameter(0.0);
        parameter[0] = .582;
        parameter[1] = .8972;

        // Create the shape function container
        NurbsSurfaceShapeFunction shape_functions(polynomial_degree_u, polynomial_degree_v, derivative_order);

        // Get the knot span where parameter parameterT lies in along u
        const int knot_span_u = NurbsUtilities::GetLowerSpan(shape_functions.PolynomialDegreeU(), knot_vector_u, parameter[0]);

        // Get the knot span where parameter parameterT lies in along v
        const int knot_span_v = NurbsUtilities::GetLowerSpan(shape_functions.PolynomialDegreeV(), knot_vector_v, parameter[1]);

        // Computer the shape functions at the given parametric location
        shape_functions.ComputeNurbsShapeFunctionValuesAtSpan(knot_vector_u, knot_vector_v, knot_span_u, knot_span_v, weights, parameter[0], parameter[1]);

        // Return the object
        return shape_functions;
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsBasisFunctions2dHemisphere, KratosCoreNurbsGeometriesFastSuite) {
        NurbsSurfaceShapeFunction shape_functions = GenerateReferenceHemisphereSurfaceShapeFunctions();

        KRATOS_CHECK_NEAR(shape_functions(0,0), 0.000595497314518, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,0), 0.000411054470482, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,0), 0.000141869470763, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,0), 0.012547339321201, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,0), 0.008661063274159, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,0), 0.002989239994181, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,0), 0.109508490651574, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,0), 0.075590524995869, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,0), 0.026088970065943, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,0), 0.395884281079147, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,0), 0.273267401151554, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,0), 0.094314268210607, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,1), -6.538093869497354e-03, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,1), 2.946556221434208e-03, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,1), 3.591537648063146e-03, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,1), -1.377599533942699e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,1), 6.208498315267830e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,1), 7.567497024159162e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,1), -1.202317414254270e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,1), 5.418545416788226e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,1), 6.604628725754473e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,1), -4.346499183204187e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,1), 1.958859029155177e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,1), 2.387640154049010e+00, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,2), -1.767126734886447e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,2), -1.219796171324385e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,2), -4.209949038186101e-03, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,2), -2.362990345977836e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,2), -1.631103485673575e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,2), -5.629516399643902e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,2), -8.750161854898419e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,2), -6.039982146362372e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,2), -2.084611972517605e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,2), 1.128986487436490e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,2), 7.793065249168387e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,2), 2.689663102863856e-01, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,3), 1.727766039820421e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,3), -6.307021943770767e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,3), 4.579255903950346e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,3), 3.640464237326134e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,3), -1.328911860816241e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,3), 9.648654370836276e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,3), 3.177261200125493e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,3), -1.159824631833007e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,3), 8.420985118204577e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,3), 1.148612092567699e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,3), -4.192883472497170e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,3), 3.044271379929471e+01, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,4), 1.940166679228831e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,4), -8.743848457074178e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,4), -1.065781833521413e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,4), 2.594378231112102e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,4), -1.169221713579426e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,4), -1.425156517532675e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,4), 9.606992036043165e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,4), -4.329632262567705e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,4), -5.277359773475458e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,4), -1.239538693507815e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,4), 5.586292460717871e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,4), 6.809094474360275e+00, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,5), 3.547481559201592e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,5), 2.448723307916013e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,5), 8.451412274685188e-02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,5), 2.047407652402728e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,5), 1.413265934037026e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,5), 4.877681779604565e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,5), -3.751909478983691e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,5), -2.589833953202036e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,5), -8.938435139132052e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,5), 1.349753670660804e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,5), 9.316956883734089e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,5), 3.215612132058969e-01, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,6), 4.701577506505209e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,6), -3.653908869085659e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,6), -1.047668637419549e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,6), 9.906390319622304e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,6), -7.698915395823892e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,6), -2.207474923798412e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,6), 8.645927426814330e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,6), -6.719325771530347e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,6), -1.926601655283984e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,6), 3.125590301958625e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,6), -2.429104297366770e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,6), -6.964860045918553e+01, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,7), -5.127112223949221e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,7), 1.871596533288967e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,7), -1.358885310894045e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,7), -6.855941030576663e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,7), 2.502686679125819e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,7), -1.817092576068153e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,7), -2.538757459898127e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,7), 9.267457885185344e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,7), -6.728700425287219e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,7), 3.275622685195290e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,7), -1.195730421764006e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,7), 8.681681532444775e+01, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,8), -3.894856764070146e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,8), 1.755315028942571e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,8), 2.139541735127575e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,8), -2.247893163273971e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,8), 1.013069514995133e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,8), 1.234823648278837e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,8), 4.119302600599678e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,8), -1.856467182643897e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,8), -2.262835417955780e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,8), -1.481923760918696e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,8), 6.678661647545073e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,8), 8.140575961641863e+00, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,9), -3.746781354842437e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,9), -2.586294158307824e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,9), -8.926218051998480e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,9), 5.796561712776662e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,9), 4.001197901940428e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,9), 1.380955249316547e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,9), 8.479541727658783e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,9), 5.853180946101294e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,9), 2.020140255696483e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,9), -2.897734530700103e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,9), -2.000221838242734e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,9), -6.903474696863472e-01, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,10), -4.907291725401470e-01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,10), 6.792515846795983e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,10), -6.301786674255835e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,10), -1.033983746451421e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,10), 1.431207145633095e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,10), -1.327808770987953e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,10), -9.024223903854219e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,10), 1.249104135274331e+03, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,10), -1.158861896235789e+03, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,10), -3.262348308535680e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,10), 4.515637916692989e+03, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,10), -4.189403085839420e+03, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,11), -1.395184009285945e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,11), 1.084290372429864e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,11), 3.108936368560812e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,11), -1.865630958844154e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,11), 1.449906015061858e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,11), 4.157249437822951e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,11), -6.908438233436947e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,11), 5.369007252939933e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,11), 1.539430980497010e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,11), 8.913587593209688e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,11), -6.927342305244773e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,11), -1.986245287964915e+02, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,12), 1.029260420735227e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,12), -3.757203180186667e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,12), 2.727942759451440e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,12), 5.940314633243176e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,12), -2.168447224987788e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,12), 1.574415761663471e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,12), -1.088572798604870e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,12), 3.973716562287830e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,12), -2.885143763682960e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,12), 3.916152932070303e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,12), -1.429549019281375e+02, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,12), 1.037933726074346e+02, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,13), 4.113672322142911e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,13), -1.853929756240066e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,13), -2.259742565902845e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,13), -6.364170530160391e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,13), 2.868173300080503e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,13), 3.495997230079885e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,13), -9.309872342682372e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,13), 4.195727809915702e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,13), 5.114144532766768e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,13), 3.181485442285721e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,13), -1.433816324832006e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,13), -1.747669117453719e+01, TOLERANCE);

        KRATOS_CHECK_NEAR(shape_functions(0,14), 4.736463414407606e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(1,14), 3.269442889665566e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(2,14), 1.128400651873466e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(3,14), -2.662005196153613e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(4,14), -1.837504737045614e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(5,14), -6.341880292990669e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(6,14), 2.621160430233409e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(7,14), 1.809310783491226e+01, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(8,14), 6.244572963750500e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(9,14), -4.328015755205562e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(10,14), -2.987503354121680e+00, TOLERANCE);
        KRATOS_CHECK_NEAR(shape_functions(11,14), -1.031093322633297e+00, TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.