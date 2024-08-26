// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "containers/model.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "includes/checks.h"
#include "geo_mechanics_fast_suite.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(Calculatecompressibility_matrix2D3NGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector n_p(3);
    n_p <<= 1.0, 2.0, 3.0;

    BoundedMatrix<double, 3, 3> compressibility_matrix  = ZeroMatrix(3, 3);
    const double                integration_coefficient = 1.0;
    const double                biot_modulus_inverse    = 0.02;
    compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix<3>(
        n_p, biot_modulus_inverse, integration_coefficient);

    BoundedMatrix<double, 3, 3> expected_compressibility_matrix;
    // clang-format off
    expected_compressibility_matrix <<= -0.02,-0.04,-0.06,
                                        -0.04,-0.08,-0.12,
                                        -0.06,-0.12,-0.18;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(compressibility_matrix, expected_compressibility_matrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(Calculatecompressibility_matrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector n_p(4);
    n_p <<= 1.0, 2.0, 3.0, 3.0;

    BoundedMatrix<double, 4, 4> compressibility_matrix = ZeroMatrix(4, 4);
    BoundedMatrix<double, 3, 3> material_compressibility_matrix;

    const double integration_coefficient = 1.0;
    const double biot_modulus_inverse    = 0.1;
    compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix<4>(
        n_p, biot_modulus_inverse, integration_coefficient);

    BoundedMatrix<double, 4, 4> expected_compressibility_matrix;
    // clang-format off
    expected_compressibility_matrix <<= -0.1,-0.2,-0.3,-0.3,
                                        -0.2,-0.4,-0.6,-0.6,
                                        -0.3,-0.6,-0.9,-0.9,
                                        -0.3,-0.6,-0.9,-0.9;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(compressibility_matrix, expected_compressibility_matrix, 1e-12)
}

} // namespace Kratos::Testing