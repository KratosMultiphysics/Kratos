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
#include "custom_utilities/ublas_utilities.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(Calculatecompressibility_matrix2D3NGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto n_p = UblasUtilities::CreateVector({1.0, 2.0, 3.0});

    BoundedMatrix<double, 3, 3> compressibility_matrix  = ZeroMatrix(3, 3);
    constexpr auto              integration_coefficient = 1.0;
    constexpr auto              biot_modulus_inverse    = 0.02;
    compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix<3>(
        n_p, biot_modulus_inverse, integration_coefficient);

    auto expected_compressibility_matrix =
        UblasUtilities::CreateMatrix({{-0.02, -0.04, -0.06}, {-0.04, -0.08, -0.12}, {-0.06, -0.12, -0.18}});

    KRATOS_CHECK_MATRIX_NEAR(compressibility_matrix, expected_compressibility_matrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(Calculatecompressibility_matrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto n_p = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 3.0});

    BoundedMatrix<double, 4, 4> compressibility_matrix = ZeroMatrix(4, 4);
    BoundedMatrix<double, 3, 3> material_compressibility_matrix;

    constexpr auto integration_coefficient = 1.0;
    constexpr auto biot_modulus_inverse    = 0.1;
    compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix<4>(
        n_p, biot_modulus_inverse, integration_coefficient);

    auto expected_compressibility_matrix = UblasUtilities::CreateMatrix(
        {{-0.1, -0.2, -0.3, -0.3}, {-0.2, -0.4, -0.6, -0.6}, {-0.3, -0.6, -0.9, -0.9}, {-0.3, -0.6, -0.9, -0.9}});

    KRATOS_CHECK_MATRIX_NEAR(compressibility_matrix, expected_compressibility_matrix, 1e-12)
}

} // namespace Kratos::Testing
