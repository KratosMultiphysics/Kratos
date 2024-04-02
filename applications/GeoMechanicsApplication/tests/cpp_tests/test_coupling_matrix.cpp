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
#include "testing/testing.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateCouplingMatrix2D3NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    constexpr SizeType voigt_size = VOIGT_SIZE_2D_PLANE_STRESS;
    Matrix             b_matrix   = ZeroMatrix(voigt_size, 2 * 3);
    // clang-format off
    b_matrix <<=  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,
                 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
                 21.0, 22.0, 23.0, 24.0, 25.0, 26.0;
    // clang-format on
    Vector voigt_vector = ZeroVector(voigt_size);
    voigt_vector <<= 1.0, 1.0, 0.0;

    Vector n_p = ZeroVector(voigt_size);
    n_p <<= 1.0, 2.0, 3.0;

    const double                integration_coefficient = 1.0;
    const double                biot_coefficient        = 0.02;
    BoundedMatrix<double, 6, 3> coupling_matrix         = ZeroMatrix(6, 3);

    coupling_matrix = GeoTransportEquationUtilities::CalculateCouplingMatrix<2, 3>(
        b_matrix, voigt_vector, n_p, biot_coefficient, integration_coefficient);

    BoundedMatrix<double, 6, 3> expected_coupling_matrix;
    // clang-format off
    expected_coupling_matrix <<= 0.24,0.48,0.72,
                                 0.28,0.56,0.84,
                                 0.32,0.64,0.96,
                                 0.36,0.72,1.08,
                                 0.40,0.80,1.20,
                                 0.44,0.88,1.32;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(coupling_matrix, expected_coupling_matrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateCouplingMatrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    constexpr SizeType voigt_size = VOIGT_SIZE_3D;
    Matrix             b_matrix   = ZeroMatrix(voigt_size, 3 * 4);
    // clang-format off
    b_matrix <<=  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 11.0, 12.0,
                 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0,
                 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0,
                 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0,
                 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0,
                 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0;

    // clang-format on
    Vector voigt_vector = ZeroVector(voigt_size);
    voigt_vector <<= 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

    Vector n_p = ZeroVector(voigt_size);
    n_p <<= 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

    const double                 integration_coefficient = 1.0;
    const double                 biot_coefficient        = 0.02;
    BoundedMatrix<double, 12, 6> coupling_matrix         = ZeroMatrix(12, 6);

    coupling_matrix = GeoTransportEquationUtilities::CalculateCouplingMatrix<3, 4>(
        b_matrix, voigt_vector, n_p, biot_coefficient, integration_coefficient);

    BoundedMatrix<double, 12, 6> expected_coupling_matrix = ZeroMatrix(12, 6);
    // clang-format off
    expected_coupling_matrix <<= 0.66,1.32,1.98,2.64,3.3,3.96,
                                 0.72,1.44,2.16,2.88,3.6,4.32,
                                 0.78,1.56,2.34,3.12,3.9,4.68,
                                 0.84,1.68,2.52,3.36,4.2,5.04,
                                 0.90,1.80,2.70,3.60,4.5,5.40,
                                 0.96,1.92,2.88,3.84,4.8,5.76,
                                 1.02,2.04,3.06,4.08,5.1,6.12,
                                 1.08,2.16,3.24,4.32,5.4,6.48,
                                 1.14,2.28,3.42,4.56,5.7,6.84,
                                 1.20,2.40,3.60,4.80,6.0,7.20,
                                 1.26,2.52,3.78,5.04,6.3,7.56,
                                 1.32,2.64,3.96,5.28,6.6,7.92;
    std::cout << expected_coupling_matrix << std::endl;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(coupling_matrix, expected_coupling_matrix, 1e-12)
}

} // namespace Kratos::Testing