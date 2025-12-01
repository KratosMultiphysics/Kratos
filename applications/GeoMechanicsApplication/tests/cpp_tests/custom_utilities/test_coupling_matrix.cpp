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
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateCouplingMatrix2D3NGivesCorrectResults)
{
    constexpr SizeType voigt_size = VOIGT_SIZE_2D_PLANE_STRESS;
    constexpr SizeType size_D     = 2;
    constexpr SizeType size_N     = 3;

    Matrix b_matrix = ZeroMatrix(voigt_size, size_D * size_N);
    // clang-format off
    b_matrix <<=  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,
                 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
                 21.0, 22.0, 23.0, 24.0, 25.0, 26.0;
    // clang-format on
    Vector voigt_vector = ZeroVector(voigt_size);
    voigt_vector <<= 1.0, 1.0, 0.0;

    Vector n_p = ZeroVector(size_N);
    n_p <<= 1.0, 2.0, 3.0;

    const double integration_coefficient = 1.0;
    const double biot_coefficient        = 0.02;
    const double bishop_coefficient      = 0.1;

    const auto coupling_matrix = GeoTransportEquationUtilities::CalculateCouplingMatrix<size_D, size_N>(
        b_matrix, voigt_vector, n_p, biot_coefficient, bishop_coefficient, integration_coefficient);

    BoundedMatrix<double, size_D * size_N, size_N> expected_coupling_matrix;
    // clang-format off
    expected_coupling_matrix <<= 0.024,0.048,0.072,
                                 0.028,0.056,0.084,
                                 0.032,0.064,0.096,
                                 0.036,0.072,0.108,
                                 0.040,0.080,0.120,
                                 0.044,0.088,0.132;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(coupling_matrix, expected_coupling_matrix, 1e-12)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateCouplingMatrix3D4NGivesCorrectResults)
{
    constexpr SizeType voigt_size = VOIGT_SIZE_3D;
    constexpr SizeType size_D     = 3;
    constexpr SizeType size_N     = 4;

    Matrix b_matrix = ZeroMatrix(voigt_size, size_D * size_N);
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

    Vector n_p = ZeroVector(size_N);
    n_p <<= 1.0, 2.0, 3.0, 4.0;

    const double integration_coefficient = 1.0;
    const double biot_coefficient        = 0.02;
    const double bishop_coefficient      = 0.1;

    BoundedMatrix<double, size_D * size_N, size_N> coupling_matrix;
    GeoTransportEquationUtilities::CalculateCouplingMatrix(coupling_matrix, b_matrix, voigt_vector,
                                                           n_p, biot_coefficient, bishop_coefficient,
                                                           integration_coefficient);

    BoundedMatrix<double, size_D * size_N, size_N> expected_coupling_matrix;
    // clang-format off
    expected_coupling_matrix <<= 0.066,0.132,0.198,0.264,
                                 0.072,0.144,0.216,0.288,
                                 0.078,0.156,0.234,0.312,
                                 0.084,0.168,0.252,0.336,
                                 0.090,0.180,0.270,0.360,
                                 0.096,0.192,0.288,0.384,
                                 0.102,0.204,0.306,0.408,
                                 0.108,0.216,0.324,0.432,
                                 0.114,0.228,0.342,0.456,
                                 0.120,0.240,0.360,0.480,
                                 0.126,0.252,0.378,0.504,
                                 0.132,0.264,0.396,0.528;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(coupling_matrix, expected_coupling_matrix, 1e-12)
}

} // namespace Kratos::Testing
