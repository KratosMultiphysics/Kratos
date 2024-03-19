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
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculatePermeabilityMatrix2D3NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    Vector Np(3);
    Np <<= 1.0, 2.0, 3.0;

    // clang-format off
    Matrix GradNpT(3, 2);
    GradNpT <<= 1.0, 2.0,
                3.0, 4.0,
                5.0, 6.0;
    // clang-format on

    BoundedMatrix<double, 3, 2> PDimMatrix              = ZeroMatrix(3, 2);
    BoundedMatrix<double, 3, 3> PMatrix                 = ZeroMatrix(3, 3);
    double                      DynamicViscosityInverse = 1;
    BoundedMatrix<double, 2, 2> PermeabilityMatrix;
    // clang-format off
    PermeabilityMatrix <<= 1.0, 2.0, 3.0, 4.0;
    // clang-format on
    double IntegrationCoefficient = 1.0;
    GeoTransportEquationUtilities::CalculatePermeabilityMatrixH<2, 3>(
        PDimMatrix, PMatrix, GradNpT, DynamicViscosityInverse, PermeabilityMatrix, IntegrationCoefficient);

    // the old original content
    BoundedMatrix<double, 3, 2> PDimMatrixOriginal;
    BoundedMatrix<double, 3, 3> PMatrixOriginal;
    noalias(PDimMatrixOriginal) = -PORE_PRESSURE_SIGN_FACTOR * prod(GradNpT, PermeabilityMatrix);

    double RelativePermeability     = 0.02;
    double PermeabilityUpdateFactor = 1.5;
    noalias(PMatrixOriginal) = DynamicViscosityInverse * RelativePermeability * PermeabilityUpdateFactor *
                               prod(PDimMatrixOriginal, trans(GradNpT)) * IntegrationCoefficient;

    PDimMatrix *= -1.0;
    KRATOS_CHECK_MATRIX_NEAR(PDimMatrix, PDimMatrixOriginal, 1e-12)

    GeoTransportEquationUtilities::PreparePermeabilityMatrixHForIntegration<3>(
        PMatrix, RelativePermeability, PermeabilityUpdateFactor);
    KRATOS_CHECK_MATRIX_NEAR(PMatrix, PMatrixOriginal, 1e-12)

    Matrix PMatrix1 = GeoTransportEquationUtilities::CalculatePermeabilityMatrixH(
        GradNpT, DynamicViscosityInverse, PermeabilityMatrix, IntegrationCoefficient);
    GeoTransportEquationUtilities::PreparePermeabilityMatrixHForIntegration(
        PMatrix1, RelativePermeability, PermeabilityUpdateFactor);
    KRATOS_CHECK_MATRIX_NEAR(PMatrix, PMatrix1, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CalculatePermeabilityMatrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    Vector Np(4);
    Np <<= 1.0, 2.0, 3.0, 4.0;

    // clang-format off
    Matrix GradNpT(4, 3);
    GradNpT <<= 1.0, 2.0, 3.0,
                3.0, 4.0, 5.0,
                5.0, 6.0, 7.0,
                7.0, 8.0, 9.0;
    // clang-format on

    BoundedMatrix<double, 4, 3> PDimMatrix              = ZeroMatrix(4, 3);
    BoundedMatrix<double, 4, 4> PMatrix                 = ZeroMatrix(4, 4);
    double                      DynamicViscosityInverse = 1;
    BoundedMatrix<double, 3, 3> PermeabilityMatrix;
    // clang-format off
    PermeabilityMatrix <<= 1.0, 2.0, 3.0,
                           4.0, 5.0, 6.0,
                           7.0, 8.0, 9.0;
    // clang-format on
    double IntegrationCoefficient = 1.0;
    GeoTransportEquationUtilities::CalculatePermeabilityMatrixH<3, 4>(
        PDimMatrix, PMatrix, GradNpT, DynamicViscosityInverse, PermeabilityMatrix, IntegrationCoefficient);

    // the old original content
    BoundedMatrix<double, 4, 3> PDimMatrixOriginal;
    BoundedMatrix<double, 4, 4> PMatrixOriginal;
    noalias(PDimMatrixOriginal) = -PORE_PRESSURE_SIGN_FACTOR * prod(GradNpT, PermeabilityMatrix);

    double RelativePermeability     = 0.1;
    double PermeabilityUpdateFactor = 2.0;
    noalias(PMatrixOriginal) = DynamicViscosityInverse * RelativePermeability * PermeabilityUpdateFactor *
                               prod(PDimMatrixOriginal, trans(GradNpT)) * IntegrationCoefficient;

    PDimMatrix *= -1.0;
    KRATOS_CHECK_MATRIX_NEAR(PDimMatrix, PDimMatrixOriginal, 1e-12)

    GeoTransportEquationUtilities::PreparePermeabilityMatrixHForIntegration<4>(
        PMatrix, RelativePermeability, PermeabilityUpdateFactor);
    KRATOS_CHECK_MATRIX_NEAR(PMatrix, PMatrixOriginal, 1e-12)

    Matrix PMatrix1 = GeoTransportEquationUtilities::CalculatePermeabilityMatrixH(
        GradNpT, DynamicViscosityInverse, PermeabilityMatrix, IntegrationCoefficient);
    GeoTransportEquationUtilities::PreparePermeabilityMatrixHForIntegration(
        PMatrix1, RelativePermeability, PermeabilityUpdateFactor);
    KRATOS_CHECK_MATRIX_NEAR(PMatrix, PMatrix1, 1e-12)
}

} // namespace Kratos::Testing