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

KRATOS_TEST_CASE_IN_SUITE(CalculateCompressibilityMatrix2D3NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    // clang-format off
    Vector Np(3);
    Np <<= 1.0, 2.0, 3.0;
    // clang-format on

    BoundedMatrix<double, 3, 3> CompressibilityMatrix  = ZeroMatrix(3, 3);
    const double                IntegrationCoefficient = 1.0;
    const double                BiotModulusInverse     = 0.02;
    CompressibilityMatrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix<3>(
        Np, BiotModulusInverse, IntegrationCoefficient);

    BoundedMatrix<double, 3, 3> CMatrix;
    // clang-format off
    CMatrix <<= -0.02,-0.04,-0.06,
                -0.04,-0.08,-0.12,
                -0.06,-0.12,-0.18;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(CompressibilityMatrix, CMatrix, 1e-12)

    const double DtPressureCoefficient = 1.5;
    CompressibilityMatrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix<3>(
        Np, BiotModulusInverse, IntegrationCoefficient, DtPressureCoefficient);
    CMatrix *= DtPressureCoefficient;

    KRATOS_CHECK_MATRIX_NEAR(CompressibilityMatrix, CMatrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateCompressibilityMatrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    // clang-format off
    Vector Np(4);
    Np <<= 1.0, 2.0, 3.0, 3.0;
    // clang-format on

    BoundedMatrix<double, 4, 4> CompressibilityMatrix   = ZeroMatrix(4, 4);
    const double                DynamicViscosityInverse = 1;
    BoundedMatrix<double, 3, 3> MaterialCompressibilityMatrix;

    const double IntegrationCoefficient = 1.0;
    const double BiotModulusInverse     = 0.1;
    CompressibilityMatrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix<4>(
        Np, BiotModulusInverse, IntegrationCoefficient);

    BoundedMatrix<double, 4, 4> CMatrix;
    // clang-format off
    CMatrix <<= -0.1,-0.2,-0.3,-0.3,
                -0.2,-0.4,-0.6,-0.6,
                -0.3,-0.6,-0.9,-0.9,
                -0.3,-0.6,-0.9,-0.9;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(CompressibilityMatrix, CMatrix, 1e-12)

    const double DtPressureCoefficient = 2.5;
    CompressibilityMatrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix<4>(
        Np, BiotModulusInverse, IntegrationCoefficient, DtPressureCoefficient);
    CMatrix *= DtPressureCoefficient;

    KRATOS_CHECK_MATRIX_NEAR(CompressibilityMatrix, CMatrix, 1e-12)
}

} // namespace Kratos::Testing