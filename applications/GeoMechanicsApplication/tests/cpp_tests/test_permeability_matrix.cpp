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

KRATOS_TEST_CASE_IN_SUITE(CalculatePermeabilityMatrix2D3NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    // clang-format off
    Matrix GradNpT(3, 2);
    GradNpT <<= 1.0, 2.0,
                3.0, 4.0,
                5.0, 6.0;
    // clang-format on

    BoundedMatrix<double, 3, 3> PermeabilityMatrix      = ZeroMatrix(3, 3);
    double                      DynamicViscosityInverse = 1;
    BoundedMatrix<double, 2, 2> MaterialPermeabilityMatrix;
    // clang-format off
    MaterialPermeabilityMatrix <<= 1.0, 2.0, 3.0, 4.0;
    // clang-format on
    double IntegrationCoefficient = 1.0;
    GeoTransportEquationUtilities::CalculatePermeabilityMatrixH<2, 3>(
        PermeabilityMatrix, GradNpT, DynamicViscosityInverse, MaterialPermeabilityMatrix, IntegrationCoefficient);

    // Kratos original calculation of the permeability matrix
    BoundedMatrix<double, 3, 2> PDimMatrix;
    BoundedMatrix<double, 3, 3> PMatrix;
    noalias(PDimMatrix) = -PORE_PRESSURE_SIGN_FACTOR * prod(GradNpT, MaterialPermeabilityMatrix);

    double RelativePermeability     = 0.02;
    double PermeabilityUpdateFactor = 1.5;
    noalias(PMatrix) = DynamicViscosityInverse * RelativePermeability * PermeabilityUpdateFactor *
                       prod(PDimMatrix, trans(GradNpT)) * IntegrationCoefficient;

    GeoTransportEquationUtilities::PreparePermeabilityMatrixHForIntegration<3>(
        PermeabilityMatrix, RelativePermeability, PermeabilityUpdateFactor);
    KRATOS_CHECK_MATRIX_NEAR(PermeabilityMatrix, PMatrix, 1e-12)

    // check of permeability matrix calculation using Matrix input parameters
    Matrix PermeabilityMatrixMatrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrixH(
        GradNpT, DynamicViscosityInverse, MaterialPermeabilityMatrix, IntegrationCoefficient);
    GeoTransportEquationUtilities::PreparePermeabilityMatrixHForIntegration(
        PermeabilityMatrixMatrix, RelativePermeability, PermeabilityUpdateFactor);
    KRATOS_CHECK_MATRIX_NEAR(PermeabilityMatrix, PermeabilityMatrixMatrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CalculatePermeabilityMatrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    // clang-format off
    Matrix GradNpT(4, 3);
    GradNpT <<= 1.0, 2.0, 3.0,
                3.0, 4.0, 5.0,
                5.0, 6.0, 7.0,
                7.0, 8.0, 9.0;
    // clang-format on

    BoundedMatrix<double, 4, 4> PermeabilityMatrix      = ZeroMatrix(4, 4);
    double                      DynamicViscosityInverse = 1;
    BoundedMatrix<double, 3, 3> MaterialPermeabilityMatrix;
    // clang-format off
    MaterialPermeabilityMatrix <<= 1.0, 2.0, 3.0,
                           4.0, 5.0, 6.0,
                           7.0, 8.0, 9.0;
    // clang-format on
    double IntegrationCoefficient = 1.0;
    GeoTransportEquationUtilities::CalculatePermeabilityMatrixH<3, 4>(
        PermeabilityMatrix, GradNpT, DynamicViscosityInverse, MaterialPermeabilityMatrix, IntegrationCoefficient);

    // Kratos original calculation of the permeability matrix
    BoundedMatrix<double, 4, 3> PDimMatrix;
    BoundedMatrix<double, 4, 4> PMatrix;
    noalias(PDimMatrix) = -PORE_PRESSURE_SIGN_FACTOR * prod(GradNpT, MaterialPermeabilityMatrix);

    double RelativePermeability     = 0.1;
    double PermeabilityUpdateFactor = 2.0;
    noalias(PMatrix) = DynamicViscosityInverse * RelativePermeability * PermeabilityUpdateFactor *
                       prod(PDimMatrix, trans(GradNpT)) * IntegrationCoefficient;

    GeoTransportEquationUtilities::PreparePermeabilityMatrixHForIntegration<4>(
        PermeabilityMatrix, RelativePermeability, PermeabilityUpdateFactor);
    KRATOS_CHECK_MATRIX_NEAR(PermeabilityMatrix, PMatrix, 1e-12)

    // check of permeability matrix calculation using Matrix input parameters
    Matrix PermeabilityMatrixMatrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrixH(
        GradNpT, DynamicViscosityInverse, MaterialPermeabilityMatrix, IntegrationCoefficient);
    GeoTransportEquationUtilities::PreparePermeabilityMatrixHForIntegration(
        PermeabilityMatrixMatrix, RelativePermeability, PermeabilityUpdateFactor);
    KRATOS_CHECK_MATRIX_NEAR(PermeabilityMatrix, PermeabilityMatrixMatrix, 1e-12)
}

} // namespace Kratos::Testing