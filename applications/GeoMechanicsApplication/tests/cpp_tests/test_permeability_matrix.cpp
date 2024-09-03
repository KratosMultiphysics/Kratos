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
#include "geo_mechanics_fast_suite.h"
#include "includes/checks.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculatePermeabilityMatrix2D3NGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // clang-format off
    Matrix GradNpT(3, 2);
    GradNpT <<= 1.0, 2.0,
                3.0, 4.0,
                5.0, 6.0;
    // clang-format on

    BoundedMatrix<double, 3, 3> PermeabilityMatrix      = ZeroMatrix(3, 3);
    const double                DynamicViscosityInverse = 1;
    BoundedMatrix<double, 2, 2> MaterialPermeabilityMatrix;
    // clang-format off
    MaterialPermeabilityMatrix <<= 1.0, 2.0,
                                   3.0, 4.0;
    // clang-format on
    const double IntegrationCoefficient   = 1.0;
    const double RelativePermeability     = 0.02;
    const double PermeabilityUpdateFactor = 1.5;
    PermeabilityMatrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrix<2, 3>(
        GradNpT, DynamicViscosityInverse, MaterialPermeabilityMatrix,
        RelativePermeability * PermeabilityUpdateFactor, IntegrationCoefficient);

    BoundedMatrix<double, 3, 3> PMatrix;
    // clang-format off
    PMatrix <<= -0.81, -1.83, -2.85,
                -1.77, -3.99, -6.21,
                -2.73, -6.15, -9.57;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(PermeabilityMatrix, PMatrix, 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CalculatePermeabilityMatrix3D4NGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // clang-format off
    Matrix GradNpT(4, 3);
    GradNpT <<= 1.0, 2.0, 3.0,
                3.0, 4.0, 5.0,
                5.0, 6.0, 7.0,
                7.0, 8.0, 9.0;
    // clang-format on

    BoundedMatrix<double, 4, 4> PermeabilityMatrix      = ZeroMatrix(4, 4);
    const double                DynamicViscosityInverse = 1;
    BoundedMatrix<double, 3, 3> MaterialPermeabilityMatrix;
    // clang-format off
    MaterialPermeabilityMatrix <<= 1.0, 2.0, 3.0,
                                   4.0, 5.0, 6.0,
                                   7.0, 8.0, 9.0;
    // clang-format on
    const double IntegrationCoefficient   = 1.0;
    const double RelativePermeability     = 0.1;
    const double PermeabilityUpdateFactor = 2.0;
    PermeabilityMatrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrix<3, 4>(
        GradNpT, DynamicViscosityInverse, MaterialPermeabilityMatrix,
        RelativePermeability * PermeabilityUpdateFactor, IntegrationCoefficient);

    BoundedMatrix<double, 4, 4> PMatrix;
    // clang-format off
    PMatrix <<= -45.6,  -88.8, -132.0, -175.2,
                -84.0, -163.2, -242.4, -321.6,
               -122.4, -237.6, -352.8, -468.0,
               -160.8, -312.0, -463.2, -614.4;
    // clang-format on

    KRATOS_CHECK_MATRIX_NEAR(PermeabilityMatrix, PMatrix, 1e-12)
}

} // namespace Kratos::Testing