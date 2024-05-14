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

#include "custom_retention/retention_law.h"
#include "includes/checks.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateFluidPressureGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    Vector N(5);
    N(0) = 1.0;
    N(1) = 2.0;
    N(2) = 3.0;
    N(3) = 4.0;
    N(4) = 5.0;

    Vector pressure_vector(5);
    pressure_vector(0) = 0.5;
    pressure_vector(1) = 0.7;
    pressure_vector(2) = 0.8;
    pressure_vector(3) = 0.9;
    pressure_vector(4) = 0.4;

    auto   fluid_pressure = RetentionLaw::Parameters::CalculateFluidPressure(N, pressure_vector);
    double expected_value = 9.9;
    KRATOS_CHECK_NEAR(fluid_pressure, expected_value, 1e-12);
}

} // namespace Kratos::Testing