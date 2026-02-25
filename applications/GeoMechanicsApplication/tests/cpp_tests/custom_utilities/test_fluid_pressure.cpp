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

#include "custom_utilities/transport_equation_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateFluidPressureGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto N               = UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0, 5.0});
    auto pressure_vector = UblasUtilities::CreateVector({0.5, 0.7, 0.8, 0.9, 0.4});

    auto fluid_pressure = GeoTransportEquationUtilities::CalculateFluidPressure(N, pressure_vector);
    double expected_value = 9.9;
    KRATOS_CHECK_NEAR(fluid_pressure, expected_value, 1e-12);
}

} // namespace Kratos::Testing
