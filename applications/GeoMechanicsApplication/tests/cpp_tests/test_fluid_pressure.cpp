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

#include "boost/numeric/ublas/assignment.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_fast_suite.h"
#include "includes/checks.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateFluidPressureGivesCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Vector N(5);
    N <<= 1.0, 2.0, 3.0, 4.0, 5.0;

    Vector pressure_vector(5);
    pressure_vector <<= 0.5, 0.7, 0.8, 0.9, 0.4;

    auto fluid_pressure = GeoTransportEquationUtilities::CalculateFluidPressure(N, pressure_vector);
    double expected_value = 9.9;
    KRATOS_CHECK_NEAR(fluid_pressure, expected_value, 1e-12);
}

} // namespace Kratos::Testing