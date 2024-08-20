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
#include "includes/checks.h"
#include "geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateSoilDensityGivesCorrectResults, KratosGeoMechanicsFastSuite)
{
    Properties properties(0);
    // Please note these are not representative values, it just ensures the values are set
    properties.SetValue(DENSITY_WATER, 1000.0);
    properties.SetValue(POROSITY, 0.3);
    properties.SetValue(DENSITY_SOLID, 2500.0);

    double degree_of_saturation = 0.0;
    auto soil_density = GeoTransportEquationUtilities::CalculateSoilDensity(degree_of_saturation, properties);
    double expected_value = 0.7 * 2500.0;
    KRATOS_CHECK_NEAR(soil_density, expected_value, 1e-12);

    degree_of_saturation = 0.1;
    soil_density = GeoTransportEquationUtilities::CalculateSoilDensity(degree_of_saturation, properties);
    expected_value = 0.1 * 0.3 * 1000.0 + 0.7 * 2500.0;
    KRATOS_CHECK_NEAR(soil_density, expected_value, 1e-12);
}

} // namespace Kratos::Testing