// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "testing/testing.h"
#include "geo_mechanics_application.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ThermalAnalysisVariablesExistAfterRegistration, KratosGeoMechanicsFastSuite)
{
    KratosGeoMechanicsApplication geo_app;

    // Note that when this unit test is run through Python, the variables have already been registered. Therefore,
    // checking whether the thermal variables don't exist prior to registration will fail. Presumably after migrating to
    // GoogleTest, we can add those checks here as well, to test the effect of `Register()`.
    geo_app.Register();

    const auto variable_names = std::vector<std::string>{"SPECIFIC_HEAT_CAPACITY_WATER",
                                                         "SPECIFIC_HEAT_CAPACITY_SOLID",
                                                         "THERMAL_CONDUCTIVITY_WATER",
                                                         "THERMAL_CONDUCTIVITY_SOLID_XX",
                                                         "THERMAL_CONDUCTIVITY_SOLID_YY",
                                                         "THERMAL_CONDUCTIVITY_SOLID_ZZ",
                                                         "THERMAL_CONDUCTIVITY_SOLID_XY",
                                                         "THERMAL_CONDUCTIVITY_SOLID_YZ",
                                                         "THERMAL_CONDUCTIVITY_SOLID_XZ",
                                                         "SOLID_COMPRESSIBILITY",
                                                         "DT_TEMPERATURE_COEFFICIENT",
                                                         "DT_TEMPERATURE",
                                                         "NORMAL_HEAT_FLUX",
                                                         "AIR_TEMPERATURE",
                                                         "SOLAR_RADIATION",
                                                         "AIR_HUMIDITY",
                                                         "PRECIPITATION",
                                                         "WIND_SPEED",
                                                         "A1_COEFFICIENT",
                                                         "A2_COEFFICIENT",
                                                         "A3_COEFFICIENT",
                                                         "ALPHA_COEFFICIENT",
                                                         "QF_COEFFICIENT",
                                                         "SMIN_COEFFICIENT",
                                                         "SMAX_COEFFICIENT"};
    for (const auto& name : variable_names) {
        KRATOS_EXPECT_TRUE(KratosComponents<VariableData>::Has(name))
    }
}

}
