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

#include "geo_mechanics_application.h"
#include "geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ThermalAnalysisVariablesExistAfterRegistration, KratosGeoMechanicsFastSuiteWithoutKernel)
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

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementsAreAvailableAfterGeoAppRegistration, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KratosGeoMechanicsApplication geo_app;
    const auto element_type_names = std::vector<std::string>{"UPw2DLineInterfaceElement2Plus2N",
                                                             "UPw2DLineInterfaceElement3Plus3N"};

    for (const auto& r_name : element_type_names) {
        KRATOS_EXPECT_FALSE(KratosComponents<Element>::Has(r_name))
    }

    geo_app.Register();

    for (const auto& r_name : element_type_names) {
        KRATOS_EXPECT_TRUE(KratosComponents<Element>::Has(r_name))
    }
}

} // namespace Kratos::Testing
