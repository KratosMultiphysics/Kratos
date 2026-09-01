// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#include "custom_retention/van_genuchten_law.h"
#include "geo_mechanics_application_variables.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "gtest/gtest.h"

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, VanGenuchtenLawReturnsCloneOfCorrectType)
{
    const auto law         = VanGenuchtenLaw();
    const auto p_law_clone = law.Clone();
    KRATOS_EXPECT_NE(&law, p_law_clone.get());
    KRATOS_EXPECT_NE(dynamic_cast<const VanGenuchtenLaw*>(p_law_clone.get()), nullptr);
}

struct TestData {
    const Variable<double>& variable;
    double                  expected_value;
    double                  input_pressure;
};

class VanGenuchtenCalculateValuesSuite : public testing::TestWithParam<TestData>
{
};

TEST_P(VanGenuchtenCalculateValuesSuite, TestValuesAreCalculatedCorrectly)
{
    const auto& test_data = GetParam();
    auto        law       = VanGenuchtenLaw();
    Properties  properties;
    properties.SetValue(SATURATED_SATURATION, 0.9);
    properties.SetValue(RESIDUAL_SATURATION, 0.1);
    properties.SetValue(MINIMUM_RELATIVE_PERMEABILITY, 0.05);
    properties.SetValue(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE, 2.5);
    properties.SetValue(VAN_GENUCHTEN_GN, 2.5);
    properties.SetValue(VAN_GENUCHTEN_GL, 1.5);
    auto retention_law_parameters = RetentionLaw::Parameters{properties};

    retention_law_parameters.SetFluidPressure(test_data.input_pressure);
    double value = 0.0;
    EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, test_data.variable, value),
                     test_data.expected_value)
        << "Incorrect value for variable = " << test_data.variable.Name()
        << " and input pressure = " << test_data.input_pressure;
}

INSTANTIATE_TEST_SUITE_P(
    SaturatedZone,
    VanGenuchtenCalculateValuesSuite,
    ::testing::Values(
        TestData{.variable = DEGREE_OF_SATURATION, .expected_value = 0.9, .input_pressure = -10.0},
        TestData{.variable = EFFECTIVE_SATURATION, .expected_value = 1.0, .input_pressure = -10.0},
        TestData{.variable = BISHOP_COEFFICIENT, .expected_value = 1.0, .input_pressure = -10.0},
        TestData{.variable = DERIVATIVE_OF_SATURATION, .expected_value = 0.0, .input_pressure = -10.0},
        TestData{.variable = RELATIVE_PERMEABILITY, .expected_value = 1.0, .input_pressure = -10.0},
        TestData{.variable = DEGREE_OF_SATURATION, .expected_value = 0.79023542376288392, .input_pressure = 1.5},
        TestData{.variable = EFFECTIVE_SATURATION, .expected_value = 0.86279427970360489, .input_pressure = 1.5},
        TestData{.variable = BISHOP_COEFFICIENT, .expected_value = 0.86279427970360489, .input_pressure = 1.5},
        TestData{.variable = DERIVATIVE_OF_SATURATION, .expected_value = -0.15050611026881838, .input_pressure = 1.5},
        TestData{.variable = RELATIVE_PERMEABILITY, .expected_value = 0.28755984470352691, .input_pressure = 1.5}));

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, VanGenuchtenLawChecksInputParameters)
{
    Properties properties;
    properties.SetId(1);
    const auto process_info = ProcessInfo{};
    auto       law          = VanGenuchtenLaw();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "SATURATED_SATURATION does not exist in the parameters of material with Id 1.");
    properties.SetValue(SATURATED_SATURATION, 1.1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "SATURATED_SATURATION in the parameters of material with Id 1 has an "
        "invalid value: 1.1 is out of the range [0, 1].");
    properties.SetValue(SATURATED_SATURATION, 0.9);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "RESIDUAL_SATURATION does not exist in the parameters of material with Id 1.");
    properties.SetValue(RESIDUAL_SATURATION, 1.1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "RESIDUAL_SATURATION in the parameters of material with Id 1 has an "
        "invalid value: 1.1 is out of the range [0, 0.9).");
    properties.SetValue(RESIDUAL_SATURATION, 0.1);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "MINIMUM_RELATIVE_PERMEABILITY does not exist in the parameters of material with Id 1.");
    properties.SetValue(MINIMUM_RELATIVE_PERMEABILITY, 1.1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, process_info),
                                      "MINIMUM_RELATIVE_PERMEABILITY in the parameters of material "
                                      "with Id 1 has an invalid value: 1.1 "
                                      "is out of the range [0, 1].");
    properties.SetValue(MINIMUM_RELATIVE_PERMEABILITY, 0.05);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_AIR_ENTRY_PRESSURE does not exist in the parameters of material with Id 1.");
    properties.SetValue(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE, -4.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, process_info),
                                      "VAN_GENUCHTEN_AIR_ENTRY_PRESSURE in the parameters of "
                                      "material with Id 1 has an invalid value: -4 "
                                      "is out of the range (0, -).");
    properties.SetValue(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE, 4.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_GN does not exist in the parameters of material with Id 1.");
    properties.SetValue(VAN_GENUCHTEN_GN, -2.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_GN in the parameters of material with Id 1 has an "
        "invalid value: -2.5 is out of the range (0, -).");
    properties.SetValue(VAN_GENUCHTEN_GN, 2.5);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_GL does not exist in the parameters of material with Id 1.");
    properties.SetValue(VAN_GENUCHTEN_GL, 1.5);

    KRATOS_EXPECT_EQ(law.Check(properties, process_info), 0);
}

} // namespace Kratos::Testing
