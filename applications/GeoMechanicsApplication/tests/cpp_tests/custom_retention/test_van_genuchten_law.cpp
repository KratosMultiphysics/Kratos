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
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(VanGenuchtenLawReturnsCloneOfCorrectType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law         = VanGenuchtenLaw();
    const auto p_law_clone = law.Clone();
    KRATOS_EXPECT_NE(&law, p_law_clone.get());
    KRATOS_EXPECT_NE(dynamic_cast<const VanGenuchtenLaw*>(p_law_clone.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(VanGenuchtenLawReturnsCalculatedValues, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law = VanGenuchtenLaw();
    Properties properties;
    properties.SetValue(SATURATED_SATURATION, 0.9);
    properties.SetValue(RESIDUAL_SATURATION, 0.1);
    properties.SetValue(MINIMUM_RELATIVE_PERMEABILITY, 0.05);
    properties.SetValue(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE, 2.5);
    properties.SetValue(VAN_GENUCHTEN_GN, 2.5);
    properties.SetValue(VAN_GENUCHTEN_GL, 1.5);

    auto retention_law_parameters = RetentionLaw::Parameters{properties};

    retention_law_parameters.SetFluidPressure(-10.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateSaturation(retention_law_parameters), 0.9);
    double value = 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, DEGREE_OF_SATURATION, value), 0.9);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateEffectiveSaturation(retention_law_parameters), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, EFFECTIVE_SATURATION, value), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateBishopCoefficient(retention_law_parameters), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, BISHOP_COEFFICIENT, value), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateDerivativeOfSaturation(retention_law_parameters), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, DERIVATIVE_OF_SATURATION, value), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateRelativePermeability(retention_law_parameters), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, RELATIVE_PERMEABILITY, value), 1.0);

    // Values below are regression values, avoiding reimplementation here of the Van Genuchten Law
    retention_law_parameters.SetFluidPressure(1.5);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateSaturation(retention_law_parameters), 0.79023542376288392);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, DEGREE_OF_SATURATION, value),
                            0.79023542376288392);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateEffectiveSaturation(retention_law_parameters), 0.86279427970360489);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, EFFECTIVE_SATURATION, value),
                            0.86279427970360489);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateBishopCoefficient(retention_law_parameters), 0.86279427970360489);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, BISHOP_COEFFICIENT, value),
                            0.86279427970360489);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateDerivativeOfSaturation(retention_law_parameters), -0.15050611026881838);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, DERIVATIVE_OF_SATURATION, value),
                            -0.15050611026881838);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateRelativePermeability(retention_law_parameters), 0.28755984470352691);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, RELATIVE_PERMEABILITY, value),
                            0.28755984470352691);
}

KRATOS_TEST_CASE_IN_SUITE(VanGenuchtenLawChecksInputParameters, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetId(1);
    const auto process_info = ProcessInfo{};
    auto       law          = VanGenuchtenLaw();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "SATURATED_SATURATION is not available in the parameters of material 1.");
    properties.SetValue(SATURATED_SATURATION, 1.1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "SATURATED_SATURATION (1.1) must be in the range [0.0, 1.0] for material 1.");
    properties.SetValue(SATURATED_SATURATION, 0.9);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "RESIDUAL_SATURATION is not available in the parameters of material 1.");
    properties.SetValue(RESIDUAL_SATURATION, 1.1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "RESIDUAL_SATURATION (1.1) must be in the range [0.0, 0.9> for material 1.");
    properties.SetValue(RESIDUAL_SATURATION, 0.1);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "MINIMUM_RELATIVE_PERMEABILITY is not available in the parameters of material 1.");
    properties.SetValue(MINIMUM_RELATIVE_PERMEABILITY, 1.1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "MINIMUM_RELATIVE_PERMEABILITY (1.1) must be in the range [0.0, 1.0] for material 1.");
    properties.SetValue(MINIMUM_RELATIVE_PERMEABILITY, 0.05);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_AIR_ENTRY_PRESSURE is not available in the parameters of material 1.");
    properties.SetValue(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE, -4.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_AIR_ENTRY_PRESSURE (-4) must be greater than 0 for material 1.");
    properties.SetValue(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE, 4.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_GN is not available in the parameters of material 1.");
    properties.SetValue(VAN_GENUCHTEN_GN, -2.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_GN (-2.5) must be greater than 0 for material 1.");
    properties.SetValue(VAN_GENUCHTEN_GN, 2.5);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_GL is not available in the parameters of material 1.");
    properties.SetValue(VAN_GENUCHTEN_GL, -1.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "VAN_GENUCHTEN_GL (-1.5) must be greater than or equal to 0 for material 1.");
    properties.SetValue(VAN_GENUCHTEN_GL, 1.5);

    KRATOS_EXPECT_EQ(law.Check(properties, process_info), 0);
}

} // namespace Kratos::Testing