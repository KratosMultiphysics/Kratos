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

#include "custom_retention/saturated_below_phreatic_level_law.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace
{

using namespace Kratos;


SaturatedBelowPhreaticLevelLaw CreateSaturatedBelowPhreaticLevelLaw()
{
    return SaturatedBelowPhreaticLevelLaw{};
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SaturatedBelowPhreaticLevelLawReturnsCloneOfCorrectType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law         = CreateSaturatedBelowPhreaticLevelLaw();
    const auto p_law_clone = law.Clone();
    KRATOS_EXPECT_NE(&law, p_law_clone.get());
    KRATOS_EXPECT_NE(dynamic_cast<const SaturatedBelowPhreaticLevelLaw*>(p_law_clone.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(SaturatedBelowPhreaticLevelLawReturnsSaturation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law         = CreateSaturatedBelowPhreaticLevelLaw();
    Properties properties;
    properties.SetValue(SATURATED_SATURATION, 0.9);
    properties.SetValue(RESIDUAL_SATURATION, 0.1);
    properties.SetValue(MINIMUM_RELATIVE_PERMEABILITY, 0.05);

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

    retention_law_parameters.SetFluidPressure(10.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateSaturation(retention_law_parameters), 0.1);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, DEGREE_OF_SATURATION, value), 0.1);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateEffectiveSaturation(retention_law_parameters), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, EFFECTIVE_SATURATION, value), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateBishopCoefficient(retention_law_parameters), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, BISHOP_COEFFICIENT, value), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateDerivativeOfSaturation(retention_law_parameters), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, DERIVATIVE_OF_SATURATION, value), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateRelativePermeability(retention_law_parameters), 0.05);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, RELATIVE_PERMEABILITY, value), 0.05);
}
} // namespace Kratos::Testing