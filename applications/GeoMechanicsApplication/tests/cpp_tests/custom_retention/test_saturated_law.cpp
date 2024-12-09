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

#include "custom_retention/saturated_law.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace
{

using namespace Kratos;

SaturatedLaw CreateSaturatedLaw() { return SaturatedLaw{}; }

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SaturatedLawReturnsCloneOfCorrectType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law         = CreateSaturatedLaw();
    const auto p_law_clone = law.Clone();
    KRATOS_EXPECT_NE(&law, p_law_clone.get());
    KRATOS_EXPECT_NE(dynamic_cast<const SaturatedLaw*>(p_law_clone.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(SaturatedLawReturnsCalculatedValues, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law = CreateSaturatedLaw();
    Properties properties;
    auto       retention_law_parameters = RetentionLaw::Parameters{properties};

    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateSaturation(retention_law_parameters), 1.0);
    double value = 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, DEGREE_OF_SATURATION, value), 1.0);

    properties.SetValue(SATURATED_SATURATION, 0.9);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateSaturation(retention_law_parameters), 0.9);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, DEGREE_OF_SATURATION, value), 0.9);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateEffectiveSaturation(retention_law_parameters), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, EFFECTIVE_SATURATION, value), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateBishopCoefficient(retention_law_parameters), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, BISHOP_COEFFICIENT, value), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateDerivativeOfSaturation(retention_law_parameters), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, DERIVATIVE_OF_SATURATION, value), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateRelativePermeability(retention_law_parameters), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(law.CalculateValue(retention_law_parameters, RELATIVE_PERMEABILITY, value), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(SaturatedLawChecksInputParameters, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetId(1);
    const auto process_info = ProcessInfo{};
    auto       law          = CreateSaturatedLaw();

    KRATOS_EXPECT_EQ(law.Check(properties, process_info), 0);

    properties.SetValue(SATURATED_SATURATION, 1.1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, process_info),
        "SATURATED_SATURATION (1.1) must be in the range [0.0, 1.0] for material 1.");

    properties.SetValue(SATURATED_SATURATION, 0.9);
    KRATOS_EXPECT_EQ(law.Check(properties, process_info), 0);
}

} // namespace Kratos::Testing