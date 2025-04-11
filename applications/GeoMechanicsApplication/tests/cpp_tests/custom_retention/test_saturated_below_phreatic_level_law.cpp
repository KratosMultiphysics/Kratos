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
#include "custom_utilities/registration_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <string>

using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SaturatedBelowPhreaticLevelLawReturnsCloneOfCorrectType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law         = SaturatedBelowPhreaticLevelLaw();
    const auto p_law_clone = law.Clone();
    KRATOS_EXPECT_NE(&law, p_law_clone.get());
    KRATOS_EXPECT_NE(dynamic_cast<const SaturatedBelowPhreaticLevelLaw*>(p_law_clone.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(SaturatedBelowPhreaticLevelLawReturnsCalculatedValues, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law = SaturatedBelowPhreaticLevelLaw();
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

KRATOS_TEST_CASE_IN_SUITE(SaturatedBelowPhreaticLevelLawChecksInputParameters, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetId(1);
    const auto process_info = ProcessInfo{};
    auto       law          = SaturatedBelowPhreaticLevelLaw();

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

    KRATOS_EXPECT_EQ(law.Check(properties, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(SaturatedBelowPhreaticLevelLaw_CanBeSavedAndLoaded, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration = ScopedSerializerRegistration{
        std::make_pair("SaturatedBelowPhreaticLevelLaw"s, SaturatedBelowPhreaticLevelLaw{})};
    const auto p_retention_law =
        std::unique_ptr<RetentionLaw>{std::make_unique<SaturatedBelowPhreaticLevelLaw>()};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_retention_law);
    auto p_loaded_retention_law = std::unique_ptr<RetentionLaw>{};
    serializer.load("test_tag"s, p_loaded_retention_law);

    // Assert
    KRATOS_EXPECT_NE(p_loaded_retention_law, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<SaturatedBelowPhreaticLevelLaw*>(p_loaded_retention_law.get()), nullptr);
}

} // namespace Kratos::Testing