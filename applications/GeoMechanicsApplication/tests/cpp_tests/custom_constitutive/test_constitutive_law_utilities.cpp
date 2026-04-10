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

#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_MakeInterfaceConstitutiveTensorReturnsConstitutiveTensor,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange & Act
    const auto constitutive_tensor =
        ConstitutiveLawUtilities::MakeInterfaceConstitutiveTensor(2.0, 1.0, 4, 2);

    // Assert
    auto expected_tensor = UblasUtilities::CreateMatrix(
        {{2.0, 0.0, 0.0, 0.0}, {0.0, 2.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}});
    KRATOS_EXPECT_MATRIX_EQ(constitutive_tensor, expected_tensor);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_MakeCotinuumConstitutiveTensorReturnsConstitutiveTensor,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange & Act
    const auto constitutive_tensor =
        ConstitutiveLawUtilities::MakeContinuumConstitutiveTensor(1.0, 0.25, 4, 2);

    // Assert
    auto expected_tensor = UblasUtilities::CreateMatrix(
        {{1.2, 0.4, 0.0, 0.0}, {0.4, 1.2, 0.0, 0.0}, {0.0, 0.0, 0.4, 0.0}, {0.0, 0.0, 0.0, 0.4}});
    KRATOS_EXPECT_MATRIX_EQ(constitutive_tensor, expected_tensor);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_CalculateK0NCFromFrictionAngleGivesK0NC,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange & Act
    const auto computed_K0NC =
        ConstitutiveLawUtilities::CalculateK0NCFromFrictionAngleInRadians(30.0 * Globals::Pi / 180.0);

    // Assert
    KRATOS_EXPECT_EQ(computed_K0NC, 0.5);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_GetUndrainedYoungsModulusGivesUndrainedYoungsModulus,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0);
    properties.SetValue(POISSON_RATIO, 0.2);

    // Act
    const auto undrained_youngs_modulus =
        ConstitutiveLawUtilities::GetUndrainedYoungsModulus(properties, 0.4);

    // Assert
    KRATOS_EXPECT_EQ(undrained_youngs_modulus, 7.0 / 6.0);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_GetSkemptonBGivesSkemptonB, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0);
    properties.SetValue(POISSON_RATIO, 0.2);
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.5);

    // Act
    auto skempton_b = ConstitutiveLawUtilities::GetSkemptonB(properties);

    // Assert
    KRATOS_EXPECT_EQ(skempton_b, 0.5);

    // Arrange
    properties.SetValue(SKEMPTON_B, 0.8);
    // Act
    skempton_b = ConstitutiveLawUtilities::GetSkemptonB(properties);

    // Assert
    KRATOS_EXPECT_EQ(skempton_b, 0.8);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_GetUndrainedPoissonsRatioGivesUndrainedPoissonsRatio,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(POISSON_RATIO, 0.2);
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.5);

    // Act
    auto undrained_poissons_ratio = ConstitutiveLawUtilities::GetUndrainedPoissonsRatio(properties);

    // Assert
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, 1.0 / 3.0, Defaults::absolute_tolerance);

    // Arrange
    properties.SetValue(SKEMPTON_B, 0.5);

    // Act
    undrained_poissons_ratio = ConstitutiveLawUtilities::GetUndrainedPoissonsRatio(properties);

    // Assert
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, 1.0 / 3.0, Defaults::absolute_tolerance);

    // Arrange
    properties.SetValue(POISSON_UNDRAINED, 0.4);

    // Act
    undrained_poissons_ratio = ConstitutiveLawUtilities::GetUndrainedPoissonsRatio(properties);

    // Assert
    KRATOS_EXPECT_NEAR(undrained_poissons_ratio, 0.4, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtitlities_UndrainedExcessPorePressureIncrementGivesDeltaPw,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties;
    properties.SetValue(POISSON_RATIO, 0.2);
    properties.SetValue(BIOT_COEFFICIENT, 1.0);
    properties.SetValue(BULK_MODULUS_FLUID, 1.E3);
    properties.SetValue(BULK_MODULUS_SOLID, 2.E3);
    properties.SetValue(POROSITY, 0.5);
}

} // namespace Kratos::Testing