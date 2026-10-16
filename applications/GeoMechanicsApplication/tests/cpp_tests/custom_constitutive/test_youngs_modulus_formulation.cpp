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

#include "custom_constitutive/youngs_modulus_formulations.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_InitializeFormulation_Constant, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto formulation =
        GeoYoungsModulusFormulations::InitializeFormulation(GeoYoungsModulusFormulations::Constant::Name);

    KRATOS_EXPECT_TRUE(std::holds_alternative<GeoYoungsModulusFormulations::Constant>(formulation))
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_InitializeFormulation_SchanzVermeer,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto formulation = GeoYoungsModulusFormulations::InitializeFormulation(
        GeoYoungsModulusFormulations::SchanzVermeer::Name);

    KRATOS_EXPECT_TRUE(std::holds_alternative<GeoYoungsModulusFormulations::SchanzVermeer>(formulation))
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_InitializeFormulation_Error, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::string invalid_name = "InvalidModel";

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GeoYoungsModulusFormulations::InitializeFormulation(invalid_name),
                                      "Unknown GEO_YOUNGS_MODULUS_FORMULATION: InvalidModel");
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulusFormulation_FromProperties_HasValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetValue(GEO_YOUNGS_MODULUS_FORMULATION, GeoYoungsModulusFormulations::SchanzVermeer::Name);

    KRATOS_EXPECT_EQ(GeoYoungsModulusFormulations::GetYoungsModulusFormulation(properties),
                     GeoYoungsModulusFormulations::SchanzVermeer::Name);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulusFormulation_FromProperties_NoValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const Properties properties;

    KRATOS_EXPECT_EQ(GeoYoungsModulusFormulations::GetYoungsModulusFormulation(properties),
                     GeoYoungsModulusFormulations::Constant::Name);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulusFormulation_FromVariant,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const GeoYoungsModulusFormulations::YoungsModulusVariant variant =
        GeoYoungsModulusFormulations::SchanzVermeer{};

    KRATOS_EXPECT_EQ(GeoYoungsModulusFormulations::GetYoungsModulusFormulation(variant),
                     GeoYoungsModulusFormulations::SchanzVermeer::Name);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CalculateMinorPrincipalEffectiveStress_ValidVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_vector_finalized = UblasUtilities::CreateVector({-10.0, 5.0, 20.0, 0.0, 0.0, 0.0});

    const auto result =
        GeoYoungsModulusFormulations::CalculateMinorPrincipalEffectiveStress(stress_vector_finalized);
    constexpr auto expected_value = -10.0;
    KRATOS_EXPECT_DOUBLE_EQ(result, expected_value);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CalculateMinorPrincipalEffectiveStress_TooSmallVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_vector_finalized = Vector(3, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoYoungsModulusFormulations::CalculateMinorPrincipalEffectiveStress(stress_vector_finalized), "Could not compute principal stresses from stress vector with size 3. Expected 3 principal stresses, got 2.");
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CalculateYoungsModulusSchanzVermeer_BaseCase,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetValue(GEO_PRESSURE_REFERENCE, 100.0);
    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 0.5);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    constexpr auto reference_youngs_modulus = 20000.0;
    const auto stress_vector_finalized = UblasUtilities::CreateVector({0.0, 0.0, 5.0, 0.0, 0.0, 0.0});

    const auto result = GeoYoungsModulusFormulations::CalculateYoungsModulusSchanzVermeer(
        properties, reference_youngs_modulus, stress_vector_finalized);
    constexpr auto expected_value = 7684.64;
    KRATOS_EXPECT_RELATIVE_NEAR(result, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CalculateYoungsModulusSchanzVermeer_NegativeBaseError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetValue(GEO_PRESSURE_REFERENCE, 100.0);
    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 0.5);
    properties.SetValue(GEO_COHESION, 1.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 5.0);
    constexpr auto reference_youngs_modulus = 20000.0;
    const auto stress_vector_finalized = UblasUtilities::CreateVector({500.0, 200.0, 500.0, 0.0, 0.0, 0.0});

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GeoYoungsModulusFormulations::CalculateYoungsModulusSchanzVermeer(
                                          properties, reference_youngs_modulus, stress_vector_finalized),
                                      "Non-positive base for std::pow");
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulus_Dispatcher_Constant,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoYoungsModulusFormulations::YoungsModulusVariant variant = GeoYoungsModulusFormulations::Constant{};
    Properties properties;
    const auto dummy_stress = Vector(6, 0.0);

    constexpr auto reference_youngs_modulus = 15000.0;
    const auto     result                   = GeoYoungsModulusFormulations::GetYoungsModulus(
        variant, properties, reference_youngs_modulus, dummy_stress);

    KRATOS_EXPECT_DOUBLE_EQ(result, reference_youngs_modulus);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulus_Dispatcher_SchanzVermeer,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoYoungsModulusFormulations::YoungsModulusVariant variant =
        GeoYoungsModulusFormulations::SchanzVermeer{};
    Properties properties;
    properties.SetValue(GEO_PRESSURE_REFERENCE, 100.0);
    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 0.5);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    const auto stress_vector_finalized = UblasUtilities::CreateVector({0.0, 0.0, 5.0, 0.0, 0.0, 0.0});

    constexpr auto reference_youngs_modulus = 20000.0;
    const auto     result                   = GeoYoungsModulusFormulations::GetYoungsModulus(
        variant, properties, reference_youngs_modulus, stress_vector_finalized);

    KRATOS_EXPECT_NE(result, reference_youngs_modulus);
    constexpr auto expected_value = 7684.64;
    KRATOS_EXPECT_RELATIVE_NEAR(result, expected_value, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CheckInputData_PassesForConstant,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetValue(GEO_YOUNGS_MODULUS_FORMULATION, GeoYoungsModulusFormulations::Constant::Name);

    EXPECT_NO_THROW(GeoYoungsModulusFormulations::CheckInputData(properties));
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CheckInputData_FailsMissingParams,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetValue(GEO_YOUNGS_MODULUS_FORMULATION, GeoYoungsModulusFormulations::SchanzVermeer::Name);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoYoungsModulusFormulations::CheckInputData(properties),
        "GEO_PRESSURE_REFERENCE does not exist in the parameters of material with Id 0.")

    properties.SetValue(GEO_PRESSURE_REFERENCE, 0.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoYoungsModulusFormulations::CheckInputData(properties),
        "GEO_PRESSURE_REFERENCE in the parameters of material with Id 0 has an "
        "invalid value: 0 is out of the range (0, -).")

    properties.SetValue(GEO_PRESSURE_REFERENCE, 50.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoYoungsModulusFormulations::CheckInputData(properties),
        "GEO_STRESS_DEPENDENCY_EXPONENT does not exist in the parameters of material with Id 0.")

    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 0.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoYoungsModulusFormulations::CheckInputData(properties),
        "GEO_STRESS_DEPENDENCY_EXPONENT in the parameters of material with Id 0 has an "
        "invalid value: 0 is out of the range (0, -).")

    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoYoungsModulusFormulations::CheckInputData(properties),
        "GEO_COHESION does not exist in the parameters of material with Id 0.")

    properties.SetValue(GEO_COHESION, 20.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        GeoYoungsModulusFormulations::CheckInputData(properties),
        "GEO_FRICTION_ANGLE does not exist in the parameters of material with Id 0.")

    properties.SetValue(GEO_FRICTION_ANGLE, 45.0);
    EXPECT_NO_THROW(GeoYoungsModulusFormulations::CheckInputData(properties));
}

} // namespace Kratos::Testing