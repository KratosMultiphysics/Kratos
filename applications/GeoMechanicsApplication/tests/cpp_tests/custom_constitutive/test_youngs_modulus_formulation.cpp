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

#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/youngs_modulus_formulations.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"
#include "testing/testing.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_InitializeFormulation_Constant, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto formulation =
        GeoYoungsModulusFormulations::InitializeFormulation(GeoYoungsModulusFormulations::Constant::Name);

    KRATOS_EXPECT_TRUE(std::holds_alternative<GeoYoungsModulusFormulations::Constant>(formulation));
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_InitializeFormulation_Eur, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto formulation =
        GeoYoungsModulusFormulations::InitializeFormulation(GeoYoungsModulusFormulations::Eur::Name);

    KRATOS_EXPECT_TRUE(std::holds_alternative<GeoYoungsModulusFormulations::Eur>(formulation));
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_InitializeFormulation_Error, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::string invalid_name = "InvalidModel";

    // Expect an exception to be thrown
    KRATOS_CHECK_EXCEPTION_IS_THROWN(GeoYoungsModulusFormulations::InitializeFormulation(invalid_name),
                                     "Unknown GEO_YOUNGS_MODULUS_FORMULATION: InvalidModel");
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulusFormulation_FromProperties_HasValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetValue(GEO_YOUNGS_MODULUS_FORMULATION, GeoYoungsModulusFormulations::Eur::Name);

    KRATOS_EXPECT_EQ(GeoYoungsModulusFormulations::GetYoungsModulusFormulation(properties),
                     GeoYoungsModulusFormulations::Eur::Name);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulusFormulation_FromProperties_NoValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties; // Empty properties

    KRATOS_EXPECT_EQ(GeoYoungsModulusFormulations::GetYoungsModulusFormulation(properties),
                     GeoYoungsModulusFormulations::Constant::Name);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulusFormulation_FromVariant,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoYoungsModulusFormulations::YoungsModulusVariant var = GeoYoungsModulusFormulations::Eur{};

    KRATOS_EXPECT_EQ(GeoYoungsModulusFormulations::GetYoungsModulusFormulation(var),
                     GeoYoungsModulusFormulations::Eur::Name);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CalculateMinorPrincipalEffectiveStress_ValidVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_vector_finalized = UblasUtilities::CreateVector({-10.0, 5.0, 20.0, 0.0, 0.0, 0.0});

    double result = GeoYoungsModulusFormulations::CalculateMinorPrincipalEffectiveStress(stress_vector_finalized);

    KRATOS_EXPECT_DOUBLE_EQ(result, -10.0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CalculateMinorPrincipalEffectiveStress_TooSmallVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto stress_vector_finalized = Vector(3, 0.0);

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        GeoYoungsModulusFormulations::CalculateMinorPrincipalEffectiveStress(stress_vector_finalized), "Could not compute principal stresses from stress vector with size 3. Expected at least 3 principal stresses, got 2");
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CalculateYoungsModulusForEur_BaseCase,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetValue(GEO_PRESSURE_REFERENCE, 100.0);
    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 0.5);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    auto stress_vector_finalized = UblasUtilities::CreateVector({0.0, 0.0, 50.0, 0.0, 0.0, 0.0});

    double youngs_input = 20000.0;
    double result       = GeoYoungsModulusFormulations::CalculateYoungsModulusForEur(
        properties, youngs_input, stress_vector_finalized);

    // Manual calculation check:
    // friction_angle = 30 -> sin=0.5, cos=sqrt(3)/2
    // stress_shift = 10 * (sqrt(3)/2) / 0.5 = 10 * sqrt(3) ~ 17.32
    // base = (17.32 - 50) / (17.32 + 100) -> This will be negative! Should throw.

    // Let's use values that yield a positive base
    stress_vector_finalized[2] = 5.0; // Lower minor stress
    result = GeoYoungsModulusFormulations::CalculateYoungsModulusForEur(properties, youngs_input,
                                                                        stress_vector_finalized);

    KRATOS_EXPECT_GT(result, 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_CalculateYoungsModulusForEur_NegativeBaseError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties;
    properties.SetValue(GEO_PRESSURE_REFERENCE, 100.0);
    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 0.5);
    properties.SetValue(GEO_COHESION, 1.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 5.0);
    double reference_youngs_modulus = 20000.0;
    const auto stress_vector_finalized = UblasUtilities::CreateVector({500.0, 200.0, 500.0, 0.0, 0.0, 0.0});

    KRATOS_CHECK_EXCEPTION_IS_THROWN(GeoYoungsModulusFormulations::CalculateYoungsModulusForEur(
                                         properties, reference_youngs_modulus, stress_vector_finalized),
                                     "Non-positive base for std::pow");
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulus_Dispatcher_Constant,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoYoungsModulusFormulations::YoungsModulusVariant variant = GeoYoungsModulusFormulations::Constant{};
    Properties properties;
    Vector     dummy_stress = ZeroVector(6);

    double reference_youngs_modulus = 15000.0;
    double result                   = GeoYoungsModulusFormulations::GetYoungsModulus(
        variant, properties, reference_youngs_modulus, dummy_stress);

    KRATOS_EXPECT_DOUBLE_EQ(result, reference_youngs_modulus);
}

KRATOS_TEST_CASE_IN_SUITE(GeoYoungsModulusFormulations_GetYoungsModulus_Dispatcher_Eur,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoYoungsModulusFormulations::YoungsModulusVariant variant = GeoYoungsModulusFormulations::Eur{};
    Properties properties;
    properties.SetValue(GEO_PRESSURE_REFERENCE, 100.0);
    properties.SetValue(GEO_STRESS_DEPENDENCY_EXPONENT, 0.5);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    Vector stress_vector_finalized(6);
    stress_vector_finalized[2] = 5.0;

    double reference_youngs_modulus = 20000.0;
    double result                   = GeoYoungsModulusFormulations::GetYoungsModulus(
        variant, properties, reference_youngs_modulus, stress_vector_finalized);

    // Just ensure it runs without throwing and returns something plausible
    KRATOS_EXPECT_NE(result, reference_youngs_modulus); // Eur should degrade modulus based on strain/pressure logic
    KRATOS_EXPECT_GT(result, 0.0);
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
    properties.SetValue(GEO_YOUNGS_MODULUS_FORMULATION, GeoYoungsModulusFormulations::Eur::Name);
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