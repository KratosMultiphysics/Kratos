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
//                   Mohamed Nabi
//

#include "custom_constitutive/mohr_coulomb_with_tension_cutoff.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_utilities/registration_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{
Vector CalculateMappedStressVector(Vector&                       rCauchyStressVector,
                                   ConstitutiveLaw::Parameters&  rParameters,
                                   MohrCoulombWithTensionCutOff& rLaw)
{
    Vector strain_vector = ZeroVector(4);
    rParameters.SetStrainVector(strain_vector);
    rParameters.SetStressVector(rCauchyStressVector);
    const auto dummy_process_info = ProcessInfo{};
    rLaw.SetValue(CAUCHY_STRESS_VECTOR, rCauchyStressVector, dummy_process_info);
    rLaw.FinalizeMaterialResponseCauchy(rParameters);
    rLaw.CalculateMaterialResponseCauchy(rParameters);

    Vector result;
    rLaw.GetValue(CAUCHY_STRESS_VECTOR, result);
    return result;
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_Clone, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto original_law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    // Act
    auto p_cloned_law = original_law.Clone();

    // Assert
    KRATOS_EXPECT_NE(p_cloned_law.get(), nullptr);
    KRATOS_EXPECT_NE(p_cloned_law.get(), &original_law);
    KRATOS_EXPECT_NE(dynamic_cast<const MohrCoulombWithTensionCutOff*>(p_cloned_law.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_Check, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto                        law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    ConstitutiveLaw::Parameters parameters;
    Properties                  properties(3);
    parameters.SetMaterialProperties(properties);
    const auto element_geometry = Geometry<Node>{};
    const auto process_info     = ProcessInfo{};

    // Act & Assert
    properties.SetValue(GEO_COHESION, 1.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "GEO_TENSILE_STRENGTH does not exist in the property with Id 3.")
    properties.SetValue(GEO_TENSILE_STRENGTH, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "GEO_TENSILE_STRENGTH in the property with Id 3 has an invalid value: -1 is out of the range [0, 1.73205].")
    properties.SetValue(GEO_TENSILE_STRENGTH, 2.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "GEO_TENSILE_STRENGTH in the property with Id 3 has an invalid value: 2 is out of the range [0, 1.73205].")
    properties.SetValue(GEO_TENSILE_STRENGTH, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "YOUNG_MODULUS does not exist in the property with Id 3.")
    properties.SetValue(YOUNG_MODULUS, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "YOUNG_MODULUS in the property with Id 3 has an invalid value: -1 is out of the range [0, -).")
    properties.SetValue(YOUNG_MODULUS, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "POISSON_RATIO does not exist in the property with Id 3.")
    properties.SetValue(POISSON_RATIO, -0.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "POISSON_RATIO in the property with Id 3 has an invalid value: -0.5 is out of the range [0, 0.5].")
    properties.SetValue(POISSON_RATIO, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), " POISSON_RATIO in the property with Id 3 has an invalid value: 1 is out of the range [0, 0.5].")
    properties.SetValue(POISSON_RATIO, 0.3);

    KRATOS_EXPECT_EQ(law.Check(properties, element_geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtElasticZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto cauchy_stress_vector          = UblasUtilities::CreateVector({6.0, 0.0, -10.0, 0.0});
    auto expected_cauchy_stress_vector = UblasUtilities::CreateVector({6.0, 0.0, -10.0, 0.0});
    KRATOS_EXPECT_VECTOR_EQ(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                            expected_cauchy_stress_vector);

    cauchy_stress_vector          = UblasUtilities::CreateVector({8.0, 6.0, 4.0, 0.0});
    expected_cauchy_stress_vector = UblasUtilities::CreateVector({8.0, 6.0, 4.0, 0.0});
    KRATOS_EXPECT_VECTOR_EQ(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                            expected_cauchy_stress_vector);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.25);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto cauchy_stress_vector = UblasUtilities::CreateVector({8.0, 0.0, -12.0, 0.0});
    auto expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({7.338673315592010089, 0.0, -11.338673315592010089, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector = UblasUtilities::CreateVector({12.0, 0.0, -16.0, 0.0});
    expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({7.338673315592010089, 0.0, -11.338673315592010089, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    // LOWEST_PRINCIPAL_STRESSES as averaging type
    cauchy_stress_vector = UblasUtilities::CreateVector({12.0, 10.0, -16.0, 0.0});
    expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({7.806379130008, 7.806379130008, -9.61275826001616129, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance * 1.0e3);

    // HIGHEST_PRINCIPAL_STRESSES as averaging type
    cauchy_stress_vector = UblasUtilities::CreateVector({12.0, -12.0, -16.0, 0.0});
    expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({7.259759295835, -11.6298796479175, -11.6298796479175, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance * 1.0e3);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCornerReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.25);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    // The following stress vector will exercise the "no averaging" type
    auto cauchy_stress_vector = UblasUtilities::CreateVector({18.0, 8.0, -2.0, 0.0});
    auto expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({10.0, 6.12052019550083, -1.5179192179966735, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    // The following stress vector will exercise the "lowest averaging" type
    cauchy_stress_vector = UblasUtilities::CreateVector({24.0, 22.0, -8.0, 0.0});
    expected_cauchy_stress_vector = UblasUtilities::CreateVector({10.0, 10.0, -1.5179192179966735, 0.0});
    constexpr double tolerance = 1.0e-10;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, tolerance);

    // The following stress vector will exercise the "highest averaging" type
    cauchy_stress_vector = UblasUtilities::CreateVector({24.0, -6.0, -8.0, 0.0});
    expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({8.48187256360842, -7.12007108043552, -7.12007108043552, 0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionCutoffReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto cauchy_stress_vector          = UblasUtilities::CreateVector({12.0, 9.0, 8.0, 0.0});
    auto expected_cauchy_stress_vector = UblasUtilities::CreateVector({10.0, 9.0, 8.0, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector          = UblasUtilities::CreateVector({14.0, 9.0, 6.0, 0.0});
    expected_cauchy_stress_vector = UblasUtilities::CreateVector({10.0, 9.0, 6.0, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector          = UblasUtilities::CreateVector({14.0, 12.0, 6.0, 0.0});
    expected_cauchy_stress_vector = UblasUtilities::CreateVector({10.0, 10.0, 6.0, 0.0});
    constexpr double tolerance    = 1.0e-10;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionApexReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);
    constexpr double tolerance = 1.0e-10;

    // Act and Assert
    auto cauchy_stress_vector          = UblasUtilities::CreateVector({19.0, 12.0, 11.0, 0.0});
    auto expected_cauchy_stress_vector = UblasUtilities::CreateVector({10.0, 10.0, 10.0, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, tolerance);

    cauchy_stress_vector          = UblasUtilities::CreateVector({11.5, 10.0, 10.5, 0.0});
    expected_cauchy_stress_vector = UblasUtilities::CreateVector({10.0, 10.0, 10.0, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyWithLargeTensileStrength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.0);
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 20.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto       cauchy_stress_vector = UblasUtilities::CreateVector({22.0, 20.0, 18.0, 0.0});
    const auto expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({14.2814800674211450, 14.2814800674211450, 14.2814800674211450, 0.0});
    constexpr double tolerance = 1.0e-10;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_Serialization, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_law = std::unique_ptr<ConstitutiveLaw>{
        std::make_unique<MohrCoulombWithTensionCutOff>(std::make_unique<PlaneStrain>())};
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    p_law->InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    auto   cauchy_stress_vector = UblasUtilities::CreateVector({6.0, 0.0, -10.0, 0.0});
    Vector strain_vector        = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    const auto dummy_process_info = ProcessInfo{};
    p_law->SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, dummy_process_info);
    p_law->InitializeMaterialResponseCauchy(parameters);
    p_law->CalculateMaterialResponseCauchy(parameters);

    Vector calculated_cauchy_stress_vector;
    p_law->GetValue(CAUCHY_STRESS_VECTOR, calculated_cauchy_stress_vector);

    const auto scoped_registration = ScopedSerializerRegistration{
        std::make_pair("PlaneStrain"s, PlaneStrain{}),
        std::make_pair("MohrCoulombWithTensionCutOff"s, MohrCoulombWithTensionCutOff{})};
    auto serializer = StreamSerializer{};

    ConstitutiveLaw::Parameters parameters_to_be_ignored;
    Vector                      any_cauchy_stress_vector = ZeroVector(4);
    Vector                      any_strain_vector        = UnitVector(4);
    parameters_to_be_ignored.SetStrainVector(any_cauchy_stress_vector);
    parameters_to_be_ignored.SetStressVector(any_strain_vector);

    // Act
    serializer.save("test_tag"s, p_law);
    auto p_loaded_law = std::unique_ptr<ConstitutiveLaw>();
    serializer.load("test_tag"s, p_loaded_law);

    // Assert
    Vector loaded_calculated_cauchy_stress_vector;
    p_loaded_law->GetValue(CAUCHY_STRESS_VECTOR, loaded_calculated_cauchy_stress_vector);
    KRATOS_EXPECT_VECTOR_EQ(loaded_calculated_cauchy_stress_vector, calculated_cauchy_stress_vector);

    p_loaded_law->InitializeMaterialResponseCauchy(parameters_to_be_ignored);
    p_loaded_law->CalculateMaterialResponseCauchy(parameters);
    p_loaded_law->GetValue(CAUCHY_STRESS_VECTOR, loaded_calculated_cauchy_stress_vector);
    KRATOS_EXPECT_VECTOR_EQ(loaded_calculated_cauchy_stress_vector, calculated_cauchy_stress_vector);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateConstitutiveMatrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0e8);
    properties.SetValue(POISSON_RATIO, 0.3);
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    parameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    Vector stress_vector = ZeroVector(4);
    parameters.SetStressVector(stress_vector);
    Matrix constitutive_matrix = ZeroMatrix(4, 4);
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act
    law.InitializeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    auto expected_constitutive_matrix = UblasUtilities::CreateMatrix({{1.35E8, 5.77E7, 5.77E7, 0.},
                                                                      {5.77E7, 1.35E8, 5.77E7, 0.},
                                                                      {5.77E7, 5.77E7, 1.35E8, 0.},
                                                                      {0., 0., 0., 3.85E7}});
    KRATOS_EXPECT_MATRIX_NEAR(constitutive_matrix, expected_constitutive_matrix, 1.E6);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_WorkingSpaceDimensionDependsOnConstitutiveLawDimension,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()}.WorkingSpaceDimension(), 2);
    KRATOS_EXPECT_EQ(
        MohrCoulombWithTensionCutOff{std::make_unique<ThreeDimensional>()}.WorkingSpaceDimension(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_StressMeasureIsAlwaysCauchy, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()}.GetStressMeasure(),
                     ConstitutiveLaw::StressMeasure_Cauchy);
    KRATOS_EXPECT_EQ(MohrCoulombWithTensionCutOff{std::make_unique<ThreeDimensional>()}.GetStressMeasure(),
                     ConstitutiveLaw::StressMeasure_Cauchy);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_StrainSizeDependsOnConstitutiveLawDimension,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()}.GetStrainSize(), 4);
    KRATOS_EXPECT_EQ(MohrCoulombWithTensionCutOff{std::make_unique<ThreeDimensional>()}.GetStrainSize(), 6);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_StrainMeasureIsAlwaysInfinitesimal,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()}.GetStrainMeasure(),
                     ConstitutiveLaw::StrainMeasure_Infinitesimal);
    KRATOS_EXPECT_EQ(MohrCoulombWithTensionCutOff{std::make_unique<ThreeDimensional>()}.GetStrainMeasure(),
                     ConstitutiveLaw::StrainMeasure_Infinitesimal);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_HasIncrementalFormulation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    EXPECT_TRUE(MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()}.IsIncremental());
    EXPECT_TRUE(MohrCoulombWithTensionCutOff{std::make_unique<ThreeDimensional>()}.IsIncremental());
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_RequiresInitializeMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    EXPECT_TRUE(MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()}.RequiresInitializeMaterialResponse());
    EXPECT_TRUE(MohrCoulombWithTensionCutOff{std::make_unique<ThreeDimensional>()}.RequiresInitializeMaterialResponse());
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_GetPlaneStrainLawFeatures, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange & act
    ConstitutiveLaw::Features features;
    MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()}.GetLawFeatures(features);

    // Assert
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::INFINITESIMAL_STRAINS))
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::ISOTROPIC))
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::PLANE_STRAIN_LAW))
    KRATOS_EXPECT_EQ(features.GetStrainMeasures().front(), ConstitutiveLaw::StrainMeasure_Infinitesimal);
    KRATOS_EXPECT_EQ(features.GetStrainSize(), 4);
    KRATOS_EXPECT_EQ(features.GetSpaceDimension(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_Get3DLawFeatures, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange & act
    ConstitutiveLaw::Features features;
    MohrCoulombWithTensionCutOff{std::make_unique<ThreeDimensional>()}.GetLawFeatures(features);

    // Assert
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::INFINITESIMAL_STRAINS))
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::ISOTROPIC))
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::THREE_DIMENSIONAL_LAW))
    KRATOS_EXPECT_EQ(features.GetStrainMeasures().front(), ConstitutiveLaw::StrainMeasure_Infinitesimal);
    KRATOS_EXPECT_EQ(features.GetStrainSize(), 6);
    KRATOS_EXPECT_EQ(features.GetSpaceDimension(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtElasticZoneInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto cauchy_stress_vector = UblasUtilities::CreateVector({-6.0, 0.0, 0.0, 8.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector = UblasUtilities::CreateVector({8.0, 0.0, 0.0, 3.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureZoneInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto       cauchy_stress_vector = UblasUtilities::CreateVector({-6.0, 0.0, 0.0, 18.0});
    const auto expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({-4.62956382113722, -1.37043617886278, 0.0, 9.77738292682331});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionCutoffReturnZoneInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto       cauchy_stress_vector = UblasUtilities::CreateVector({18.0, 0.0, 0.0, -2.0});
    const auto expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({9.87832130144553, -0.0978657587384199, 0.0, -1.10846522890933});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_TrialStressInDegeneratedTensileApexReturnZoneIsReturnedToApex,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto           law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties     properties;
    constexpr auto phi_in_degrees = 30.0;
    properties.SetValue(GEO_FRICTION_ANGLE, phi_in_degrees);
    constexpr auto cohesion = 10.0;
    properties.SetValue(GEO_COHESION, cohesion);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    const auto tensile_strength = cohesion / std::tan(MathUtils<>::DegreesToRadians(phi_in_degrees));
    properties.SetValue(GEO_TENSILE_STRENGTH, tensile_strength);
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.15);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto cauchy_stress_vector = UblasUtilities::CreateVector(
        {tensile_strength + 20.0, tensile_strength + 10.0, tensile_strength, 0.0});
    const auto expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({tensile_strength, tensile_strength, tensile_strength, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCornerReturnZoneWithShearComponent,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.25);
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    // The following stress vector will exercise the "highest averaging" type
    auto       cauchy_stress_vector          = UblasUtilities::CreateVector({18.0, 0.0, 0.0, -8.0});
    const auto expected_cauchy_stress_vector = UblasUtilities::CreateVector(
        {8.54534046868983, -0.0632596866864898, -1.51791921799666, -3.82604451350059});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionApexReturnZoneInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, -1.0);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto cauchy_stress_vector = UblasUtilities::CreateVector({10.0, 0.0, 0.0, -2.0});
    const auto expected_cauchy_stress_vector = UblasUtilities::CreateVector({-1.0, -1.0, -1.0, 0.0});
    constexpr double tolerance = 1.0e-10;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureZoneWithLinearHardening,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_COULOMB_HARDENING_TYPE, "Linear");
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.0);
    properties.SetValue(GEO_MAX_PLASTIC_ITERATIONS, 100);
    properties.SetValue(GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({0.0}));
    properties.SetValue(GEO_COHESION_FUNCTION_COEFFICIENTS,
                        UblasUtilities::CreateVector({(2.5 - std::sqrt(3.0)) * 1.0e6}));
    properties.SetValue(GEO_DILATANCY_ANGLE_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({0.0}));

    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto cauchy_stress_vector          = UblasUtilities::CreateVector({10.0, 0.0, -40.0, 0.0});
    auto expected_cauchy_stress_vector = UblasUtilities::CreateVector({5.0, 0.0, -35.0, 0.0});
    constexpr auto tolerance           = 1.0e-8;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, tolerance);

    // Arrange
    properties.SetValue(GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({1.0e5}));
    properties.SetValue(GEO_COHESION_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({0.0}));
    parameters.SetMaterialProperties(properties);
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    cauchy_stress_vector = UblasUtilities::CreateVector({10.0, 0.0, -40.0, 0.0});
    expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({2.780289249892, 2.780289249892, -35.5605784998, 0}); // regression values
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, tolerance);

    // Arrange
    properties.SetValue(GEO_FRICTION_ANGLE_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({0.0}));
    properties.SetValue(GEO_DILATANCY_ANGLE_FUNCTION_COEFFICIENTS, UblasUtilities::CreateVector({1.0}));
    parameters.SetMaterialProperties(properties);
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    cauchy_stress_vector = UblasUtilities::CreateVector({10.0, 0.0, -40.0, 0.0});
    expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({1.16025403784, 0.0, -31.16025403784, 0.0}); // regression values
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, tolerance);
}

Vector ComputeStressVectorUsingCPhiReductionTestData(double  Cohesion,
                                                     double  FrictionAngle,
                                                     Vector& rStrainVectorFinalized,
                                                     Vector& rStressVectorFinalized,
                                                     Vector& rStrainVector)
{
    Properties properties;
    properties.SetValue(GEO_COHESION, Cohesion);
    properties.SetValue(GEO_FRICTION_ANGLE, FrictionAngle);
    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 1000.0);
    properties.SetValue(YOUNG_MODULUS, 30.0e6);
    properties.SetValue(POISSON_RATIO, 0.2);

    auto law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    parameters.SetStrainVector(rStrainVectorFinalized);
    parameters.SetStressVector(rStressVectorFinalized);
    const auto dummy_process_info = ProcessInfo{};
    law.SetValue(CAUCHY_STRESS_VECTOR, rStressVectorFinalized, dummy_process_info);
    law.FinalizeMaterialResponseCauchy(parameters);

    parameters.SetStrainVector(rStrainVector);
    law.CalculateMaterialResponseCauchy(parameters);

    Vector result;
    law.GetValue(CAUCHY_STRESS_VECTOR, result);
    return result;
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureRegion_CPhiVariation,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange, data are taken from C-Phi reduction test for GPoint 0 of Element 8
    constexpr auto cohesion                = 800.0;
    constexpr auto friction_angle          = 24.79128089714489;
    auto           strain_vector_finalized = UblasUtilities::CreateVector(
        {6.9595905431284102e-04, -5.1669244190907055e-04, 0.0000000000000000e+00, 3.6495572640421296e-05});
    auto stress_vector_finalized = UblasUtilities::CreateVector(
        {-1.8485198418257996e+03, -7.9720205330738399e+03, -3.6748820224188021e+03, 7.4955824230365451e+01});
    auto strain_vector = UblasUtilities::CreateVector({6.9595905431284102e-04, -5.1669244190907055e-04,
                                                       0.0000000000000000e+00, 3.6495572640421296e-05});
    // Act
    const auto actual_cauchy_stress_vector = ComputeStressVectorUsingCPhiReductionTestData(
        cohesion, friction_angle, strain_vector_finalized, stress_vector_finalized, strain_vector);
    // Assert
    const auto expected_cauchy_stress_vector = UblasUtilities::CreateVector(
        {-2.1258867028468462e+03, -7.6946536720527947e+03, -3.6748820224188021e+03, 6.8165505185660834e+01});
    KRATOS_EXPECT_VECTOR_NEAR(actual_cauchy_stress_vector, expected_cauchy_stress_vector,
                              Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureRegion,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange, data are taken from C-Phi reduction test for GPoint 0 of Element 8
    constexpr auto cohesion                = 900.0;
    constexpr auto friction_angle          = 27.457076095938259;
    auto           strain_vector_finalized = Vector(4, 0.0);
    auto           stress_vector_finalized = UblasUtilities::CreateVector(
        {-5.1687704591168895e+03, -1.2121212099273185e+04, -5.1687704591168895e+03, 0.0000000000000000e+00});
    auto strain_vector = UblasUtilities::CreateVector({2.0849072892302034e-04, -5.2137806109330782e-05,
                                                       0.0000000000000000e+00, 1.3759667604308981e-08});
    // Act
    const auto actual_cauchy_stress_vector = ComputeStressVectorUsingCPhiReductionTestData(
        cohesion, friction_angle, strain_vector_finalized, stress_vector_finalized, strain_vector);
    // Assert
    const auto expected_cauchy_stress_vector = UblasUtilities::CreateVector(
        {-2.1048640257824104e+03, -8.6704134153705982e+03, -3.8658294356694760e+03, 8.3845724538103905e-02});
    KRATOS_EXPECT_VECTOR_NEAR(actual_cauchy_stress_vector, expected_cauchy_stress_vector,
                              Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCornerWithTensionCutoff_1,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange, data are taken from C-Phi reduction test for GPoint 1 of Element 8
    constexpr auto cohesion                = 900.0;
    constexpr auto friction_angle          = 27.457076095938259;
    auto           strain_vector_finalized = Vector(4, 0.0);
    auto           stress_vector_finalized = UblasUtilities::CreateVector(
        {-1.2921926147786228e+03, -3.0303030248168902e+03, -1.2921926147786228e+03, 0.0000000000000000e+00});
    auto strain_vector = UblasUtilities::CreateVector({2.1223519849231651e-04, -5.3061831008594524e-05,
                                                       0.0000000000000000e+00, 8.1947498012640424e-09});
    // Act
    const auto actual_cauchy_stress_vector = ComputeStressVectorUsingCPhiReductionTestData(
        cohesion, friction_angle, strain_vector_finalized, stress_vector_finalized, strain_vector);
    // Assert
    const auto expected_cauchy_stress_vector = UblasUtilities::CreateVector(
        {999.999999812408, -252.651138668248, -252.65113885584, 0.0153293087909138});
    KRATOS_EXPECT_VECTOR_NEAR(actual_cauchy_stress_vector, expected_cauchy_stress_vector,
                              Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCornerWithTensionCutoff_2,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange, data are taken from C-Phi reduction test for GPoint 2 of Element 8
    constexpr auto cohesion                = 900.0;
    constexpr auto friction_angle          = 27.457076095938259;
    auto           strain_vector_finalized = Vector(4, 0.0);
    auto           stress_vector_finalized = UblasUtilities::CreateVector(
        {-1.2921926147789015e+03, -3.0303030248175442e+03, -1.2921926147789015e+03, 0.0000000000000000e+00});
    auto strain_vector = UblasUtilities::CreateVector({2.1221334908075760e-04, -5.3056189639024483e-05,
                                                       0.0000000000000000e+00, 3.4835553978606573e-08});
    // Act
    const auto actual_cauchy_stress_vector = ComputeStressVectorUsingCPhiReductionTestData(
        cohesion, friction_angle, strain_vector_finalized, stress_vector_finalized, strain_vector);
    // Assert
    const auto expected_cauchy_stress_vector = UblasUtilities::CreateVector(
        {999.999996609523, -252.651135465368, -252.651138855841, 0.065169629648695});
    KRATOS_EXPECT_VECTOR_NEAR(actual_cauchy_stress_vector, expected_cauchy_stress_vector,
                              Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCompressionCapZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);
    properties.SetValue(GEO_COHESION, 5.0);
    properties.SetValue(YOUNG_MODULUS, 1.0);
    properties.SetValue(POISSON_RATIO, 0.0);
    properties.SetValue(GEO_ENABLE_COMPRESSION_CAP, true);
    properties.SetValue(GEO_COMPRESSION_CAP_SIZE, 1.0);
    properties.SetValue(GEO_PRECONSOLIDATION_STRESS, 20.0);

    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);

    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    auto cauchy_stress_vector =
        UblasUtilities::CreateVector({-19.834494905518027796261932699488, -26.30497103308104685498440364009,
                                      -26.30497103308104685498440364009, 0.0});
    auto expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({-15.867595924414422237, 0.0, -11.338673315592010089, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCapCornerZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);
    properties.SetValue(GEO_COHESION, 5.0);
    properties.SetValue(YOUNG_MODULUS, 1.0e6);
    properties.SetValue(POISSON_RATIO, 0.0);
    properties.SetValue(GEO_ENABLE_COMPRESSION_CAP, true);
    properties.SetValue(GEO_COMPRESSION_CAP_SIZE, 1.0);
    properties.SetValue(GEO_PRECONSOLIDATION_STRESS, 20.0);

    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);

    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    // auto cauchy_stress_vector =
    //    UblasUtilities::CreateVector({8, -11.4501655647293, -26.5498344352707, 0.0});
    auto cauchy_stress_vector = UblasUtilities::CreateVector({8.5, -16, -30, 0.0});
    auto expected_cauchy_stress_vector =
        UblasUtilities::CreateVector({-15.867595924414422237, 0.0, -11.338673315592010089, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
