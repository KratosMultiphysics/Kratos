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
#include "custom_utilities/registration_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

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
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: GEO_COHESION is not defined for property 3")
    properties.SetValue(GEO_COHESION, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: value of GEO_COHESION for property 3 is out of range: -1 is not in [0.0, ->")
    properties.SetValue(GEO_COHESION, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: GEO_FRICTION_ANGLE is not defined for property 3")
    properties.SetValue(GEO_FRICTION_ANGLE, -30.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: value of GEO_FRICTION_ANGLE for property 3 is out of range: -30 is not in [0.0, ->")
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: GEO_DILATANCY_ANGLE is not defined for property 3")
    properties.SetValue(GEO_DILATANCY_ANGLE, -30.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: value of GEO_DILATANCY_ANGLE for property 3 is out "
        "of range: -30 is not in [0.0, 30.000000]")
    properties.SetValue(GEO_DILATANCY_ANGLE, 40.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: value of GEO_DILATANCY_ANGLE for property 3 is out "
        "of range: 40 is not in [0.0, 30.000000]")
    properties.SetValue(GEO_DILATANCY_ANGLE, 30.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: GEO_TENSILE_STRENGTH is not defined for property 3")
    properties.SetValue(GEO_TENSILE_STRENGTH, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: value of GEO_TENSILE_STRENGTH for property 3 is out "
        "of range: -1 is not in [0.0, 1.732051]")
    properties.SetValue(GEO_TENSILE_STRENGTH, 2.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: value of GEO_TENSILE_STRENGTH for property 3 is out "
        "of range: 2 is not in [0.0, 1.732051]")
    properties.SetValue(GEO_TENSILE_STRENGTH, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "YOUNG_MODULUS is not defined for property 3")
    properties.SetValue(YOUNG_MODULUS, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: value of YOUNG_MODULUS for property 3 is out of range: -1 is not in [0.0, ->")
    properties.SetValue(YOUNG_MODULUS, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "POISSON_RATIO is not defined for property 3")
    properties.SetValue(POISSON_RATIO, -0.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: value of POISSON_RATIO for property 3 is out of "
        "range: -0.5 is not in [0.0, 0.500000]")
    properties.SetValue(POISSON_RATIO, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: value of POISSON_RATIO for property 3 is out of "
        "range: 1 is not in [0.0, 0.500000]")
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
    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 6.0, 0.0, -10.0, 0.0;
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 6.0, 0.0, -10.0, 0.0;
    KRATOS_EXPECT_VECTOR_EQ(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                            expected_cauchy_stress_vector);

    cauchy_stress_vector <<= 8.0, 6.0, 4.0, 0.0;
    expected_cauchy_stress_vector <<= 8.0, 6.0, 4.0, 0.0;
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
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 8.0, 0.0, -12.0, 0.0;
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 7.338673315592010089, 0.0, -11.338673315592010089, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector <<= 12.0, 0.0, -16.0, 0.0;
    expected_cauchy_stress_vector <<= 7.338673315592010089, 0.0, -11.338673315592010089, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector <<= 12.0, 10.0, -16.0, 0.0;
    expected_cauchy_stress_vector <<= 7.338673315592010089, 7.90609951999166506, -9.24477283558367515, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
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
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 18.0, 8.0, -2.0, 0.0;
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 10.0, 8.0, -1.5179192179966735, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector <<= 24.0, 8.0, -8.0, 0.0;
    expected_cauchy_stress_vector <<= 10.0, 8.0, -1.5179192179966735, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector <<= 24.0, 22.0, -8.0, 0.0;
    expected_cauchy_stress_vector <<= 10.0, 10.0, -1.5179192179966735, 0.0;
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
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 12.0, 9.0, 8.0, 0.0;
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 10.0, 9.0, 8.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector <<= 14.0, 9.0, 6.0, 0.0;
    expected_cauchy_stress_vector <<= 10.0, 9.0, 6.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector <<= 14.0, 12.0, 6.0, 0.0;
    expected_cauchy_stress_vector <<= 10.0, 10.0, 6.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
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
    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(properties);
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act and Assert
    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 19.0, 12.0, 11.0, 0.0;
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 10.0, 10.0, 10.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    cauchy_stress_vector <<= 11.5, 10.0, 10.5, 0.0;
    expected_cauchy_stress_vector <<= 10.0, 10.0, 10.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyWithLargeTensileStrength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());
    Properties properties;
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
    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 22.0, 20.0, 18.0, 0.0;
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 14.2814800674211450, 14.2814800674211450, 14.2814800674211450, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedStressVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
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

    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 6.0, 0.0, -10.0, 0.0;
    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    const auto dummy_process_info = ProcessInfo{};
    p_law->SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, dummy_process_info);
    p_law->InitializeMaterialResponseCauchy(parameters);
    p_law->CalculateMaterialResponseCauchy(parameters);

    Vector calculated_cauchy_stress_vector;
    p_law->GetValue(CAUCHY_STRESS_VECTOR, calculated_cauchy_stress_vector);

    const auto scoped_registration_dimension = ScopedSerializerRegistration{"PlaneStrain"s, PlaneStrain{}};
    const auto scoped_registration_law =
        ScopedSerializerRegistration{"MohrCoulombWithTensionCutOff"s, MohrCoulombWithTensionCutOff{}};
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
    Matrix expected_constitutive_matrix(4, 4);
    // clang-format off
    expected_constitutive_matrix <<= 1.35E8, 5.77E7, 5.77E7, 0.,
                                     5.77E7, 1.35E8, 5.77E7, 0.,
                                     5.77E7, 5.77E7, 1.35E8, 0.,
                                     0.,     0.,     0.,     3.85E7;
    // clang-format on
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

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_GetLawFeatures, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto mc_law_plane_strain = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    // Act
    ConstitutiveLaw::Features features;
    mc_law_plane_strain.GetLawFeatures(features);

    // Assert
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::INFINITESIMAL_STRAINS))
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::ISOTROPIC))
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::PLANE_STRAIN_LAW))
    KRATOS_EXPECT_EQ(features.GetStrainMeasures().front(), ConstitutiveLaw::StrainMeasure_Infinitesimal);
    KRATOS_EXPECT_EQ(features.GetStrainSize(), 4);
    KRATOS_EXPECT_EQ(features.GetSpaceDimension(), 2);

    // Arrange
    auto mc_law_3d = MohrCoulombWithTensionCutOff(std::make_unique<ThreeDimensional>());

    // Act
    mc_law_3d.GetLawFeatures(features);

    // Assert
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::INFINITESIMAL_STRAINS))
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::ISOTROPIC))
    KRATOS_EXPECT_TRUE(features.GetOptions().Is(ConstitutiveLaw::THREE_DIMENSIONAL_LAW))
    KRATOS_EXPECT_EQ(features.GetStrainMeasures().front(), ConstitutiveLaw::StrainMeasure_Infinitesimal);
    KRATOS_EXPECT_EQ(features.GetStrainSize(), 6);
    KRATOS_EXPECT_EQ(features.GetSpaceDimension(), 3);
}

} // namespace Kratos::Testing