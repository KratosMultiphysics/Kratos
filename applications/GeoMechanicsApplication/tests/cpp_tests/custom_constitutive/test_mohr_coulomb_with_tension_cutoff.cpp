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
#include "geo_mechanics_application_variables.h"

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

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

} // namespace Kratos::Testing