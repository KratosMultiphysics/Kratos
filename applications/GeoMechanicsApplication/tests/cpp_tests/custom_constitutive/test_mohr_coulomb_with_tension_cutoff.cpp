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
//

#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_constitutive/mohr_coulomb_with_tension_cutoff.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include "utilities/math_utils.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <sstream>
#include <string>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

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
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_COHESION is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_COHESION, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_FRICTION_ANGLE is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_FRICTION_ANGLE, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_DILATANCY_ANGLE is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_DILATANCY_ANGLE, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_TENSILE_STRENGTH is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_TENSILE_STRENGTH, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: YOUNG_MODULUS is not defined or has an invalid value for property: 3")
    properties.SetValue(YOUNG_MODULUS, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: POISSON_RATIO is not defined or has an invalid value for property: 3")
    properties.SetValue(POISSON_RATIO, 1.0);
    KRATOS_EXPECT_EQ(law.Check(properties, element_geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtElasticZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    ConstitutiveLaw::Parameters parameters;
    Properties                  properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    parameters.SetMaterialProperties(properties);
    ProcessInfo process;

    Geometry<Node> dummyGeometry;
    Vector         dummyVector;
    law.InitializeMaterial(properties, dummyGeometry, dummyVector);

    // Act
    Vector cauchy_stress_vector = ZeroVector(4);
    cauchy_stress_vector <<= 6.0, 0.0, -10.0, 0.0;
    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    Vector mapped_stress_vector;
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 6.0, 0.0, -10.0, 0.0;
    KRATOS_EXPECT_VECTOR_EQ(mapped_stress_vector, expected_cauchy_stress_vector);

    // Act
    cauchy_stress_vector <<= 8.0, 6.0, 4.0, 0.0;
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    expected_cauchy_stress_vector <<= 8.0, 6.0, 4.0, 0.0;
    KRATOS_EXPECT_VECTOR_EQ(mapped_stress_vector, expected_cauchy_stress_vector);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    ConstitutiveLaw::Parameters parameters;
    Properties                  properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    parameters.SetMaterialProperties(properties);
    ProcessInfo process;

    Geometry<Node> dummyGeometry;
    Vector         dummyVector;
    law.InitializeMaterial(properties, dummyGeometry, dummyVector);

    // Act
    Vector cauchy_stress_vector = ZeroVector(4);
    cauchy_stress_vector <<= 8.0, 0.0, -12.0, 0.0;
    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    Vector mapped_stress_vector;
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 7.338673315592010089, 0.0, -11.338673315592010089, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    // Act
    cauchy_stress_vector <<= 12.0, 0.0, -16.0, 0.0;
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    expected_cauchy_stress_vector <<= 7.338673315592010089, 0.0, -11.338673315592010089, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    // Act
    cauchy_stress_vector <<= 12.0, 10.0, -16.0, 0.0;
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    expected_cauchy_stress_vector <<= 7.338673315592010089, 7.90609951999166506, -9.24477283558367515, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCornerReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    ConstitutiveLaw::Parameters parameters;
    Properties                  properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    parameters.SetMaterialProperties(properties);
    ProcessInfo process;

    Geometry<Node> dummyGeometry;
    Vector         dummyVector;
    law.InitializeMaterial(properties, dummyGeometry, dummyVector);

    // Act
    Vector cauchy_stress_vector = ZeroVector(4);
    cauchy_stress_vector <<= 14.0, 8.0, 2.0, 0.0;
    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    Vector mapped_stress_vector;
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 10.0, 8.0, -1.5179192179966735, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    // Act
    cauchy_stress_vector <<= 24.0, 8.0, -8.0, 0.0;
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    expected_cauchy_stress_vector <<= 10.0, 8.0, -1.5179192179966735, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    // Assert
    expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 10.0, 8.0, -1.5179192179966735, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    // Act
    cauchy_stress_vector <<= 24.0, 22.0, -8.0, 0.0;
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    expected_cauchy_stress_vector <<= 10.0, 10.0, -1.5179192179966735, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtAxialTensionZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    ConstitutiveLaw::Parameters parameters;
    Properties                  properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    parameters.SetMaterialProperties(properties);
    ProcessInfo process;

    Geometry<Node> dummyGeometry;
    Vector         dummyVector;
    law.InitializeMaterial(properties, dummyGeometry, dummyVector);

    // Act
    Vector cauchy_stress_vector = ZeroVector(4);
    cauchy_stress_vector <<= 12.0, 9.0, 8.0, 0.0;
    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    Vector mapped_stress_vector;
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 10.0, 9.0, 6.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    // Act
    cauchy_stress_vector <<= 14.0, 9.0, 6.0, 0.0;
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    expected_cauchy_stress_vector <<= 10.0, 9.0, 2.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);

    // Act
    cauchy_stress_vector <<= 14.0, 12.0, 6.0, 0.0;
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    expected_cauchy_stress_vector <<= 10.0, 10.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyWithLargeTensialStrength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = MohrCoulombWithTensionCutOff(std::make_unique<PlaneStrain>());

    ConstitutiveLaw::Parameters parameters;
    Properties                  properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 20.0);
    parameters.SetMaterialProperties(properties);
    ProcessInfo process;

    Geometry<Node> dummyGeometry;
    Vector         dummyVector;
    law.InitializeMaterial(properties, dummyGeometry, dummyVector);

    // Act
    Vector cauchy_stress_vector = ZeroVector(4);
    cauchy_stress_vector <<= 22.0, 20.0, 18.0, 0.0;
    Vector strain_vector = ZeroVector(4);
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(cauchy_stress_vector);
    law.SetValue(CAUCHY_STRESS_VECTOR, cauchy_stress_vector, process);
    law.FinalizeMaterialResponseCauchy(parameters);
    law.CalculateMaterialResponseCauchy(parameters);
    Vector mapped_stress_vector;
    law.GetValue(CAUCHY_STRESS_VECTOR, mapped_stress_vector);

    // Assert
    Vector expected_cauchy_stress_vector(4);
    expected_cauchy_stress_vector <<= 14.2814800674211450, 14.2814800674211450, 14.2814800674211450, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(mapped_stress_vector, expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing