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

#include "custom_constitutive/interface_coulomb_with_tension_cut_off.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_utilities/registration_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;
using namespace std::string_literals;

namespace
{

using namespace Kratos;

Vector CalculateMappedTractionVector(Vector&                                rTractionVector,
                                     ConstitutiveLaw::Parameters&           rParameters,
                                     InterfaceMohrCoulombWithTensionCutOff& rLaw)
{
    Vector strain_vector = ZeroVector(2);
    rParameters.SetStrainVector(strain_vector);
    rParameters.SetStressVector(rTractionVector);
    const auto dummy_process_info = ProcessInfo{};
    rLaw.SetValue(CAUCHY_STRESS_VECTOR, rTractionVector, dummy_process_info);
    rLaw.FinalizeMaterialResponseCauchy(rParameters);
    rLaw.CalculateMaterialResponseCauchy(rParameters);

    Vector result;
    rLaw.GetValue(CAUCHY_STRESS_VECTOR, result);
    return result;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceMohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtElasticZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: set elastic compressive state
    Properties properties;
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);

    auto p_initial_state         = make_intrusive<InitialState>();
    auto initial_traction_vector = Vector{2};
    initial_traction_vector <<= -6.0, 8.0;
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});
    auto law = InterfaceMohrCoulombWithTensionCutOff{};
    law.SetInitialState(p_initial_state);

    auto traction_output_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_output_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    law.InitializeMaterial(properties, dummy_element_geometry, dummy_shape_function_values);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    auto expected_traction_vector = Vector{2};
    expected_traction_vector <<= -6.0, 8.0;
    KRATOS_EXPECT_VECTOR_NEAR(traction_output_vector, expected_traction_vector, Defaults::absolute_tolerance);

    // Arrange: set elastic tensile state
    traction_output_vector <<= 5.0, 4.0;
    const auto dummy_process_info = ProcessInfo{};
    law.SetValue(CAUCHY_STRESS_VECTOR, traction_output_vector, dummy_process_info);
    law.FinalizeMaterialResponseCauchy(parameters);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    expected_traction_vector <<= 5.0, 4.0;
    KRATOS_EXPECT_VECTOR_NEAR(traction_output_vector, expected_traction_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceMohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionApexReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = InterfaceMohrCoulombWithTensionCutOff{};
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
    Vector cauchy_stress_vector(2);
    cauchy_stress_vector <<= 15.0, 2.0;
    Vector expected_cauchy_stress_vector(2);
    expected_cauchy_stress_vector <<= 10.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedTractionVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceMohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionCutoffReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = InterfaceMohrCoulombWithTensionCutOff{};
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
    Vector cauchy_stress_vector(2);
    cauchy_stress_vector <<= 10.0, 4.0;
    Vector expected_cauchy_stress_vector(2);
    expected_cauchy_stress_vector <<= 8.0, 2.0;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedTractionVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceMohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCornerReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = InterfaceMohrCoulombWithTensionCutOff{};
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
    Vector cauchy_stress_vector(2);
    cauchy_stress_vector <<= 10.0, 14.0;
    Vector expected_cauchy_stress_vector(2);
    expected_cauchy_stress_vector <<= 4.2410403910016600, 5.7589596089983400;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedTractionVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceMohrCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law = InterfaceMohrCoulombWithTensionCutOff{};
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
    Vector cauchy_stress_vector(2);
    cauchy_stress_vector <<= -5.0, 18.0;
    Vector expected_cauchy_stress_vector(2);
    expected_cauchy_stress_vector <<= -5.0, 11.059402624645148377;
    KRATOS_EXPECT_VECTOR_NEAR(CalculateMappedTractionVector(cauchy_stress_vector, parameters, law),
                              expected_cauchy_stress_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing