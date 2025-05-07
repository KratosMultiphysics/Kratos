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
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;
using namespace std::string_literals;

namespace
{

using namespace Kratos;

void InitializeLawMaterial(ConstitutiveLaw& rLaw, const Properties& rProperties)
{
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    rLaw.InitializeMaterial(rProperties, dummy_element_geometry, dummy_shape_function_values);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_Clone, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto original_law = InterfaceCoulombWithTensionCutOff{};

    // Act
    auto p_cloned_law = original_law.Clone();

    // Assert
    KRATOS_EXPECT_NE(p_cloned_law.get(), nullptr);
    KRATOS_EXPECT_NE(p_cloned_law.get(), &original_law);
    KRATOS_EXPECT_NE(dynamic_cast<const InterfaceCoulombWithTensionCutOff*>(p_cloned_law.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtElasticZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: set elastic compressive state
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto p_initial_state         = make_intrusive<InitialState>();
    auto initial_traction_vector = Vector{2};
    initial_traction_vector <<= -6.0, 8.0;
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombWithTensionCutOff{};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    auto expected_traction_vector = Vector{2};
    expected_traction_vector <<= -6.0, 8.0;
    KRATOS_EXPECT_VECTOR_NEAR(parameters.GetStressVector(), expected_traction_vector, Defaults::absolute_tolerance);

    // Arrange: set elastic tensile state
    traction_vector <<= 5.0, 4.0;
    const auto dummy_process_info = ProcessInfo{};
    law.SetValue(CAUCHY_STRESS_VECTOR, traction_vector, dummy_process_info);
    law.FinalizeMaterialResponseCauchy(parameters);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    expected_traction_vector <<= 5.0, 4.0;
    KRATOS_EXPECT_VECTOR_NEAR(parameters.GetStressVector(), expected_traction_vector, Defaults::absolute_tolerance);

    // Arrange: set elastic tensile state (reverse shear)
    traction_vector <<= 5.0, -2.0;
    law.SetValue(CAUCHY_STRESS_VECTOR, traction_vector, dummy_process_info);
    law.FinalizeMaterialResponseCauchy(parameters);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    expected_traction_vector <<= 5.0, -2.0;
    KRATOS_EXPECT_VECTOR_NEAR(parameters.GetStressVector(), expected_traction_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionApexReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto p_initial_state         = make_intrusive<InitialState>();
    auto initial_traction_vector = Vector{2};
    initial_traction_vector <<= 15.0, 2.0;
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombWithTensionCutOff{};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    Vector expected_cauchy_stress_vector(2);
    expected_cauchy_stress_vector <<= 10.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(parameters.GetStressVector(), expected_cauchy_stress_vector,
                              Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionCutoffReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto p_initial_state         = make_intrusive<InitialState>();
    auto initial_traction_vector = Vector{2};
    initial_traction_vector <<= 10.0, 4.0;
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombWithTensionCutOff{};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    Vector expected_cauchy_stress_vector(2);
    expected_cauchy_stress_vector <<= 8.0, 2.0;
    KRATOS_EXPECT_VECTOR_NEAR(parameters.GetStressVector(), expected_cauchy_stress_vector,
                              Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCornerReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto p_initial_state         = make_intrusive<InitialState>();
    auto initial_traction_vector = Vector{2};
    initial_traction_vector <<= 10.0, 14.0;
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombWithTensionCutOff{};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    Vector expected_cauchy_stress_vector(2);
    expected_cauchy_stress_vector <<= 4.2410403910016600, 5.7589596089983400;
    KRATOS_EXPECT_VECTOR_NEAR(parameters.GetStressVector(), expected_cauchy_stress_vector,
                              Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto p_initial_state         = make_intrusive<InitialState>();
    auto initial_traction_vector = Vector{2};
    initial_traction_vector <<= -5.0, 18.0;
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombWithTensionCutOff{};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    Vector expected_cauchy_stress_vector(2);
    expected_cauchy_stress_vector <<= -5.0, 11.059402624645148377;
    KRATOS_EXPECT_VECTOR_NEAR(parameters.GetStressVector(), expected_cauchy_stress_vector,
                              Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_Serialization, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto p_initial_state         = make_intrusive<InitialState>();
    auto initial_traction_vector = Vector{2};
    initial_traction_vector <<= -2.0, 8.0;
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto p_law = std::unique_ptr<ConstitutiveLaw>{std::make_unique<InterfaceCoulombWithTensionCutOff>()};
    p_law->SetInitialState(p_initial_state);
    InitializeLawMaterial(*p_law, properties);

    p_law->CalculateMaterialResponseCauchy(parameters);
    const auto calculated_traction_vector = parameters.GetStressVector();

    const auto scoped_registration_law = ScopedSerializerRegistration{
        "InterfaceCoulombWithTensionCutOff"s, InterfaceCoulombWithTensionCutOff{}};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_law);
    auto p_loaded_law = std::unique_ptr<ConstitutiveLaw>();
    serializer.load("test_tag"s, p_loaded_law);

    ASSERT_NE(p_loaded_law.get(), nullptr);

    auto loaded_calculated_traction_vector = Vector{};
    p_loaded_law->GetValue(CAUCHY_STRESS_VECTOR, loaded_calculated_traction_vector);
    KRATOS_EXPECT_VECTOR_EQ(loaded_calculated_traction_vector, calculated_traction_vector);

    // Check whether the finalized traction and relative displacement have been restored properly
    p_loaded_law->CalculateMaterialResponseCauchy(parameters);
    KRATOS_EXPECT_VECTOR_EQ(parameters.GetStressVector(), calculated_traction_vector);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_WorkingSpaceDimensionIsAlways2D,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(InterfaceCoulombWithTensionCutOff{}.WorkingSpaceDimension(), N_DIM_2D);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_StressMeasureIsAlwaysCauchy,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(InterfaceCoulombWithTensionCutOff{}.GetStressMeasure(), ConstitutiveLaw::StressMeasure_Cauchy);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_StrainSizeIsAlways2, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(InterfaceCoulombWithTensionCutOff{}.GetStrainSize(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_StrainMeasureIsAlwaysInfinitesimal,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(InterfaceCoulombWithTensionCutOff{}.GetStrainMeasure(),
                     ConstitutiveLaw::StrainMeasure_Infinitesimal);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_HasIncrementalFormulation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    EXPECT_TRUE(InterfaceCoulombWithTensionCutOff{}.IsIncremental());
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_RequiresInitializeMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    EXPECT_TRUE(InterfaceCoulombWithTensionCutOff{}.RequiresInitializeMaterialResponse());
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_Check, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto                        law = InterfaceCoulombWithTensionCutOff{};
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
        "Error: INTERFACE_NORMAL_STIFFNESS is not defined for property 3")
    properties.SetValue(INTERFACE_NORMAL_STIFFNESS, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "Error: value of INTERFACE_NORMAL_STIFFNESS for property 3 is out of range: -1 is not in [0.0, ->")
    properties.SetValue(INTERFACE_NORMAL_STIFFNESS, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "Error: INTERFACE_SHEAR_STIFFNESS is not defined for property 3")
    properties.SetValue(INTERFACE_SHEAR_STIFFNESS, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "Error: value of INTERFACE_SHEAR_STIFFNESS for property 3 is out of range: -1 is not in [0.0, ->")
    properties.SetValue(INTERFACE_SHEAR_STIFFNESS, 1.0);

    KRATOS_EXPECT_EQ(law.Check(properties, element_geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateConstitutiveMatrix,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(INTERFACE_NORMAL_STIFFNESS, 1.0E8);
    properties.SetValue(INTERFACE_SHEAR_STIFFNESS, 5.0E7);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);

    auto law = InterfaceCoulombWithTensionCutOff{};
    InitializeLawMaterial(law, properties);

    // Act
    auto constitutive_matrix = Matrix{ZeroMatrix{2, 2}};
    law.CalculateValue(parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);

    // Assert
    auto expected_constitutive_matrix = Matrix{2, 2};
    // clang-format off
    expected_constitutive_matrix <<= 1.0E8, 0.0,
                                     0.0,   5.0E7;
    // clang-format on
    KRATOS_EXPECT_MATRIX_NEAR(constitutive_matrix, expected_constitutive_matrix, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing