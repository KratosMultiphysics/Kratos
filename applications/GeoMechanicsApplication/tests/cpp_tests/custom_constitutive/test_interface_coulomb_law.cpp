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

#include "custom_constitutive/interface_coulomb_law.h"
#include "custom_constitutive/interface_plane_strain.h"
#include "custom_utilities/registration_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

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

// Integrates the interface Coulomb law along a sequence of cumulative relative-displacement
// states and returns the final traction vector.
Vector IntegrateInterfaceCoulombDisplacementPath(const Properties& rProperties,
                                                 const std::vector<Vector>& rCumulativeRelativeDisplacementStates)
{
    auto law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};
    InitializeLawMaterial(law, rProperties);

    ConstitutiveLaw::Parameters parameters;
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    parameters.SetMaterialProperties(rProperties);

    Vector zero_state = ZeroVector(2);
    parameters.SetStrainVector(zero_state);
    parameters.SetStressVector(zero_state);
    law.InitializeMaterialResponseCauchy(parameters);

    for (const auto& r_relative_displacement_state : rCumulativeRelativeDisplacementStates) {
        Vector relative_displacement = r_relative_displacement_state;
        parameters.SetStrainVector(relative_displacement);
        law.CalculateMaterialResponseCauchy(parameters);
        law.FinalizeMaterialResponseCauchy(parameters);
    }

    Vector result;
    law.GetValue(GEO_EFFECTIVE_TRACTION_VECTOR, result);
    return result;
}

// Builds a list of cumulative relative-displacement states that follow the given vertices.
// Between each pair of consecutive vertices, NIncrementsPerLeg equally spaced points are added.
std::vector<Vector> BuildLinearlyInterpolatedDisplacementPath(const std::vector<Vector>& rPathVertices,
                                                              const std::size_t NIncrementsPerLeg)
{
    KRATOS_ERROR_IF(NIncrementsPerLeg == 0)
        << "BuildLinearlyInterpolatedDisplacementPath requires NIncrementsPerLeg > 0" << std::endl;

    std::vector<Vector> relative_displacement_states;
    for (std::size_t leg = 0; leg + 1 < rPathVertices.size(); ++leg) {
        const Vector& r_start = rPathVertices[leg];
        const Vector& r_end   = rPathVertices[leg + 1];
        for (std::size_t increment = 1; increment <= NIncrementsPerLeg; ++increment) {
            const auto fraction = static_cast<double>(increment) / static_cast<double>(NIncrementsPerLeg);
            relative_displacement_states.emplace_back(r_start + fraction * (r_end - r_start));
        }
    }
    return relative_displacement_states;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulomb_Clone, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto original_law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};

    // Act
    auto p_cloned_law = original_law.Clone();

    // Assert
    KRATOS_EXPECT_NE(p_cloned_law.get(), nullptr);
    KRATOS_EXPECT_NE(p_cloned_law.get(), &original_law);
    KRATOS_EXPECT_NE(dynamic_cast<const InterfaceCoulombLaw*>(p_cloned_law.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtElasticZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: set elastic compressive state
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_ENABLE_TENSION_CUT_OFF, true);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto       p_initial_state         = make_intrusive<InitialState>();
    const auto initial_traction_vector = UblasUtilities::CreateVector({-6.0, 8.0});
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);
    auto resulting_traction = parameters.GetStressVector();
    int  plasticity_status;
    law.GetValue(GEO_PLASTICITY_STATUS, plasticity_status);

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(resulting_traction, initial_traction_vector, Defaults::absolute_tolerance);
    KRATOS_EXPECT_EQ(plasticity_status, static_cast<int>(PlasticityStatus::ELASTIC));

    // Arrange: set elastic tensile state
    traction_vector = UblasUtilities::CreateVector({5.0, 4.0});
    law.SetValue(GEO_EFFECTIVE_TRACTION_VECTOR, traction_vector, ProcessInfo{});
    law.FinalizeMaterialResponseCauchy(parameters);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);
    resulting_traction = parameters.GetStressVector();
    law.GetValue(GEO_PLASTICITY_STATUS, plasticity_status);

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(resulting_traction, traction_vector, Defaults::absolute_tolerance);
    KRATOS_EXPECT_EQ(plasticity_status, static_cast<int>(PlasticityStatus::ELASTIC));

    // Arrange: set elastic tensile state (reverse shear)
    traction_vector = UblasUtilities::CreateVector({5.0, -2.0});
    law.SetValue(GEO_EFFECTIVE_TRACTION_VECTOR, traction_vector, ProcessInfo{});
    law.FinalizeMaterialResponseCauchy(parameters);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);
    resulting_traction = parameters.GetStressVector();
    law.GetValue(GEO_PLASTICITY_STATUS, plasticity_status);

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(resulting_traction, traction_vector, Defaults::absolute_tolerance);
    KRATOS_EXPECT_EQ(plasticity_status, static_cast<int>(PlasticityStatus::ELASTIC));
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionApexReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_ENABLE_TENSION_CUT_OFF, true);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto       p_initial_state         = make_intrusive<InitialState>();
    const auto initial_traction_vector = UblasUtilities::CreateVector({15.0, 2.0});
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);
    const auto& r_resulting_traction = parameters.GetStressVector();
    int         plasticity_status;
    law.GetValue(GEO_PLASTICITY_STATUS, plasticity_status);

    // Assert
    const auto expected_traction = UblasUtilities::CreateVector({10.0, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(r_resulting_traction, expected_traction, Defaults::absolute_tolerance);
    KRATOS_EXPECT_EQ(plasticity_status, static_cast<int>(PlasticityStatus::TENSION_APEX));
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtTensionCutoffReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_ENABLE_TENSION_CUT_OFF, true);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(INTERFACE_NORMAL_STIFFNESS, 25.0);
    properties.SetValue(INTERFACE_SHEAR_STIFFNESS, 12.5);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto       p_initial_state         = make_intrusive<InitialState>();
    const auto initial_traction_vector = UblasUtilities::CreateVector({10.0, 4.0});
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);
    const auto& r_resulting_traction = parameters.GetStressVector();
    int         plasticity_status;
    law.GetValue(GEO_PLASTICITY_STATUS, plasticity_status);

    // Assert
    const auto expected_traction = UblasUtilities::CreateVector({22.0 / 3.0, 8.0 / 3.0});
    KRATOS_EXPECT_VECTOR_NEAR(r_resulting_traction, expected_traction, Defaults::absolute_tolerance);
    KRATOS_EXPECT_EQ(plasticity_status, static_cast<int>(PlasticityStatus::TENSION_CUT_OFF));
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtCornerReturnZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_ENABLE_TENSION_CUT_OFF, true);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto       p_initial_state         = make_intrusive<InitialState>();
    const auto initial_traction_vector = UblasUtilities::CreateVector({10.0, 14.0});
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);
    const auto& r_resulting_traction = parameters.GetStressVector();
    int         plasticity_status;
    law.GetValue(GEO_PLASTICITY_STATUS, plasticity_status);

    // Assert
    const auto expected_traction = UblasUtilities::CreateVector({4.2410403910016600, 5.7589596089983400});
    KRATOS_EXPECT_VECTOR_NEAR(r_resulting_traction, expected_traction, Defaults::absolute_tolerance);
    KRATOS_EXPECT_EQ(plasticity_status, static_cast<int>(PlasticityStatus::TENSION_MOHR_COULOMB_CORNER));
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateMaterialResponseCauchyAtRegularFailureZone,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    properties.SetValue(GEO_ENABLE_TENSION_CUT_OFF, true);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(INTERFACE_NORMAL_STIFFNESS, 25.0);
    properties.SetValue(INTERFACE_SHEAR_STIFFNESS, 12.5);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto       p_initial_state         = make_intrusive<InitialState>();
    const auto initial_traction_vector = UblasUtilities::CreateVector({-5.0, 18.0});
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};
    law.SetInitialState(p_initial_state);
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);
    const auto& r_resulting_traction = parameters.GetStressVector();
    int         plasticity_status;
    law.GetValue(GEO_PLASTICITY_STATUS, plasticity_status);

    // Assert
    const auto expected_traction = UblasUtilities::CreateVector({-5.0, 11.059402624645148377});
    KRATOS_EXPECT_VECTOR_NEAR(r_resulting_traction, expected_traction, Defaults::absolute_tolerance);
    KRATOS_EXPECT_EQ(plasticity_status, static_cast<int>(PlasticityStatus::MOHR_COULOMB_FAILURE));
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_DoesNotCalculateStressWhenComputeStressOptionIsNotSet,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    properties.SetValue(GEO_ENABLE_TENSION_CUT_OFF, true);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties.SetValue(INTERFACE_NORMAL_STIFFNESS, 25.0);
    properties.SetValue(INTERFACE_SHEAR_STIFFNESS, 12.5);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{2, 1.0};
    parameters.SetStrainVector(relative_displacement_vector);

    auto law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};
    InitializeLawMaterial(law, properties);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);
    const auto& r_resulting_traction = parameters.GetStressVector();
    int         plasticity_status;
    law.GetValue(GEO_PLASTICITY_STATUS, plasticity_status);

    // Assert
    const auto expected_traction = UblasUtilities::CreateVector({0.0, 0.0});
    KRATOS_EXPECT_VECTOR_NEAR(r_resulting_traction, expected_traction, Defaults::absolute_tolerance);
    KRATOS_EXPECT_EQ(plasticity_status, static_cast<int>(PlasticityStatus::ELASTIC));
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_Serialization, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("InterfacePlaneStrain"s, InterfacePlaneStrain{})};
    auto properties = Properties{};
    properties.SetValue(GEO_FRICTION_ANGLE, 35.0);
    properties.SetValue(GEO_COHESION, 10.0);
    properties.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties.SetValue(GEO_ENABLE_TENSION_CUT_OFF, true);
    properties.SetValue(GEO_TENSILE_STRENGTH, 10.0);

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto traction_vector = Vector{ZeroVector{2}};
    parameters.SetStressVector(traction_vector);
    auto relative_displacement_vector = Vector{ZeroVector{2}};
    parameters.SetStrainVector(relative_displacement_vector);

    auto       p_initial_state         = make_intrusive<InitialState>();
    const auto initial_traction_vector = UblasUtilities::CreateVector({-2.0, 8.0});
    p_initial_state->SetInitialStressVector(initial_traction_vector);
    p_initial_state->SetInitialStrainVector(Vector{ZeroVector{2}});

    auto p_law = std::unique_ptr<ConstitutiveLaw>{
        std::make_unique<InterfaceCoulombLaw>(std::make_unique<InterfacePlaneStrain>())};
    p_law->SetInitialState(p_initial_state);
    InitializeLawMaterial(*p_law, properties);

    p_law->CalculateMaterialResponseCauchy(parameters);
    const auto calculated_traction_vector = parameters.GetStressVector();

    const auto scoped_registration_law =
        ScopedSerializerRegistration{std::make_pair("InterfaceCoulombLaw"s, InterfaceCoulombLaw{})};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_law);
    auto p_loaded_law = std::unique_ptr<ConstitutiveLaw>();
    serializer.load("test_tag"s, p_loaded_law);

    ASSERT_NE(p_loaded_law.get(), nullptr);

    auto loaded_calculated_traction_vector = Vector{};
    p_loaded_law->GetValue(GEO_EFFECTIVE_TRACTION_VECTOR, loaded_calculated_traction_vector);
    KRATOS_EXPECT_VECTOR_EQ(loaded_calculated_traction_vector, calculated_traction_vector);

    // Check whether the finalized traction and relative displacement have been restored properly
    p_loaded_law->CalculateMaterialResponseCauchy(parameters);
    const auto& r_resulting_traction = parameters.GetStressVector();
    KRATOS_EXPECT_VECTOR_EQ(r_resulting_traction, calculated_traction_vector);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_Has2DWorkingSpaceDimension,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(
        InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()}.WorkingSpaceDimension(), N_DIM_2D);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_HasCauchyStressMeasure, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()}.GetStressMeasure(),
                     ConstitutiveLaw::StressMeasure_Cauchy);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_StrainSizeEqualsTwo, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()}.GetStrainSize(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_HasInfinitesimalStrainMeasure,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()}.GetStrainMeasure(),
                     ConstitutiveLaw::StrainMeasure_Infinitesimal);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_HasIncrementalFormulation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    EXPECT_TRUE(InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()}.IsIncremental());
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_RequiresInitializeMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    EXPECT_TRUE(InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()}.RequiresInitializeMaterialResponse());
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_Check, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto                        law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};
    ConstitutiveLaw::Parameters parameters;
    Properties                  properties(3);
    parameters.SetMaterialProperties(properties);
    const auto element_geometry = Geometry<Node>{};
    const auto process_info     = ProcessInfo{};

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "GEO_COHESION does not exist in the property with Id 3.")
    properties.SetValue(GEO_COHESION, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "GEO_COHESION in the property with Id 3 has an invalid value: -1 is out of the range [0, -).")
    properties.SetValue(GEO_COHESION, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "GEO_FRICTION_ANGLE does not exist in the property with Id 3.")
    properties.SetValue(GEO_FRICTION_ANGLE, -30.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "GEO_FRICTION_ANGLE in the property with Id 3 has an invalid value: -30 is out of the range (0, 90).")
    properties.SetValue(GEO_FRICTION_ANGLE, 30.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "GEO_DILATANCY_ANGLE does not exist in the property with Id 3.")
    properties.SetValue(GEO_DILATANCY_ANGLE, -30.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "GEO_DILATANCY_ANGLE in the property with Id 3 has an invalid value: -30 is out of the range [0, 30].")
    properties.SetValue(GEO_DILATANCY_ANGLE, 40.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "GEO_DILATANCY_ANGLE in the property with Id 3 has an invalid value: 40 is out of the range [0, 30].")
    properties.SetValue(GEO_DILATANCY_ANGLE, 30.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "INTERFACE_NORMAL_STIFFNESS does not exist in the property with Id 3.")
    properties.SetValue(INTERFACE_NORMAL_STIFFNESS, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "INTERFACE_NORMAL_STIFFNESS in the property with Id 3 has an invalid value: -1 is out of the range [0, -).")
    properties.SetValue(INTERFACE_NORMAL_STIFFNESS, 1.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info),
        "INTERFACE_SHEAR_STIFFNESS does not exist in the property with Id 3.")
    properties.SetValue(INTERFACE_SHEAR_STIFFNESS, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = law.Check(properties, element_geometry, process_info), "INTERFACE_SHEAR_STIFFNESS in the property with Id 3 has an invalid value: -1 is out of the range [0, -).")
    properties.SetValue(INTERFACE_SHEAR_STIFFNESS, 1.0);
    KRATOS_EXPECT_EQ(law.Check(properties, element_geometry, process_info), 0);

    properties.SetValue(GEO_ENABLE_TENSION_CUT_OFF, true);
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

    KRATOS_EXPECT_EQ(law.Check(properties, element_geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_CalculateConstitutiveMatrix,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties                        = Properties{};
    properties[INTERFACE_NORMAL_STIFFNESS] = 1.0E8;
    properties[INTERFACE_SHEAR_STIFFNESS]  = 5.0E7;
    properties[GEO_FRICTION_ANGLE]         = 0.0;
    properties[GEO_COHESION]               = 0.0;
    properties[GEO_DILATANCY_ANGLE]        = 0.0;

    auto parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);

    auto law = InterfaceCoulombLaw{std::make_unique<InterfacePlaneStrain>()};
    InitializeLawMaterial(law, properties);

    // Act
    auto constitutive_matrix = Matrix{ZeroMatrix{2, 2}};
    law.CalculateValue(parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);

    // Assert
    const auto expected_constitutive_matrix = UblasUtilities::CreateMatrix({{1.0E8, 0.0}, {0.0, 5.0E7}});
    KRATOS_EXPECT_MATRIX_NEAR(constitutive_matrix, expected_constitutive_matrix, Defaults::absolute_tolerance);
}

/// <summary>
///  This test checks if substepping improves the accuracy of the Interface coulomb integration
///  along a strain path that changes yield condition.
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulomb_SubStepping, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: perfect plasticity (no hardening) with non-associated dilatancy.
    Properties properties_no_sub_stepping;
    properties_no_sub_stepping.SetValue(GEO_COULOMB_HARDENING_TYPE, "None");
    properties_no_sub_stepping.SetValue(GEO_FRICTION_ANGLE, 30.0);
    properties_no_sub_stepping.SetValue(GEO_COHESION, 10.0);
    properties_no_sub_stepping.SetValue(GEO_DILATANCY_ANGLE, 20.0);
    properties_no_sub_stepping.SetValue(GEO_ENABLE_TENSION_CUT_OFF, true);
    properties_no_sub_stepping.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    properties_no_sub_stepping.SetValue(INTERFACE_NORMAL_STIFFNESS, 25.0);
    properties_no_sub_stepping.SetValue(INTERFACE_SHEAR_STIFFNESS, 12.5);

    Properties properties_with_sub_stepping(properties_no_sub_stepping);
    properties_with_sub_stepping.SetValue(GEO_MAX_RELATIVE_OVERSHOOT, 0.1);
    properties_with_sub_stepping.SetValue(GEO_MAX_NUMBER_OF_SUB_STEPS, 100);

    // Non-proportional relative-displacement history (normal, shear):
    //   part 1: add normal compression
    //   part 2: add shear and extension
    const auto zero_state        = Vector{ZeroVector(2)};
    const auto compression_state = UblasUtilities::CreateVector({-22.0, 0.0});
    const auto end_state         = UblasUtilities::CreateVector({2.0, 8.0});

    const std::vector<Vector> displacement_vertices{zero_state, compression_state, end_state};

    // Reference solution: the SAME path, integrated with a large number of steps.
    constexpr std::size_t reference_increments_per_leg = 2000;
    const auto            reference_traction           = IntegrateInterfaceCoulombDisplacementPath(
        properties_no_sub_stepping,
        BuildLinearlyInterpolatedDisplacementPath(displacement_vertices, reference_increments_per_leg));

    // Calculate the traction without sub-stepping, using only the vertices as the path.
    const auto traction_no_sub_stepping = IntegrateInterfaceCoulombDisplacementPath(
        properties_no_sub_stepping, BuildLinearlyInterpolatedDisplacementPath(displacement_vertices, 1));

    // Calculate the traction with sub-stepping, using only the vertices as the path.
    const auto traction_with_sub_stepping = IntegrateInterfaceCoulombDisplacementPath(
        properties_with_sub_stepping, BuildLinearlyInterpolatedDisplacementPath(displacement_vertices, 1));

    // Calculate the errors without and with sub-stepping, using the reference as the "exact" solution.
    const auto error_no_sub_step = norm_2(Vector{traction_no_sub_stepping - reference_traction});
    const auto error_with_sub_step = norm_2(Vector{traction_with_sub_stepping - reference_traction});

    // Check that the error without sub-stepping is at least 10x larger than the error with sub-stepping.
    KRATOS_EXPECT_GT(error_no_sub_step, error_with_sub_step * 10);
}

} // namespace Kratos::Testing
