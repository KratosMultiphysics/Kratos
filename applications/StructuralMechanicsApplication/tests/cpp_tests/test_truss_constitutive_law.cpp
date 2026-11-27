// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Richard Faasse
//

#include "structural_mechanics_fast_suite.h"
#include "custom_constitutive/truss_constitutive_law.h"
#include "custom_constitutive/thermal_truss_constitutive_law.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TrussConstitutiveLaw_CalculatesLinearElasticStress, KratosStructuralMechanicsFastSuite) {
    TrussConstitutiveLaw law;
    ConstitutiveLaw::Parameters parameters;
    constexpr auto induced_strain = 0.005;
    Vector temp_strain = ScalarVector(1, induced_strain);
    Vector temp_stress = ZeroVector(1);
    parameters.SetStrainVector(temp_strain);
    parameters.SetStressVector(temp_stress);

    Properties properties;
    constexpr auto youngs_modulus = 1.0e6;
    properties[YOUNG_MODULUS] = youngs_modulus;
    parameters.SetMaterialProperties(properties);
    law.CalculateMaterialResponsePK2(parameters);

    constexpr double expected_stress = induced_strain * youngs_modulus;
    KRATOS_EXPECT_EQ(expected_stress, parameters.GetStressVector()[0]);
}

KRATOS_TEST_CASE_IN_SUITE(TrussConstitutiveLaw_CalculatesLinearElasticStress_WithInitialState, KratosStructuralMechanicsFastSuite) {
    TrussConstitutiveLaw law;
    constexpr auto induced_strain = 0.005;
    Vector temp_strain = ScalarVector(1, induced_strain);
    Vector temp_stress = ZeroVector(1);

    ConstitutiveLaw::Parameters parameters;
    parameters.SetStrainVector(temp_strain);
    parameters.SetStressVector(temp_stress);

    Properties properties;
    constexpr auto youngs_modulus = 1.0e6;
    properties[YOUNG_MODULUS] = youngs_modulus;
    parameters.SetMaterialProperties(properties);

    constexpr auto initial_stress = 5000.0;
    constexpr auto initial_strain = 0.01;
    const auto p_initial_state = Kratos::make_intrusive<InitialState>(
        ScalarVector(1, initial_strain), ScalarVector(1, initial_stress));
    law.SetInitialState(p_initial_state);

    law.CalculateMaterialResponsePK2(parameters);

    constexpr double expected_stress =
        (induced_strain - initial_strain) * youngs_modulus + initial_stress;
    KRATOS_EXPECT_EQ(expected_stress, parameters.GetStressVector()[0]);
}

KRATOS_TEST_CASE_IN_SUITE(ThermalTrussConstitutiveLaw_CalculatesLinearElasticStress, KratosStructuralMechanicsFastSuite) {
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    p_node_1->FastGetSolutionStepValue(TEMPERATURE) = 30.0;
    p_node_2->FastGetSolutionStepValue(TEMPERATURE) = 50.0;
    auto p_geometry = r_model_part.CreateNewGeometry("Line3D2", 1, std::vector<ModelPart::IndexType>{1, 2});

    Properties properties;
    constexpr auto youngs_modulus = 1.0e6;
    properties[YOUNG_MODULUS] = youngs_modulus;
    properties[THERMAL_EXPANSION_COEFFICIENT] = 1.0e-5;
    properties[REFERENCE_TEMPERATURE] = 20.0;

    constexpr auto induced_strain = 0.005;
    Vector strain = ScalarVector(1, induced_strain);
    Vector stress = ZeroVector(1);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetElementGeometry(*p_geometry);
    parameters.SetMaterialProperties(properties);
    parameters.SetStrainVector(strain);
    parameters.SetStressVector(stress);

    ThermalTrussConstitutiveLaw law;
    law.InitializeMaterial(properties, *p_geometry, Vector());
    law.CalculateMaterialResponsePK2(parameters);

    constexpr auto thermal_strain = 1.0e-5 * (40.0 - 20.0);
    constexpr auto expected_stress = youngs_modulus * (induced_strain - thermal_strain);
    KRATOS_EXPECT_NEAR(parameters.GetStressVector()[0], expected_stress, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(ThermalTrussConstitutiveLaw_RejectsMultipleIntegrationPoints, KratosStructuralMechanicsFastSuite) {
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_geometry = r_model_part.CreateNewGeometry("Line3D2", 1, std::vector<ModelPart::IndexType>{1, 2});

    Properties properties;
    properties[YOUNG_MODULUS] = 1.0e6;
    properties[DENSITY] = 0.0;
    properties[THERMAL_EXPANSION_COEFFICIENT] = 1.0e-5;
    properties[REFERENCE_TEMPERATURE] = 20.0;
    properties[INTEGRATION_ORDER] = 2;
    
    ThermalTrussConstitutiveLaw law;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, *p_geometry, r_model_part.GetProcessInfo()),
        "ThermalTrussConstitutiveLaw requires one integration point");
}

KRATOS_TEST_CASE_IN_SUITE(ThermalTrussConstitutiveLaw_CalculatesLinearElasticStress_WithInitialState, KratosStructuralMechanicsFastSuite) {
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    p_node_1->FastGetSolutionStepValue(TEMPERATURE) = 30.0;
    p_node_2->FastGetSolutionStepValue(TEMPERATURE) = 50.0;
    auto p_geometry = r_model_part.CreateNewGeometry("Line3D2", 1, std::vector<ModelPart::IndexType>{1, 2});

    Properties properties;
    constexpr auto youngs_modulus = 1.0e6;
    properties[YOUNG_MODULUS] = youngs_modulus;
    properties[THERMAL_EXPANSION_COEFFICIENT] = 1.0e-5;
    properties[REFERENCE_TEMPERATURE] = 20.0;

    constexpr auto induced_strain = 0.005;
    constexpr auto initial_strain = 0.01;
    constexpr auto initial_stress = 5000.0;
    Vector strain = ScalarVector(1, induced_strain);
    Vector stress = ZeroVector(1);
    ConstitutiveLaw::Parameters parameters;
    parameters.SetElementGeometry(*p_geometry);
    parameters.SetMaterialProperties(properties);
    parameters.SetStrainVector(strain);
    parameters.SetStressVector(stress);

    ThermalTrussConstitutiveLaw law;
    law.InitializeMaterial(properties, *p_geometry, Vector());
    const auto p_initial_state = Kratos::make_intrusive<InitialState>(
        ScalarVector(1, initial_strain), ScalarVector(1, initial_stress));
    law.SetInitialState(p_initial_state);
    law.CalculateMaterialResponsePK2(parameters);

    constexpr auto thermal_strain = 1.0e-5 * (40.0 - 20.0);
    constexpr auto expected_stress =
        (induced_strain - initial_strain - thermal_strain) * youngs_modulus + initial_stress;
    KRATOS_EXPECT_NEAR(parameters.GetStressVector()[0], expected_stress, 1.0e-12);
}

}
