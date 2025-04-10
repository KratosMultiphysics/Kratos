// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_constitutive_law.h"
#include "tests/cpp_tests/test_utilities.h"
#include "tests/cpp_tests/test_utilities/element_setup_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <custom_constitutive/incremental_linear_elastic_law.h>
#include <custom_constitutive/plane_strain.h>
#include <custom_elements/plane_strain_stress_state.h>
#include <custom_utilities/registration_utilities.h>

namespace
{

using namespace Kratos;

ModelPart& CreateModelPartWithUPwSolutionStepVariables(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.SetBufferSize(2);
    r_result.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_result.AddNodalSolutionStepVariable(VELOCITY);
    r_result.AddNodalSolutionStepVariable(ACCELERATION);
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_result.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);

    return r_result;
}

std::shared_ptr<Properties> CreateProperties()
{
    const auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(CONSTITUTIVE_LAW, std::make_shared<GeoIncrementalLinearElasticLaw>(
                                                 std::make_unique<PlaneStrain>()));
    p_properties->SetValue(YOUNG_MODULUS, 1.000000e+07);
    p_properties->SetValue(POISSON_RATIO, 0.000000e+00);
    p_properties->SetValue(DENSITY_SOLID, 2.650000e+03);
    p_properties->SetValue(DENSITY_WATER, 1.000000e+03);
    p_properties->SetValue(POROSITY, 1.000000e-01);
    p_properties->SetValue(BULK_MODULUS_SOLID, 1.000000e+12);
    p_properties->SetValue(BULK_MODULUS_FLUID, 200.0); // small to get a significant value for the compressibility term
    p_properties->SetValue(PERMEABILITY_XX, 9.084000e-06);
    p_properties->SetValue(PERMEABILITY_YY, 9.084000e-06);
    p_properties->SetValue(PERMEABILITY_XY, 0.000000e+00);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    // Biot alpha = 0, no coupling
    p_properties->SetValue(BIOT_COEFFICIENT, 0.000000e+00);
    p_properties->SetValue(RETENTION_LAW, "SaturatedLaw");
    p_properties->SetValue(SATURATED_SATURATION, 1.000000e+00);
    p_properties->SetValue(IGNORE_UNDRAINED, false);

    return p_properties;
}

void SetSolutionStepValuesForFluidFluxCheck(const Element::Pointer& rElement)
{
    const auto zero_values = array_1d<double, 3>{0.0, 0.0, 0.0};
    for (auto& r_node : rElement->GetGeometry()) {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) = zero_values;
        r_node.FastGetSolutionStepValue(VELOCITY)     = zero_values;
        // Zero acceleration -> no Fluid Body Flow
        r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = zero_values;
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)      = 1.0E4;
        r_node.FastGetSolutionStepValue(DT_WATER_PRESSURE)   = 0.0;
    }
    rElement->GetGeometry()[2].FastGetSolutionStepValue(WATER_PRESSURE) = 2.0E4;
}

void SetSolutionStepValuesForGeneralCheck(const Element::Pointer& rElement)
{
    const auto zero_values = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto gravity     = array_1d<double, 3>{0.0, -10.0, 0.0};

    for (auto& r_node : rElement->GetGeometry()) {
        r_node.FastGetSolutionStepValue(VELOCITY)            = zero_values;
        r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity;
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)      = 1.0E4;
        r_node.FastGetSolutionStepValue(DT_WATER_PRESSURE)   = 0.0;
    }
    rElement->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{-0.015, 0.0, 0.0};
    rElement->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.015, 0.00, 0.0};
    rElement->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.0, 0.015, 0.0};
}

Element::Pointer CreateSmallStrainUPwDiffOrderElementWithUPwDofs(const Properties::Pointer& rProperties,
                                                                 const Geometry<Node>::Pointer& rGeometry)
{
    auto p_result = make_intrusive<SmallStrainUPwDiffOrderElement>(
        1, rGeometry, rProperties, std::make_unique<PlaneStrainStressState>());
    for (auto& r_node : p_result->GetGeometry()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        r_node.AddDof(WATER_PRESSURE);
        r_node.AddDof(DT_WATER_PRESSURE);
    }

    return p_result;
}

auto CreateSmallStrainUPwDiffOrderElementWithUPwDofs(Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithUPwSolutionStepVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, -1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 0.0, -0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.5, -0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(6, 0.5, 0.05, 0.0));
    const auto p_geometry = std::make_shared<Triangle2D6<Node>>(nodes);
    return CreateSmallStrainUPwDiffOrderElementWithUPwDofs(rProperties, p_geometry);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateShearCapacity, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_element = ElementSetupUtilities::Create2D6NDiffOrderElement();

    auto& p_properties = p_element->GetProperties();
    p_properties.SetValue(CONSTITUTIVE_LAW, std::make_shared<StubConstitutiveLaw>());
    p_properties.SetValue(GEO_COHESION, 2.0);
    p_properties.SetValue(GEO_FRICTION_ANGLE, 0.0);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    auto stress_vector = Vector{4};
    stress_vector <<= -1.5, 0.0, 1.5, 0.0;
    p_element->SetValuesOnIntegrationPoints(
        CAUCHY_STRESS_VECTOR, std::vector<Vector>{3, stress_vector}, dummy_process_info);

    // Act
    auto actual_shear_capacity_values = std::vector<double>{};
    p_element->CalculateOnIntegrationPoints(GEO_SHEAR_CAPACITY, actual_shear_capacity_values, dummy_process_info);

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(actual_shear_capacity_values, (Vector{ScalarVector{3, 0.75}}),
                              Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateLHS, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  p_element = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(model, CreateProperties());

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    // Act
    auto actual_lhs_values = Matrix{};
    p_element->CalculateLeftHandSide(actual_lhs_values, dummy_process_info);

    EXPECT_EQ(actual_lhs_values.size1(), 15);
    EXPECT_EQ(actual_lhs_values.size2(), 15);
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateLHS_WithSaveAndLoad, KratosGeoMechanicsFastSuite)
{
    ScopedSerializerRegistration registration("SaturatedLaw", SaturatedLaw{});
    ScopedSerializerRegistration registration2("PlaneStrain", PlaneStrain{});
    ScopedSerializerRegistration registration3("PlaneStrainStressState", PlaneStrainStressState{});
    // Arrange
    Model model;
    auto  p_element = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(model, CreateProperties());

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    auto serializer = StreamSerializer{};
    serializer.save("test_tag", p_element);

    // Act
    auto p_loaded_element = make_intrusive<SmallStrainUPwDiffOrderElement>();
    serializer.load("test_tag", p_loaded_element);
    // Act
    auto actual_lhs_values = Matrix{};
    p_loaded_element->CalculateLeftHandSide(actual_lhs_values, dummy_process_info);

    EXPECT_EQ(actual_lhs_values.size1(), 15);
    EXPECT_EQ(actual_lhs_values.size2(), 15);
}

} // namespace Kratos::Testing
