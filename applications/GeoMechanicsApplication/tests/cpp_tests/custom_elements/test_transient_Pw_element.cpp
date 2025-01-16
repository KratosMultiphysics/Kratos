// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Gennady Markelov
//

#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/three_dimensional_stress_state.h"
#include "custom_elements/transient_Pw_element.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>

namespace
{
using namespace Kratos;

PointerVector<Node> CreateThreeNodes()
{
    PointerVector<Node> result;
    result.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    return result;
}

PointerVector<Node> CreateThreeNodesOnModelPart(ModelPart& rModelPart)
{
    PointerVector<Node> result;
    result.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    return result;
}

PointerVector<Node> CreateThetrahedralOnModelPart(ModelPart& rModelPart)
{
    PointerVector<Node> result;
    result.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(4, 1.0, 1.0, 1.0));
    return result;
}

PointerVector<Node> CreateCoincidentNodes()
{
    PointerVector<Node> result;
    result.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    return result;
}

ModelPart& CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    return r_result;
}

void RemoveThreeNodes(ModelPart& rModelPart)
{
    rModelPart.RemoveNodeFromAllLevels(1);
    rModelPart.RemoveNodeFromAllLevels(2);
    rModelPart.RemoveNodeFromAllLevels(3);
}

Element::IndexType NextElementNumber(const ModelPart& rModelPart)
{
    return rModelPart.NumberOfElements() + 1;
}

intrusive_ptr<TransientPwElement<2, 3>> CreateTransientPwElementWithPWDofs(const ModelPart& rModelPart,
                                                                           const Properties::Pointer& rProperties,
                                                                           const Geometry<Node>::Pointer& rGeometry)
{
    auto p_result = make_intrusive<TransientPwElement<2, 3>>(
        NextElementNumber(rModelPart), rGeometry, rProperties, std::make_unique<PlaneStrainStressState>());
    for (auto& node : p_result->GetGeometry()) {
        node.AddDof(WATER_PRESSURE);
    }

    return p_result;
}

intrusive_ptr<TransientPwElement<3, 4>> CreateTransientPwElement3D4NWithPWDofs(
    const ModelPart& rModelPart, const Properties::Pointer& rProperties, const Geometry<Node>::Pointer& rGeometry)
{
    auto p_result =
        make_intrusive<TransientPwElement<3, 4>>(NextElementNumber(rModelPart), rGeometry, rProperties,
                                                 std::make_unique<ThreeDimensionalStressState>());
    for (auto& node : p_result->GetGeometry()) {
        node.AddDof(WATER_PRESSURE);
    }

    return p_result;
}

intrusive_ptr<TransientPwElement<2, 3>> CreateTriangleTransientPwElementWithPWDofs(ModelPart& rModelPart,
                                                                                   const Properties::Pointer& rProperties)
{
    const auto p_geometry = std::make_shared<Triangle2D3<Node>>(CreateThreeNodesOnModelPart(rModelPart));
    auto p_element = CreateTransientPwElementWithPWDofs(rModelPart, rProperties, p_geometry);

    rModelPart.AddElement(p_element);
    return p_element;
}

intrusive_ptr<TransientPwElement<3, 4>> CreateThreeDTransientPwElement3D4NWithPWDofs(ModelPart& rModelPart,
                                                                                     const Properties::Pointer& rProperties)
{
    const auto p_geometry = std::make_shared<Tetrahedra3D4<Node>>(CreateThetrahedralOnModelPart(rModelPart));
    auto p_element = CreateTransientPwElement3D4NWithPWDofs(rModelPart, rProperties, p_geometry);

    rModelPart.AddElement(p_element);
    return p_element;
}

intrusive_ptr<TransientPwElement<2, 3>> CreateTriangleTransientPwElementWithoutPWDofs(ModelPart& rModelPart,
                                                                                      const Properties::Pointer& rProperties)
{
    const auto p_geometry = std::make_shared<Triangle2D3<Node>>(CreateThreeNodesOnModelPart(rModelPart));
    auto p_element = make_intrusive<TransientPwElement<2, 3>>(
        NextElementNumber(rModelPart), p_geometry, rProperties, std::make_unique<PlaneStrainStressState>());

    rModelPart.AddElement(p_element);
    return p_element;
}

template <unsigned int TDim, unsigned int TNumNodes>
void SetBasicPropertiesAndVariables(intrusive_ptr<TransientPwElement<TDim, TNumNodes>> rElement)
{
    rElement->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    rElement->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    rElement->GetProperties().SetValue(PERMEABILITY_XX, 1.0);
    rElement->GetProperties().SetValue(PERMEABILITY_YY, 1.0);
    rElement->GetProperties().SetValue(PERMEABILITY_XY, 1.0);
    if constexpr (TDim == 3) {
        rElement->GetProperties().SetValue(PERMEABILITY_ZZ, 1.0);
        rElement->GetProperties().SetValue(PERMEABILITY_YZ, 1.0);
        rElement->GetProperties().SetValue(PERMEABILITY_ZX, 1.0);
    }
    const auto gravity_acceleration = array_1d<double, 3>{0.0, -10.0, 0.0};
    rElement->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    rElement->GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    rElement->GetGeometry()[2].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    if constexpr (TNumNodes == 4) {
        rElement->GetGeometry()[2].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    }
    rElement->GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    rElement->GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    rElement->GetGeometry()[2].FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    if constexpr (TNumNodes == 4) {
        rElement->GetGeometry()[3].FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    }
}

} // namespace

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_CreateInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_geometry   = std::make_shared<Triangle2D3<Node>>(CreateThreeNodes());
    const auto p_properties = std::make_shared<Properties>();
    const TransientPwElement<2, 3> element(0, p_geometry, p_properties,
                                           std::make_unique<PlaneStrainStressState>());

    // Act
    const auto p_created_element = element.Create(1, p_geometry, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_CreateInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    const auto p_geometry   = std::make_shared<Triangle2D3<Node>>(CreateThreeNodes());
    const TransientPwElement<2, 3> element(0, p_geometry, p_properties,
                                           std::make_unique<PlaneStrainStressState>());

    // Act
    const auto p_created_element = element.Create(1, CreateThreeNodes(), p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_DoFList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model      model;
    auto&      r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    const auto p_element =
        CreateTriangleTransientPwElementWithPWDofs(r_model_part, std::make_shared<Properties>());

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    p_element->GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(degrees_of_freedom.size(), 3);
    for (auto p_dof : degrees_of_freedom) {
        KRATOS_EXPECT_EQ(p_dof->GetVariable(), WATER_PRESSURE);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_EquationIdVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateTriangleTransientPwElementWithPWDofs(r_model_part, std::make_shared<Properties>());

    unsigned int i = 0;
    for (const auto& node : p_element->GetGeometry()) {
        ++i;
        node.pGetDof(WATER_PRESSURE)->SetEquationId(i);
    }

    // Act
    const auto                    dummy_process_info = ProcessInfo{};
    Element::EquationIdVectorType equation_id_vector;
    p_element->EquationIdVector(equation_id_vector, dummy_process_info);

    // Assert
    const Element::EquationIdVectorType expected_ids = {1, 2, 3};
    KRATOS_EXPECT_VECTOR_EQ(equation_id_vector, expected_ids)
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_IntegrationMethod, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const TransientPwElement<2, 3> element(
        0, std::make_shared<Triangle2D3<Node>>(CreateCoincidentNodes()),
        std::make_shared<Properties>(), std::make_unique<PlaneStrainStressState>());

    // Act
    const auto p_integration_method = element.GetIntegrationMethod();

    // Assert
    constexpr auto expected_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    KRATOS_EXPECT_EQ(p_integration_method, expected_integration_method);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_CheckThrowsOnFaultyInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto                     p_properties = std::make_shared<Properties>();
    const TransientPwElement<2, 3> element_with_coincident_nodes(
        1, std::make_shared<Triangle2D3<Node>>(CreateCoincidentNodes()), p_properties,
        std::make_unique<PlaneStrainStressState>());

    // Act and Assert
    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element_with_coincident_nodes.Check(dummy_process_info),
                                      "Error: DomainSize < 1.0e-15 for the element 1")

    const TransientPwElement<2, 3> element_with_correct_domain_size(
        1, std::make_shared<Triangle2D3<Node>>(CreateThreeNodes()), p_properties,
        std::make_unique<PlaneStrainStressState>());
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element_with_correct_domain_size.Check(dummy_process_info),
                                      "Error: Missing variable WATER_PRESSURE on node 1")

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_element = CreateTriangleTransientPwElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: Missing variable DT_WATER_PRESSURE on node 1")

    RemoveThreeNodes(model_part);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    p_element = CreateTriangleTransientPwElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable VOLUME_ACCELERATION on node 1")

    RemoveThreeNodes(model_part);
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    p_element = CreateTriangleTransientPwElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable WATER_PRESSURE on node 1")

    RemoveThreeNodes(model_part);
    p_element = CreateTriangleTransientPwElementWithPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "DENSITY_WATER does not exist in the material properties or "
                                      "has an invalid value at element 4")

    p_element->GetProperties().SetValue(DENSITY_WATER, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "DENSITY_WATER does not exist in the material properties or "
                                      "has an invalid value at element 4")

    p_element->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: BULK_MODULUS_SOLID does not exist in the material "
                                      "properties or has an invalid value at element 4")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, -1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: BULK_MODULUS_SOLID does not exist in the material "
                                      "properties or has an invalid value at element 4")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: POROSITY does not exist in the material properties "
                                      "or has an invalid value at element 4")

    p_element->GetProperties().SetValue(POROSITY, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: POROSITY does not exist in the material properties "
                                      "or has an invalid value at element 4")

    p_element->GetProperties().SetValue(POROSITY, 2.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: POROSITY does not exist in the material properties "
                                      "or has an invalid value at element 4")

    p_element->GetProperties().SetValue(POROSITY, 0.5);

    p_element->GetGeometry().begin()->Z() += 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error:  Node with non-zero Z coordinate found. Id: 1")
    p_element->GetGeometry().begin()->Z() = 0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: BULK_MODULUS_FLUID does not exist in the material "
                                      "properties or has an invalid value at element 4")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, -1.0e6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: BULK_MODULUS_FLUID does not exist in the material "
                                      "properties or has an invalid value at element 4")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0e6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: DYNAMIC_VISCOSITY does not exist in the material "
                                      "properties or has an invalid value at element 4")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Error: DYNAMIC_VISCOSITY does not exist in the material "
                                      "properties or has an invalid value at element 4")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "PERMEABILITY_XX does not exist in the material properties "
                                      "or has an invalid value at element 4")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "PERMEABILITY_XX does not exist in the material properties "
                                      "or has an invalid value at element 4")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "PERMEABILITY_YY does not exist in the material properties "
                                      "or has an invalid value at element 4")

    p_element->GetProperties().SetValue(PERMEABILITY_YY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "PERMEABILITY_YY does not exist in the material properties "
                                      "or has an invalid value at element 4")

    p_element->GetProperties().SetValue(PERMEABILITY_YY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "PERMEABILITY_XY does not exist in the material properties "
                                      "or has an invalid value at element 4")

    p_element->GetProperties().SetValue(PERMEABILITY_XY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "PERMEABILITY_XY does not exist in the material properties "
                                      "or has an invalid value at element 4")

    p_element->GetProperties().SetValue(PERMEABILITY_XY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "Error: BIOT_COEFFICIENT does not exist in the material properties in element 4")

    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.0E-2);

    // No exceptions on correct input for 2D element
    KRATOS_EXPECT_EQ(p_element->Check(dummy_process_info), 0);

    auto p_3D_element = CreateThreeDTransientPwElement3D4NWithPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_3D_element->Check(dummy_process_info),
                                      "PERMEABILITY_ZZ does not exist in the material properties "
                                      "or has an invalid value at element 5")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZZ, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_3D_element->Check(dummy_process_info),
                                      "PERMEABILITY_ZZ does not exist in the material properties "
                                      "or has an invalid value at element 5")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZZ, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_3D_element->Check(dummy_process_info),
                                      "PERMEABILITY_YZ does not exist in the material properties "
                                      "or has an invalid value at element 5")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_YZ, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_3D_element->Check(dummy_process_info),
                                      "PERMEABILITY_YZ does not exist in the material properties "
                                      "or has an invalid value at element 5")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_YZ, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_3D_element->Check(dummy_process_info),
                                      "PERMEABILITY_ZX does not exist in the material properties "
                                      "or has an invalid value at element 5")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZX, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_3D_element->Check(dummy_process_info),
                                      "PERMEABILITY_ZX does not exist in the material properties "
                                      "or has an invalid value at element 5")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZX, 1.0E-2);

    // No exceptions on correct input for 3D element
    KRATOS_EXPECT_EQ(p_3D_element->Check(dummy_process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_Initialize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    TransientPwElement<2, 3> element(0, std::make_shared<Triangle2D3<Node>>(CreateCoincidentNodes()),
                                     std::make_shared<Properties>(),
                                     std::make_unique<PlaneStrainStressState>());
    const auto dummy_process_info = ProcessInfo{};

    // Act
    element.Initialize(dummy_process_info);

    // Assert
    const auto number_of_integration_points =
        element.GetGeometry().IntegrationPointsNumber(element.GetIntegrationMethod());
    const auto& r_constitutive_law_vector = element.mConstitutiveLawVector;
    KRATOS_EXPECT_EQ(r_constitutive_law_vector.size(), number_of_integration_points);
    for (const auto& constitutive_law : r_constitutive_law_vector) {
        KRATOS_EXPECT_EQ(constitutive_law, nullptr);
    }

    const auto& r_retention_law_vector = element.mRetentionLawVector;
    KRATOS_EXPECT_EQ(r_retention_law_vector.size(), number_of_integration_points);
    for (const auto& retention_law : r_retention_law_vector) {
        KRATOS_EXPECT_NE(dynamic_cast<SaturatedLaw*>(retention_law.get()), nullptr);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_InitializeSolution, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTriangleTransientPwElementWithPWDofs(r_model_part, std::make_shared<Properties>());
    const auto dummy_process_info = ProcessInfo{};

    // Act
    p_element->InitializeSolutionStep(dummy_process_info);

    // Assert
    for (auto& r_node : p_element->GetGeometry()) {
        KRATOS_EXPECT_EQ(r_node.FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 0.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_FinalizeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTriangleTransientPwElementWithPWDofs(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    const auto dummy_process_info = ProcessInfo{};
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    p_element->FinalizeSolutionStep(dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(p_element->GetGeometry()[0].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 500000);
    KRATOS_EXPECT_EQ(p_element->GetGeometry()[1].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 0);
    KRATOS_EXPECT_EQ(p_element->GetGeometry()[2].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), -500000);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_CalculateOnIntegrationPoints_Vector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTriangleTransientPwElementWithPWDofs(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    const auto dummy_process_info = ProcessInfo{};
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    std::vector<double> results{};
    p_element->CalculateOnIntegrationPoints(DEGREE_OF_SATURATION, results, dummy_process_info);
    auto const number_of_integration_points =
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod());

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    Vector expected_results(3);
    expected_results <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(EFFECTIVE_SATURATION, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(BISHOP_COEFFICIENT, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(DERIVATIVE_OF_SATURATION, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(RELATIVE_PERMEABILITY, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(HYDRAULIC_HEAD, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0.166667, 0.166667, 0.666667;
    KRATOS_EXPECT_VECTOR_NEAR(results, expected_results, Defaults::relative_tolerance);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(DT_WATER_PRESSURE, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0, 0, 0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_CalculateOnIntegrationPoints_1DArray, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTriangleTransientPwElementWithPWDofs(r_model_part, p_properties);
    SetBasicPropertiesAndVariables(p_element);
    const auto dummy_process_info = ProcessInfo{};
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    std::vector<array_1d<double, 3>> results{};
    p_element->CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, results, dummy_process_info);
    auto const number_of_integration_points =
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod());

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    array_1d<double, 3> expected_nonzero_component{-1e+06, -1e+06, 0};
    for (const auto& component : results) {
        KRATOS_EXPECT_VECTOR_EQ(component, expected_nonzero_component);
    }

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(LOCAL_FLUID_FLUX_VECTOR, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    array_1d<double, 3> expected_zero_component{0, 0, 0};
    for (const auto& component : results) {
        KRATOS_EXPECT_VECTOR_EQ(component, expected_zero_component);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement_CalculateOnIntegrationPoints_Matrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTriangleTransientPwElementWithPWDofs(r_model_part, p_properties);
    SetBasicPropertiesAndVariables(p_element);
    const auto dummy_process_info = ProcessInfo{};
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    std::vector<Matrix> results{};
    p_element->CalculateOnIntegrationPoints(PERMEABILITY_MATRIX, results, dummy_process_info);
    auto const number_of_integration_points =
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod());

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    Matrix expected_nonzero_component(2, 2);
    expected_nonzero_component <<= 1, 1, 1, 1;
    for (const auto& component : results) {
        KRATOS_EXPECT_MATRIX_EQ(component, expected_nonzero_component);
    }

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(LOCAL_PERMEABILITY_MATRIX, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    Matrix expected_zero_component = ZeroMatrix(2, 2);
    for (const auto& component : results) {
        KRATOS_EXPECT_MATRIX_EQ(component, expected_zero_component);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement2D3N_CalculateLocalSystem, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTriangleTransientPwElementWithPWDofs(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 0.5);
    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0E6);
    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    p_element->GetProperties().SetValue(POROSITY, 0.1);
    p_element->GetProperties().SetValue(IGNORE_UNDRAINED, false);
    const auto dummy_process_info = ProcessInfo{};
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    p_element->CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, dummy_process_info);

    // Assert
    Matrix expected_left_hand_side(3, 3);
    // clang-format off
    expected_left_hand_side <<= -49.999999999999993,0,49.999999999999993,
                                 0,0,0,
                                 49.999999999999993,0,-49.999999999999993;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    Vector expected_right_hand_side(3);
    expected_right_hand_side <<= 500000, 0, -500000;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwElement3D4N_CalculateLocalSystem, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateThreeDTransientPwElement3D4NWithPWDofs(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 0.5);
    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0E6);
    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    p_element->GetProperties().SetValue(POROSITY, 0.1);
    p_element->GetProperties().SetValue(IGNORE_UNDRAINED, false);
    const auto dummy_process_info = ProcessInfo{};
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    p_element->CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, dummy_process_info);

    // Assert
    Matrix expected_left_hand_side(4, 4);
    // clang-format off
    expected_left_hand_side <<= -16.666666666666664,0,0,16.666666666666664,
                                 0,0,0,0,
                                 0,0,0,0,
                                 16.666666666666664,0,0,-16.666666666666664;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    Vector expected_right_hand_side(4);
    expected_right_hand_side <<= 125000, 0, 0, -125000;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

} // namespace Kratos::Testing