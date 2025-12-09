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

#include "containers/model.h"
#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/three_dimensional_stress_state.h"
#include "custom_elements/transient_Pw_element.hpp"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "includes/cfd_variables.h"
#include "includes/expect.h"
#include "includes/model_part.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
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

PointerVector<Node> CreateThreeCoincidentNodes()
{
    PointerVector<Node> result;
    for (unsigned int id = 1; id <= 3; id++) {
        result.push_back(make_intrusive<Node>(id, 0.0, 0.0, 0.0));
    }
    return result;
}

template <unsigned int TNumNodes>
PointerVector<Node> CreateNodesOnModelPart(ModelPart& rModelPart)
{
    PointerVector<Node> result;
    result.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    if constexpr (TNumNodes == 4) {
        result.push_back(rModelPart.CreateNewNode(4, 1.0, 1.0, 1.0));
    }
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

template <unsigned int TDim, unsigned int TNumNodes>
intrusive_ptr<TransientPwElement<TDim, TNumNodes>> CreateTransientPwElementWithPWDofs(ModelPart& rModelPart,
                                                                                      const Properties::Pointer& rProperties)
{
    intrusive_ptr<TransientPwElement<TDim, TNumNodes>> p_element;
    if constexpr (TDim == 2) {
        p_element = make_intrusive<TransientPwElement<TDim, TNumNodes>>(
            NextElementNumber(rModelPart),
            std::make_shared<Triangle2D3<Node>>(CreateNodesOnModelPart<TNumNodes>(rModelPart)),
            rProperties, std::make_unique<PlaneStrainStressState>());
    } else {
        p_element = make_intrusive<TransientPwElement<TDim, TNumNodes>>(
            NextElementNumber(rModelPart),
            std::make_shared<Tetrahedra3D4<Node>>(CreateNodesOnModelPart<TNumNodes>(rModelPart)),
            rProperties, std::make_unique<ThreeDimensionalStressState>());
    }
    for (auto& r_node : p_element->GetGeometry()) {
        r_node.AddDof(WATER_PRESSURE);
    }
    rModelPart.AddElement(p_element);
    return p_element;
}

intrusive_ptr<TransientPwElement<2, 3>> CreateTriangleTransientPwElementWithoutPWDofs(ModelPart& rModelPart,
                                                                                      const Properties::Pointer& rProperties)
{
    auto p_element = make_intrusive<TransientPwElement<2, 3>>(
        NextElementNumber(rModelPart),
        std::make_shared<Triangle2D3<Node>>(CreateNodesOnModelPart<3>(rModelPart)), rProperties,
        std::make_unique<PlaneStrainStressState>(), nullptr);

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
    for (auto& r_node : rElement->GetGeometry()) {
        r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)      = 0.0;
    }
}

} // namespace

namespace Kratos::Testing
{

using namespace Kratos;

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_CreateInstanceWithGeometryInput)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_CreateInstanceWithNodeInput)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    const TransientPwElement<2, 3> element(0, std::make_shared<Triangle2D3<Node>>(CreateThreeNodes()),
                                           p_properties, std::make_unique<PlaneStrainStressState>());

    // Act
    const auto p_created_element = element.Create(1, CreateThreeNodes(), p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_DoFList)
{
    // Arrange
    Model      model;
    auto&      r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    const auto p_element =
        CreateTransientPwElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    p_element->GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    EXPECT_EQ(degrees_of_freedom.size(), 3);
    KRATOS_EXPECT_TRUE(std::all_of(degrees_of_freedom.begin(), degrees_of_freedom.end(),
                                   [](auto p_dof) { return p_dof->GetVariable() == WATER_PRESSURE; }))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_EquationIdVector)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto p_element = CreateTransientPwElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());

    unsigned int i = 0;
    for (const auto& r_node : p_element->GetGeometry()) {
        ++i;
        r_node.pGetDof(WATER_PRESSURE)->SetEquationId(i);
    }

    // Act
    const auto                    dummy_process_info = ProcessInfo{};
    Element::EquationIdVectorType equation_id_vector;
    p_element->EquationIdVector(equation_id_vector, dummy_process_info);

    // Assert
    const Element::EquationIdVectorType expected_ids = {1, 2, 3};
    KRATOS_EXPECT_VECTOR_EQ(equation_id_vector, expected_ids)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_IntegrationMethod)
{
    // Arrange
    const TransientPwElement<2, 3> element(
        0, std::make_shared<Triangle2D3<Node>>(CreateThreeCoincidentNodes()),
        std::make_shared<Properties>(), std::make_unique<PlaneStrainStressState>(), nullptr);

    // Act
    const auto p_integration_method = element.GetIntegrationMethod();

    // Assert
    constexpr auto expected_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    EXPECT_EQ(p_integration_method, expected_integration_method);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_CheckThrowsOnFaultyInput)
{
    // Arrange
    const auto                     p_properties = std::make_shared<Properties>();
    const TransientPwElement<2, 3> element_with_coincident_nodes(
        1, std::make_shared<Triangle2D3<Node>>(CreateThreeCoincidentNodes()), p_properties,
        std::make_unique<PlaneStrainStressState>(), nullptr);

    // Act and Assert
    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element_with_coincident_nodes.Check(dummy_process_info),
                                      "DomainSize (0) is smaller than 1e-15 for element 1")

    const TransientPwElement<2, 3> element_with_correct_domain_size(
        1, std::make_shared<Triangle2D3<Node>>(CreateThreeNodes()), p_properties,
        std::make_unique<PlaneStrainStressState>(), nullptr);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element_with_correct_domain_size.Check(dummy_process_info),
                                      "Missing variable WATER_PRESSURE on nodes 1 2 3")

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_element = CreateTriangleTransientPwElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable DT_WATER_PRESSURE on nodes 1 2 3")

    RemoveThreeNodes(model_part);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    p_element = CreateTriangleTransientPwElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable VOLUME_ACCELERATION on nodes 1 2 3")

    RemoveThreeNodes(model_part);
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    p_element = CreateTriangleTransientPwElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "Missing the DoF for the variable WATER_PRESSURE on nodes 1 2 3")

    RemoveThreeNodes(model_part);
    p_element = CreateTransientPwElementWithPWDofs<2, 3>(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_WATER does not exist in the material properties with Id 0 at element with Id 4.")

    p_element->GetProperties().SetValue(DENSITY_WATER, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_WATER in the material properties with Id 0 at element with Id 4 has an invalid "
        "value: -1000 is out of the range [0, -).")

    p_element->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "BULK_MODULUS_SOLID does not exist in the material "
                                      "properties with Id 0 at element with Id 4.")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, -1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_SOLID in the material properties with Id 0 at element with Id 4 has an "
        "invalid value: -1e+06 is out of the range [0, -).")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "POROSITY does not exist in the material properties with Id 0 at element with Id 4.")

    p_element->GetProperties().SetValue(POROSITY, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "POROSITY in the material properties with Id 0 at element with Id 4 has an invalid value: "
        "-1 is out of the range [0, 1].")

    p_element->GetProperties().SetValue(POROSITY, 2.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "POROSITY in the material properties with Id 0 at element with Id 4 has an invalid value: "
        "2 is out of the range [0, 1].")

    p_element->GetProperties().SetValue(POROSITY, 0.5);

    p_element->GetGeometry().begin()->Z() += 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Node with Id: 1 has non-zero Z coordinate.")
    p_element->GetGeometry().begin()->Z() = 0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "BULK_MODULUS_FLUID does not exist in the material "
                                      "properties with Id 0 at element with Id 4.")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, -1.0e6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_FLUID in the material properties with Id 0 at element with Id 4 has an "
        "invalid value: -1e+06 is out of the range (0, -).")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0e6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "DYNAMIC_VISCOSITY does not exist in the material properties "
                                      "with Id 0 at element with Id 4.")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DYNAMIC_VISCOSITY in the material properties with Id 0 at element with Id 4 has an "
        "invalid value: -0.01 is out of the range (0, -).")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "BIOT_COEFFICIENT does not exist in the material properties "
                                      "with Id 0 at element with Id 4.")

    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XX does not exist in the material properties with Id 0 at element with Id 4.")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XX in the material properties with Id 0 at element with Id 4 has an invalid "
        "value: -0.01 is out of the range [0, -).")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_YY does not exist in the material properties with Id 0 at element with Id 4.")

    p_element->GetProperties().SetValue(PERMEABILITY_YY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_YY in the material properties with Id 0 at element with Id 4 has an invalid "
        "value: -0.01 is out of the range [0, -).")

    p_element->GetProperties().SetValue(PERMEABILITY_YY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XY does not exist in the material properties with Id 0 at element with Id 4.")

    p_element->GetProperties().SetValue(PERMEABILITY_XY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XY in the material properties with Id 0 at element with Id 4 has an invalid "
        "value: -0.01 is out of the range [0, -).")

    p_element->GetProperties().SetValue(PERMEABILITY_XY, 1.0E-2);

    // No exceptions on correct input for 2D element when retention law vector is initialized
    p_element->Initialize(dummy_process_info);

    EXPECT_EQ(p_element->Check(dummy_process_info), 0);

    auto p_3D_element = CreateTransientPwElementWithPWDofs<3, 4>(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_ZZ does not exist in the material properties with Id 0 at element with Id 5.")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZZ, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_ZZ in the material properties with Id 0 at element with Id 5 has an invalid "
        "value: -0.01 is out of the range [0, -).")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZZ, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_YZ does not exist in the material properties with Id 0 at element with Id 5.")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_YZ, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_YZ in the material properties with Id 0 at element with Id 5 has an invalid "
        "value: -0.01 is out of the range [0, -).")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_YZ, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_ZX does not exist in the material properties with Id 0 at element with Id 5.")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZX, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_ZX in the material properties with Id 0 at element with Id 5 has an invalid "
        "value: -0.01 is out of the range [0, -).")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZX, 1.0E-2);

    // to enable a call of RetentionLaw check
    p_3D_element->Initialize(dummy_process_info);

    // No exceptions on correct input for 3D element
    EXPECT_EQ(p_3D_element->Check(dummy_process_info), 0);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_Initialize)
{
    // Arrange
    TransientPwElement<2, 3> element(
        0, std::make_shared<Triangle2D3<Node>>(CreateThreeCoincidentNodes()),
        std::make_shared<Properties>(), std::make_unique<PlaneStrainStressState>(), nullptr);
    const auto dummy_process_info = ProcessInfo{};

    // Act
    element.Initialize(dummy_process_info);

    // Assert
    const auto number_of_integration_points =
        element.GetGeometry().IntegrationPointsNumber(element.GetIntegrationMethod());
    const auto& r_constitutive_law_vector = element.mConstitutiveLawVector;
    EXPECT_EQ(r_constitutive_law_vector.size(), number_of_integration_points);
    for (const auto& constitutive_law : r_constitutive_law_vector) {
        EXPECT_EQ(constitutive_law, nullptr);
    }

    const auto& r_retention_law_vector = element.mRetentionLawVector;
    EXPECT_EQ(r_retention_law_vector.size(), number_of_integration_points);
    KRATOS_EXPECT_TRUE(std::none_of(r_retention_law_vector.begin(), r_retention_law_vector.end(), [](auto p_retention_law) {
        return dynamic_cast<SaturatedLaw*>(p_retention_law.get()) == nullptr;
    }))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_InitializeSolution)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTransientPwElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    // Act
    p_element->InitializeSolutionStep(dummy_process_info);

    // Assert
    KRATOS_EXPECT_TRUE(std::all_of(
        p_element->GetGeometry().begin(), p_element->GetGeometry().end(),
        [](auto& r_node) { return r_node.FastGetSolutionStepValue(HYDRAULIC_DISCHARGE) == 0.0; }))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_FinalizeSolutionStep)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTransientPwElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    p_element->FinalizeSolutionStep(dummy_process_info);

    // Assert
    EXPECT_EQ(p_element->GetGeometry()[0].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 500000);
    EXPECT_EQ(p_element->GetGeometry()[1].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 0);
    EXPECT_EQ(p_element->GetGeometry()[2].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), -500000);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_CalculateOnIntegrationPoints_Vector)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTransientPwElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    std::vector<double> results{};
    p_element->CalculateOnIntegrationPoints(DEGREE_OF_SATURATION, results, dummy_process_info);
    auto const number_of_integration_points =
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod());

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    Vector expected_results(3);
    expected_results <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(EFFECTIVE_SATURATION, results, dummy_process_info);

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(BISHOP_COEFFICIENT, results, dummy_process_info);

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(DERIVATIVE_OF_SATURATION, results, dummy_process_info);

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(RELATIVE_PERMEABILITY, results, dummy_process_info);

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(HYDRAULIC_HEAD, results, dummy_process_info);

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0.166667, 0.166667, 0.666667;
    KRATOS_EXPECT_VECTOR_NEAR(results, expected_results, Defaults::relative_tolerance);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(DT_WATER_PRESSURE, results, dummy_process_info);

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0, 0, 0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_CalculateOnIntegrationPoints_1DArray)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTransientPwElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    std::vector<array_1d<double, 3>> results{};
    p_element->CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, results, dummy_process_info);
    auto const number_of_integration_points =
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod());

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    array_1d<double, 3> expected_nonzero_component{-1e+06, -1e+06, 0};
    for (const auto& component : results) {
        KRATOS_EXPECT_VECTOR_EQ(component, expected_nonzero_component);
    }

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(LOCAL_FLUID_FLUX_VECTOR, results, dummy_process_info);

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    array_1d<double, 3> expected_zero_component{0, 0, 0};
    for (const auto& component : results) {
        KRATOS_EXPECT_VECTOR_EQ(component, expected_zero_component);
    }
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_CalculateOnIntegrationPoints_Matrix)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTransientPwElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    std::vector<Matrix> results{};
    p_element->CalculateOnIntegrationPoints(PERMEABILITY_MATRIX, results, dummy_process_info);
    auto const number_of_integration_points =
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod());

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    Matrix expected_nonzero_component(2, 2);
    expected_nonzero_component <<= 1, 1, 1, 1;
    for (const auto& component : results) {
        KRATOS_EXPECT_MATRIX_EQ(component, expected_nonzero_component);
    }

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(LOCAL_PERMEABILITY_MATRIX, results, dummy_process_info);

    // Assert
    EXPECT_EQ(results.size(), number_of_integration_points);
    Matrix expected_zero_component = ZeroMatrix(2, 2);
    for (const auto& component : results) {
        KRATOS_EXPECT_MATRIX_EQ(component, expected_zero_component);
    }
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement2D3N_CalculateLocalSystem)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTransientPwElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 0.5);
    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0E6);
    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    p_element->GetProperties().SetValue(POROSITY, 0.1);
    p_element->GetProperties().SetValue(IGNORE_UNDRAINED, false);
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement3D4N_CalculateLocalSystem)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTransientPwElementWithPWDofs<3, 4>(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 0.5);
    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0E6);
    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    p_element->GetProperties().SetValue(POROSITY, 0.1);
    p_element->GetProperties().SetValue(IGNORE_UNDRAINED, false);
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
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
    expected_right_hand_side <<= 166666.666, 0, 0, -166666.666;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TransientPwElement_ZeroReturnFunctions)
{
    // Arrange
    TransientPwElement<2, 3> element(
        0, std::make_shared<Triangle2D3<Node>>(CreateThreeCoincidentNodes()),
        std::make_shared<Properties>(), std::make_unique<PlaneStrainStressState>(), nullptr);
    const auto   dummy_process_info = ProcessInfo{};
    const auto   n_DoF              = 3;
    const Matrix expected_matrix    = ZeroMatrix(n_DoF, n_DoF);
    const Vector expected_vector    = ZeroVector(n_DoF);

    // Act
    Matrix actual_matrix;
    element.CalculateMassMatrix(actual_matrix, dummy_process_info);

    // Assert
    EXPECT_EQ(actual_matrix.size1(), n_DoF);
    EXPECT_EQ(actual_matrix.size2(), n_DoF);
    KRATOS_EXPECT_MATRIX_EQ(actual_matrix, expected_matrix);

    // Act
    element.CalculateDampingMatrix(actual_matrix, dummy_process_info);

    // Assert
    EXPECT_EQ(actual_matrix.size1(), n_DoF);
    EXPECT_EQ(actual_matrix.size2(), n_DoF);
    KRATOS_EXPECT_MATRIX_EQ(actual_matrix, expected_matrix);

    // Act
    Vector actual_vector;
    element.GetValuesVector(actual_vector, dummy_process_info);

    // Assert
    EXPECT_EQ(actual_vector.size(), n_DoF);
    KRATOS_EXPECT_VECTOR_EQ(actual_vector, expected_vector);

    // Act
    element.GetFirstDerivativesVector(actual_vector, dummy_process_info);

    // Assert
    EXPECT_EQ(actual_vector.size(), n_DoF);
    KRATOS_EXPECT_VECTOR_EQ(actual_vector, expected_vector);

    // Act
    element.GetSecondDerivativesVector(actual_vector, dummy_process_info);

    // Assert
    EXPECT_EQ(actual_vector.size(), n_DoF);
    KRATOS_EXPECT_VECTOR_EQ(actual_vector, expected_vector);
}

} // namespace Kratos::Testing