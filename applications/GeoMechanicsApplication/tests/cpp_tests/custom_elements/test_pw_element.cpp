// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_elements/Pw_element.h"
#include "custom_elements/calculation_contribution.h"
#include "custom_elements/integration_coefficient_modifier_for_line_element.h"
#include "geometries/line_2d_4.h"
#include "geometries/line_2d_5.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include "custom_utilities/registration_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;
using namespace std::string_literals;

ModelPart& CreateModelPartWithSolutionStepVariables(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    return r_result;
}

PwElement<2, 2> TransientPwLineElementWithoutPWDofs(const Properties::Pointer&     rProperties,
                                                    const Geometry<Node>::Pointer& rGeometry)
{
    using enum CalculationContribution;
    return PwElement<2, 2>{1,
                           rGeometry,
                           rProperties,
                           {Permeability, Compressibility, FluidBodyFlow},
                           std::make_unique<IntegrationCoefficientModifierForLineElement>()};
}

intrusive_ptr<Element> CreatePwLineElementWithoutPWDofs(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(1, 1.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    using enum CalculationContribution;
    const std::vector contributions = {Permeability, Compressibility, FluidBodyFlow};
    auto              element       = make_intrusive<PwElement<2, 2>>(
        rModelPart.NumberOfElements() + 1, p_geometry, rProperties, contributions,
        std::make_unique<IntegrationCoefficientModifierForLineElement>());
    rModelPart.AddElement(element);
    return element;
}

void RemoveTwoNodes(ModelPart& rModelPart)
{
    rModelPart.RemoveNodeFromAllLevels(0);
    rModelPart.RemoveNodeFromAllLevels(1);
}

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

PointerVector<Node> CreateInclinedLine(ModelPart& rModelPart)
{
    PointerVector<Node> result;
    result.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(2, 1.0, -1.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(3, 0.5, -0.5, 0.0));

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

void DeleteThreeNodes(ModelPart& rModelPart)
{
    rModelPart.RemoveNodeFromAllLevels(1);
    rModelPart.RemoveNodeFromAllLevels(2);
    rModelPart.RemoveNodeFromAllLevels(3);
}

Element::IndexType GetNextElementNumber(const ModelPart& rModelPart)
{
    return rModelPart.NumberOfElements() + 1;
}

template <unsigned int TDim, unsigned int TNumNodes>
intrusive_ptr<PwElement<TDim, TNumNodes>> CreateTransientPwLineElementWithPWDofs(ModelPart& rModelPart,
                                                                                 const Properties::Pointer& rProperties)
{
    intrusive_ptr<PwElement<TDim, TNumNodes>> p_element;
    using enum CalculationContribution;
    const std::vector contributions = {Permeability, Compressibility, FluidBodyFlow};
    if constexpr (TDim == 2) {
        p_element = make_intrusive<PwElement<TDim, TNumNodes>>(
            GetNextElementNumber(rModelPart),
            std::make_shared<Triangle2D3<Node>>(CreateNodesOnModelPart<TNumNodes>(rModelPart)),
            rProperties, contributions, nullptr);
    } else {
        p_element = make_intrusive<PwElement<TDim, TNumNodes>>(
            GetNextElementNumber(rModelPart),
            std::make_shared<Tetrahedra3D4<Node>>(CreateNodesOnModelPart<TNumNodes>(rModelPart)),
            rProperties, contributions, nullptr);
    }
    for (auto& r_node : p_element->GetGeometry()) {
        r_node.AddDof(WATER_PRESSURE);
    }
    rModelPart.AddElement(p_element);
    return p_element;
}

intrusive_ptr<PwElement<2, 3>> CreateTriangleTransientPwLineElementWithoutPWDofs(ModelPart& rModelPart,
                                                                                 const Properties::Pointer& rProperties)
{
    using enum CalculationContribution;
    const std::vector contributions = {Permeability, Compressibility, FluidBodyFlow};
    auto              p_element     = make_intrusive<PwElement<2, 3>>(
        GetNextElementNumber(rModelPart),
        std::make_shared<Triangle2D3<Node>>(CreateNodesOnModelPart<3>(rModelPart)), rProperties,
        contributions, nullptr);

    rModelPart.AddElement(p_element);
    return p_element;
}

template <unsigned int TDim, unsigned int TNumNodes>
void SetBasicPropertiesAndVariables(intrusive_ptr<PwElement<TDim, TNumNodes>> rElement)
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

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_ReturnsTheExpectedLeftHandSideAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(YOUNG_MODULUS, 1.000000e+07);
    p_properties->SetValue(POISSON_RATIO, 0.000000e+00);
    p_properties->SetValue(DENSITY_SOLID, 2.650000e+03);
    p_properties->SetValue(DENSITY_WATER, 1.000000e+03);
    p_properties->SetValue(POROSITY, 1.000000e-01);
    p_properties->SetValue(BULK_MODULUS_SOLID, 1.000000e+12);
    p_properties->SetValue(BULK_MODULUS_FLUID, 200.0); // small to get a significant value for the compressibility term
    p_properties->SetValue(PERMEABILITY_XX, 9.084000e-06);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    p_properties->SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    p_properties->SetValue(RETENTION_LAW, "SaturatedLaw");
    p_properties->SetValue(SATURATED_SATURATION, 1.000000e+00);
    p_properties->SetValue(CROSS_AREA, 1.0);
    p_properties->SetValue(IGNORE_UNDRAINED, false);

    auto process_info                     = ProcessInfo{};
    process_info[DT_PRESSURE_COEFFICIENT] = 1.5;

    Model model;
    auto& r_model_part = CreateModelPartWithSolutionStepVariables(model);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 1.0, 1.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    using enum CalculationContribution;
    auto element = PwElement<2, 2>{1,
                                   p_geometry,
                                   p_properties,
                                   {Permeability, Compressibility, FluidBodyFlow},
                                   std::make_unique<IntegrationCoefficientModifierForLineElement>()};
    for (auto& node : element.GetGeometry()) {
        node.AddDof(WATER_PRESSURE);
        node.AddDof(DT_WATER_PRESSURE);
    }
    element.GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    element.GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    element.GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE)    = 10.0;
    element.GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE)    = 0.0;
    element.GetGeometry()[0].FastGetSolutionStepValue(DT_WATER_PRESSURE) = 4.0;
    element.GetGeometry()[1].FastGetSolutionStepValue(DT_WATER_PRESSURE) = 5.0;

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.Initialize(process_info);
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, process_info);

    Vector actual_isolated_right_hand_side;
    element.CalculateRightHandSide(actual_isolated_right_hand_side, process_info);
    Matrix actual_isolated_left_hand_side;
    element.CalculateLeftHandSide(actual_isolated_left_hand_side, process_info);

    // Assert
    // clang-format off
    auto expected_left_hand_side = Matrix{2, 2};
    expected_left_hand_side <<= -0.00099588919125952972,  0.00046555910441502474,
                                 0.00046555910441502474, -0.00099588919125952972;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_isolated_left_hand_side, expected_left_hand_side,
                                       Defaults::relative_tolerance)

    auto expected_right_hand_side = Vector{2};
    expected_right_hand_side <<= 6.43131, -6.42813;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_isolated_right_hand_side, expected_right_hand_side,
                                       Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_CheckThrowsOnFaultyInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto          p_properties = std::make_shared<Properties>();
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(0, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    const auto p_geometry_with_coincide_nodes = std::make_shared<Line2D2<Node>>(nodes);
    const auto element_with_coincident_nodes =
        TransientPwLineElementWithoutPWDofs(p_properties, p_geometry_with_coincide_nodes);

    // Act and Assert
    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element_with_coincident_nodes.Check(dummy_process_info),
                                      "Error: Length (0) is smaller than 1e-15 for element 1")
    nodes.erase(nodes.begin() + 1);
    nodes.push_back(make_intrusive<Node>(1, 1.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    const auto element_with_correct_domain_size = TransientPwLineElementWithoutPWDofs(p_properties, p_geometry);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element_with_correct_domain_size.Check(dummy_process_info),
                                      "Missing variable WATER_PRESSURE on nodes 0 1")

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_element = CreatePwLineElementWithoutPWDofs(model_part, p_properties);

    RemoveTwoNodes(model_part);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    p_element = CreatePwLineElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable VOLUME_ACCELERATION on nodes 0 1")

    RemoveTwoNodes(model_part);
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    p_element = CreatePwLineElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "Missing the DoF for the variable WATER_PRESSURE on nodes 0 1")

    RemoveTwoNodes(model_part);
    p_element = CreatePwLineElementWithoutPWDofs(model_part, p_properties);
    for (auto& node : p_element->GetGeometry()) {
        node.AddDof(WATER_PRESSURE);
    }
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_WATER does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(DENSITY_WATER, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "DENSITY_WATER in the material properties at element with Id "
                                      "4 has an invalid value: -1000 is out of the range [0; -).")

    p_element->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_SOLID does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(DENSITY_SOLID, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "DENSITY_SOLID in the material properties at element with Id "
                                      "4 has an invalid value: -1000 is out of the range [0; -).")

    p_element->GetProperties().SetValue(DENSITY_SOLID, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "POROSITY does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(POROSITY, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "POROSITY in the material properties at element with Id 4 "
                                      "has an invalid value: -1 is out of the range [0; 1.000000].")

    p_element->GetProperties().SetValue(POROSITY, 2.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "POROSITY in the material properties at element with Id 4 "
                                      "has an invalid value: 2 is out of the range [0; 1.000000].")

    p_element->GetProperties().SetValue(POROSITY, 0.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_SOLID does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, -1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_SOLID in the material properties at element with Id 4 has an invalid value: "
        "-1e+06 is out of the range [0; -).")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_FLUID does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, -1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_FLUID in the material properties at element with Id 4 has an invalid value: "
        "-1e+06 is out of the range [0; -).")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DYNAMIC_VISCOSITY does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, -1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DYNAMIC_VISCOSITY in the material properties at element with Id 4 has an invalid value: "
        "-1e+06 is out of the range [0; -).")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BIOT_COEFFICIENT does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "BIOT_COEFFICIENT in the material properties at element with "
                                      "Id 4 has an invalid value: -1 is out of the range [0; -).")

    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XX does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "PERMEABILITY_XX in the material properties at element with "
                                      "Id 4 has an invalid value: -1 is out of the range [0; -).")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, 1.0);
    p_element->GetGeometry().begin()->Z() += 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Node with Id: 0 has non-zero Z coordinate.")
    p_element->GetGeometry().begin()->Z() = 0;

    // to enable a call of RetentionLaw check
    p_element->Initialize(dummy_process_info);

    // No exceptions on correct input
    KRATOS_EXPECT_EQ(p_element->Check(dummy_process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_CreateInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_geometry = std::make_shared<Triangle2D3<Node>>(CreateThreeNodes());
    using enum CalculationContribution;
    const std::vector     contributions = {Permeability, Compressibility, FluidBodyFlow};
    const auto            p_properties  = std::make_shared<Properties>();
    const PwElement<2, 3> element(0, p_geometry, p_properties, contributions, nullptr);

    // Act
    const auto p_created_element = element.Create(1, p_geometry, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_CreateInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    using enum CalculationContribution;
    const std::vector     contributions = {Permeability, Compressibility, FluidBodyFlow};
    const PwElement<2, 3> element(0, std::make_shared<Triangle2D3<Node>>(CreateThreeNodes()),
                                  p_properties, contributions, nullptr);

    // Act
    const auto p_created_element = element.Create(1, CreateThreeNodes(), p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_DoFList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model      model;
    auto&      r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    const auto p_element =
        CreateTransientPwLineElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    p_element->GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(degrees_of_freedom.size(), 3);
    KRATOS_EXPECT_TRUE(std::all_of(degrees_of_freedom.begin(), degrees_of_freedom.end(),
                                   [](auto p_dof) { return p_dof->GetVariable() == WATER_PRESSURE; }))
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_EquationIdVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateTransientPwLineElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());

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

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_IntegrationMethod, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    using enum CalculationContribution;
    const std::vector contributions = {Permeability, Compressibility, FluidBodyFlow};
    const PwElement<2, 3> element(0, std::make_shared<Triangle2D3<Node>>(CreateThreeCoincidentNodes()),
                                  std::make_shared<Properties>(), contributions, nullptr);

    // Act
    const auto p_integration_method = element.GetIntegrationMethod();

    // Assert
    constexpr auto expected_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    KRATOS_EXPECT_EQ(p_integration_method, expected_integration_method);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_CheckThrowsOnFaultyInput2, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    using enum CalculationContribution;
    const std::vector     contributions = {Permeability, Compressibility, FluidBodyFlow};
    const auto            p_properties  = std::make_shared<Properties>();
    const PwElement<2, 3> element_with_coincident_nodes(
        1, std::make_shared<Triangle2D3<Node>>(CreateThreeCoincidentNodes()), p_properties, contributions, nullptr);

    // Act and Assert
    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element_with_coincident_nodes.Check(dummy_process_info),
                                      "Error: DomainSize (0) is smaller than 1e-15 for element 1")

    const PwElement<2, 3> element_with_correct_domain_size(
        1, std::make_shared<Triangle2D3<Node>>(CreateThreeNodes()), p_properties, contributions, nullptr);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element_with_correct_domain_size.Check(dummy_process_info),
                                      "Missing variable WATER_PRESSURE on nodes 1 2 3")

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_element = CreateTriangleTransientPwLineElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable DT_WATER_PRESSURE on nodes 1 2 3")

    DeleteThreeNodes(model_part);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    p_element = CreateTriangleTransientPwLineElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable VOLUME_ACCELERATION on nodes 1 2 3")

    DeleteThreeNodes(model_part);
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    p_element = CreateTriangleTransientPwLineElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "Missing the DoF for the variable WATER_PRESSURE on nodes 1 2 3")

    DeleteThreeNodes(model_part);
    p_element = CreateTransientPwLineElementWithPWDofs<2, 3>(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_WATER does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(DENSITY_WATER, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "DENSITY_WATER in the material properties at element with Id "
                                      "4 has an invalid value: -1000 is out of the range [0; -).")

    p_element->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_SOLID does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(DENSITY_SOLID, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "DENSITY_SOLID in the material properties at element with Id "
                                      "4 has an invalid value: -1000 is out of the range [0; -).")

    p_element->GetProperties().SetValue(DENSITY_SOLID, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "POROSITY does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(POROSITY, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "POROSITY in the material properties at element with Id 4 "
                                      "has an invalid value: -1 is out of the range [0; 1.000000].")

    p_element->GetProperties().SetValue(POROSITY, 2.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "POROSITY in the material properties at element with Id 4 "
                                      "has an invalid value: 2 is out of the range [0; 1.000000].")

    p_element->GetProperties().SetValue(POROSITY, 0.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_SOLID does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, -1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_SOLID in the material properties at element with Id 4 has an invalid value: "
        "-1e+06 is out of the range [0; -).")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_FLUID does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, -1.0e6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_FLUID in the material properties at element with Id 4 has an invalid value: "
        "-1e+06 is out of the range [0; -).")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0e6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DYNAMIC_VISCOSITY does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DYNAMIC_VISCOSITY in the material properties at element with Id 4 has an invalid value: "
        "-0.01 is out of the range [0; -).")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BIOT_COEFFICIENT does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XX does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XX in the material properties at element with Id 4 has an invalid value: "
        "-0.01 is out of the range [0; -).")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_YY does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(PERMEABILITY_YY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_YY in the material properties at element with Id 4 has an invalid value: "
        "-0.01 is out of the range [0; -).")

    p_element->GetProperties().SetValue(PERMEABILITY_YY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XY does not exist in the material properties at element with Id 4.")

    p_element->GetProperties().SetValue(PERMEABILITY_XY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XY in the material properties at element with Id 4 has an invalid value: "
        "-0.01 is out of the range [0; -).")

    p_element->GetProperties().SetValue(PERMEABILITY_XY, 1.0E-2);

    p_element->GetGeometry().begin()->Z() += 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Node with Id: 1 has non-zero Z coordinate.")
    p_element->GetGeometry().begin()->Z() = 0;

    p_element->Initialize(dummy_process_info);
    // No exceptions on correct input for 2D element
    KRATOS_EXPECT_EQ(p_element->Check(dummy_process_info), 0);

    auto p_3D_element = CreateTransientPwLineElementWithPWDofs<3, 4>(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_ZZ does not exist in the material properties at element with Id 5.")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZZ, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_ZZ in the material properties at element with Id 5 has an invalid value: "
        "-0.01 is out of the range [0; -).")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZZ, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_YZ does not exist in the material properties at element with Id 5.")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_YZ, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_YZ in the material properties at element with Id 5 has an invalid value: "
        "-0.01 is out of the range [0; -).")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_YZ, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_ZX does not exist in the material properties at element with Id 5.")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZX, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_3D_element->Check(dummy_process_info),
        "PERMEABILITY_ZX in the material properties at element with Id 5 has an invalid value: "
        "-0.01 is out of the range [0; -).")

    p_3D_element->GetProperties().SetValue(PERMEABILITY_ZX, 1.0E-2);

    // to enable a call of RetentionLaw check
    p_3D_element->Initialize(dummy_process_info);

    // No exceptions on correct input for 3D element
    KRATOS_EXPECT_EQ(p_3D_element->Check(dummy_process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_InitializeSolution, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTransientPwLineElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    // Act
    p_element->InitializeSolutionStep(dummy_process_info);

    // Assert
    KRATOS_EXPECT_TRUE(std::all_of(
        p_element->GetGeometry().begin(), p_element->GetGeometry().end(),
        [](auto& r_node) { return r_node.FastGetSolutionStepValue(HYDRAULIC_DISCHARGE) == 0.0; }))
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_FinalizeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTransientPwLineElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    p_element->InitializeSolutionStep(dummy_process_info);

    // Act
    p_element->FinalizeSolutionStep(dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(p_element->GetGeometry()[0].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 500000);
    KRATOS_EXPECT_EQ(p_element->GetGeometry()[1].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 0);
    KRATOS_EXPECT_EQ(p_element->GetGeometry()[2].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), -500000);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_CalculateOnIntegrationPoints_Vector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTransientPwLineElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
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
    p_element->CalculateOnIntegrationPoints(DYNAMIC_VISCOSITY, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0.01, 0.01, 0.01;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // Act
    results.clear();
    p_element->CalculateOnIntegrationPoints(DT_WATER_PRESSURE, results, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0, 0, 0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_CalculateOnIntegrationPoints_1DArray, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTransientPwLineElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
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

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_CalculateOnIntegrationPoints_Matrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTransientPwLineElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
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

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement2D3N_CalculateLocalSystem, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTransientPwLineElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
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

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement2D3N_Case_A1_2D3N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    PointerVector<Node> nodes;
    // Node 6, 5, 8  from test_transient_groundwater_flow.py: case A1_2D3N
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0100000000, 1.9600000000, 0.0000000000));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0000000000, 1.9600000000, 0.0000000000));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0100000000, 1.9400000000, 0.0000000000));

    using enum CalculationContribution;
    const std::vector contributions = {Permeability, Compressibility, FluidBodyFlow};
    auto              properties    = std::make_shared<Properties>();
    properties->SetValue(IGNORE_UNDRAINED, false);
    properties->SetValue(YOUNG_MODULUS, 10000);
    properties->SetValue(POISSON_RATIO, 0.2);
    properties->SetValue(DENSITY_SOLID, 2.65);
    properties->SetValue(DENSITY_WATER, 1.0);
    properties->SetValue(POROSITY, 0.36);
    properties->SetValue(BULK_MODULUS_SOLID, 1.0e9);
    properties->SetValue(BULK_MODULUS_FLUID, 175500);
    properties->SetValue(PERMEABILITY_XX, 0.1521);
    properties->SetValue(PERMEABILITY_YY, 0.1521);
    properties->SetValue(PERMEABILITY_XY, 0.0);
    properties->SetValue(DYNAMIC_VISCOSITY, 9.81);
    properties->SetValue(THICKNESS, 1.0);
    properties->SetValue(BIOT_COEFFICIENT, 1.0);
    properties->SetValue(RETENTION_LAW, "VanGenuchtenLaw");
    properties->SetValue(SATURATED_SATURATION, 1.0);
    properties->SetValue(RESIDUAL_SATURATION, 0.06203);
    properties->SetValue(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE, 4.379464286);
    properties->SetValue(VAN_GENUCHTEN_GN, 2.286);
    properties->SetValue(VAN_GENUCHTEN_GL, 0);
    properties->SetValue(MINIMUM_RELATIVE_PERMEABILITY, 0.0001);
    auto process_info                     = ProcessInfo{};
    process_info[DT_PRESSURE_COEFFICIENT] = 250.0;

    auto element =
        PwElement<2, 3>(1, std::make_shared<Triangle2D3<Node>>(nodes), properties, contributions, nullptr);
    element.GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -9.81, 0.0};
    element.GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -9.81, 0.0};
    element.GetGeometry()[2].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -9.81, 0.0};
    element.GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE)    = 10.9671575;
    element.GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE)    = 10.2635811;
    element.GetGeometry()[2].FastGetSolutionStepValue(WATER_PRESSURE)    = 9.736040251;
    element.GetGeometry()[0].FastGetSolutionStepValue(DT_WATER_PRESSURE) = 289.2893748;
    element.GetGeometry()[1].FastGetSolutionStepValue(DT_WATER_PRESSURE) = 113.395276;
    element.GetGeometry()[2].FastGetSolutionStepValue(DT_WATER_PRESSURE) = -18.48993726;
    element.Initialize(process_info);
    element.InitializeSolutionStep(process_info);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, process_info);

    // Assert
    Matrix expected_left_hand_side(3, 3);
    // clang-format off
    expected_left_hand_side <<= -0.000144736222,5.634212252e-05,-3.774083572e-06,
                                 5.634212252e-05,-0.000127227034,-2.416421349e-05,
                                 -3.774083572e-06,-2.416421349e-05,-6.942855172e-05;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    Vector expected_right_hand_side(3);
    expected_right_hand_side <<= 0.0001404394389, -9.276110236e-06, 1.140103478e-05;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)

    // Copy of TransientPwLineElement_CalculateOnIntegrationPoints_Vector
    // Act
    std::vector<double> results{};
    element.CalculateOnIntegrationPoints(DEGREE_OF_SATURATION, results, process_info);
    auto const number_of_integration_points =
        element.GetGeometry().IntegrationPointsNumber(element.GetIntegrationMethod());

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    Vector expected_results(3);
    expected_results <<= 0.341298621521, 0.352123307571, 0.360698051689;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(results, expected_results, Defaults::relative_tolerance);

    // Act
    results.clear();
    element.CalculateOnIntegrationPoints(EFFECTIVE_SATURATION, results, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0.297737264007, 0.309277810133, 0.318419620765;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(results, expected_results, Defaults::relative_tolerance);

    // Act
    results.clear();
    element.CalculateOnIntegrationPoints(BISHOP_COEFFICIENT, results, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0.297737264007, 0.309277810133, 0.318419620765;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(results, expected_results, Defaults::relative_tolerance);

    // Act
    results.clear();
    element.CalculateOnIntegrationPoints(DERIVATIVE_OF_SATURATION, results, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= -0.029823048151, -0.031743614419, -0.033288714765;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(results, expected_results, Defaults::relative_tolerance);

    // Act
    results.clear();
    element.CalculateOnIntegrationPoints(RELATIVE_PERMEABILITY, results, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0.004495279802, 0.005166152553, 0.005748179834;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(results, expected_results, Defaults::relative_tolerance);

    // Act
    results.clear();
    element.CalculateOnIntegrationPoints(HYDRAULIC_HEAD, results, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0.871579, 0.907439, 0.924327;
    KRATOS_EXPECT_VECTOR_NEAR(results, expected_results, Defaults::relative_tolerance);

    // Act
    results.clear();
    element.CalculateOnIntegrationPoints(DT_WATER_PRESSURE, results, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    expected_results <<= 0, 0, 0;
    KRATOS_EXPECT_VECTOR_EQ(results, expected_results);

    // copy of CalculateOnIntegrationPoints_1DArray
    // Act
    std::vector<array_1d<double, 3>> results_1d{};
    element.CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, results_1d, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    array_1d<double, 3> expected_nonzero_component_1d{0.004903748623, 0.003606555048, 0.000000000000};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(results_1d[0], expected_nonzero_component_1d, Defaults::relative_tolerance);
    expected_nonzero_component_1d[0] = 0.005635581005;
    expected_nonzero_component_1d[1] = 0.004144795072;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(results_1d[1], expected_nonzero_component_1d, Defaults::relative_tolerance);
    expected_nonzero_component_1d[0] = 0.006270494871;
    expected_nonzero_component_1d[1] = 0.004611754531;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(results_1d[2], expected_nonzero_component_1d, Defaults::relative_tolerance);

    // Act
    results_1d.clear();
    element.CalculateOnIntegrationPoints(LOCAL_FLUID_FLUX_VECTOR, results_1d, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results.size(), number_of_integration_points);
    array_1d<double, 3> expected_zero_component_1d{0, 0, 0};
    for (const auto& component : results_1d) {
        KRATOS_EXPECT_VECTOR_EQ(component, expected_zero_component_1d);
    }

    // copy of CalculateOnIntegrationPoints_Matrix
    // Act
    std::vector<Matrix> results_m{};
    element.CalculateOnIntegrationPoints(PERMEABILITY_MATRIX, results_m, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results_m.size(), number_of_integration_points);
    Matrix expected_nonzero_component_m(2, 2);
    expected_nonzero_component_m <<= 0.1521, 0.0, 0.0, 0.1521;
    for (const auto& component : results_m) {
        KRATOS_EXPECT_MATRIX_EQ(component, expected_nonzero_component_m);
    }

    // Act
    results_m.clear();
    element.CalculateOnIntegrationPoints(LOCAL_PERMEABILITY_MATRIX, results_m, process_info);

    // Assert
    KRATOS_EXPECT_EQ(results_m.size(), number_of_integration_points);
    Matrix expected_zero_component_m = ZeroMatrix(2, 2);
    for (const auto& component : results_m) {
        KRATOS_EXPECT_MATRIX_EQ(component, expected_zero_component_m);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement3D4N_CalculateLocalSystem, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTransientPwLineElementWithPWDofs<3, 4>(r_model_part, std::make_shared<Properties>());
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

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_GetIntegrationMethodForAllRegisteredElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::vector<CalculationContribution> contributions;
    PointerVector<Node>                        nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));

    // Act and Assert
    auto p_transient_pw_line_element_2D2N = make_intrusive<PwElement<2, 2>>(
        1, std::make_shared<Line2D2<Node>>(nodes), nullptr, contributions,
        std::make_unique<IntegrationCoefficientModifierForLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_2D2N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    auto p_transient_pw_line_element_3D2N = make_intrusive<PwElement<3, 2>>(
        1, std::make_shared<Line3D2<Node>>(nodes), nullptr, contributions,
        std::make_unique<IntegrationCoefficientModifierForLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_3D2N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    auto p_transient_pw_line_element_2D3N = make_intrusive<PwElement<2, 3>>(
        1, std::make_shared<Line2D3<Node>>(nodes), nullptr, contributions,
        std::make_unique<IntegrationCoefficientModifierForLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_2D3N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    auto p_transient_pw_line_element_3D3N = make_intrusive<PwElement<3, 3>>(
        1, std::make_shared<Line3D3<Node>>(nodes), nullptr, contributions,
        std::make_unique<IntegrationCoefficientModifierForLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_3D3N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(4, 0.5, 0.0, 0.0));
    auto p_transient_pw_line_element_2D4N = make_intrusive<PwElement<2, 4>>(
        1, std::make_shared<Line2D4<Node>>(nodes), nullptr, contributions,
        std::make_unique<IntegrationCoefficientModifierForLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_2D4N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_3);

    nodes.push_back(make_intrusive<Node>(5, 1.0, 0.5, 0.0));
    auto p_transient_pw_line_element_2D5N = make_intrusive<PwElement<2, 5>>(
        1, std::make_shared<Line2D5<Node>>(nodes), nullptr, contributions,
        std::make_unique<IntegrationCoefficientModifierForLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_2D5N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_5);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement2D3N_CalculateLocalSystemWithPressureFilterLaw,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    using enum CalculationContribution;
    const std::vector contributions = {Permeability, Compressibility, FluidBodyFlow};
    auto              p_element     = make_intrusive<PwElement<2, 3>>(
        GetNextElementNumber(r_model_part), std::make_shared<Line2D3<Node>>(CreateInclinedLine(r_model_part)),
        std::make_shared<Properties>(), contributions, nullptr);
    SetBasicPropertiesAndVariables(p_element);
    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0E20);
    p_element->GetProperties().SetValue(RETENTION_LAW, "PressureFilterLaw");
    p_element->GetProperties().SetValue(CROSS_AREA, 1.0);
    p_element->GetProperties().SetValue(FILTER_LENGTH, 3.0);
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
    expected_left_hand_side <<= -6.5648064719788994, -0.93782949599698517, 7.5026359679758832,
                                -0.93782949599698517,-6.5648064719788994,  7.5026359679758832,
                                 7.5026359679758832,  7.5026359679758832,-15.005271935951766;
    // clang-format on
    KRATOS_EXPECT_EQ(actual_left_hand_side.size1(), expected_left_hand_side.size1());
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    Vector expected_right_hand_side(3);
    expected_right_hand_side <<= -28134.9, 28134.9, 0;
    KRATOS_EXPECT_EQ(actual_right_hand_side.size(), expected_right_hand_side.size());
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement2D3N_SaveLoad, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element =
        CreateTransientPwLineElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 0.5);
    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0E6);
    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    p_element->GetProperties().SetValue(POROSITY, 0.1);
    p_element->GetProperties().SetValue(IGNORE_UNDRAINED, false);
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    p_element->InitializeSolutionStep(dummy_process_info);
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    p_element->CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, dummy_process_info);

    const auto triangle_2D3N       = static_cast<Triangle2D3<Node>>(p_element->GetGeometry());
    const auto scoped_registration = ScopedSerializerRegistration{
        std::make_pair(std::string{p_element->GetGeometry().Name()}, triangle_2D3N),
        std::make_pair("SaturatedLaw"s, SaturatedLaw())};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_element);
    auto p_loaded_element = make_intrusive<PwElement<2, 3>>();
    // the following Add actions are needed to load the saved element
    KratosComponents<VariableData>::Add("WATER_PRESSURE", WATER_PRESSURE);
    KratosComponents<VariableData>::Add("DT_WATER_PRESSURE", DT_WATER_PRESSURE);
    KratosComponents<VariableData>::Add("VOLUME_ACCELERATION", VOLUME_ACCELERATION);
    KratosComponents<VariableData>::Add("HYDRAULIC_DISCHARGE", HYDRAULIC_DISCHARGE);
    KratosComponents<VariableData>::Add("DENSITY_WATER", DENSITY_WATER);
    KratosComponents<VariableData>::Add("DYNAMIC_VISCOSITY", DYNAMIC_VISCOSITY);
    KratosComponents<VariableData>::Add("PERMEABILITY_XX", PERMEABILITY_XX);
    KratosComponents<VariableData>::Add("PERMEABILITY_YY", PERMEABILITY_YY);
    KratosComponents<VariableData>::Add("PERMEABILITY_XY", PERMEABILITY_XY);
    KratosComponents<VariableData>::Add("PERMEABILITY_ZZ", PERMEABILITY_ZZ);
    KratosComponents<VariableData>::Add("PERMEABILITY_YZ", PERMEABILITY_YZ);
    KratosComponents<VariableData>::Add("PERMEABILITY_ZX", PERMEABILITY_ZX);
    KratosComponents<VariableData>::Add("BIOT_COEFFICIENT", BIOT_COEFFICIENT);
    KratosComponents<VariableData>::Add("BULK_MODULUS_FLUID", BULK_MODULUS_FLUID);
    KratosComponents<VariableData>::Add("BULK_MODULUS_SOLID", BULK_MODULUS_SOLID);
    KratosComponents<VariableData>::Add("POROSITY", POROSITY);
    KratosComponents<VariableData>::Add("IGNORE_UNDRAINED", IGNORE_UNDRAINED);
    KratosComponents<VariableData>::Add("RETENTION_LAW", RETENTION_LAW);

    serializer.load("test_tag"s, p_loaded_element);

    // Assert
    Vector loaded_right_hand_side;
    Matrix loaded_left_hand_side;
    p_loaded_element->CalculateLocalSystem(loaded_left_hand_side, loaded_right_hand_side, dummy_process_info);

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, loaded_left_hand_side, Defaults::relative_tolerance)

    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, loaded_right_hand_side, Defaults::relative_tolerance)
}
} // namespace Kratos::Testing