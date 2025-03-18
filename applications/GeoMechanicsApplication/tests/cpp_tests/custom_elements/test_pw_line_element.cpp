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

#include "custom_elements/calculation_contribution.h"
#include "custom_elements/transient_Pw_line_element.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

ModelPart& CreateModelPartWithSolutionStepVariables(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    return r_result;
}

TransientPwLineElement<2, 2> TransientPwLineElementWithoutPWDofs(const Properties::Pointer& rProperties,
                                                                 const Geometry<Node>::Pointer& rGeometry)
{
    auto result = TransientPwLineElement<2, 2>{
        1, rGeometry, rProperties, {CalculationContribution::Permeability, CalculationContribution::Compressibility, CalculationContribution::FluidBodyFlow}};

    return result;
}

TransientPwLineElement<2, 2> TransientPwLineElementWithPWDofs(const Properties::Pointer& rProperties,
                                                              const Geometry<Node>::Pointer& rGeometry)
{
    auto result = TransientPwLineElementWithoutPWDofs(rProperties, rGeometry);
    for (auto& node : result.GetGeometry()) {
        node.AddDof(WATER_PRESSURE);
        node.AddDof(DT_WATER_PRESSURE);
    }

    return result;
}

TransientPwLineElement<2, 2> TransientPwLineElementWithPWDofs(Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithSolutionStepVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 1.0, 1.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    return TransientPwLineElementWithPWDofs(rProperties, p_geometry);
}

intrusive_ptr<Element> CreatePwLineElementWithoutPWDofs(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(1, 1.0, 1.0, 0.0));
    const auto                                 p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    const std::vector<CalculationContribution> contirbutions = {
        CalculationContribution::Permeability, CalculationContribution::Compressibility, CalculationContribution::FluidBodyFlow};
    auto element = make_intrusive<TransientPwLineElement<2, 2>>(
        rModelPart.NumberOfElements() + 1, p_geometry, rProperties, contirbutions);
    rModelPart.AddElement(element);
    return element;
}

void RemoveTwoNodes(ModelPart& rModelPart)
{
    rModelPart.RemoveNodeFromAllLevels(0);
    rModelPart.RemoveNodeFromAllLevels(1);
}

} // namespace

namespace Kratos::Testing
{
using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElementReturnsTheExpectedLeftHandSideAndRightHandSide,
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
    auto  element = TransientPwLineElementWithPWDofs(model, p_properties);
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

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElementCheckThrowsOnFaultyInput, KratosGeoMechanicsFastSuiteWithoutKernel)
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
                                      "Error: Length smaller than 1e-15 for element 1")
    nodes.erase(nodes.begin() + 1);
    nodes.push_back(make_intrusive<Node>(1, 1.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    const auto element_with_correct_domain_size = TransientPwLineElementWithoutPWDofs(p_properties, p_geometry);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element_with_correct_domain_size.Check(dummy_process_info),
                                      "Error: Missing variable WATER_PRESSURE on node 0")

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_element = CreatePwLineElementWithoutPWDofs(model_part, p_properties);

    RemoveTwoNodes(model_part);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    p_element = CreatePwLineElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable VOLUME_ACCELERATION on node 0")

    RemoveTwoNodes(model_part);
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    p_element = CreatePwLineElementWithoutPWDofs(model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing degree of freedom for WATER_PRESSURE on node 0")

    RemoveTwoNodes(model_part);
    p_element = CreatePwLineElementWithoutPWDofs(model_part, p_properties);
    for (auto& node : p_element->GetGeometry()) {
        node.AddDof(WATER_PRESSURE);
    }
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_WATER does not exist in the material properties (Id = 0) at element 4")

    p_element->GetProperties().SetValue(DENSITY_WATER, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_WATER of material Id = 0 at element 4 has an invalid value -1000 which is below "
        "the minimum allowed value of 0")

    p_element->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_SOLID does not exist in the material properties (Id = 0) at element 4")

    p_element->GetProperties().SetValue(DENSITY_SOLID, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_SOLID of material Id = 0 at element 4 has an invalid value -1000 which is below "
        "the minimum allowed value of 0")

    p_element->GetProperties().SetValue(DENSITY_SOLID, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "POROSITY does not exist in the material properties (Id = 0) at element 4")

    p_element->GetProperties().SetValue(POROSITY, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "POROSITY of material Id = 0 at element 4 has an invalid "
                                      "value -1 which is outside of the range [ 0, 1]")

    p_element->GetProperties().SetValue(POROSITY, 2.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "POROSITY of material Id = 0 at element 4 has an invalid "
                                      "value 2 which is outside of the range [ 0, 1]")

    p_element->GetProperties().SetValue(POROSITY, 0.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_SOLID does not exist in the material properties (Id = 0) at element 4")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, -1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_SOLID of material Id = 0 at element 4 has an invalid value -1e+06 which is "
        "below the minimum allowed value of 0")

    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_FLUID does not exist in the material properties (Id = 0) at element 4")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, -1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BULK_MODULUS_FLUID of material Id = 0 at element 4 has an invalid value -1e+06 which is "
        "below the minimum allowed value of 0")

    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DYNAMIC_VISCOSITY does not exist in the material properties (Id = 0) at element 4")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, -1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DYNAMIC_VISCOSITY of material Id = 0 at element 4 has an invalid value -1e+06 which is "
        "below the minimum allowed value of 0")

    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E6);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BIOT_COEFFICIENT does not exist in the material properties (Id = 0) at element 4")

    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "BIOT_COEFFICIENT of material Id = 0 at element 4 has an invalid value -1 which is below "
        "the minimum allowed value of 0")

    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XX does not exist in the material properties (Id = 0) at element 4")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XX of material Id = 0 at element 4 has an invalid value -1 which is below "
        "the minimum allowed value of 0")

    p_element->GetProperties().SetValue(PERMEABILITY_XX, 1.0);
    p_element->GetGeometry().begin()->Z() += 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Node with non-zero Z coordinate found. Id: 0")
    p_element->GetGeometry().begin()->Z() = 0;

    // to enable a call of RetentionLaw check
    p_element->Initialize(dummy_process_info);

    // No exceptions on correct input
    KRATOS_EXPECT_EQ(p_element->Check(dummy_process_info), 0);
}

} // namespace Kratos::Testing