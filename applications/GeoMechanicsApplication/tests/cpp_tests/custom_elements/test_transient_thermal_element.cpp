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

#include "containers/model.h"
#include "custom_elements/transient_thermal_element.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_15.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

void SetProperties(Element::Pointer p_element)
{
    Properties::Pointer p_properties = p_element->pGetProperties();

    // Please note these are not representative values, it just ensures the values are set
    p_properties->SetValue(DENSITY_WATER, 1.0);
    p_properties->SetValue(POROSITY, 0.0);
    p_properties->SetValue(SATURATED_SATURATION, 1.0);
    p_properties->SetValue(RETENTION_LAW, "SaturatedLaw");
    p_properties->SetValue(DENSITY_SOLID, 1.0);
    p_properties->SetValue(SPECIFIC_HEAT_CAPACITY_WATER, 1.0);
    p_properties->SetValue(SPECIFIC_HEAT_CAPACITY_SOLID, 1.0);
    p_properties->SetValue(THERMAL_CONDUCTIVITY_WATER, 1.0);

    // This ensures the lhs matrix has differences for x/y
    p_properties->SetValue(THERMAL_CONDUCTIVITY_SOLID_XX, 10.0);
    p_properties->SetValue(THERMAL_CONDUCTIVITY_SOLID_YY, 1.0);
    p_properties->SetValue(THERMAL_CONDUCTIVITY_SOLID_XY, 0.0);
}

void SetupElement(ModelPart& rModelPart)
{
    Element::Pointer        p_element;
    Element::DofsVectorType elemental_dofs;

    p_element                                 = rModelPart.pGetElement(1);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    auto               element_geometry       = p_element->GetGeometry();
    for (unsigned int i = 0; i < element_geometry.size(); i++) {
        element_geometry[i].AddDof(TEMPERATURE);
        // Set temperature to some values to get a non-zero right hand side
        element_geometry[i].GetSolutionStepValue(TEMPERATURE) = i * 1.0;
    }

    p_element->GetDofList(elemental_dofs, r_current_process_info);
    for (unsigned int i = 0; i < elemental_dofs.size(); i++) {
        elemental_dofs[i]->SetEquationId(i);
    }
    SetProperties(p_element);
}

void ValidateThermalElement(ModelPart& rModelPart)
{
    auto               p_element              = rModelPart.pGetElement(1);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    Element::DofsVectorType elemental_dofs;
    p_element->GetDofList(elemental_dofs, r_current_process_info);

    KRATOS_EXPECT_EQ(elemental_dofs.size(), p_element->GetGeometry().size());
    for (const auto& element_dof : elemental_dofs) {
        // Only the TEMPERATURE dofs should be returned by the element
        KRATOS_EXPECT_EQ(element_dof->GetVariable(), TEMPERATURE);
    }

    Element::EquationIdVectorType equation_ids;
    p_element->EquationIdVector(equation_ids, r_current_process_info);
    for (unsigned int i = 0; i < equation_ids.size(); i++) {
        KRATOS_EXPECT_EQ(equation_ids[i], i);
    }
}

Element::Pointer GenerateTransientThermalElementWithZeroDomainSize()
{
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    auto transient_thermal_element =
        make_intrusive<TransientThermalElement<2, 3>>(1, Kratos::make_shared<Triangle2D3<Node>>(nodes));

    return transient_thermal_element;
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenDomainSizeIsInvalid, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Element::Pointer p_element = GenerateTransientThermalElementWithZeroDomainSize();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(ProcessInfo()),
                                      "DomainSize (0) is smaller than 1e-15 for element 1")
}

void GenerateTransientThermalElement2D3N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    auto transient_thermal_element = make_intrusive<TransientThermalElement<2, 3>>(
        1, Kratos::make_shared<Triangle2D3<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenTemperatureIsMissing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    GenerateTransientThermalElement2D3N(model_part);

    Element::Pointer   p_element              = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(r_current_process_info),
                                      "Missing variable TEMPERATURE on nodes 1 2 3")
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenDtTemperatureIsMissing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    GenerateTransientThermalElement2D3N(model_part);

    Element::Pointer   p_element              = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(r_current_process_info),
                                      "Missing variable DT_TEMPERATURE on nodes 1 2 3")
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenTemperatureDofIsMissing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement2D3N(model_part);

    Element::Pointer   p_element              = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(r_current_process_info),
                                      "Missing the DoF for the variable TEMPERATURE on nodes 1 2 3")
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenPropertyIsMissing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement2D3N(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    for (auto& node : p_element->GetGeometry()) {
        node.AddDof(TEMPERATURE);
    }
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(r_current_process_info),
        "DENSITY_WATER does not exist in the properties with Id 0 at element with Id 1.")
}

void GenerateTransientThermalElement2D3NWithNonZeroZ(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 1.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, -1.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    auto transient_thermal_element = make_intrusive<TransientThermalElement<2, 3>>(
        1, Kratos::make_shared<Triangle2D3<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_When2DElementHasNonZeroZValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement2D3NWithNonZeroZ(model_part);
    SetupElement(model_part);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(model_part.GetElement(1).Check(r_current_process_info),
                                      "Node with Id: 1 has non-zero Z coordinate.")
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Returns0_When2DElementIsCorrectlySet, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement2D3N(model_part);
    SetupElement(model_part);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    KRATOS_EXPECT_EQ(model_part.GetElement(1).Check(r_current_process_info), 0);
}

void ExpectDoubleMatrixEqual(const Matrix& rExpectedMatrix, const Matrix& rActualMatrix)
{
    KRATOS_EXPECT_EQ(rActualMatrix.size1(), rExpectedMatrix.size1());
    KRATOS_EXPECT_EQ(rActualMatrix.size2(), rExpectedMatrix.size2());
    for (std::size_t i = 0; i < rActualMatrix.size1(); i++) {
        for (std::size_t j = 0; j < rActualMatrix.size2(); j++) {
            KRATOS_EXPECT_DOUBLE_EQ(rActualMatrix(i, j), rExpectedMatrix(i, j));
        }
    }
}

void ExpectDoubleVectorEqual(const Vector& rExpectedVector, const Vector& rActualVector)
{
    KRATOS_EXPECT_EQ(rActualVector.size(), rExpectedVector.size());
    for (std::size_t i = 0; i < rActualVector.size(); i++) {
        KRATOS_EXPECT_DOUBLE_EQ(rActualVector[i], rExpectedVector[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ThermalElement_ReturnsExpectedMatrixAndVector_WhenCalculateLocalSystem,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement2D3N(model_part);
    SetupElement(model_part);

    Element::Pointer   p_element              = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    Matrix left_hand_side_matrix;
    Vector right_hand_side;
    p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side, r_current_process_info);

    Matrix expected_matrix(3, 3);
    // clang-format off
    expected_matrix <<=  5, -5,    0,
                        -5,  5.5, -0.5,
                         0, -0.5,  0.5;
    // clang-format on
    ExpectDoubleMatrixEqual(expected_matrix, left_hand_side_matrix);

    // Calculated by hand (matrix multiplication between the
    // lhs and the temperature vector, which is 0.0, 1.0, 2.0)
    Vector expected_vector(3);
    expected_vector <<= 5, -4.5, -0.5;
    ExpectDoubleVectorEqual(expected_vector, right_hand_side);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D3N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D3N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement2D4N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.6666, 0.3333, 0.0));
    auto transient_thermal_element = make_intrusive<TransientThermalElement<2, 4>>(
        1, Kratos::make_shared<Quadrilateral2D4<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D4N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D4N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement2D6N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<2, 6>>(
        1, Kratos::make_shared<Triangle2D6<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D6N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D6N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement2D8N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<2, 8>>(
        1, Kratos::make_shared<Quadrilateral2D8<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D8N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D8N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement2D9N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<2, 9>>(
        1, Kratos::make_shared<Quadrilateral2D9<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D9N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D9N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement2D10N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<2, 10>>(
        1, Kratos::make_shared<Triangle2D10<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D10N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D10N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement2D15N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(11, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(12, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(13, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(14, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(15, 0.5, 0.5, 0.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<2, 15>>(
        1, Kratos::make_shared<Triangle2D15<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D15N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D15N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement3D4N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    // The z coordinate is non-zero to ensure the element has a non-zero domain size
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 1.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<3, 4>>(
        1, Kratos::make_shared<Tetrahedra3D4<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_When3DPropertyHasInvalidValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement3D4N(model_part);
    SetupElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    p_element->GetProperties().SetValue(THERMAL_CONDUCTIVITY_SOLID_ZZ, -5.0);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(r_current_process_info),
        "THERMAL_CONDUCTIVITY_SOLID_ZZ in the properties with Id 0 at element with Id 1 has an "
        "invalid value: -5 is out of the range [0, -).")
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D4N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D4N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement3D8N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<3, 8>>(
        1, Kratos::make_shared<Hexahedra3D8<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D8N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D8N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement3D10N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<3, 10>>(
        1, Kratos::make_shared<Tetrahedra3D10<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D10N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D10N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement3D20N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(11, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(12, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(13, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(14, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(15, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(16, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(17, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(18, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(19, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(20, 0.5, 0.5, 0.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<3, 20>>(
        1, Kratos::make_shared<Hexahedra3D20<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D20N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D20N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

void GenerateTransientThermalElement3D27N(ModelPart& rModelPart)
{
    // Geometry creation
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(11, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(12, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(13, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(14, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(15, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(16, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(17, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(18, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(19, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(20, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(21, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(22, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(23, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(24, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(25, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(26, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(27, 0.5, 0.5, 0.0));

    auto transient_thermal_element = make_intrusive<TransientThermalElement<3, 27>>(
        1, Kratos::make_shared<Hexahedra3D27<Node>>(nodes), rModelPart.CreateNewProperties(0));
    rModelPart.AddElement(transient_thermal_element);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D27N, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D27N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

KRATOS_TEST_CASE_IN_SUITE(TransientThermalElement_GetIntegrationMethodForAllRegisteredElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto          p_properties = std::make_shared<Properties>();
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));

    // Act and Assert
    auto p_transient_thermal_line_element_2D2N = make_intrusive<TransientThermalElement<2, 2>>(
        1, std::make_shared<Line2D2<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_line_element_2D2N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    auto p_transient_thermal_line_element_3D2N = make_intrusive<TransientThermalElement<3, 2>>(
        1, std::make_shared<Line3D2<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_line_element_3D2N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    auto p_transient_thermal_line_element_2D3N = make_intrusive<TransientThermalElement<2, 3>>(
        1, std::make_shared<Line2D3<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_line_element_2D3N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    auto p_transient_thermal_line_element_3D3N = make_intrusive<TransientThermalElement<3, 3>>(
        1, std::make_shared<Line3D3<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_line_element_3D3N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    auto p_transient_thermal_element_2D3N = make_intrusive<TransientThermalElement<2, 3>>(
        1, std::make_shared<Triangle2D3<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_2D3N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(4, 0.5, 0.0, 0.0));
    auto p_transient_thermal_line_element_2D4N = make_intrusive<TransientThermalElement<2, 4>>(
        1, std::make_shared<Line2D4<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_line_element_2D4N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_3);

    auto p_transient_thermal_element_2D4N = make_intrusive<TransientThermalElement<2, 4>>(
        1, std::make_shared<Quadrilateral2D4<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_2D4N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    auto p_transient_thermal_element_3D4N = make_intrusive<TransientThermalElement<3, 4>>(
        1, std::make_shared<Tetrahedra3D4<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_3D4N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(5, 1.0, 0.5, 0.0));
    auto p_transient_thermal_line_element_2D5N = make_intrusive<TransientThermalElement<2, 5>>(
        1, std::make_shared<Line2D5<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_line_element_2D5N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_5);

    nodes.push_back(make_intrusive<Node>(6, 0.5, 0.5, 0.0));
    auto p_transient_thermal_element_2D6N = make_intrusive<TransientThermalElement<2, 6>>(
        1, std::make_shared<Triangle2D6<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_2D6N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(7, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(8, 0.5, 0.5, 0.0));
    auto p_transient_thermal_element_2D8N = make_intrusive<TransientThermalElement<2, 8>>(
        1, std::make_shared<Quadrilateral2D8<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_2D8N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    auto p_transient_thermal_element_3D8N = make_intrusive<TransientThermalElement<3, 8>>(
        1, std::make_shared<Hexahedra3D8<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_3D8N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(9, 0.5, 0.5, 0.0));
    auto p_transient_thermal_element_2D9N = make_intrusive<TransientThermalElement<2, 9>>(
        1, std::make_shared<Quadrilateral2D9<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_2D9N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(10, 0.5, 0.5, 0.0));
    auto p_transient_thermal_element_2D10N = make_intrusive<TransientThermalElement<2, 10>>(
        1, std::make_shared<Triangle2D10<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_2D10N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_4);

    auto p_transient_thermal_element_3D10N = make_intrusive<TransientThermalElement<3, 10>>(
        1, std::make_shared<Tetrahedra3D10<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_3D10N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(11, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(12, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(13, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(14, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(15, 0.5, 0.5, 0.0));

    auto p_transient_thermal_element_2D15N = make_intrusive<TransientThermalElement<2, 15>>(
        1, std::make_shared<Triangle2D15<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_2D15N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_5);

    nodes.push_back(make_intrusive<Node>(16, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(17, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(18, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(19, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(20, 0.5, 0.5, 0.0));

    auto p_transient_thermal_element_3D20N = make_intrusive<TransientThermalElement<3, 20>>(
        1, std::make_shared<Hexahedra3D20<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_3D20N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(21, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(22, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(23, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(24, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(25, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(26, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(27, 0.5, 0.5, 0.0));

    auto p_transient_thermal_element_3D27N = make_intrusive<TransientThermalElement<3, 27>>(
        1, std::make_shared<Hexahedra3D27<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_transient_thermal_element_3D27N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);
}

} // namespace Kratos::Testing
