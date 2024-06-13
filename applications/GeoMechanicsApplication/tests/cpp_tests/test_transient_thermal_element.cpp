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
#include "testing/testing.h"
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
    Element::Pointer p_element;
    Element::DofsVectorType elemental_dofs;

    p_element = rModelPart.pGetElement(1);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    auto element_geometry = p_element->GetGeometry();
    for (unsigned int i = 0; i < element_geometry.size(); i++)
    {
        element_geometry[i].AddDof(TEMPERATURE);
        // Set temperature to some values to get a non-zero right hand side
        element_geometry[i].GetSolutionStepValue(TEMPERATURE) = i * 1.0;
    }

    p_element->GetDofList(elemental_dofs, r_current_process_info);
    for (unsigned int i = 0; i < elemental_dofs.size(); i++)
    {
        elemental_dofs[i]->SetEquationId(i);
    }
    SetProperties(p_element);
}

void ValidateThermalElement(ModelPart& rModelPart)
{
    auto p_element = rModelPart.pGetElement(1);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    Element::DofsVectorType elemental_dofs;
    p_element->GetDofList(elemental_dofs, r_current_process_info);

    KRATOS_EXPECT_EQ(elemental_dofs.size(), p_element->GetGeometry().size());
    for (const auto& element_dof : elemental_dofs)
    {
        // Only the TEMPERATURE dofs should be returned by the element
        KRATOS_EXPECT_EQ(element_dof->GetVariable(), TEMPERATURE);
    }

    Element::EquationIdVectorType equation_ids;
    p_element->EquationIdVector(equation_ids, r_current_process_info);
    for (unsigned int i = 0; i < equation_ids.size(); i++)
    {
        KRATOS_EXPECT_EQ(equation_ids[i], i);
    }
}

void GenerateTransientThermalElementWithZeroDomainSize(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 0.0, 0.0);
    std::vector<ModelPart::IndexType> node_ids{1, 2, 3};
    rModelPart.CreateNewElement("GeoTransientThermalElement2D3N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenDomainSizeIsInvalid, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransientThermalElementWithZeroDomainSize(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(r_current_process_info),
        "DomainSize smaller than 1e-15 for element 1")
}

void GenerateTransientThermalElement2D3N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> node_ids{1, 2, 3};
    rModelPart.CreateNewElement("GeoTransientThermalElement2D3N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenTemperatureIsMissing, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    GenerateTransientThermalElement2D3N(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(r_current_process_info),
                                      "Missing variable TEMPERATURE on node 1")
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenDtTemperatureIsMissing, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    GenerateTransientThermalElement2D3N(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(r_current_process_info),
        "Missing variable DT_TEMPERATURE on node 1")
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenTemperatureDofIsMissing, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement2D3N(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(r_current_process_info),
        "Missing degree of freedom for TEMPERATURE on node 1")
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_WhenPropertyIsMissing, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement2D3N(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    for (auto& node : p_element->GetGeometry())
    {
        node.AddDof(TEMPERATURE);
    }
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(r_current_process_info),
        "DENSITY_WATER does not exist in the thermal element's properties")
}

void GenerateTransientThermalElement2D3NWithNonZeroZ(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 1.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, -1.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> node_ids{1, 2, 3};
    rModelPart.CreateNewElement("GeoTransientThermalElement2D3N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_When2DElementHasNonZeroZValue, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement2D3NWithNonZeroZ(model_part);
    SetupElement(model_part);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        model_part.GetElement(1).Check(r_current_process_info),
        "Node with non-zero Z coordinate found. Id: 1")
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Returns0_When2DElementIsCorrectlySet, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    for (std::size_t i = 0; i < rActualMatrix.size1(); i++)
    {
        for (std::size_t j = 0; j < rActualMatrix.size2(); j++)
        {
            KRATOS_EXPECT_DOUBLE_EQ(rActualMatrix(i, j), rExpectedMatrix(i, j));
        }
    }
}

void ExpectDoubleVectorEqual(const Vector& rExpectedVector, const Vector& rActualVector)
{
    KRATOS_EXPECT_EQ(rActualVector.size(), rExpectedVector.size());
    for (std::size_t i = 0; i < rActualVector.size(); i++)
    {
        KRATOS_EXPECT_DOUBLE_EQ(rActualVector[i], rExpectedVector[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ThermalElement_ReturnsExpectedMatrixAndVector_WhenCalculateLocalSystem,
                          KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    GenerateTransientThermalElement2D3N(model_part);
    SetupElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    Matrix left_hand_side_matrix;
    Vector right_hand_side;
    p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side,
                                    r_current_process_info);

    Matrix expected_matrix(3, 3);
    expected_matrix <<=  5,   -5,    0,
                        -5,  5.5, -0.5,
                         0, -0.5,  0.5;
    ExpectDoubleMatrixEqual(expected_matrix, left_hand_side_matrix);

    // Calculated by hand (matrix multiplication between the
    // lhs and the temperature vector, which is 0.0, 1.0, 2.0)
    Vector expected_vector(3);
    expected_vector <<= 5, -4.5, -0.5;
    ExpectDoubleVectorEqual(expected_vector, right_hand_side);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D3N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.6666, 0.3333, 0.0);
    std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4};
    rModelPart.CreateNewElement("GeoTransientThermalElement2D4N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D4N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6};
    rModelPart.CreateNewElement("GeoTransientThermalElement2D6N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D6N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6, 7, 8};
    rModelPart.CreateNewElement("GeoTransientThermalElement2D8N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D8N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6, 7, 8, 9};
    rModelPart.CreateNewElement("GeoTransientThermalElement2D9N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D9N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    rModelPart.CreateNewElement("GeoTransientThermalElement2D10N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D10N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(11, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(12, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(13, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(14, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(15, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2,  3,  4,  5,  6,  7, 8,
                                               9, 10, 11, 12, 13, 14, 15};
    rModelPart.CreateNewElement("GeoTransientThermalElement2D15N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D15N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    // The z coordinate is non-zero to ensure the element has a non-zero domain size
    rModelPart.CreateNewNode(4, 0.5, 0.0, 1.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4};
    rModelPart.CreateNewElement("GeoTransientThermalElement3D4N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(CheckElement_Throws_When3DPropertyHasInvalidValue, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
        "THERMAL_CONDUCTIVITY_SOLID_ZZ has an invalid value at element 1")
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D4N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6, 7, 8};
    rModelPart.CreateNewElement("GeoTransientThermalElement3D8N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D8N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    rModelPart.CreateNewElement("GeoTransientThermalElement3D10N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D10N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(11, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(12, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(13, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(14, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(15, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(16, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(17, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(18, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(19, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(20, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> node_ids{
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    rModelPart.CreateNewElement("GeoTransientThermalElement3D20N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D20N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
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
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(7, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(8, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(9, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(11, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(12, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(13, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(14, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(15, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(16, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(17, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(18, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(19, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(20, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(21, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(22, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(23, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(24, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(25, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(26, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(27, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> node_ids{
        1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14,
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
    rModelPart.CreateNewElement("GeoTransientThermalElement3D27N", 1, node_ids,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D27N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D27N(model_part);
    SetupElement(model_part);

    ValidateThermalElement(model_part);
}

} // namespace Kratos::Testing
