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

using namespace Kratos;

namespace Kratos::Testing {

void TestThermalElement(ModelPart& rModelPart)
{
    Element::Pointer p_element = rModelPart.pGetElement(1);

    for (unsigned int i = 0; i < rModelPart.NumberOfNodes(); i++) {
        p_element->GetGeometry()[i].AddDof(TEMPERATURE);
    }

    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    Element::DofsVectorType ElementalDofList;
    p_element->GetDofList(ElementalDofList, r_current_process_info);

    for (unsigned int i = 0; i < ElementalDofList.size(); i++) {
        ElementalDofList[i]->SetEquationId(i);
    }

    Element::EquationIdVectorType EquationIdVector;
    p_element->EquationIdVector(EquationIdVector, r_current_process_info);

    for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_EXPECT_EQ(EquationIdVector[i], i);
    }
}

void GenerateTransientThermalElement2D3N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
    rModelPart.CreateNewElement("TransientThermalElement2D3N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D3N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D3N(model_part);

    TestThermalElement(model_part);
}

void GenerateTransientThermalElement2D4N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.6666, 0.3333, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4};
    rModelPart.CreateNewElement("TransientThermalElement2D4N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D4N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D4N(model_part);

    TestThermalElement(model_part);
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

    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4, 5, 6};
    rModelPart.CreateNewElement("TransientThermalElement2D6N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D6N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D6N(model_part);

    TestThermalElement(model_part);
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

    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4, 5, 6, 7, 8};
    rModelPart.CreateNewElement("TransientThermalElement2D8N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D8N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D8N(model_part);

    TestThermalElement(model_part);
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

    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4, 5, 6, 7, 8, 9};
    rModelPart.CreateNewElement("TransientThermalElement2D9N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D9N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D9N(model_part);

    TestThermalElement(model_part);
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

    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    rModelPart.CreateNewElement("TransientThermalElement2D10N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D10N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D10N(model_part);

    TestThermalElement(model_part);
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

    std::vector<ModelPart::IndexType> elemNodes{1, 2,  3,  4,  5,  6,  7, 8,
                                                9, 10, 11, 12, 13, 14, 15};
    rModelPart.CreateNewElement("TransientThermalElement2D15N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D15N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement2D15N(model_part);

    TestThermalElement(model_part);
}

void GenerateTransientThermalElement3D4N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);

    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4};
    rModelPart.CreateNewElement("TransientThermalElement3D4N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D4N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D4N(model_part);

    TestThermalElement(model_part);
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

    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4, 5, 6, 7, 8};
    rModelPart.CreateNewElement("TransientThermalElement3D8N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D8N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D8N(model_part);

    TestThermalElement(model_part);
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

    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    rModelPart.CreateNewElement("TransientThermalElement3D10N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D10N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D10N(model_part);

    TestThermalElement(model_part);
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

    std::vector<ModelPart::IndexType> elemNodes{
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    rModelPart.CreateNewElement("TransientThermalElement3D20N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D20N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D20N(model_part);

    TestThermalElement(model_part);
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

    std::vector<ModelPart::IndexType> elemNodes{
        1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14,
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
    rModelPart.CreateNewElement("TransientThermalElement3D27N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement3D27N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTransientThermalElement3D27N(model_part);

    TestThermalElement(model_part);
}

} // namespace Kratos::Testing
