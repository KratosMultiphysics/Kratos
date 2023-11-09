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
#include "custom_conditions/T_normal_flux_condition.h"
#include "geo_mechanics_application_variables.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing {

void TestTnormalFluxCondition(ModelPart& rModelPart)
{
    Element::Pointer p_element = rModelPart.pGetElement(1);

    for (unsigned int i = 0; i < rModelPart.NumberOfNodes(); i++) {
        p_element->GetGeometry()[i].AddDof(TEMPERATURE);
    }

    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    // p_element->Check(r_current_process_info);
    // Check requires a solution?

    Element::DofsVectorType ElementalDofList;
    p_element->GetDofList(ElementalDofList, r_current_process_info);

    for (unsigned int i = 0; i < ElementalDofList.size(); i++) {
        ElementalDofList[i]->SetEquationId(i);
    }

    Element::EquationIdVectorType EquationIdVector;
    p_element->EquationIdVector(EquationIdVector, r_current_process_info);

    // Check the EquationIdVector values
    for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_EXPECT_TRUE(EquationIdVector[i] == i);
    }
}

void GenerateTnormalFluxCondition2D2N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2};
    rModelPart.CreateNewElement("mTNormalFluxCondition2D2N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTnormalFluxCondition2D2N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTnormalFluxCondition2D2N(model_part);

    TestTnormalFluxCondition(model_part);
}

void GenerateTnormalFluxCondition2D3N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
    rModelPart.CreateNewElement("mTNormalFluxCondition2D3N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTnormalFluxCondition2D3N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTnormalFluxCondition2D3N(model_part);

    TestTnormalFluxCondition(model_part);
}

void GenerateTnormalFluxCondition2D4N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.6666, 0.3333, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4};
    rModelPart.CreateNewElement("mTNormalFluxCondition2D4N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectormTNormalFluxCondition2D4N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTnormalFluxCondition2D4N(model_part);

    TestTnormalFluxCondition(model_part);
}

void GenerateTnormalFluxCondition2D5N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);

    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4, 5, 6};
    rModelPart.CreateNewElement("mTNormalFluxCondition2D5N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectormTNormalFluxCondition2D5N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTnormalFluxCondition2D5N(model_part);

    TestTnormalFluxCondition(model_part);
}

void GenerateTnormalFluxCondition3D4N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);

    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4};
    rModelPart.CreateNewElement("mTNormalFluxCondition3D4N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectormTNormalFluxCondition3D4N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTnormalFluxCondition3D4N(model_part);

    TestTnormalFluxCondition(model_part);
}

void GenerateTnormalFluxCondition3D6N(ModelPart& rModelPart)
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
    rModelPart.CreateNewElement("mTNormalFluxCondition3D6N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectormTNormalFluxCondition3D6N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTnormalFluxCondition3D6N(model_part);

    TestTnormalFluxCondition(model_part);
}

void GenerateTnormalFluxCondition3D8N(ModelPart& rModelPart)
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
    rModelPart.CreateNewElement("mTNormalFluxCondition3D8N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectormTNormalFluxCondition3D8N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTnormalFluxCondition3D8N(model_part);

    TestTnormalFluxCondition(model_part);
}

void GenerateTnormalFluxCondition3D9N(ModelPart& rModelPart)
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
    rModelPart.CreateNewElement("mTNormalFluxCondition3D9N", 1, elemNodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectormTNormalFluxCondition3D9N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTnormalFluxCondition3D9N(model_part);

    TestTnormalFluxCondition(model_part);
}

} // namespace Kratos::Testing
