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
    Condition::Pointer p_condition = rModelPart.pGetCondition(1);

    for (unsigned int i = 0; i < rModelPart.NumberOfNodes(); i++) {
        p_condition->GetGeometry()[i].AddDof(TEMPERATURE);
    }

    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
 
    Condition::DofsVectorType dof_list;
    p_condition->GetDofList(dof_list, r_current_process_info);

    for (unsigned int i = 0; i < dof_list.size(); i++) {
        dof_list[i]->SetEquationId(i);
    }

    Condition::EquationIdVectorType equation_id_vector;
    p_condition->EquationIdVector(equation_id_vector, r_current_process_info);

    // Check the EquationIdVector values
    for (unsigned int i = 0; i < equation_id_vector.size(); i++) {
        KRATOS_EXPECT_EQ(equation_id_vector[i], i);
    }
}

void GenerateTnormalFluxCondition2D2N(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo()[DOMAIN_SIZE] = 2;

    Properties::Pointer p_cond_prop = rModelPart.CreateNewProperties(0);
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    std::vector<ModelPart::IndexType> cond_nodes{1, 2};
    rModelPart.CreateNewCondition("TNormalFluxCondition2D2N", 1, cond_nodes,
                                p_cond_prop);
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
    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3};
    rModelPart.CreateNewCondition("TNormalFluxCondition2D3N", 1, cond_nodes,
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
    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4};
    rModelPart.CreateNewCondition("TNormalFluxCondition2D4N", 1, cond_nodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTNormalFluxCondition2D4N, KratosGeoMechanicsFastSuite)
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

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5};
    rModelPart.CreateNewCondition("TNormalFluxCondition2D5N", 1, cond_nodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTNormalFluxCondition2D5N, KratosGeoMechanicsFastSuite)
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

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4};
    rModelPart.CreateNewCondition("TNormalFluxCondition3D4N", 1, cond_nodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTNormalFluxCondition3D4N, KratosGeoMechanicsFastSuite)
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

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5, 6};
    rModelPart.CreateNewCondition("TNormalFluxCondition3D6N", 1, cond_nodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTNormalFluxCondition3D6N, KratosGeoMechanicsFastSuite)
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

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5, 6, 7, 8};
    rModelPart.CreateNewCondition("TNormalFluxCondition3D8N", 1, cond_nodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTNormalFluxCondition3D8N, KratosGeoMechanicsFastSuite)
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

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9};
    rModelPart.CreateNewCondition("TNormalFluxCondition3D9N", 1, cond_nodes,
                                rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTNormalFluxCondition3D9N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    GenerateTnormalFluxCondition3D9N(model_part);

    TestTnormalFluxCondition(model_part);
}

} // namespace Kratos::Testing
