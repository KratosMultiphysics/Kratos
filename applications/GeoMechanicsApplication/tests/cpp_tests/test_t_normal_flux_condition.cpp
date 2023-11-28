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
#include "includes/condition.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing {

void TestGeoTnormalFluxCondition(ModelPart& rModelPart, const std::vector<double> & rExpectedRightHandSide)
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

    Matrix left_hand_side_matrix = ZeroMatrix(1, 1);
    Vector right_hand_side_vector = ZeroVector(1);

    p_condition->CalculateLocalSystem(
        left_hand_side_matrix, right_hand_side_vector, r_current_process_info);

    for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
        KRATOS_EXPECT_NEAR(right_hand_side_vector[i], rExpectedRightHandSide[i], 1.0e-4);
    }
}

void AssignNormalHeatFlux(ModelPart& rModelPart, double value)
{
    rModelPart.GetProcessInfo()[DELTA_TIME] = 1.0;

    for (unsigned int i = 0; i < rModelPart.NumberOfNodes(); i++) {
        auto p_node = rModelPart.pGetNode(i + 1);
        p_node->FastGetSolutionStepValue(NORMAL_HEAT_FLUX) = value;
        p_node->FastGetSolutionStepValue(NORMAL_HEAT_FLUX, 1) = 0;
    }
}

void GenerateGeoTnormalFluxCondition2D2N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);

    AssignNormalHeatFlux(rModelPart, 10.0);

    std::vector<ModelPart::IndexType> cond_nodes{1, 2};
    rModelPart.CreateNewCondition("GeoTNormalFluxCondition2D2N", 1, cond_nodes,
                                  rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(GeoTNormalFluxCondition2D2N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(NORMAL_HEAT_FLUX);

    GenerateGeoTnormalFluxCondition2D2N(model_part);

    const std::vector<double> expected_right_hand_side{5.0, 5.0};
    TestGeoTnormalFluxCondition(model_part, expected_right_hand_side);
}

void GenerateGeoTnormalFluxCondition2D3N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);

    AssignNormalHeatFlux(rModelPart, 10.0);

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3};
    rModelPart.CreateNewCondition("GeoTNormalFluxCondition2D3N", 1, cond_nodes,
                                  rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(GeoTNormalFluxCondition2D3N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(NORMAL_HEAT_FLUX);

    GenerateGeoTnormalFluxCondition2D3N(model_part);

    const std::vector<double> expected_right_hand_side{5.77898, 3.3428, 18.24365};
    TestGeoTnormalFluxCondition(model_part, expected_right_hand_side);
}

void GenerateGeoTnormalFluxCondition2D4N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.6666, 0.3333, 0.0);

    AssignNormalHeatFlux(rModelPart, 10.0);

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4};
    rModelPart.CreateNewCondition("GeoTNormalFluxCondition2D4N", 1, cond_nodes,
                                  rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(GeoTNormalFluxCondition2D4N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(NORMAL_HEAT_FLUX);

    GenerateGeoTnormalFluxCondition2D4N(model_part);

    const std::vector<double> expected_right_hand_side{7.33331, 2.60131, 16.89195, 5.27702};
    TestGeoTnormalFluxCondition(model_part, expected_right_hand_side);
}

void GenerateGeoTnormalFluxCondition2D5N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);

    AssignNormalHeatFlux(rModelPart, 10.0);

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5};
    rModelPart.CreateNewCondition("GeoTNormalFluxCondition2D5N", 1, cond_nodes,
                                  rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(GeoTNormalFluxCondition2D5N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(NORMAL_HEAT_FLUX);

    GenerateGeoTnormalFluxCondition2D5N(model_part);

    const std::vector<double> expected_right_hand_side{9.11795, 3.92042, 27.99155, 1.3485, 18.8636};
    TestGeoTnormalFluxCondition(model_part, expected_right_hand_side);
}

void GenerateGeoTnormalFluxCondition3D4N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);

    AssignNormalHeatFlux(rModelPart, 10.0);

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4};
    rModelPart.CreateNewCondition("GeoTNormalFluxCondition3D4N", 1, cond_nodes,
                                  rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(GeoTNormalFluxCondition3D4N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(NORMAL_HEAT_FLUX);

    GenerateGeoTnormalFluxCondition3D4N(model_part);

    const std::vector<double> expected_right_hand_side{0.569177, 1.29087, 0.985844, 0.569177};
    TestGeoTnormalFluxCondition(model_part, expected_right_hand_side);
}

void GenerateGeoTnormalFluxCondition3D6N(ModelPart& rModelPart)
{
    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 0.5, 0.0);

    AssignNormalHeatFlux(rModelPart, 10.0);

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5, 6};
    rModelPart.CreateNewCondition("GeoTNormalFluxCondition3D6N", 1, cond_nodes,
                                  rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(GeoTNormalFluxCondition3D6N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(NORMAL_HEAT_FLUX);

    GenerateGeoTnormalFluxCondition3D6N(model_part);

    const std::vector<double> expected_right_hand_side{3.60822e-16, -2.77556e-16, -1.66533e-16,
                            1.66667,     1.66667,      1.66667};
    TestGeoTnormalFluxCondition(model_part, expected_right_hand_side);
}

void GenerateGeoTnormalFluxCondition3D8N(ModelPart& rModelPart)
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

    AssignNormalHeatFlux(rModelPart, 10.0);

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5, 6, 7, 8};
    rModelPart.CreateNewCondition("GeoTNormalFluxCondition3D8N", 1, cond_nodes,
                                  rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(GeoTNormalFluxCondition3D8N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(NORMAL_HEAT_FLUX);

    GenerateGeoTnormalFluxCondition3D8N(model_part);

    const std::vector<double> expected_right_hand_side{0.133919, 0.635973, -0.403995, -0.117107,
                            1.86947,  1.53697,  0.721923,  1.09723};
    TestGeoTnormalFluxCondition(model_part, expected_right_hand_side);
}

void GenerateGeoTnormalFluxCondition3D9N(ModelPart& rModelPart)
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

    AssignNormalHeatFlux(rModelPart, 10.0);

    std::vector<ModelPart::IndexType> cond_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9};
    rModelPart.CreateNewCondition("GeoTNormalFluxCondition3D9N", 1, cond_nodes,
                                  rModelPart.CreateNewProperties(0));
}

KRATOS_TEST_CASE_IN_SUITE(GeoTNormalFluxCondition3D9N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);
    // Variables addition
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(NORMAL_HEAT_FLUX);

    GenerateGeoTnormalFluxCondition3D9N(model_part);

    const std::vector<double> expected_right_hand_side{0.325617, 0.741605, -0.0401643, 0.160657, 0.785555,
                            0.467627, 0.228477, 0.324183,   0.676021};
    TestGeoTnormalFluxCondition(model_part, expected_right_hand_side);
}

} // namespace Kratos::Testing
