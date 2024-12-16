//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/kratos_flags.h"
#include "includes/mapping_variables.h"
#include "utilities/normal_calculation_utils.h"
// #include "includes/gid_io.h"

/* Processes */
#include "processes/simple_mortar_mapper_process.h"

namespace Kratos::Testing
{

// void GiDIODebugMapper(ModelPart& ThisModelPart)
// {
//     GidIO<> gid_io("TEST_MAPPER", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
//     const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//     const double label = static_cast<double>(nl_iter);

//     gid_io.InitializeMesh(label);
//     gid_io.WriteMesh(ThisModelPart.GetMesh());
//     gid_io.FinalizeMesh();
//     gid_io.InitializeResults(label, ThisModelPart.GetMesh());
//     gid_io.WriteNodalResults(TEMPERATURE, ThisModelPart.Nodes(), label, 0);
//     gid_io.WriteNodalResults(NORMAL, ThisModelPart.Nodes(), label, 0);
// }

/**
* Checks the correct work of the simple mortar mapper
* Test 1
*/
KRATOS_TEST_CASE_IN_SUITE(SimpleMortarMapper1, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);
    this_model_part.CreateSubModelPart("SlaveModelPart");
    ModelPart& slave_model_part = this_model_part.GetSubModelPart("SlaveModelPart");
    this_model_part.CreateSubModelPart("MasterModelPart");
    ModelPart& master_model_part = this_model_part.GetSubModelPart("MasterModelPart");

    this_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    this_model_part.AddNodalSolutionStepVariable(NORMAL);

    Properties::Pointer p_cond_prop = this_model_part.CreateNewProperties(0);

    auto& process_info = this_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    // First we create the nodes
    Node::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
    Node::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
    Node::Pointer p_node_3 = this_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);

    Node::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
    Node::Pointer p_node_5 = this_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
    Node::Pointer p_node_6 = this_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);

    // Now we create the "conditions"
    std::vector<Node::Pointer> condition_nodes_0 (3);
    condition_nodes_0[0] = p_node_3;
    condition_nodes_0[1] = p_node_2;
    condition_nodes_0[2] = p_node_1;
    Triangle3D3 <Node> triangle_0( PointerVector<Node>{condition_nodes_0} );

    std::vector<Node::Pointer> condition_nodes_1 (3);
    condition_nodes_1[0] = p_node_4;
    condition_nodes_1[1] = p_node_5;
    condition_nodes_1[2] = p_node_6;
    Triangle3D3 <Node> triangle_1( PointerVector<Node>{condition_nodes_1} );

    Condition::Pointer p_cond_0 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_0, p_cond_prop);
    Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, triangle_1, p_cond_prop);

    // Adding map
    IndexSet this_set;
    this_set.AddId(2);
    p_cond_0->SetValue(INDEX_SET, Kratos::make_shared<IndexSet>(this_set));

    // SLAVE
    slave_model_part.AddNode(p_node_1);
    slave_model_part.AddNode(p_node_2);
    slave_model_part.AddNode(p_node_3);
    slave_model_part.AddCondition(p_cond_0);
    // MASTER
    master_model_part.AddNode(p_node_4);
    master_model_part.AddNode(p_node_5);
    master_model_part.AddNode(p_node_6);
    master_model_part.AddCondition(p_cond_1);

    // We compute the normals
    NormalCalculationUtils().CalculateUnitNormals<ModelPart::ConditionsContainerType>(this_model_part, true);

    p_node_4->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_4->X(), 2) + std::pow(p_node_4->Y(), 2);
    p_node_5->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_5->X(), 2) + std::pow(p_node_5->Y(), 2);
    p_node_6->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_6->X(), 2) + std::pow(p_node_6->Y(), 2);

    const double tolerance = 1.0e-4;
    auto int_util = ExactMortarIntegrationUtility<3, 3>();
    double area;
    int_util.GetExactAreaIntegration(p_cond_0->GetGeometry(), p_cond_0->GetValue(NORMAL), p_cond_1->GetGeometry(), p_cond_1->GetValue(NORMAL), area);
    KRATOS_EXPECT_LE((area - 0.499925)/0.499925, tolerance);

    typedef SimpleMortarMapperProcess<3, 3, Variable<double>> MapperType;
    MapperType process = MapperType(master_model_part, slave_model_part, TEMPERATURE);
    process.Execute();

    // // DEBUG
    // GiDIODebugMapper(this_model_part);

    KRATOS_EXPECT_LE(std::abs(p_node_1->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_1->X(), 2) + std::pow(p_node_1->Y(), 2))), tolerance);
    KRATOS_EXPECT_LE(std::abs(p_node_2->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_2->X(), 2) + std::pow(p_node_2->Y(), 2))), tolerance);
    KRATOS_EXPECT_LE(std::abs(p_node_3->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_3->X(), 2) + std::pow(p_node_3->Y(), 2))), tolerance);
}

/**
* Checks the correct work of the simple mortar mapper
* Test 2
*/
KRATOS_TEST_CASE_IN_SUITE(SimpleMortarMapper2, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);
    this_model_part.CreateSubModelPart("SlaveModelPart");
    ModelPart& slave_model_part = this_model_part.GetSubModelPart("SlaveModelPart");
    this_model_part.CreateSubModelPart("MasterModelPart");
    ModelPart& master_model_part = this_model_part.GetSubModelPart("MasterModelPart");

    this_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    this_model_part.AddNodalSolutionStepVariable(NORMAL);

    Properties::Pointer p_cond_prop = this_model_part.CreateNewProperties(0);

    auto& process_info = this_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    // First we create the nodes
    Node::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
    Node::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
    Node::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.01);
    Node::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.01);

    Node::Pointer p_node_5 = this_model_part.CreateNewNode(5, 0.0 , 0.0 , 0.01);
    Node::Pointer p_node_6 = this_model_part.CreateNewNode(6, 1.0 , 0.0 , 0.01);
    Node::Pointer p_node_7 = this_model_part.CreateNewNode(7, 1.0 , 1.0 , 0.02);
    Node::Pointer p_node_8 = this_model_part.CreateNewNode(8, 0.0 , 1.0 , 0.02);

    // Now we create the "conditions"
    std::vector<Node::Pointer> condition_nodes_0 (4);
    condition_nodes_0[0] = p_node_1;
    condition_nodes_0[1] = p_node_2;
    condition_nodes_0[2] = p_node_3;
    condition_nodes_0[3] = p_node_4;
    Quadrilateral3D4 <Node> quad_0( PointerVector<Node>{condition_nodes_0} );

    std::vector<Node::Pointer> condition_nodes_1 (4);
    condition_nodes_1[0] = p_node_8;
    condition_nodes_1[1] = p_node_7;
    condition_nodes_1[2] = p_node_6;
    condition_nodes_1[3] = p_node_5;
    Quadrilateral3D4 <Node> quad_1( PointerVector<Node>{condition_nodes_1} );

    Condition::Pointer p_cond_0 = this_model_part.CreateNewCondition("SurfaceCondition3D4N", 1, quad_0, p_cond_prop);
    Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("SurfaceCondition3D4N", 2, quad_1, p_cond_prop);

    // Adding map
    IndexSet this_set;
    this_set.AddId(2);
    p_cond_0->SetValue(INDEX_SET, Kratos::make_shared<IndexSet>(this_set));

    // SLAVE
    slave_model_part.AddNode(p_node_1);
    slave_model_part.AddNode(p_node_2);
    slave_model_part.AddNode(p_node_3);
    slave_model_part.AddNode(p_node_4);
    slave_model_part.AddCondition(p_cond_0);
    // MASTER
    master_model_part.AddNode(p_node_5);
    master_model_part.AddNode(p_node_6);
    master_model_part.AddNode(p_node_7);
    master_model_part.AddNode(p_node_8);
    master_model_part.AddCondition(p_cond_1);

    // We compute the normals
    NormalCalculationUtils().CalculateUnitNormals<ModelPart::ConditionsContainerType>(this_model_part, true);

    p_node_5->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_5->X(), 2) + std::pow(p_node_5->Y(), 2);
    p_node_6->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_6->X(), 2) + std::pow(p_node_6->Y(), 2);
    p_node_7->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_7->X(), 2) + std::pow(p_node_7->Y(), 2);
    p_node_8->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_8->X(), 2) + std::pow(p_node_8->Y(), 2);

    typedef SimpleMortarMapperProcess<3, 4, Variable<double>> MapperType;
    MapperType process = MapperType(master_model_part, slave_model_part, TEMPERATURE);
    process.Execute();

    // // DEBUG
    // GiDIODebugMapper(this_model_part);

    const double tolerance = 1.0e-4;
    KRATOS_EXPECT_LE(std::abs(p_node_1->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_1->X(), 2) + std::pow(p_node_1->Y(), 2))), tolerance);
    KRATOS_EXPECT_LE(std::abs(p_node_2->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_2->X(), 2) + std::pow(p_node_2->Y(), 2))), tolerance);
    KRATOS_EXPECT_LE(std::abs(p_node_3->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_3->X(), 2) + std::pow(p_node_3->Y(), 2))), tolerance);
    KRATOS_EXPECT_LE(std::abs(p_node_4->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_4->X(), 2) + std::pow(p_node_4->Y(), 2))), tolerance);
}

/**
* Checks the correct work of the simple mortar mapper
* Test 3
*/
KRATOS_TEST_CASE_IN_SUITE(SimpleMortarMapper3, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);
    this_model_part.CreateSubModelPart("SlaveModelPart");
    ModelPart& slave_model_part = this_model_part.GetSubModelPart("SlaveModelPart");
    this_model_part.CreateSubModelPart("MasterModelPart");
    ModelPart& master_model_part = this_model_part.GetSubModelPart("MasterModelPart");

    this_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    this_model_part.AddNodalSolutionStepVariable(NORMAL);

    Properties::Pointer p_cond_prop = this_model_part.CreateNewProperties(0);

    auto& process_info = this_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    // First we create the nodes
    Node::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.000, 0.000, 0.000);
    Node::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.000, 1.000, 1.000);
    Node::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.000, 0.000, 1.000);
    Node::Pointer p_node_4 = this_model_part.CreateNewNode(4, 1.000, 1.000, 0.000);

    Node::Pointer p_node_5 = this_model_part.CreateNewNode(5, 1.001, 0.000, 1.000);
    Node::Pointer p_node_6 = this_model_part.CreateNewNode(6, 1.001, 1.000, 0.000);
    Node::Pointer p_node_7 = this_model_part.CreateNewNode(7, 0.001, 0.000, 0.000);
    Node::Pointer p_node_8 = this_model_part.CreateNewNode(8, 1.001, 1.000, 1.000);

    // Now we create the "conditions"
    std::vector<Node::Pointer> condition_nodes_0 (3);
    std::vector<Node::Pointer> condition_nodes_1 (3);
    condition_nodes_0[0] = p_node_4;
    condition_nodes_0[1] = p_node_3;
    condition_nodes_0[2] = p_node_1;
    Triangle3D3 <Node> triangle_0( PointerVector<Node>{condition_nodes_0} );
    condition_nodes_1[0] = p_node_4;
    condition_nodes_1[1] = p_node_2;
    condition_nodes_1[2] = p_node_3;
    Triangle3D3 <Node> triangle_1( PointerVector<Node>{condition_nodes_1} );

    std::vector<Node::Pointer> condition_nodes_2 (3);
    std::vector<Node::Pointer> condition_nodes_3 (3);
    condition_nodes_2[0] = p_node_7;
    condition_nodes_2[1] = p_node_5;
    condition_nodes_2[2] = p_node_6;
    Triangle3D3 <Node> triangle_3( PointerVector<Node>{condition_nodes_2} );
    condition_nodes_3[0] = p_node_5;
    condition_nodes_3[1] = p_node_8;
    condition_nodes_3[2] = p_node_6;
    Triangle3D3 <Node> triangle_4( PointerVector<Node>{condition_nodes_3} );

    Condition::Pointer p_cond_0 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_0, p_cond_prop);
    Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, triangle_1, p_cond_prop);
    Condition::Pointer p_cond_2 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 3, triangle_3, p_cond_prop);
    Condition::Pointer p_cond_3 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 4, triangle_4, p_cond_prop);

    // Adding map
    IndexSet this_set0, this_set1;
    this_set0.AddId(3);
    this_set0.AddId(4);
    this_set1.AddId(3);
    this_set1.AddId(4);
    p_cond_0->SetValue(INDEX_SET, Kratos::make_shared<IndexSet>(this_set0));
    p_cond_1->SetValue(INDEX_SET, Kratos::make_shared<IndexSet>(this_set1));

    // SLAVE
    slave_model_part.AddNode(p_node_1);
    slave_model_part.AddNode(p_node_2);
    slave_model_part.AddNode(p_node_3);
    slave_model_part.AddNode(p_node_4);
    slave_model_part.AddCondition(p_cond_0);
    slave_model_part.AddCondition(p_cond_1);
    // MASTER
    master_model_part.AddNode(p_node_5);
    master_model_part.AddNode(p_node_6);
    master_model_part.AddNode(p_node_7);
    master_model_part.AddNode(p_node_8);
    master_model_part.AddCondition(p_cond_2);
    master_model_part.AddCondition(p_cond_3);

    // We compute the normals
    NormalCalculationUtils().CalculateUnitNormals<ModelPart::ConditionsContainerType>(this_model_part, true);

    p_node_5->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_5->Z(), 2) + std::pow(p_node_5->Y(), 2);
    p_node_6->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_6->Z(), 2) + std::pow(p_node_6->Y(), 2);
    p_node_7->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_7->Z(), 2) + std::pow(p_node_7->Y(), 2);
    p_node_8->FastGetSolutionStepValue(TEMPERATURE) = std::pow(p_node_8->Z(), 2) + std::pow(p_node_8->Y(), 2);

    typedef SimpleMortarMapperProcess<3, 3, Variable<double>> MapperType;
    MapperType process = MapperType(master_model_part, slave_model_part, TEMPERATURE);
    process.Execute();

    // // DEBUG
    // GiDIODebugMapper(this_model_part);

    const double tolerance = 1.0e-3;
    KRATOS_EXPECT_LE(std::abs(p_node_1->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_1->Z(), 2) + std::pow(p_node_1->Y(), 2))), tolerance);
    KRATOS_EXPECT_LE(std::abs(p_node_2->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_2->Z(), 2) + std::pow(p_node_2->Y(), 2))), tolerance);
    KRATOS_EXPECT_LE(std::abs(p_node_3->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_3->Z(), 2) + std::pow(p_node_3->Y(), 2))), tolerance);
    KRATOS_EXPECT_LE(std::abs(p_node_4->FastGetSolutionStepValue(TEMPERATURE) - (std::pow(p_node_4->Z(), 2) + std::pow(p_node_4->Y(), 2))), tolerance);
}
}  // namespace Kratos::Testing.