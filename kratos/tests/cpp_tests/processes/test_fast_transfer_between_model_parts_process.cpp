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
#include "geometries/triangle_3d_3.h"

/* Processes */
#include "processes/fast_transfer_between_model_parts_process.h"

namespace Kratos::Testing
{
/**
* Checks the correct work of the fast_transfer_between_model_parts_process
* Test 1
*/
KRATOS_TEST_CASE_IN_SUITE(FastTransferBetweenModelPartsProcess1, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& r_origin_model_part = current_model.CreateModelPart("Origin");
    ModelPart& r_destination_model_part = current_model.CreateModelPart("Destination");

    Properties::Pointer p_prop = r_origin_model_part.CreateNewProperties(0);

    auto& r_process_info = r_origin_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;

    // First we create the nodes
    Node::Pointer p_node_1 = r_origin_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
    Node::Pointer p_node_2 = r_origin_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
    Node::Pointer p_node_3 = r_origin_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);
    std::vector<Node::Pointer> nodes_0 = {p_node_3, p_node_2, p_node_1};

    Node::Pointer p_node_4 = r_origin_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
    Node::Pointer p_node_5 = r_origin_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
    Node::Pointer p_node_6 = r_origin_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);
    std::vector<Node::Pointer> nodes_1 = { p_node_4, p_node_5, p_node_6};

    // Now we create the "geometries"
    Triangle3D3<Node> triangle_0( PointerVector<Node>{nodes_0} );
    Triangle3D3<Node> triangle_1( PointerVector<Node>{nodes_1} );
    r_origin_model_part.CreateNewGeometry("Triangle3D3", 1, triangle_0);
    r_origin_model_part.CreateNewGeometry("Triangle3D3", 2, triangle_1);

    // Now we create the "elements"
    r_origin_model_part.CreateNewElement("Element3D3N", 1, triangle_0, p_prop);
    r_origin_model_part.CreateNewElement("Element3D3N", 2, triangle_1, p_prop);

    // Now we create the "conditions"
    r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_0, p_prop);
    r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, triangle_1, p_prop);

    // This will copy all
    FastTransferBetweenModelPartsProcess process = FastTransferBetweenModelPartsProcess(r_destination_model_part, r_origin_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL);
    process.Execute();

    std::size_t count = 0;
    KRATOS_EXPECT_EQ(r_origin_model_part.NumberOfNodes(), r_destination_model_part.NumberOfNodes());
    for (auto& r_node : r_origin_model_part.Nodes()) {
        ++count;
        KRATOS_EXPECT_EQ(r_node.Id(), r_destination_model_part.GetNode(count).Id());
    }

    // Number of geometries (not ordered)
    KRATOS_EXPECT_EQ(r_origin_model_part.NumberOfGeometries(), r_destination_model_part.NumberOfGeometries());

    count = 0;
    KRATOS_EXPECT_EQ(r_origin_model_part.NumberOfElements(), r_destination_model_part.NumberOfElements());
    for (auto& r_elem : r_origin_model_part.Elements()) {
        ++count;
        KRATOS_EXPECT_EQ(r_elem.Id(), r_destination_model_part.GetElement(count).Id());
    }

    count = 0;
    KRATOS_EXPECT_EQ(r_origin_model_part.NumberOfConditions(), r_destination_model_part.NumberOfConditions());
    for (auto& r_cond : r_origin_model_part.Conditions()) {
        ++count;
        KRATOS_EXPECT_EQ(r_cond.Id(), r_destination_model_part.GetCondition(count).Id());
    }
}

/**
* Checks the correct work of the fast_transfer_between_model_parts_process
* Test 2 (with flags)
*/
KRATOS_TEST_CASE_IN_SUITE(FastTransferBetweenModelPartsProcess2, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& r_origin_model_part = current_model.CreateModelPart("Origin");
    ModelPart& r_destination_model_part = current_model.CreateModelPart("Destination");

    Properties::Pointer p_cond_prop = r_origin_model_part.CreateNewProperties(0);

    auto& r_process_info = r_origin_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;

    // First we create the nodes
    Node::Pointer p_node_1 = r_origin_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
    Node::Pointer p_node_2 = r_origin_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
    Node::Pointer p_node_3 = r_origin_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);
    std::vector<Node::Pointer> nodes_0 = {p_node_3, p_node_2, p_node_1};

    Node::Pointer p_node_4 = r_origin_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
    Node::Pointer p_node_5 = r_origin_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
    Node::Pointer p_node_6 = r_origin_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);
    std::vector<Node::Pointer> nodes_1 = { p_node_4, p_node_5, p_node_6};

    // Now we create the "conditions"
    Triangle3D3<Node> triangle_0( PointerVector<Node>{nodes_0} );
    Triangle3D3<Node> triangle_1( PointerVector<Node>{nodes_1} );

    Condition::Pointer p_cond_0 = r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_0, p_cond_prop);
    Condition::Pointer p_cond_1 = r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, triangle_1, p_cond_prop);

    // Setting flags
    // SLAVE
    p_node_1->Set(SLAVE, true);
    p_node_1->Set(MASTER, false);
    p_node_2->Set(SLAVE, true);
    p_node_2->Set(MASTER, false);
    p_node_3->Set(SLAVE, true);
    p_node_3->Set(MASTER, false);
    p_cond_0->Set(SLAVE, true);
    p_cond_0->Set(MASTER, false);
    // MASTER
    p_node_4->Set(SLAVE, false);
    p_node_4->Set(MASTER, true);
    p_node_5->Set(SLAVE, false);
    p_node_5->Set(MASTER, true);
    p_node_6->Set(SLAVE, false);
    p_node_6->Set(MASTER, true);
    p_cond_1->Set(SLAVE, false);
    p_cond_1->Set(MASTER, true);

    // This will copy all
    FastTransferBetweenModelPartsProcess process = FastTransferBetweenModelPartsProcess(r_destination_model_part, r_origin_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, MASTER);
    process.Execute();

    std::size_t count = 0;
    for (auto& r_node : r_origin_model_part.Nodes()) {
        ++count;
        if (r_node.Is(MASTER))
            KRATOS_EXPECT_EQ(r_node.Id(), r_destination_model_part.GetNode(count).Id());
    }

    count = 0;
    for (auto& cond : r_origin_model_part.Conditions()) {
        ++count;
        if (cond.Is(MASTER))
            KRATOS_EXPECT_EQ(cond.Id(), r_destination_model_part.GetCondition(count).Id());
    }
}

/**
* Checks the correct work of the fast_transfer_between_model_parts_process
* Test 3 (clone/replicate)
*/
KRATOS_TEST_CASE_IN_SUITE(FastTransferBetweenModelPartsProcess3, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& r_origin_model_part = current_model.CreateModelPart("Origin");
    ModelPart& r_destination_model_part = current_model.CreateModelPart("Destination");

    Properties::Pointer p_cond_prop = r_origin_model_part.CreateNewProperties(0);

    auto& r_process_info = r_origin_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;

    // First we create the nodes
    Node::Pointer p_node_1 = r_origin_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
    Node::Pointer p_node_2 = r_origin_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
    Node::Pointer p_node_3 = r_origin_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);
    std::vector<Node::Pointer> nodes_0 = {p_node_3, p_node_2, p_node_1};

    Node::Pointer p_node_4 = r_origin_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
    Node::Pointer p_node_5 = r_origin_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
    Node::Pointer p_node_6 = r_origin_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);
    std::vector<Node::Pointer> nodes_1 = { p_node_4, p_node_5, p_node_6};

    // Now we create the "conditions"
    Triangle3D3<Node> triangle_0( PointerVector<Node>{nodes_0} );
    Triangle3D3<Node> triangle_1( PointerVector<Node>{nodes_1} );

    Condition::Pointer p_cond_0 = r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_0, p_cond_prop);
    Condition::Pointer p_cond_1 = r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, triangle_1, p_cond_prop);

    // This will copy all
    FastTransferBetweenModelPartsProcess process = FastTransferBetweenModelPartsProcess(r_destination_model_part, r_origin_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, Flags(), true);
    process.Execute();

    std::size_t count = 0;
    for (auto& r_node : r_destination_model_part.Nodes()) {
        ++count;
        KRATOS_EXPECT_EQ(r_node.Id(), r_origin_model_part.GetNode(count).Id() + 6);
    }

    count = 0;
    for (auto& r_cond : r_destination_model_part.Conditions()) {
        ++count;
        KRATOS_EXPECT_EQ(r_cond.Id(), r_origin_model_part.GetCondition(count).Id() + 2);
    }
}

/**
* Checks the correct work of the fast_transfer_between_model_parts_process
* Test 4 (clone/replicate with flags)
*/
KRATOS_TEST_CASE_IN_SUITE(FastTransferBetweenModelPartsProcess4, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& r_origin_model_part = current_model.CreateModelPart("Origin");
    ModelPart& r_destination_model_part = current_model.CreateModelPart("Destination");

    Properties::Pointer p_cond_prop = r_origin_model_part.CreateNewProperties(0);

    auto& r_process_info = r_origin_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;

    // First we create the nodes
    Node::Pointer p_node_1 = r_origin_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
    Node::Pointer p_node_2 = r_origin_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
    Node::Pointer p_node_3 = r_origin_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);
    std::vector<Node::Pointer> nodes_0 = {p_node_3, p_node_2, p_node_1};

    Node::Pointer p_node_4 = r_origin_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
    Node::Pointer p_node_5 = r_origin_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
    Node::Pointer p_node_6 = r_origin_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);
    std::vector<Node::Pointer> nodes_1 = { p_node_4, p_node_5, p_node_6};

    // Now we create the "conditions"
    Triangle3D3 <Node> triangle_0( PointerVector<Node>{nodes_0} );
    Triangle3D3 <Node> triangle_1( PointerVector<Node>{nodes_1} );

    Condition::Pointer p_cond_0 = r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_0, p_cond_prop);
    Condition::Pointer p_cond_1 = r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, triangle_1, p_cond_prop);

    // Setting flags
    // SLAVE
    p_node_1->Set(SLAVE, true);
    p_node_1->Set(MASTER, false);
    p_node_2->Set(SLAVE, true);
    p_node_2->Set(MASTER, false);
    p_node_3->Set(SLAVE, true);
    p_node_3->Set(MASTER, false);
    p_cond_0->Set(SLAVE, true);
    p_cond_0->Set(MASTER, false);
    // MASTER
    p_node_4->Set(SLAVE, false);
    p_node_4->Set(MASTER, true);
    p_node_5->Set(SLAVE, false);
    p_node_5->Set(MASTER, true);
    p_node_6->Set(SLAVE, false);
    p_node_6->Set(MASTER, true);
    p_cond_1->Set(SLAVE, false);
    p_cond_1->Set(MASTER, true);

    // This will copy all
    FastTransferBetweenModelPartsProcess process = FastTransferBetweenModelPartsProcess(r_destination_model_part, r_origin_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, MASTER, true);
    process.Execute();

    std::size_t count = 0;
    for (auto& r_node : r_origin_model_part.Nodes()) {
        ++count;
        if (r_node.Is(MASTER))
            KRATOS_EXPECT_EQ(r_node.Id() + 6, r_destination_model_part.GetNode(count + 6).Id());
    }

    count = 0;
    for (auto& r_cond : r_origin_model_part.Conditions()) {
        ++count;
        if (r_cond.Is(MASTER))
            KRATOS_EXPECT_EQ(r_cond.Id() + 2, r_destination_model_part.GetCondition(count  + 2).Id());
    }
}

}  // namespace Kratos::Testing.
