//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Philipp Bucher
//                   Vicente Mataix Ferrandiz
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "utilities/auxiliar_model_part_utilities.h"

namespace Kratos {
    namespace Testing {

    typedef Node<3> NodeType;

    void GenerateGenericModelPart(ModelPart& rModelPart)
    {
        Properties::Pointer p_elem_prop = rModelPart.pGetProperties(0);

        // First we create the nodes
        NodeType::Pointer p_node_1 = rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_2 = rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_3 = rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
        NodeType::Pointer p_node_4 = rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
        NodeType::Pointer p_node_5 = rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
        NodeType::Pointer p_node_6 = rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

        // Now we create the "conditions"
        std::vector<NodeType::Pointer> condition_nodes_0 (2);
        condition_nodes_0[0] = p_node_1;
        condition_nodes_0[1] = p_node_2;

        std::vector<NodeType::Pointer> condition_nodes_1 (2);
        condition_nodes_1[0] = p_node_1;
        condition_nodes_1[1] = p_node_4;

        std::vector<NodeType::Pointer> condition_nodes_2 (2);
        condition_nodes_2[0] = p_node_2;
        condition_nodes_2[1] = p_node_5;

        std::vector<NodeType::Pointer> condition_nodes_3 (2);
        condition_nodes_3[0] = p_node_5;
        condition_nodes_3[1] = p_node_6;

        Condition::Pointer p_cond_0 = rModelPart.CreateNewCondition("Condition2D2N", 1, PointerVector<NodeType>{condition_nodes_0}, p_elem_prop);
        Condition::Pointer p_cond_1 = rModelPart.CreateNewCondition("Condition2D2N", 2, PointerVector<NodeType>{condition_nodes_1}, p_elem_prop);
        Condition::Pointer p_cond_2 = rModelPart.CreateNewCondition("Condition2D2N", 3, PointerVector<NodeType>{condition_nodes_2}, p_elem_prop);
        Condition::Pointer p_cond_3 = rModelPart.CreateNewCondition("Condition2D2N", 4, PointerVector<NodeType>{condition_nodes_3}, p_elem_prop);

        // Now we create the "elements"
        std::vector<NodeType::Pointer> element_nodes_0 (3);
        element_nodes_0[0] = p_node_1;
        element_nodes_0[1] = p_node_2;
        element_nodes_0[2] = p_node_3;

        std::vector<NodeType::Pointer> element_nodes_1 (3);
        element_nodes_1[0] = p_node_1;
        element_nodes_1[1] = p_node_3;
        element_nodes_1[2] = p_node_4;

        std::vector<NodeType::Pointer> element_nodes_2 (3);
        element_nodes_2[0] = p_node_2;
        element_nodes_2[1] = p_node_5;
        element_nodes_2[2] = p_node_3;

        std::vector<NodeType::Pointer> element_nodes_3 (3);
        element_nodes_3[0] = p_node_5;
        element_nodes_3[1] = p_node_6;
        element_nodes_3[2] = p_node_3;

        Element::Pointer p_elem_0 = rModelPart.CreateNewElement("Element2D3N", 1, PointerVector<NodeType>{element_nodes_0}, p_elem_prop);
        Element::Pointer p_elem_1 = rModelPart.CreateNewElement("Element2D3N", 2, PointerVector<NodeType>{element_nodes_1}, p_elem_prop);
        Element::Pointer p_elem_2 = rModelPart.CreateNewElement("Element2D3N", 3, PointerVector<NodeType>{element_nodes_2}, p_elem_prop);
        Element::Pointer p_elem_3 = rModelPart.CreateNewElement("Element2D3N", 4, PointerVector<NodeType>{element_nodes_3}, p_elem_prop);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartSubModelPartsIterator, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& model_part = current_model.CreateModelPart("Main");

        model_part.CreateSubModelPart("Inlet1");
        model_part.CreateSubModelPart("Inlet2");
        model_part.CreateSubModelPart("Outlet");
        model_part.CreateSubModelPart("AnotherOutlet");

        std::size_t id = 1;
        for(auto i_SubModelPart = model_part.SubModelPartsBegin() ; i_SubModelPart != model_part.SubModelPartsEnd() ; i_SubModelPart++){
            i_SubModelPart->CreateNewNode(id++, 0.00,0.00,0.00);
        }

        KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), 4);
        KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("Inlet1").NumberOfNodes(), 1);
        KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("Outlet").NumberOfNodes(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartAddNodalSolutionStepVariable, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& model_part = current_model.CreateModelPart("Main");

        model_part.AddNodalSolutionStepVariable(VELOCITY);

        model_part.CreateNewNode(123, 0.00,0.00,0.00);

        KRATOS_CHECK_EXCEPTION_IS_THROWN(model_part.AddNodalSolutionStepVariable(PRESSURE),
            "Error: Attempting to add the variable \"PRESSURE\" to the model part with name \"Main\" which is not empty");

        model_part.AddNodalSolutionStepVariable(VELOCITY); // Adding the same Variable twice is fine bcs it wont do anything
    }

		KRATOS_TEST_CASE_IN_SUITE(ModelPartHasNodalSolutionStepVariable, KratosCoreFastSuite)
		{
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("Main");


        model_part.AddNodalSolutionStepVariable(VELOCITY);

        KRATOS_CHECK(model_part.HasNodalSolutionStepVariable(VELOCITY));
        KRATOS_CHECK_IS_FALSE(model_part.HasNodalSolutionStepVariable(PRESSURE));
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveElements, KratosCoreFastSuite)
    {
        ModelPart model_part("Main");

        // Fill model part
        GenerateGenericModelPart(model_part);

        // Set flags
        model_part.pGetElement(1)->Set(TO_ERASE, true);
        model_part.pGetElement(2)->Set(TO_ERASE, true);

        // Call method
        model_part.RemoveElements(TO_ERASE);

        // Check results
        KRATOS_CHECK(model_part.NumberOfNodes() == 6);
        KRATOS_CHECK(model_part.NumberOfElements() == 2);
        KRATOS_CHECK(model_part.NumberOfConditions() == 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveElementsAndBelongings, KratosCoreFastSuite)
    {
        ModelPart model_part("Main");

        // Fill model part
        GenerateGenericModelPart(model_part);

        // Set flags
        model_part.pGetElement(1)->Set(TO_ERASE, true);
        model_part.pGetElement(2)->Set(TO_ERASE, true);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(model_part);
        aux_util.RemoveElementsAndBelongings(TO_ERASE);

        // Check results
        KRATOS_CHECK(model_part.NumberOfNodes() == 4);
        KRATOS_CHECK(model_part.NumberOfElements() == 2);
        KRATOS_CHECK(model_part.NumberOfConditions() == 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveConditions, KratosCoreFastSuite)
    {
        ModelPart model_part("Main");

        // Fill model part
        GenerateGenericModelPart(model_part);

        // Set flags
        model_part.pGetCondition(1)->Set(TO_ERASE, true);
        model_part.pGetCondition(2)->Set(TO_ERASE, true);

        // Call method
        model_part.RemoveConditions(TO_ERASE);

        // Check results
        KRATOS_CHECK(model_part.NumberOfNodes() == 6);
        KRATOS_CHECK(model_part.NumberOfConditions() == 2);
        KRATOS_CHECK(model_part.NumberOfElements() == 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveConditionsAndBelongings, KratosCoreFastSuite)
    {
        ModelPart model_part("Main");

        // Fill model part
        GenerateGenericModelPart(model_part);

        // Set flags
        model_part.pGetCondition(1)->Set(TO_ERASE, true);
        model_part.pGetCondition(2)->Set(TO_ERASE, true);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(model_part);
        aux_util.RemoveConditionsAndBelongings(TO_ERASE);

        // Check results
        KRATOS_CHECK(model_part.NumberOfNodes() == 3);
        KRATOS_CHECK(model_part.NumberOfConditions() == 2);
        KRATOS_CHECK(model_part.NumberOfElements() == 2);
    }
}  // namespace Testing.
}  // namespace Kratos.
