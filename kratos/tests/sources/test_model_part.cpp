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

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        r_model_part.CreateSubModelPart("Inlet1");
        r_model_part.CreateSubModelPart("Inlet2");
        r_model_part.CreateSubModelPart("Outlet");
        r_model_part.CreateSubModelPart("AnotherOutlet");

        std::size_t id = 1;
        for(auto i_SubModelPart = r_model_part.SubModelPartsBegin() ; i_SubModelPart != r_model_part.SubModelPartsEnd() ; i_SubModelPart++){
            i_SubModelPart->CreateNewNode(id++, 0.00,0.00,0.00);
        }

        KRATOS_CHECK_EQUAL(r_model_part.NumberOfNodes(), 4);
        KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("Inlet1").NumberOfNodes(), 1);
        KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("Outlet").NumberOfNodes(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartAddNodalSolutionStepVariable, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        r_model_part.AddNodalSolutionStepVariable(VELOCITY);

        r_model_part.CreateNewNode(123, 0.00,0.00,0.00);

        KRATOS_CHECK_EXCEPTION_IS_THROWN(r_model_part.AddNodalSolutionStepVariable(PRESSURE),
            "Error: Attempting to add the variable \"PRESSURE\" to the model part with name \"Main\" which is not empty");

        r_model_part.AddNodalSolutionStepVariable(VELOCITY); // Adding the same Variable twice is fine bcs it wont do anything
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartHasNodalSolutionStepVariable, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");


        r_model_part.AddNodalSolutionStepVariable(VELOCITY);

        KRATOS_CHECK(r_model_part.HasNodalSolutionStepVariable(VELOCITY));
        KRATOS_CHECK_IS_FALSE(r_model_part.HasNodalSolutionStepVariable(PRESSURE));
    }


    KRATOS_TEST_CASE_IN_SUITE(ModelPartEmptyName, KratosCoreFastSuite)
    {
        Model current_model;

        // Constructor with name
        KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart(""),
            "Error: Please don't use empty names (\"\") when creating a ModelPart");

        // Constructor with name and bufferSize
        KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart("", 2),
            "Error: Please don't use empty names (\"\") when creating a ModelPart");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartNameContainingPoint, KratosCoreFastSuite)
    {
        Model current_model;

        // Constructor with name
        KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart("name.other"),
            "Error: Please don't use names containing (\".\") when creating a ModelPart");

        // Constructor with name and bufferSize
        KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart("name.other", 2),
            "Error: Please don't use names containing (\".\") when creating a ModelPart");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveElements, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetElement(1)->Set(TO_ERASE, true);
        r_model_part.pGetElement(2)->Set(TO_ERASE, true);

        // Call method
        r_model_part.RemoveElements(TO_ERASE);

        // Check results
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 6);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 2);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveElementsAndBelongings, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetElement(1)->Set(TO_ERASE, true);
        r_model_part.pGetElement(2)->Set(TO_ERASE, true);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_model_part);
        aux_util.RemoveElementsAndBelongings(TO_ERASE);

        // Check results
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 4);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 2);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveConditions, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetCondition(1)->Set(TO_ERASE, true);
        r_model_part.pGetCondition(2)->Set(TO_ERASE, true);

        // Call method
        r_model_part.RemoveConditions(TO_ERASE);

        // Check results
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 6);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveConditionsAndBelongings, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetCondition(1)->Set(TO_ERASE, true);
        r_model_part.pGetCondition(2)->Set(TO_ERASE, true);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_model_part);
        aux_util.RemoveConditionsAndBelongings(TO_ERASE);

        // Check results
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 3);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 2);
    }
  }  // namespace Testing.
}  // namespace Kratos.
