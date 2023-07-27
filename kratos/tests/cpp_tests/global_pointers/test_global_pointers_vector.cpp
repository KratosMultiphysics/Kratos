//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/global_pointer.h"
#include "includes/global_pointer_variables.h"
#include "containers/global_pointers_vector.h"
#include "containers/model.h"
#include "includes/model_part.h" 
#include "includes/mpi_serializer.h" 

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(GlobalPointersContainerTest, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& mp = current_model.CreateModelPart("test");
    mp.AddNodalSolutionStepVariable(TEMPERATURE); //not to have an empty var list

    mp.CreateNewNode(1,1.0,2.0,3.0);
    mp.CreateNewNode(2,1.0,2.0,3.0);
    mp.CreateNewNode(3,1.0,2.0,3.0);

    GlobalPointersVector<Node> global_pointers_container;
    global_pointers_container.FillFromContainer(mp.Nodes());

    MpiSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);

    serializer.save("global_pointers_container", global_pointers_container);

    GlobalPointersVector<Node> new_global_pointers;
    serializer.load("global_pointers_container",new_global_pointers);

    for(std::size_t i=0; i<global_pointers_container.size(); ++i)
    {
        KRATOS_CHECK_EQUAL(&new_global_pointers[i], &global_pointers_container[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointersContainerInVariableTest, KratosCoreFastSuite)
{
    Model current_model;
    Model loaded_model;
    ModelPart& mp = current_model.CreateModelPart("test");
    mp.AddNodalSolutionStepVariable(TEMPERATURE); //not to have an empty var list

    const auto& node_1 = mp.CreateNewNode(1,1.0,2.0,3.0);
    const auto& node_2 = mp.CreateNewNode(2,1.0,2.0,3.0);
    const auto& node_3 = mp.CreateNewNode(3,1.0,2.0,3.0);

    auto target_node = mp.CreateNewNode(4,1.0,2.0,4.0);
    auto& global_pointers_container = target_node->GetValue(NEIGHBOUR_NODES);

    global_pointers_container.push_back(GlobalPointer<Node>(&*node_1));
    global_pointers_container.push_back(GlobalPointer<Node>(&*node_2));
    global_pointers_container.push_back(GlobalPointer<Node>(&*node_3));

    MpiSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);

    serializer.save("model", current_model);
    serializer.load("model", loaded_model);

    auto& new_global_pointers = loaded_model.GetModelPart("test").pGetNode(4)->GetValue(NEIGHBOUR_NODES);

    for(std::size_t i = 0; i < global_pointers_container.size(); i++)
    {
        KRATOS_CHECK_EQUAL(&new_global_pointers[i], &global_pointers_container[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointersContainerInVariableWithRecursion, KratosCoreFastSuite)
{
    Model current_model;
    Model loaded_model;
    ModelPart& mp = current_model.CreateModelPart("test");
    mp.AddNodalSolutionStepVariable(TEMPERATURE); //not to have an empty var list

    const auto& node_1 = mp.CreateNewNode(1,1.0,2.0,3.0);
    const auto& node_2 = mp.CreateNewNode(2,1.0,2.0,3.0);
    const auto& node_3 = mp.CreateNewNode(3,1.0,2.0,3.0);

    node_1->GetValue(NEIGHBOUR_NODES).push_back(GlobalPointer<Node>(&*node_2));
    node_2->GetValue(NEIGHBOUR_NODES).push_back(GlobalPointer<Node>(&*node_3));
    node_3->GetValue(NEIGHBOUR_NODES).push_back(GlobalPointer<Node>(&*node_1));

    MpiSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);

    serializer.save("model", current_model);
    serializer.load("model", loaded_model);

    KRATOS_CHECK_EQUAL(current_model.GetModelPart("test").NumberOfNodes(), loaded_model.GetModelPart("test").NumberOfNodes());

    for(std::size_t i = 1; i <= loaded_model.GetModelPart("test").NumberOfNodes(); i++) {
        auto& old_global_pointers = current_model.GetModelPart("test").pGetNode(i)->GetValue(NEIGHBOUR_NODES);
        auto& new_global_pointers = loaded_model.GetModelPart("test").pGetNode(i)->GetValue(NEIGHBOUR_NODES);

        KRATOS_CHECK_EQUAL(old_global_pointers.size(), new_global_pointers.size());

        for(std::size_t j = 0; j < new_global_pointers.size(); j++) {
            KRATOS_CHECK_EQUAL(&old_global_pointers[j], &new_global_pointers[j]);
        }
    }
}

} // namespace Kratos::Testing
