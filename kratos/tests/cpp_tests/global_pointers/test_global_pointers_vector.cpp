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

#include "testing/testing.h"
#include "includes/global_pointer.h"
#include "includes/global_pointer_variables.h"
#include "containers/global_pointers_vector.h"
#include "containers/model.h"
#include "includes/model_part.h" 
#include "includes/mpi_serializer.h" 


namespace Kratos {
namespace Testing {


KRATOS_TEST_CASE_IN_SUITE(GlobalPointersContainerTest, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& mp = current_model.CreateModelPart("test");
    mp.AddNodalSolutionStepVariable(TEMPERATURE); //not to have an empty var list

    mp.CreateNewNode(1,1.0,2.0,3.0);
    mp.CreateNewNode(2,1.0,2.0,3.0);
    mp.CreateNewNode(3,1.0,2.0,3.0);

    GlobalPointersVector<Node<3>> global_pointers_container;
    global_pointers_container.FillFromContainer(mp.Nodes());

    MpiSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);

    serializer.save("global_pointers_container", global_pointers_container);

    GlobalPointersVector<Node<3>> new_global_pointers;
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

    global_pointers_container.push_back(GlobalPointer<Node<3>>(&*node_1));
    global_pointers_container.push_back(GlobalPointer<Node<3>>(&*node_2));
    global_pointers_container.push_back(GlobalPointer<Node<3>>(&*node_3));

    MpiSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);

    serializer.save("model", current_model);
    serializer.load("model", loaded_model);

    auto& new_global_pointers = loaded_model.GetModelPart("test").pGetNode(4)->GetValue(NEIGHBOUR_NODES);

    for(std::size_t i=0; i<global_pointers_container.size(); ++i)
    {
        KRATOS_CHECK_EQUAL(&new_global_pointers[i], &global_pointers_container[i]);
    }
}

} // namespace Testing
} // namespace Kratos
