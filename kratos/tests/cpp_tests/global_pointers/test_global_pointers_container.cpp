//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/global_pointer.h"
#include "containers/global_pointers_vector.h"
#include "containers/global_pointers_unordered_map.h"
#include "containers/model.h"
#include "includes/model_part.h" 
#include "includes/mpi_serializer.h" 

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(GlobalPointersVectorTest, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& mp = current_model.CreateModelPart("test");
    mp.AddNodalSolutionStepVariable(TEMPERATURE); //not to have an empty var list

    mp.CreateNewNode(1,1.0,2.0,3.0);
    mp.CreateNewNode(2,2.0,2.0,3.0);
    mp.CreateNewNode(3,3.0,2.0,3.0);

    GlobalPointersVector<Node> global_pointers_container;
    global_pointers_container.FillFromContainer(mp.Nodes());

    MpiSerializer serializer;
    serializer.save("global_pointers_container", global_pointers_container);

    GlobalPointersVector<Node> new_global_pointers;
    serializer.load("global_pointers_container",new_global_pointers);

    for(std::size_t i=0; i<global_pointers_container.size(); ++i)
    {
        KRATOS_EXPECT_EQ(&new_global_pointers[i], &global_pointers_container[i]);
    }
    

};

KRATOS_TEST_CASE_IN_SUITE(GlobalPointersUnorderedMapTest , KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& mp = current_model.CreateModelPart("test");
    mp.AddNodalSolutionStepVariable(TEMPERATURE); //not to have an empty var list

    GlobalPointer<Node> gp1( mp.CreateNewNode(1,1.0,2.0,3.0).get());
    GlobalPointer<Node> gp2( mp.CreateNewNode(2,2.0,2.0,3.0).get());
    GlobalPointer<Node> gp3( mp.CreateNewNode(3,3.0,2.0,3.0).get());

    GlobalPointersUnorderedMap<Node, GlobalPointersVector<Node>> global_pointers_map;
    global_pointers_map[ gp1 ] = {gp1};
    global_pointers_map[ gp2 ] = {gp1,gp2};
    global_pointers_map[ gp3 ] = {gp1,gp2,gp3};

    MpiSerializer serializer;
    serializer.save("global_pointers_map", global_pointers_map);
    global_pointers_map.clear();  

    GlobalPointersUnorderedMap<Node, GlobalPointersVector<Node>> new_global_pointers;
    serializer.load("global_pointers_map",new_global_pointers);

    KRATOS_EXPECT_EQ(&new_global_pointers[gp1][0], &*gp1);   
    KRATOS_EXPECT_EQ(&new_global_pointers[gp2][0], &*gp1);
    KRATOS_EXPECT_EQ(&new_global_pointers[gp2][1], &*gp2);
    KRATOS_EXPECT_EQ(&new_global_pointers[gp3][0], &*gp1);
    KRATOS_EXPECT_EQ(&new_global_pointers[gp3][1], &*gp2);
    KRATOS_EXPECT_EQ(&new_global_pointers[gp3][2], &*gp3);

};

KRATOS_TEST_CASE_IN_SUITE(GlobalPointersVectorDeepSerializationTest , KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& mp = current_model.CreateModelPart("test");
    mp.AddNodalSolutionStepVariable(TEMPERATURE); //not to have an empty var list

    mp.CreateNewNode(1,1.0,2.0,3.0)->GetSolutionStepValue(TEMPERATURE) = 100.0;
    mp.CreateNewNode(2,2.0,2.0,3.0)->GetSolutionStepValue(TEMPERATURE) = 200.0;
    mp.CreateNewNode(3,3.0,2.0,3.0)->GetSolutionStepValue(TEMPERATURE) = 300.0;

    GlobalPointersVector<Node> global_pointers_container;
    global_pointers_container.FillFromContainer(mp.Nodes());

    StreamSerializer serializer; //NOTE: here StreamSerializer is used, hence global pointers are deep copied
    serializer.save("global_pointers_container", global_pointers_container);

    current_model.Reset();
    global_pointers_container.clear();

    GlobalPointersVector<Node> new_global_pointers;
    serializer.load("global_pointers_container",new_global_pointers);

    for(std::size_t i=0; i<new_global_pointers.size(); ++i)
    {
        KRATOS_EXPECT_EQ(new_global_pointers[i].FastGetSolutionStepValue(TEMPERATURE), (i+1)*100.0);
    }
};

} // namespace Kratos::Testing
