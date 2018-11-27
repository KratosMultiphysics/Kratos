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
#include "containers/global_pointers_vector.h"
#include "containers/global_pointers_unordered_map.h"
#include "containers/model.h"
#include "includes/model_part.h" 
#include "includes/stream_serializer.h" 


namespace Kratos {
namespace Testing {


KRATOS_TEST_CASE_IN_SUITE(GlobalPointersVectorTest, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& mp = current_model.CreateModelPart("test");

    mp.CreateNewNode(1,1.0,2.0,3.0);
    mp.CreateNewNode(2,1.0,2.0,3.0);
    mp.CreateNewNode(3,1.0,2.0,3.0);

    GlobalPointersVector<Node<3>> global_pointers_container;
    global_pointers_container.FillFromContainer(mp.Nodes());

    StreamSerializer serializer;
    serializer.save("global_pointers_container", global_pointers_container);

    GlobalPointersVector<Node<3>> new_global_pointers;
    serializer.load("global_pointers_container",new_global_pointers);

    for(std::size_t i=0; i<global_pointers_container.size(); ++i)
    {
        KRATOS_CHECK_EQUAL(&*new_global_pointers[i], &*global_pointers_container[i]);
    }
    

};

KRATOS_TEST_CASE_IN_SUITE(GlobalPointersUnorderedMapTest , KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& mp = current_model.CreateModelPart("test");

    auto pn1 = mp.CreateNewNode(1,1.0,2.0,3.0);
    auto pn2 = mp.CreateNewNode(2,1.0,2.0,3.0);
    auto pn3 = mp.CreateNewNode(3,1.0,2.0,3.0);

    GlobalPointersUnorderedMap<Node<3>, double> global_pointers_map;
    global_pointers_map[ GlobalPointer<Node<3>>(pn1) ] = 10.0;
    global_pointers_map[ GlobalPointer<Node<3>>(pn2) ] = 20.0;
    global_pointers_map[ GlobalPointer<Node<3>>(pn3) ] = 30.0;

    StreamSerializer serializer;
    serializer.save("global_pointers_map", global_pointers_map);

    GlobalPointersVector<Node<3>> new_global_pointers;
    serializer.load("global_pointers_map",global_pointers_map);

    KRATOS_CHECK_EQUAL(global_pointers_map[pn1], 10.0);
    

};


} // namespace Testing
} // namespace Kratos
