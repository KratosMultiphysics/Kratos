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
#include "containers/model.h"
#include "includes/model_part.h" 
#include "includes/stream_serializer.h" 


namespace Kratos {
namespace Testing {


KRATOS_TEST_CASE_IN_SUITE(GlobalPointersContainerTest, KratosCoreFastSuite)
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
    

}

} // namespace Testing
} // namespace Kratos
