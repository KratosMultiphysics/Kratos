//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "mapping_application_variables.h"
#include "custom_utilities/interface_communicator.h"

namespace Kratos {
namespace Testing {

void CreateNodesForMapping(ModelPart& rModelPart, const int NumNodes)
{
    for (int i=0; i< NumNodes; ++i)
        rModelPart.CreateNewNode(i+1, i*0.1, i*0.2, i*0.3);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCommunicatorInterfaceEquationIds, KratosMappingApplicationGeneralTestSuite)
{
    const int num_nodes_1 = 11;
    const int num_nodes_2 = 23;
    ModelPart model_part("ForTest");

    ModelPart::Pointer p_model_part = Kratos::make_shared<ModelPart>("pForTest");

    CreateNodesForMapping(model_part, num_nodes_1);
    CreateNodesForMapping(*p_model_part, num_nodes_2);

    // std::cout << model_part << std::endl;
    // std::cout << *p_model_part << std::endl;

    InterfaceCommunicator interface_comm(model_part, p_model_part);

    int idx = 0;
    for (const auto& r_node : model_part.Nodes())
    {
        // std::cout << "Idx = " << idx << " ; INTERFACE_EQUATION_ID = "
        //     << r_node.GetValue(INTERFACE_EQUATION_ID) << std::endl; // Print for debugging

        KRATOS_CHECK_EQUAL(idx, r_node.GetValue(INTERFACE_EQUATION_ID));
        idx += 1;
    }

    idx = 0;
    for (const auto& r_node : p_model_part->Nodes())
    {
        // std::cout << "Idx = " << idx << " ; 222 INTERFACE_EQUATION_ID = "
        //     << r_node.GetValue(INTERFACE_EQUATION_ID) << std::endl; // Print for debugging

        KRATOS_CHECK_EQUAL(idx, r_node.GetValue(INTERFACE_EQUATION_ID));
        idx += 1;
    }


    // KRATOS_CHECK(false);

}

}  // namespace Testing
}  // namespace Kratos