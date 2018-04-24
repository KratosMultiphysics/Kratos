//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "mapping_application_variables.h"
#include "interface_communicator.h"

namespace Kratos
{
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/
    InterfaceCommunicator::InterfaceCommunicator(ModelPart& rModelPartOrigin,
                                                 ModelPart::Pointer mpInterfaceModelPart)
        : mrModelPartOrigin(rModelPartOrigin),
          mpInterfaceModelPart(mpInterfaceModelPart)
    {
        PrepareInterface();

    }

    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/

    void InterfaceCommunicator::PrepareInterface()
    {
        AssignInterfaceEquationIds();
    }

    void InterfaceCommunicator::AssignInterfaceEquationIds()
    {
        AssignInterfaceEquationIds(mrModelPartOrigin.GetCommunicator());
        AssignInterfaceEquationIds(mpInterfaceModelPart->GetCommunicator());
    }

    void InterfaceCommunicator::AssignInterfaceEquationIds(Communicator& rModelPartCommunicator)
    {
        const int num_nodes_local = rModelPartCommunicator.LocalMesh().NumberOfNodes();

        int num_nodes_accumulated;

        rModelPartCommunicator.ScanSum(num_nodes_local, num_nodes_accumulated);

        const int start_equation_id = num_nodes_accumulated - num_nodes_local;

        const auto nodes_begin = rModelPartCommunicator.LocalMesh().NodesBegin();

        #pragma omp parallel for
        for (int i=0; i<num_nodes_local; ++i)
        {
            // std::cout << "INTERFACE_EQUATION_ID = " << start_equation_id + i << std::endl;
            ( nodes_begin + i )->SetValue(INTERFACE_EQUATION_ID, start_equation_id + i);
        }

        rModelPartCommunicator.SynchronizeVariable(INTERFACE_EQUATION_ID);
    }

    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/


}  // namespace Kratos.
