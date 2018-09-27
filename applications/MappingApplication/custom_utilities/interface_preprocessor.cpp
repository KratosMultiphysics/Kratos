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
#include "interface_preprocessor.h"


namespace Kratos
{
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/
    void InterfacePreprocessor::CreateMapperLocalSystems(const MapperLocalSystemPointer& rpLocalSystem)
    {
        mpMapperLocalSystems->clear();

        const bool use_nodes = rpLocalSystem->UseNodesAsBasis();

        if (use_nodes) {
            CreateMapperLocalSystemsFromNodes(rpLocalSystem);
        }
        else {
            CreateMapperLocalSystemsFromGeometries(rpLocalSystem);
        }

        int num_local_systems = mpMapperLocalSystems->size(); // int bcs of MPI
        mrModelPartDestination.GetCommunicator().SumAll(num_local_systems);

        KRATOS_ERROR_IF_NOT(num_local_systems > 0)
            << "No mapper local systems were created in Destination-ModelPart \""
            << mrModelPartDestination.Name() << "\"!" << std::endl;
    }
    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/

    void InterfacePreprocessor::CreateMapperLocalSystemsFromNodes(const MapperLocalSystemPointer& rpLocalSystem)
    {
        const std::size_t num_nodes = mrModelPartDestination.GetCommunicator().LocalMesh().NumberOfNodes();
        const auto nodes_ptr_begin = mrModelPartDestination.GetCommunicator().LocalMesh().Nodes().ptr_begin();

        if (mpMapperLocalSystems->size() != num_nodes) {
            mpMapperLocalSystems->resize(num_nodes);
        }

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(num_nodes); ++i) {
            auto it_node = nodes_ptr_begin + i;
            (*mpMapperLocalSystems)[i] = rpLocalSystem->Create(*(it_node));
        }
    }

    void InterfacePreprocessor::CreateMapperLocalSystemsFromGeometries(const MapperLocalSystemPointer& rpLocalSystem)
    {
        KRATOS_ERROR << "This function is not implemented yet" << std::endl;
    }


}  // namespace Kratos.
