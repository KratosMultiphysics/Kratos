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
    InterfacePreprocessor::InterfacePreprocessor(ModelPart& rModelPartDestination,
                                                 MapperLocalSystemPointerVectorPointer pMapperLocalSystems)
        : mrModelPartDestination(rModelPartDestination),
          mpMapperLocalSystems(pMapperLocalSystems)
    {
    }

    void InterfacePreprocessor::GenerateInterfaceModelPart(const MapperLocalSystemPointer& rpLocalSystem)
    {
        mpMapperLocalSystems->clear();

        const bool use_nodes = rpLocalSystem->UseNodesAsBasis();

        if (use_nodes)
            CreateMapperLocalSystemsFromNodes(rpLocalSystem);
        else
            CreateMapperLocalSystemsFromGeometries(rpLocalSystem);
    }
    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/

    void InterfacePreprocessor::CreateMapperLocalSystemsFromNodes(const MapperLocalSystemPointer& rpLocalSystem)
    {
        const std::size_t num_nodes = mrModelPartDestination.NumberOfNodes();

        mpMapperLocalSystems->resize(num_nodes);

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(num_nodes); ++i)
        {
            auto it_node = mrModelPartDestination.Nodes().begin() + i;
            (*mpMapperLocalSystems)[i] = rpLocalSystem->Create(*(it_node));
        }
    }

    void InterfacePreprocessor::CreateMapperLocalSystemsFromGeometries(const MapperLocalSystemPointer& rpLocalSystem)
    {
        KRATOS_ERROR << "This function is not implemented yet" << std::endl;
    }


}  // namespace Kratos.
