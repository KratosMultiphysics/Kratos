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
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{

FindNodalNeighboursProcess::FindNodalNeighboursProcess(ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
    mpNodeNeighboursCalculator = Kratos::make_unique<FindGlobalNodalNeighboursProcess>(mrModelPart);
    mpElemNeighboursCalculator = Kratos::make_unique<FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>>(mrModelPart);

    KRATOS_INFO("FindNodalNeighboursProcess") <<
        R"(please call separetely FindGlobalNodalNeighboursProcess
        and FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>.
        The two calculations are currently independent,
            hence memory savings can be achieved)" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

FindNodalNeighboursProcess::FindNodalNeighboursProcess(
    ModelPart& rModelPart,
    const SizeType AverageElements,
    const SizeType AverageNodes
    ) : FindNodalNeighboursProcess(rModelPart)
{
    KRATOS_WARNING("FindNodalNeighboursProcess") << "Parameters AverageElements and AverageNodes are currently ignored. This constructor will be removed on the 2 of April 2020" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void FindNodalNeighboursProcess::Execute()
{
    mpNodeNeighboursCalculator->Execute();
    mpElemNeighboursCalculator->Execute();
}

/***********************************************************************************/
/***********************************************************************************/

void FindNodalNeighboursProcess::ClearNeighbours()
{
    mpNodeNeighboursCalculator->ClearNeighbours();
    mpElemNeighboursCalculator->Clear();
}

}  // namespace Kratos.


