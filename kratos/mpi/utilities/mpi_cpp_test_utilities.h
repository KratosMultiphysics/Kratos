//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "utilities/pointer_communicator.h"
#include "utilities/global_pointer_utilities.h"
#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos
{

class KRATOS_API(KRATOS_MPI_CORE) MPICppTestUtilities
{
public:
/**
* @brief It generates a truss structure with an expected solution
* @param rModelPart The model part to be filled
* @param rDataCommunicator The data communicator
*/
static void GenerateDistributedBarStructure(
    ModelPart& rModelPart,
    const DataCommunicator& rDataCommunicator
    );

/**
 * @brief Synchronizes local pointers to global pointers of the given ModelPart using the provided
 * @details DataCommunicator and returns a shared pointer to a GlobalPointerCommunicator<Node>.
 * @param rModelPart the ModelPart to be synchronized
 * @param rGlobalPointers the GlobalPointersVector to be filled with the synchronized pointers
 * @param rDataCommunicator the DataCommunicator to be used for synchronization
 * @return A shared pointer to a GlobalPointerCommunicator<Node> containing synchronized pointers
 */
static GlobalPointerCommunicator<Node>::Pointer SynchronizeNodes(
    ModelPart& rModelPart,
    GlobalPointersVector<Node>& rGlobalPointers,
    const DataCommunicator& rDataCommunicator
    );

};

} // namespace Kratos