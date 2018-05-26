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
#include "interface_search_structure_mpi.h"

namespace Kratos
{
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/
    void InterfaceSearchStructureMPI::PrepareSearching(const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo,
                                                    InterfaceObject::ConstructionType InterfaceObjectTypeDestination)
    {
        // TODO pass mappingoptions!
        // TODO print info saying that ORIGIN_ONLY has no effect in MPI, the destination also has to be updated
        // => since the origin changes also the destination might be sent to other partitions!

        /*
        1. Check which partitions have part of Interface (or send BBoxes directly...? => Don't think this is a good solution)
        2. Compute Graph for DataExchange
        3. Send BBoxes to partitions that have part of Interface
        4. Compute CandidatePartitions
        5. Send stuff to CandidatePartitions
        6. Create Objects in the CandidatePartitions
        Afterwards do local search ...
        */
    }

        void InterfaceSearchStructureMPI::FinalizeSearching()
    {
        /*
        1. Check with which Partitions I have to communicate
        2. Exchange this info
        3. Compute CommunicationGraph
        4. Searialize InterfaceInfos and exchange them
        5. Deserialize InterfaceInfos and assign them to the LocalSystems
        */

    }

    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/


}  // namespace Kratos.
