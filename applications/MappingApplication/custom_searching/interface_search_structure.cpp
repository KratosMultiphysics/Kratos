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
#include "includes/model_part.h"
#include "interface_search_structure.h"
#include "custom_utilities/mapper_flags.h"

namespace Kratos
{
    using SizeType = std::size_t;
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/
    void InterfaceSearchStructure::PrepareSearching(const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo,
                                                    InterfaceObject::ConstructionType InterfaceObjectTypeDestination)
    {
        if (mpInterfaceObjectsDestination != nullptr)
        {
            mpInterfaceObjectsDestination = Kratos::make_unique<InterfaceObjectContainerType>();
            // mpInterfaceInfos = Kratos::make_unique<InterfaceObjectContainerType>();



        }
    }

    void InterfaceSearchStructure::FinalizeSearching()
    {
        /*
        TODO change to OMP
        This is threadsafe bcs there will never be more than one InterfaceInfo per LocalSystem (per partition, but the mpi stuff is handled differently anyway!
        => will also make it easier to assign the shared_ptrs to the LocalSystem
        for (const auto& r_info : Infos)
        {
            if (info.GetLocalSearchWasSuccessful() == true) // Attention, do NOT do this check in MPI, the InterfaceInfo would not have been sent to a partition if it didn't have a successful local search!
            {
                IndexType local_sys_idx =
                mpMapperLocalSystems[local_sys_idx].AddInterfaceInfo(r_info);
            }
        }

        */

    }

    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/


}  // namespace Kratos.
