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
using IndexType = std::size_t;
/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/


/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/
void InterfaceSearchStructure::InitializeSearchIteration(const Kratos::Flags& rOptions,
                                                const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    // Creating the MapperInterfaceInfos
    (*mpMapperInterfaceInfosContainer)[0].clear();

    auto& r_mapper_interface_infos = (*mpMapperInterfaceInfosContainer)[0];
    r_mapper_interface_infos.reserve(mpMapperLocalSystems->size());

    IndexType local_sys_idx = 0;
    for (const auto& r_local_sys : (*mpMapperLocalSystems))
    {
        if (!r_local_sys->HasInterfaceInfo()) // Only the local_systems that have not received an InterfaceInfo create a new one
        {
            const auto& r_coords = r_local_sys->Coordinates();
            r_mapper_interface_infos.push_back(rpRefInterfaceInfo->Create(r_coords, local_sys_idx));
        }
        ++local_sys_idx;
    }
}

void InterfaceSearchStructure::FinalizeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo)
{
    InterfaceSearchStructureBase::FilterInterfaceInfosSuccessfulSearch();
    InterfaceSearchStructureBase::AssignInterfaceInfos();
}

/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/

}  // namespace Kratos.
