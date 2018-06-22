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
    void InterfaceSearchStructure::PrepareSearchIteration(const Kratos::Flags& rOptions,
                                                    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                                    InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
    {
        // Creating the MapperInterfaceInfos
        mpMapperInterfaceInfosContainer->clear();
        mpMapperInterfaceInfosContainer->resize(1);

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
        const auto& r_mapper_interface_infos = (*mpMapperInterfaceInfosContainer)[0];
        const int num_interface_infos = r_mapper_interface_infos.size();

        // This is threadsafe bcs there will never be more than one InterfaceInfo per LocalSystem (only true in serial!)
        #pragma omp parallel for
        for (int i = 0; i<num_interface_infos; ++i)
        {
            const auto& rp_interface_info = r_mapper_interface_infos[i];

            if (rp_interface_info->GetLocalSearchWasSuccessful())
            {
                const IndexType local_sys_idx = rp_interface_info->GetLocalSystemIndex();
                KRATOS_DEBUG_ERROR_IF_NOT(local_sys_idx == static_cast<IndexType>(i))
                    << "Index mismatch!" << std::endl; // This has to be ensured in serial!
                (*mpMapperLocalSystems)[local_sys_idx]->AddInterfaceInfo(rp_interface_info);
            }
        }
    }

    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/


}  // namespace Kratos.
