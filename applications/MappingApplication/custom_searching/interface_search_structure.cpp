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
    void InterfaceSearchStructure::PrepareSearch(const Kratos::Flags& rOptions,
                                        const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                        InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
    {
        if (mpInterfaceObjectsOrigin == nullptr || rOptions.Is(MapperFlags::REMESHED))
            CreateInterfaceObjectsOrigin(InterfaceObjectTypeOrigin);
        else
            UpdateInterfaceObjectsOrigin();

        if (mpLocalBinStructure == nullptr || !rOptions.Is(MapperFlags::DESTINATION_ONLY))
            InitializeBinsSearchStructure(); // This cannot be updated, has to be recreated
    }

    void InterfaceSearchStructure::FinalizeSearch()
    {
        mpMapperInterfaceInfos->clear();
    }


    void InterfaceSearchStructure::PrepareSearchIteration(const Kratos::Flags& rOptions,
                                                    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                                    InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
    {
        CreateInterfaceInfos(rpRefInterfaceInfo);
    }

    void InterfaceSearchStructure::FinalizeSearchIteration()
    {
        const int num_interface_infos = mpMapperInterfaceInfos->size();

        // This is threadsafe bcs there will never be more than one InterfaceInfo per LocalSystem (only true in serial!)
        #pragma omp parallel for
        for (int i = 0; i<num_interface_infos; ++i)
        {
            const auto& r_interface_info = (*mpMapperInterfaceInfos)[i];

            if (r_interface_info->GetLocalSearchWasSuccessful())
            {
                const IndexType local_sys_idx = r_interface_info->GetLocalSystemIndex();
                KRATOS_DEBUG_ERROR_IF_NOT(local_sys_idx == static_cast<IndexType>(i))
                    << "Index mismatch!" << std::endl; // This has to be ensured in serial!
                (*mpMapperLocalSystems)[local_sys_idx]->AddInterfaceInfo(r_interface_info);
            }
        }
    }

    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/
    void InterfaceSearchStructure::CreateInterfaceInfos(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
    {
        mpMapperInterfaceInfos->clear();
        mpMapperInterfaceInfos->reserve(mpMapperLocalSystems->size());

        IndexType local_sys_idx = 0;
        for (const auto& r_local_sys : (*mpMapperLocalSystems))
        {
            if (!r_local_sys->HasInterfaceInfo())
            {
                const auto& r_coords = r_local_sys->GetCoordinates();
                (*mpMapperInterfaceInfos).push_back(rpRefInterfaceInfo->Create(r_coords, local_sys_idx));
            }
            ++local_sys_idx;
        }
    }


}  // namespace Kratos.
