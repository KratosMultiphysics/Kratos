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
    void InterfaceSearchStructure::PrepareSearching(const Kratos::Flags& rOptions,
                                                    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                                    InterfaceObject::ConstructionType InterfaceObjectTypeOrigin,
                                                    InterfaceObject::ConstructionType InterfaceObjectTypeDestination)
    {
        if (mpInterfaceObjectsOrigin == nullptr || rOptions.Is(MapperFlags::REMESHED))
            CreateInterfaceObjectsOrigin(InterfaceObjectTypeOrigin);
        else
            UpdateInterfaceObjectsOrigin();

        if (mpLocalBinStructure == nullptr || !rOptions.Is(MapperFlags::DESTINATION_ONLY))
            InitializeBinsSearchStructure(); // This cannot be updated, has to be recreated

        if (mpInterfaceObjectsDestination == nullptr || rOptions.Is(MapperFlags::REMESHED))
            CreateInterfaceObjectsDestination(rpRefInterfaceInfo);
        else
            UpdateInterfaceObjectsDestination();
    }

    void InterfaceSearchStructure::FinalizeSearching()
    {
        const int num_interface_infos = mpMapperInterfaceInfos->size();

        // This is threadsafe bcs there will never be more than one InterfaceInfo per LocalSystem
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
    void InterfaceSearchStructure::CreateInterfaceObjectsDestination(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
    {
        const int num_objects = mpMapperLocalSystems->size();

        mpInterfaceObjectsDestination = Kratos::make_unique<InterfaceObjectContainerType>(num_objects);
        mpMapperInterfaceInfos = Kratos::make_unique<MapperInterfaceInfoPointerVectorType>(num_objects);

        #pragma omp parallel for
        for (int i = 0; i<num_objects; ++i)
        {
            const auto& r_coords = (*mpMapperLocalSystems)[i]->GetCoordinates();

            (*mpInterfaceObjectsDestination)[i] = Kratos::make_shared<InterfaceObject>(r_coords);
            (*mpMapperInterfaceInfos)[i] = rpRefInterfaceInfo->Create(i);
        }

        // Making sure that the data-structure was correctly initialized
        KRATOS_ERROR_IF_NOT(mpInterfaceObjectsDestination->size() > 0)
            << "No interface objects were created in Destination-ModelPart!" << std::endl;
    }

    void InterfaceSearchStructure::UpdateInterfaceObjectsDestination()
    {
        const int num_objects = mpMapperLocalSystems->size();

        #pragma omp parallel for
        for (int i = 0; i<num_objects; ++i)
        {
            const auto& r_coords = (*mpMapperLocalSystems)[i]->GetCoordinates();
            (*mpInterfaceObjectsDestination)[i]->UpdateCoordinates(r_coords);
        }
    }


}  // namespace Kratos.
