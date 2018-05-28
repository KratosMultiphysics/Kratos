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
            CreateInterfaceObjectsDestination(InterfaceObjectTypeDestination);
        else
            UpdateInterfaceObjectsDestination();
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
    void InterfaceSearchStructure::CreateInterfaceObjectsDestination(InterfaceObject::ConstructionType InterfaceObjectTypeDestination)
    {
        const int num_objects = mpMapperLocalSystems->size();

        mpInterfaceObjectsDestination = Kratos::make_unique<InterfaceObjectContainerType>(num_objects);
        mpMapperInterfaceInfos = Kratos::make_unique<MapperInterfaceInfoPointerVectorType>(num_objects);


        #pragma omp parallel for
        for (int i = 0; i< num_objects; ++i)
        {
            const auto& r_coords = (*mpMapperLocalSystems)[i]->GetCoordinates();

            // mpInterfaceObjectsDestination = Kratos::make_shared<InterfaceObject>(r_coords)

            // auto it_node = nodes_begin + i;
            // (*mpInterfaceObjectsDestination)[i] = Kratos::make_unique<InterfaceNode>(*(it_node));
        }


        //         rNumObjects = CoordinateListSize / 3;

        // for (int i = 0; i < rNumObjects; ++i)   // create InterfaceObjects
        // {
        //     rRemotePointList[i] = InterfaceObject::Pointer(new InterfaceObject(
        //                               pCoordinateList[(i * 3) + 0], pCoordinateList[(i * 3) + 1], pCoordinateList[(i * 3) + 2]));
        // }



        // Making sure that the data-structure was correctly initialized
        KRATOS_ERROR_IF_NOT(mpInterfaceObjectsDestination->size() > 0)
            << "No interface objects were created in Destination-ModelPart!" << std::endl;
    }

    void InterfaceSearchStructure::UpdateInterfaceObjectsDestination()
    {
        const auto begin = mpInterfaceObjectsDestination->begin();

        // TODO uncomment
        // Commented bcs it would cause a
        // #pragma omp parallel for a segfault due to mpInterfaceObjectsDestination being uninitialized
        // for (int i = 0; i< static_cast<int>(mpInterfaceObjectsDestination->size()); ++i)
        //     (*(begin + i))->UpdateCoordinates();
    }


}  // namespace Kratos.
