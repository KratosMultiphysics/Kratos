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
#include "interface_search_structure.h"
#include "custom_searching/interface_node.h"
#include "custom_searching/interface_geometry_object.h"

namespace Kratos
{
    using SizeType = std::size_t;
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/
    void InterfaceSearchStructure::PrepareSearching(InterfaceObject::ConstructionType InterfaceObjectTypeDestination)
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
    void InterfaceSearchStructure::CreateInterfaceObjectsOrigin(InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
    {
        mpInterfaceObjectsOrigin = Kratos::make_unique<InterfaceObjectContainerType>();

        if (InterfaceObjectTypeOrigin == InterfaceObject::Node_Coords)
        {
            const SizeType num_nodes = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes();
            const auto nodes_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Nodes().begin();

            mpInterfaceObjectsOrigin->resize(num_nodes);

            #pragma omp parallel for
            for (int i = 0; i< static_cast<int>(num_nodes); ++i)
            {
                auto it_node = nodes_begin + i;
                (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceNode>(*(it_node), 0);
            }
        }
        else if (InterfaceObjectTypeOrigin == InterfaceObject::Geometry_Center)
        {
            const SizeType num_elements = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfElements();
            const SizeType num_conditions = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfConditions();

            const auto elements_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Elements().begin();
            const auto conditions_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Conditions().begin();

            mpInterfaceObjectsOrigin->resize(num_elements+num_conditions); // one of them has to be zero!!!

            #pragma omp parallel for
            for (int i = 0; i< static_cast<int>(num_elements); ++i)
            {
                auto it_elem = elements_begin + i;
                (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceGeometryObject>(it_elem->GetGeometry(),
                    0.0,0,0);
            }
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int>(num_conditions); ++i)
            {
                auto it_cond = conditions_begin + i;
                (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceGeometryObject>(it_cond->GetGeometry(),
                    0.0,0,0);
            }
        }
        else
        {
            KRATOS_ERROR << "Type of interface object construction not implemented" << std::endl;
        }

        // Making sure that the data-structure was correctly initialized
        int num_interface_objects = mpInterfaceObjectsOrigin->size(); // int bcs of MPI
        mrModelPartOrigin.GetCommunicator().SumAll(num_interface_objects);

        KRATOS_ERROR_IF_NOT(num_interface_objects > 0)
            << "No interface objects were created in Origin-ModelPart \""
            << mrModelPartOrigin.Name() << "\"!" << std::endl;
    }

    void InterfaceSearchStructure::InitializeBinsSearchStructure()
    {
        if (mpInterfaceObjectsOrigin->size() > 0)   // only construct the bins if the partition has a part of the interface
        {
            mpLocalBinStructure = Kratos::make_unique<BinsObjectDynamic<InterfaceObjectConfigure>>(
                mpInterfaceObjectsOrigin->begin(), mpInterfaceObjectsOrigin->end());
        }
    }

    void InterfaceSearchStructure::ConductLocalSearch()
    {
        // This function finds neighbors of the InterfaceObjects in rInterfaceObjects in bin_structure
        // It must be executable by serial and parallel version!
        // InterfaceObjectsSize must be passed bcs rInterfaceObjects might contain old entries (it has
        // the max receive buffer size as size)!

        SizeType num_interface_obj_bin = mpInterfaceObjectsOrigin->size();

        if (num_interface_obj_bin > 0)   // this partition has a bin structure
        {
            // InterfaceObjectConfigure::ResultContainerType neighbor_results(num_interface_obj_bin);
            // std::vector<double> neighbor_distances(num_interface_obj_bin);

            // InterfaceObjectConfigure::IteratorType interface_object_itr;
            // InterfaceObjectConfigure::ResultIteratorType results_itr;
            // std::vector<double>::iterator distance_itr;

            // //   Searching the neighbors
            // for (int i = 0; i < InterfaceObjectsSize; ++i)
            // {
            //     interface_object_itr = rInterfaceObjects.begin() + i;
            //     double search_radius = mSearchRadius; // reset search radius

            //     results_itr = neighbor_results.begin();
            //     distance_itr = neighbor_distances.begin();

            //     SizeType number_of_results = mpLocalBinStructure->SearchObjectsInRadius(
            //                                         *interface_object_itr, search_radius, results_itr,
            //                                         distance_itr, num_interface_obj_bin);

            //     if (number_of_results > 0)   // neighbors were found
            //     {
            //         SelectBestResult(interface_object_itr, neighbor_results,
            //                          neighbor_distances, number_of_results,
            //                          rInterfaceObjectResults[i], rMinDistances[i],
            //                          rShapeFunctionValues[i], rPairingIndices[i]);
            //     }
            //     else
            //     {
            //         rMinDistances[i] = -1.0f; // indicates that the search was not succesful
            //         rInterfaceObjectResults[i].reset(); // Release an old pointer, that is probably existing from a previous search
            //     }
            // }
        }
        else     // this partition has no part of the point receiving interface, i.e. the origin of the mapped values
        {
            // Actually this should never happen, since the bounding box shouldn't allow sending to this partition if it doesn't
            // have a part of the Interface



            // for (int i = 0; i < InterfaceObjectsSize; ++i)   // no results in this partition
            // {
            //     rMinDistances[i] = -1.0f; // indicates that the search was not succesful
            //     rInterfaceObjectResults[i].reset(); // Release an old pointer, that is probably existing from a previous search
            // }
        }
    }



}  // namespace Kratos.
