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
#include "interface_search_structure_base.h"

namespace Kratos
{
    using SizeType = std::size_t;
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/
    void InterfaceSearchStructureBase::ExchangeInterfaceData(const Communicator& rComm,
                                const Kratos::Flags& rOptions,
                                const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo,
                                InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
    {

        // mSearchRadius = SearchRadius;
        // mMaxSearchIterations = MaxSearchIterations;
        const int increase_factor = 4;
        int num_iteration = 1;
        bool last_iteration = (mMaxSearchIterations == 1) ? true : false; // true in case only one search iteration is conducted

        PrepareSearch(rOptions, rpInterfaceInfo, InterfaceObjectTypeOrigin);

        // First Iteration is done outside the search loop bcs it has
        // to be done in any case
        // one search iteration should be enough in most cases (if the search
        // radius was either computed or specified properly)
        // only if some points did not find a neighbor or dont have a valid
        // projection, more search iterations are necessary
        ConductSearchIteration(rOptions, rpInterfaceInfo, InterfaceObjectTypeOrigin);

        while (++num_iteration < mMaxSearchIterations && !AllNeighborsFound(rComm))
        {
            mSearchRadius *= increase_factor;

            if (num_iteration == mMaxSearchIterations) last_iteration = true;

            if (mEchoLevel >= 2 && mCommRank == 0)
            {
                std::cout << "MAPPER WARNING, search radius was increased, "
                          << "another search iteration is conducted, "
                          << "search iteration " << num_iteration << " / "
                          << mMaxSearchIterations << ", search radius "
                          << mSearchRadius << std::endl;
            }

            ConductSearchIteration(rOptions, rpInterfaceInfo, InterfaceObjectTypeOrigin);
        }

        FinalizeSearch();
    }

    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/
    void InterfaceSearchStructureBase::CreateInterfaceObjectsOrigin(InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
    {
        mpInterfaceObjectsOrigin = Kratos::make_unique<InterfaceObjectContainerType>();

        if (InterfaceObjectTypeOrigin == InterfaceObject::Node_Coords)
        {
            const SizeType num_nodes = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes();
            const auto nodes_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Nodes().ptr_begin();

            mpInterfaceObjectsOrigin->resize(num_nodes);

            #pragma omp parallel for
            for (int i = 0; i< static_cast<int>(num_nodes); ++i)
            {
                auto it_node = nodes_begin + i;
                (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceNode>(*(it_node));
            }
        }
        else if (InterfaceObjectTypeOrigin == InterfaceObject::Geometry_Center)
        {
            const SizeType num_elements = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfElements();
            const SizeType num_conditions = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfConditions();

            const auto elements_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Elements().ptr_begin();
            const auto conditions_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Conditions().ptr_begin();

            mpInterfaceObjectsOrigin->resize(num_elements+num_conditions); // one of them has to be zero!!!

            #pragma omp parallel for
            for (int i = 0; i< static_cast<int>(num_elements); ++i)
            {
                auto it_elem = elements_begin + i;
                (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceGeometryObject>((*it_elem)->pGetGeometry());
            }
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int>(num_conditions); ++i)
            {
                auto it_cond = conditions_begin + i;
                (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceGeometryObject>((*it_cond)->pGetGeometry());
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

    void InterfaceSearchStructureBase::UpdateInterfaceObjectsOrigin()
    {
        const auto begin = mpInterfaceObjectsOrigin->begin();

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(mpInterfaceObjectsOrigin->size()); ++i)
            (*(begin + i))->UpdateCoordinates();
    }

    void InterfaceSearchStructureBase::InitializeBinsSearchStructure()
    {
        if (mpInterfaceObjectsOrigin->size() > 0)   // only construct the bins if the partition has a part of the interface
        {
            mpLocalBinStructure = Kratos::make_unique<BinsObjectDynamic<InterfaceObjectConfigure>>(
                mpInterfaceObjectsOrigin->begin(), mpInterfaceObjectsOrigin->end());
        }
    }

    void InterfaceSearchStructureBase::ConductLocalSearch()
    {
        SizeType num_interface_obj_bin = mpInterfaceObjectsOrigin->size();

        if (num_interface_obj_bin > 0)   // this partition has a bin structure
        {
            InterfaceObjectConfigure::ResultContainerType neighbor_results(num_interface_obj_bin);
            std::vector<double> neighbor_distances(num_interface_obj_bin);
            auto interface_obj = Kratos::make_shared<InterfaceObject>(array_1d<double, 3>(0.0));

            // #pragma omp parallel for // TODO this requires to make some things thread-local!
            for (SizeType i = 0; i < mpMapperInterfaceInfos->size(); ++i)
            {
                const auto& r_interface_info = (*mpMapperInterfaceInfos)[i];

                interface_obj->UpdateCoordinates(r_interface_info->GetCoordinates());
                double search_radius = mSearchRadius; // reset search radius // TODO check this

                // reset the containers
                auto results_itr = neighbor_results.begin();
                auto distance_itr = neighbor_distances.begin();

                const SizeType number_of_results = mpLocalBinStructure->SearchObjectsInRadius(
                    interface_obj, search_radius, results_itr,
                    distance_itr, num_interface_obj_bin);

                for (SizeType j=0; j<number_of_results; ++j)
                    r_interface_info->ProcessSearchResult(neighbor_results[j], neighbor_distances[j]);

                // If the search did not result in a "valid" results we try to compute an approximation
                if (!r_interface_info->GetLocalSearchWasSuccessful())
                {
                    for (SizeType j=0; j<number_of_results; ++j)
                    {
                        r_interface_info->ProcessSearchResultForApproximation(
                            neighbor_results[j], neighbor_distances[j]);
                    }
                }
            }
        }
    }

    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/
    void InterfaceSearchStructureBase::ConductSearchIteration(const Kratos::Flags& rOptions,
                               const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo,
                               InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
    {
        PrepareSearchIteration(rOptions, rpInterfaceInfo, InterfaceObjectTypeOrigin);

        ConductLocalSearch();

        FinalizeSearchIteration();
    }

    bool InterfaceSearchStructureBase::AllNeighborsFound(const Communicator& rComm) const
    {
        int all_neighbors_found = 1; // set to "1" aka "true" by default in case
        // this partition doesn't have a part of the interface!

        for (const auto& local_sys : (*mpMapperLocalSystems))
        {
            if (!local_sys->HasInterfaceInfo())
            {
                all_neighbors_found = 0;
                break;
            }
        }

        // This is necessary bcs not all partitions would start a new search iteration!
        rComm.MinAll(all_neighbors_found);

        return all_neighbors_found > 0;
   }



}  // namespace Kratos.
