//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

// Project includes
#include "mpm_mpi_search.h"
#include "mpi_utilities.h"
#include "custom_utilities/mpm_search_element_utility.h"
#include "utilities/builtin_timer.h"

namespace Kratos {

void MPM_MPI_SEARCH::SearchElementMPI(ModelPart& rBackgroundGridModelPart, ModelPart& rMPMModelPart,
                                      const std::size_t MaxNumberOfResults, const double Tolerance)
{
    BuiltinTimer timer;
    MPMSearchElementUtility::ResetElementsAndNodes(rBackgroundGridModelPart);
    // Attention: Run search on submodel_parts, to ensure that elements/conditions are added to
    //            correct model_partp after mpi-send.
    // TODO: This should be bottom up style!! Here only one submodel_part hierachy of rMPMModelPart is considered
    //       However, in mpm_particle_generator also only one submodel_part hierachy level is considered.
    for( auto& r_submodel_part : rMPMModelPart.SubModelParts() ){
        std::vector<typename Element::Pointer> missing_elements;
        std::vector<typename Condition::Pointer> missing_conditions;

        MPM_MPI_Utilities::SetMPICommunicator(rMPMModelPart, r_submodel_part);

        MPMSearchElementUtility::NeighbourSearchElements(r_submodel_part, rBackgroundGridModelPart, missing_elements, Tolerance);
        MPMSearchElementUtility::NeighbourSearchConditions(r_submodel_part, rBackgroundGridModelPart, missing_conditions, Tolerance);

        std::vector<int> element_search_result_dummy(missing_elements.size());
        std::vector<int> condition_search_result_dummy(missing_conditions.size());
        if (missing_conditions.size() > 0 || missing_elements.size() > 0)
            BinBasedSearchElementsAndConditionsMPI(r_submodel_part,
                rBackgroundGridModelPart, missing_elements, element_search_result_dummy, missing_conditions,
                condition_search_result_dummy, MaxNumberOfResults, Tolerance);

        SearchElementsInOtherPartitions(r_submodel_part, rBackgroundGridModelPart, missing_elements,
            missing_conditions, MaxNumberOfResults, Tolerance);

    }
    MPM_MPI_Utilities::SynchronizeNodalDisplacementAtInterface(rMPMModelPart);
    rMPMModelPart.GetCommunicator().GetDataCommunicator().Barrier();
}

void MPM_MPI_SEARCH::SearchElementsInOtherPartitions(ModelPart& rMPMSubModelPart, ModelPart& rBackgroundGridModelPart,
                                                     std::vector<Element::Pointer>& rMissingElements,
                                                     std::vector<Condition::Pointer>& rMissingConditions,
                                                     const std::size_t MaxNumberOfResults, const double Tolerance)
{

    const unsigned size = rMPMSubModelPart.GetCommunicator().TotalProcesses();
    const unsigned rank = rMPMSubModelPart.GetCommunicator().MyPID();

    const unsigned int number_of_missing_elements = rMissingElements.size();
    const unsigned int number_of_missing_conditions = rMissingConditions.size();

    // Get local and gloabl number of elements/conditions
    // unsigned int local_number_of_elements = rMPMSubModelPart.NumberOfElements();
    // const unsigned int global_number_of_elements_begin = rBackgroundGridModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_number_of_elements);
    // unsigned int local_number_of_conditions = rMPMSubModelPart.NumberOfConditions();
    // const unsigned int global_number_of_conditions_begin = rBackgroundGridModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_number_of_conditions);

    // Construct element containers
    std::vector<ElementsArrayType> send_elements_container(size);
    ElementsArrayType send_elements;
    // Constrict condition containers
    std::vector<ConditionsArrayType> send_conditions_container(size);
    ConditionsArrayType send_conditions;

    // Fill element and condition containers
    send_elements.insert(rMissingElements.begin(), rMissingElements.end());
    send_conditions.insert(rMissingConditions.begin(), rMissingConditions.end());
    for( int i = 0; i < size; ++i){
        if( i != rank){
            send_elements_container[i] = send_elements;
            send_conditions_container[i] = send_conditions;
        }
    }

    // Construct recieving containers
    std::vector<ElementsArrayType> recv_elements_container(size);
    std::vector<ConditionsArrayType> recv_conditions_container(size);

    // Transfer elements and conditions
    MPM_MPI_Utilities::TransferElements(rMPMSubModelPart, send_elements_container, recv_elements_container);
    MPM_MPI_Utilities::TransferConditions(rMPMSubModelPart, send_conditions_container, recv_conditions_container);

    // Loop over all available ranks
    for( int i = 0; i < size; ++i){
        // Consttruct missing element/condition containers
        std::vector<Element::Pointer> missing_elements;
        std::vector<Condition::Pointer> missing_conditions;
        // Initialize search results containers
        // local search results: [0, 1, 0] -> [not found, found, not found]
        // global search results concatenates all local results: [ local_results_1, local_results_1, ..]
        std::vector<int> local_element_search_results ;
        std::vector<int> global_element_search_results;
        std::vector<int> local_condition_search_results;
        std::vector<int> global_condition_search_results;
        unsigned int number_sent_elements;
        unsigned int number_sent_conditions;
        if( rank == i){ // Current sender
            number_sent_elements = number_of_missing_elements;
            number_sent_conditions = number_of_missing_conditions;
        }
        else{ // Current recievers
            number_sent_elements = recv_elements_container[i].size();
            number_sent_conditions = recv_conditions_container[i].size();
        }
        if( number_sent_elements > 0 || number_sent_conditions > 0){
            // Elements
            local_element_search_results.resize(number_sent_elements);
            std::fill(local_element_search_results.begin(), local_element_search_results.end(), 0);
            global_element_search_results.resize(number_sent_elements*size);
            // Conditions
            local_condition_search_results.resize(number_sent_conditions);
            std::fill(local_condition_search_results.begin(), local_condition_search_results.end(), 0);
            global_condition_search_results.resize(number_sent_conditions*size);

            if( i != rank){ //If reciever perform search
                missing_elements.reserve(number_sent_elements);
                missing_conditions.reserve(number_sent_conditions);

                // Prepare element and conditions containers
                // TODO: Switch SetGeometryParent call inside BinBasedSearch to avoid this loop here!
                for( auto it = recv_elements_container[i].begin(); it != recv_elements_container[i].end(); ++it){
                    it->GetGeometry().SetGeometryParent( &rBackgroundGridModelPart.ElementsBegin()->GetGeometry());
                    missing_elements.push_back(*it.base());
                }

                for( auto it = recv_conditions_container[i].begin(); it != recv_conditions_container[i].end(); ++it){
                    it->GetGeometry() = rBackgroundGridModelPart.ElementsBegin()->GetGeometry();
                    missing_conditions.push_back(*it.base());
                }
                // Run search and fill local_search_results
                BuiltinTimer timer_2;
                BinBasedSearchElementsAndConditionsMPI(rMPMSubModelPart,
                            rBackgroundGridModelPart, missing_elements, local_element_search_results, missing_conditions,
                            local_condition_search_results, MaxNumberOfResults, Tolerance);

            }

            // Gather all local_element_search_results
            rMPMSubModelPart.GetCommunicator().GetDataCommunicator().AllGather(local_element_search_results, global_element_search_results);
            // Gather all local_condition_search_results
            rMPMSubModelPart.GetCommunicator().GetDataCommunicator().AllGather(local_condition_search_results, global_condition_search_results);


            // Make sure every element is only found ones!
            std::vector<int> global_element_search_results_combined(number_sent_elements);
            std::fill(global_element_search_results_combined.begin(), global_element_search_results_combined.end(), 0);
            for(int j = 0; j < size; ++j){
                for(int k = 0; k < number_sent_elements; ++k){
                    if( global_element_search_results_combined[k] >= 1 && global_element_search_results[k + j * number_sent_elements] >= 1 ){
                        global_element_search_results[k + j * number_sent_elements] = 0;
                    }
                    global_element_search_results_combined[k] += global_element_search_results[k + j * number_sent_elements];
                }
            }
            // Add found elements to current ModelPart
            int j = rank*number_sent_elements;
            for( auto it = recv_elements_container[i].begin(); it != recv_elements_container[i].end(); ++it){
                if( global_element_search_results[j] == 1){
                    it->Set(ACTIVE);
                    rMPMSubModelPart.AddElement(*it.base());
                    rMPMSubModelPart.GetCommunicator().LocalMesh().AddElement(*it.base());
                }
                j++;
            }

            // Make sure every condition is only found ones!
            // global_search_results_combined -> [local_search_results_1 + local_search_results_2 + ...]
            // global_search_results_combined must in the end only contain 1's
            std::vector<int> global_condition_search_results_combined(number_sent_conditions);
            std::fill(global_condition_search_results_combined.begin(), global_condition_search_results_combined.end(), 0);
            for(int j = 0; j < size; ++j){
                for(int k = 0; k < number_sent_conditions; ++k){
                    if( global_condition_search_results_combined[k] >= 1 && global_condition_search_results[k + j * number_sent_conditions] >= 1 ){
                        global_condition_search_results[k + j * number_sent_conditions] = 0;
                    }
                    global_condition_search_results_combined[k] += global_condition_search_results[k + j * number_sent_conditions];
                }
            }
            // Add found condition to current ModelPart
            j = rank*number_sent_conditions;
            for( auto it = recv_conditions_container[i].begin(); it != recv_conditions_container[i].end(); ++it){
                if( global_condition_search_results[j] == 1){
                    it->Set(ACTIVE);
                    rMPMSubModelPart.AddCondition(*it.base());
                    rMPMSubModelPart.GetCommunicator().LocalMesh().AddCondition(*it.base());
                }
                j++;
            }
        }
    }

    // TODO: If everything works, change the comparison "!=" to "<". This will than indicate that not all elements/conditions were found
    //       and hence left the global domain (global backgroundgrid)
    // local_number_of_conditions = rMPMSubModelPart.NumberOfConditions();
    // const unsigned int global_number_of_conditions_end = rBackgroundGridModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_number_of_conditions);
    // KRATOS_ERROR_IF(global_number_of_conditions_begin != global_number_of_conditions_end) << "MPM_MPI_SEARCH::SearchElements: "
    //     << "Global number of conditions before (" << global_number_of_conditions_begin << ") and after (" << global_number_of_conditions_end
    //     << ") mpi-send/recv operations are not identical!" << std::endl;

    // local_number_of_elements = rMPMSubModelPart.NumberOfElements();
    // const unsigned int global_number_of_elements_end = rBackgroundGridModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_number_of_elements);
    // KRATOS_ERROR_IF(global_number_of_elements_begin != global_number_of_elements_end) << "MPM_MPI_SEARCH::SearchElements: "
    //     << "Global number of elements before (" << global_number_of_elements_begin << ") and after (" << global_number_of_elements_end
    //     << ") mpi-send/recv operations are not identical!" << std::endl;


}

void MPM_MPI_SEARCH::BinBasedSearchElementsAndConditionsMPI(ModelPart& rMPMModelPart,
        ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        std::vector<int>& element_search_results,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        std::vector<int>& condition_search_results,
        const std::size_t MaxNumberOfResults, const double Tolerance)
{
    constexpr int TDimension = 2;
    const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();
    bool is_pqmpm = (r_process_info.Has(IS_PQMPM))
        ? r_process_info.GetValue(IS_PQMPM) : false;

    // Search background grid and make element active
    Vector N;
    const int max_result = 1000;

    // Temporary copy of missing elements and conditions
    std::vector<Element::Pointer> missing_elements;
    missing_elements.reserve(rMissingElements.size());
    std::vector<Condition::Pointer> missing_conditions;
    missing_conditions.reserve(rMissingConditions.size());


    #pragma omp parallel
    {

        BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
        SearchStructure.UpdateSearchDatabase();

        typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(max_result);

        // Element search and assign background grid
        #pragma omp for
        for (int i = 0; i < static_cast<int>(rMissingElements.size()); ++i) {
            auto element_itr = *(rMissingElements.begin() + i);
            std::vector<array_1d<double, 3>> xg;
            element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());
            typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
            Element::Pointer pelem;

            // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
            bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

            if (is_found == true) {
                if (MPMSearchElementUtility::IsFixExplicitAndOnElementEdge(N, r_process_info) && !is_pqmpm) {
                    KRATOS_ERROR << "MPM_MPI_SEARCH: IsFixExplicitAndOnElementEdge is not yet implemented/tested for MPI" << std::endl;
                    // MP is exactly on the edge. Now we give it a little 'nudge'
                    // array_1d<double, 3> xg_nudged = array_1d<double, 3>(xg[0]);
                    // std::vector<array_1d<double, 3>> mp_vel;
                    // element_itr->CalculateOnIntegrationPoints(MP_VELOCITY, mp_vel, rMPMModelPart.GetProcessInfo());
                    // xg_nudged += r_process_info[DELTA_TIME] / 1000.0 * mp_vel[0];
                    // if (SearchStructure.FindPointOnMesh(xg_nudged, N, pelem, result_begin, MaxNumberOfResults, Tolerance)) {
                    //     element_itr->SetValuesOnIntegrationPoints(MP_COORD, { xg_nudged }, rMPMModelPart.GetProcessInfo());
                    //     KRATOS_INFO("MPMSearchElementUtility") << "WARNING: To prevent spurious explicit stresses, Material Point "
                    //         << element_itr->Id() << " was nudged." << std::endl;
                    // } else {
                    //     is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                    //     KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Material Point " << element_itr->Id()
                    //         << " lies exactly on an element edge and may give spurious results." << std::endl;
                    //}
                }
                pelem->Set(ACTIVE);

                const bool is_pqmpm = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_PQMPM))
                    ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_PQMPM) : false;
                if (is_pqmpm)
                {
                    KRATOS_ERROR << "MPM_MPI_SEARCH: PQMPM is not yet implemented/tested for MPI" << std::endl;
                    // Updates the quadrature point geometry.
                    // UpdatePartitionedQuadraturePoint(rBackgroundGridModelPart, xg[0],
                    //     *element_itr, pelem->pGetGeometry(), Tolerance);
                }
                else
                {
                    auto p_quadrature_point_geometry = element_itr->pGetGeometry();
                    array_1d<double, 3> local_coordinates;

                    //Shouldn't this be also in the normal search!!! BUG???
                    pelem->GetGeometry().PointLocalCoordinates(local_coordinates, xg[0]);

                    CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                        p_quadrature_point_geometry, local_coordinates,
                        p_quadrature_point_geometry->IntegrationPoints()[0].Weight(), pelem->GetGeometry());

                }
                element_itr->Set(ACTIVE);
                auto& r_geometry = element_itr->GetGeometry();
                for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j){
                    r_geometry[j].Set(ACTIVE);
                }
                element_search_results[i] = 1;
            }
            else {

                //TODO: Add that it still clears conditions if only one mpi node is available
                #pragma omp critical
                missing_elements.push_back(&*element_itr);
                element_search_results[i] = 0;
            }
        }

        // Condition search and assign background grid
        #pragma omp for
        for (int i = 0; i < static_cast<int>(rMissingConditions.size()); ++i) {
            auto condition_itr = *(rMissingConditions.begin() + i);
            std::vector<array_1d<double, 3>> xg;
            condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

            if (xg.size() > 0) {
                // Only search for particle based BCs!
                // Grid BCs are still applied on MP_model_part but we don't want to search for them.
                typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                Element::Pointer pelem;

                // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                if (is_found == true) {
                    condition_itr->GetGeometry() = pelem->GetGeometry();
                    auto& r_geometry = condition_itr->GetGeometry();

                    for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                        r_geometry[j].Set(ACTIVE);

                    condition_search_results[i] = 1;
                } else {
                    //TODO: Add that it still clears conditions if only one mpi node is available
                    #pragma omp critical
                    missing_conditions.push_back(&*condition_itr);

                    condition_search_results[i] = 0;
                }
            }
        }
        }
        rMissingElements.clear();
        rMissingElements.insert(rMissingElements.begin(),missing_elements.begin(), missing_elements.end());

        rMissingConditions.clear();
        rMissingConditions.insert(rMissingConditions.begin(),missing_conditions.begin(), missing_conditions.end());
    }

} // namespace Kratos