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
#include "interface_communicator.h"

namespace Kratos
{

typedef std::size_t IndexType;
typedef std::size_t SizeType;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
void InterfaceCommunicator::ExchangeInterfaceData(const Communicator& rComm,
                            const Kratos::Flags& rOptions,
                            const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo,
                            InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
{

    mSearchRadius = mSearchSettings["search_radius"].GetDouble();
    const IndexType max_search_iterations = mSearchSettings["search_iterations"].GetDouble();

    const double increase_factor = 4.0;
    IndexType num_iteration = 1;
    bool last_iteration = (max_search_iterations == 1) ? true : false; // true in case only one search iteration is conducted // TODO needed???

    InitializeSearch(rOptions, rpInterfaceInfo, InterfaceObjectTypeOrigin);

    // First Iteration is done outside the search loop bcs it has
    // to be done in any case
    // one search iteration should be enough in most cases (if the search
    // radius was either computed or specified properly)
    // only if some points did not find a neighbor or dont have a valid
    // projection, more search iterations are necessary
    ConductSearchIteration(rOptions, rpInterfaceInfo);

    while (++num_iteration <= max_search_iterations && !AllNeighborsFound(rComm)) {
        mSearchRadius *= increase_factor;

        if (num_iteration == max_search_iterations) last_iteration = true; // TODO test if this works...

        KRATOS_WARNING_IF("Mapper", mEchoLevel >= 2 && rComm.MyPID() == 0)
            << "search radius was increased, another search iteration is conducted | "
            << "search iteration " << num_iteration << " / "<< max_search_iterations << " | "
            << "search radius " << mSearchRadius << std::endl;

        ConductSearchIteration(rOptions, rpInterfaceInfo);
    }

    FinalizeSearch();

    rComm.Barrier();
}

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/
void InterfaceCommunicator::InitializeSearch(const Kratos::Flags& rOptions,
                                    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                    InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
{
    if (mpInterfaceObjectsOrigin == nullptr || rOptions.Is(MapperFlags::REMESHED)) {
        CreateInterfaceObjectsOrigin(InterfaceObjectTypeOrigin);
    }
    else {
        UpdateInterfaceObjectsOrigin();
    }

    if (mpLocalBinStructure == nullptr || !rOptions.Is(MapperFlags::DESTINATION_ONLY)) {
        InitializeBinsSearchStructure(); // This cannot be updated, has to be recreated
    }
}

void InterfaceCommunicator::FinalizeSearch()
{
    for (auto& r_interface_infos_rank : (*mpMapperInterfaceInfosContainer)) {
        r_interface_infos_rank.clear();
    }
}

void InterfaceCommunicator::InitializeSearchIteration(const Kratos::Flags& rOptions,
                                                const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    // Creating the MapperInterfaceInfos
    (*mpMapperInterfaceInfosContainer)[0].clear();

    auto& r_mapper_interface_infos = (*mpMapperInterfaceInfosContainer)[0];
    r_mapper_interface_infos.reserve(mpMapperLocalSystems->size());

    IndexType local_sys_idx = 0;
    for (const auto& r_local_sys : (*mpMapperLocalSystems)) {
        if (!r_local_sys->HasInterfaceInfo()) { // Only the local_systems that have not received an InterfaceInfo create a new one
            const auto& r_coords = r_local_sys->Coordinates();
            r_mapper_interface_infos.push_back(rpRefInterfaceInfo->Create(r_coords, local_sys_idx));
        }
        ++local_sys_idx;
    }
}

void InterfaceCommunicator::FinalizeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo)
{
    FilterInterfaceInfosSuccessfulSearch();
    AssignInterfaceInfos();
}

void InterfaceCommunicator::FilterInterfaceInfosSuccessfulSearch()
{
    // Erasing all the MapperLocalSystems that don't have a successful search
    // using the erase–remove idiom
    for (IndexType i_rank=0; i_rank<mpMapperInterfaceInfosContainer->size(); ++i_rank) {
        auto& r_interface_infos_rank = (*mpMapperInterfaceInfosContainer)[i_rank];

        auto new_end = std::remove_if(
            r_interface_infos_rank.begin(),
            r_interface_infos_rank.end(),
            [](const auto& rp_interface_info)
            { return !((*rp_interface_info).GetLocalSearchWasSuccessful()); });

        r_interface_infos_rank.erase(new_end, r_interface_infos_rank.end());
    }
}

void InterfaceCommunicator::AssignInterfaceInfos()
{
    // NOTE: mpMapperInterfaceInfosContainer must contain only the ones that are a successfuls search!
    const SizeType comm_size = mpMapperInterfaceInfosContainer->size();

    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        for (const auto& rp_interface_info : (*mpMapperInterfaceInfosContainer)[i_rank]) {
            (*mpMapperLocalSystems)[rp_interface_info->GetLocalSystemIndex()]
                ->AddInterfaceInfo(rp_interface_info);
        }
    }
}

/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/
void InterfaceCommunicator::CreateInterfaceObjectsOrigin(InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
{
    mpInterfaceObjectsOrigin = Kratos::make_unique<InterfaceObjectContainerType>();

    if (InterfaceObjectTypeOrigin == InterfaceObject::Node_Coords) {
        const SizeType num_nodes = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes();
        const auto nodes_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Nodes().ptr_begin();

        mpInterfaceObjectsOrigin->resize(num_nodes);

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(num_nodes); ++i) {
            auto it_node = nodes_begin + i;
            (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceNode>(*(it_node));
        }
    }

    else if (InterfaceObjectTypeOrigin == InterfaceObject::Geometry_Center) {
        const SizeType num_elements = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfElements();
        const SizeType num_conditions = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfConditions();

        const auto elements_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Elements().ptr_begin();
        const auto conditions_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Conditions().ptr_begin();

        // TODO sum those accross ranks!
        KRATOS_ERROR_IF( (num_elements > 0 && num_conditions != 0) ||
                         (num_elements != 0 && num_conditions > 0) )
            << "Both Elements and Conditions are present which is not allowed!\n"
            << "Name of ModelPart: " << mrModelPartOrigin.Name()
            << "\nNumber of Elements: " << num_elements
            << "; Number of Condition: " << num_conditions << std::endl;

        mpInterfaceObjectsOrigin->resize(num_elements+num_conditions); // one of them has to be zero!!!

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(num_elements); ++i) {
            auto it_elem = elements_begin + i;
            (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceGeometryObject>((*it_elem)->pGetGeometry());
        }
        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(num_conditions); ++i) {
            auto it_cond = conditions_begin + i;
            (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceGeometryObject>((*it_cond)->pGetGeometry());
        }
    }
    else {
        KRATOS_ERROR << "Type of interface object construction not implemented" << std::endl;
    }

    // Making sure that the data-structure was correctly initialized
    int num_interface_objects = mpInterfaceObjectsOrigin->size(); // int bcs of MPI
    mrModelPartOrigin.GetCommunicator().SumAll(num_interface_objects);

    KRATOS_ERROR_IF_NOT(num_interface_objects > 0)
        << "No interface objects were created in Origin-ModelPart \""
        << mrModelPartOrigin.Name() << "\"!" << std::endl;
}

void InterfaceCommunicator::UpdateInterfaceObjectsOrigin()
{
    const auto begin = mpInterfaceObjectsOrigin->begin();

    #pragma omp parallel for
    for (int i = 0; i< static_cast<int>(mpInterfaceObjectsOrigin->size()); ++i) {
        (*(begin + i))->UpdateCoordinates();
    }
}

void InterfaceCommunicator::InitializeBinsSearchStructure()
{
    if (mpInterfaceObjectsOrigin->size() > 0) { // only construct the bins if the partition has a part of the interface
        mpLocalBinStructure = Kratos::make_unique<BinsObjectDynamic<InterfaceObjectConfigure>>(
            mpInterfaceObjectsOrigin->begin(), mpInterfaceObjectsOrigin->end());
    }
}

void InterfaceCommunicator::ConductLocalSearch()
{
    SizeType num_interface_obj_bin = mpInterfaceObjectsOrigin->size();

    if (num_interface_obj_bin > 0) { // this partition has a bin structure
        InterfaceObjectConfigure::ResultContainerType neighbor_results(num_interface_obj_bin);
        std::vector<double> neighbor_distances(num_interface_obj_bin);
        auto interface_obj = Kratos::make_shared<InterfaceObject>(array_1d<double, 3>(0.0));

        for (auto& r_interface_infos_rank : (*mpMapperInterfaceInfosContainer)) { // loop the ranks
            // #pragma omp parallel for // TODO this requires to make some things thread-local!
            // it makes more sense to omp this loop even though it is not the outermost one ...
            for (IndexType i=0; i<r_interface_infos_rank.size(); ++i) {
                auto& r_interface_info = r_interface_infos_rank[i];

                interface_obj->UpdateCoordinates(r_interface_info->Coordinates());
                double search_radius = mSearchRadius; // reset search radius // TODO check this

                // reset the containers
                auto results_itr = neighbor_results.begin();
                auto distance_itr = neighbor_distances.begin();

                const SizeType number_of_results = mpLocalBinStructure->SearchObjectsInRadius(
                    interface_obj, search_radius, results_itr,
                    distance_itr, num_interface_obj_bin);

                for (IndexType j=0; j<number_of_results; ++j) {
                    r_interface_info->ProcessSearchResult(neighbor_results[j], neighbor_distances[j]);
                }

                // If the search did not result in a "valid" result (e.g. the projection fails)
                // we try to compute an approximation
                if (!r_interface_info->GetLocalSearchWasSuccessful()) {
                    for (IndexType j=0; j<number_of_results; ++j) {
                        r_interface_info->ProcessSearchResultForApproximation(
                            neighbor_results[j], neighbor_distances[j]);
                    }
                }
            }
        }
    }
}

void InterfaceCommunicator::ConductSearchIteration(const Kratos::Flags& rOptions,
                            const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo)
{
    InitializeSearchIteration(rOptions, rpInterfaceInfo);

    ConductLocalSearch();

    FinalizeSearchIteration(rpInterfaceInfo);
}

bool InterfaceCommunicator::AllNeighborsFound(const Communicator& rComm) const
{
    int all_neighbors_found = 1; // set to "1" aka "true" by default in case
    // this partition doesn't have a part of the interface!

    for (const auto& local_sys : (*mpMapperLocalSystems)) {
        if (!local_sys->HasInterfaceInfo()) {
            all_neighbors_found = 0;
            break;
        }
    }

    // This is necessary bcs not all partitions would start a new search iteration!
    rComm.MinAll(all_neighbors_found);

    return all_neighbors_found > 0;
}

}  // namespace Kratos.
