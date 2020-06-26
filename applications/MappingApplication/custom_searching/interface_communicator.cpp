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

namespace Kratos {

namespace {

bool SearchNotSuccessful(const InterfaceCommunicator::MapperInterfaceInfoPointerType& rpInterfaceInfo)
{
    return !(rpInterfaceInfo->GetLocalSearchWasSuccessful());
}

}

typedef std::size_t IndexType;
typedef std::size_t SizeType;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
void InterfaceCommunicator::ExchangeInterfaceData(const Communicator& rComm,
                            const Kratos::Flags& rOptions,
                            const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo)
{

    mSearchRadius = mSearchSettings["search_radius"].GetDouble();
    const int max_search_iterations = mSearchSettings["search_iterations"].GetInt();

    const double increase_factor = 4.0;
    int num_iteration = 1;
    InitializeSearch(rOptions, rpInterfaceInfo);

    // First Iteration is done outside the search loop bcs it has
    // to be done in any case
    // one search iteration should be enough in most cases (if the search
    // radius was either computed or specified properly)
    // only if some points did not find a neighbor or dont have a valid
    // projection, more search iterations are necessary
    ConductSearchIteration(rOptions, rpInterfaceInfo);

    while (++num_iteration <= max_search_iterations && !AllNeighborsFound(rComm)) {
        mSearchRadius *= increase_factor;

        KRATOS_WARNING_IF("Mapper", mEchoLevel >= 1 && rComm.MyPID() == 0)
            << "search radius was increased, another search iteration is conducted\n"
            << "search iteration " << num_iteration << " / "<< max_search_iterations << " | "
            << "search radius " << mSearchRadius << std::endl;

        ConductSearchIteration(rOptions, rpInterfaceInfo);
    }

    FinalizeSearch();

    rComm.GetDataCommunicator().Barrier();
}

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/
void InterfaceCommunicator::InitializeSearch(const Kratos::Flags& rOptions,
                                             const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    if (mpInterfaceObjectsOrigin == nullptr || rOptions.Is(MapperFlags::REMESHED)) {
        CreateInterfaceObjectsOrigin(rpRefInterfaceInfo);
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
    for (auto& r_interface_infos_rank : mMapperInterfaceInfosContainer) {
        r_interface_infos_rank.clear();
    }
}

void InterfaceCommunicator::InitializeSearchIteration(const Kratos::Flags& rOptions,
                                                      const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    // Creating the MapperInterfaceInfos
    mMapperInterfaceInfosContainer[0].clear();

    auto& r_mapper_interface_infos = mMapperInterfaceInfosContainer[0];
    r_mapper_interface_infos.reserve(mrMapperLocalSystems.size());

    IndexType local_sys_idx = 0;
    for (const auto& r_local_sys : mrMapperLocalSystems) {
        if (!r_local_sys->HasInterfaceInfo()) { // Only the local_systems that have not received an InterfaceInfo create a new one
            const auto& r_coords = r_local_sys->Coordinates();
            r_mapper_interface_infos.push_back(rpRefInterfaceInfo->Create(r_coords, local_sys_idx, 0)); // dummy-rank of 0
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
    for (IndexType i_rank=0; i_rank<mMapperInterfaceInfosContainer.size(); ++i_rank) {
        auto& r_interface_infos_rank = mMapperInterfaceInfosContainer[i_rank];

        auto new_end = std::remove_if(
            r_interface_infos_rank.begin(),
            r_interface_infos_rank.end(),
            SearchNotSuccessful); // cannot use lambda here bcs it is not supported by some older compilers

        r_interface_infos_rank.erase(new_end, r_interface_infos_rank.end());
    }
}

void InterfaceCommunicator::AssignInterfaceInfos()
{
    // NOTE: mMapperInterfaceInfosContainer must contain only the ones that are a successfuls search!
    const SizeType comm_size = mMapperInterfaceInfosContainer.size();

    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        for (const auto& rp_interface_info : mMapperInterfaceInfosContainer[i_rank]) {
            mrMapperLocalSystems[rp_interface_info->GetLocalSystemIndex()]
                ->AddInterfaceInfo(rp_interface_info);
        }
    }
}

/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/
void InterfaceCommunicator::CreateInterfaceObjectsOrigin(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    mpInterfaceObjectsOrigin = Kratos::make_unique<InterfaceObjectContainerType>();

    const auto interface_obj_type = rpRefInterfaceInfo->GetInterfaceObjectType();

    if (interface_obj_type == InterfaceObject::ConstructionType::Node_Coords) {
        const SizeType num_nodes = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes();
        const auto nodes_begin = mrModelPartOrigin.GetCommunicator().LocalMesh().Nodes().ptr_begin();

        mpInterfaceObjectsOrigin->resize(num_nodes);

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(num_nodes); ++i) {
            auto it_node = nodes_begin + i;
            (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceNode>((*it_node).get());
        }
    }

    else if (interface_obj_type == InterfaceObject::ConstructionType::Geometry_Center) {
        Communicator& r_comm = mrModelPartOrigin.GetCommunicator();
        const SizeType num_elements = r_comm.LocalMesh().NumberOfElements();
        const SizeType num_conditions = r_comm.LocalMesh().NumberOfConditions();

        const auto elements_begin = r_comm.LocalMesh().Elements().ptr_begin();
        const auto conditions_begin = r_comm.LocalMesh().Conditions().ptr_begin();

        int num_elements_global = r_comm.GetDataCommunicator().SumAll(static_cast<int>(num_elements));
        int num_conditions_global = r_comm.GetDataCommunicator().SumAll(static_cast<int>(num_conditions));

        KRATOS_ERROR_IF(num_elements_global > 0 && num_conditions_global > 0)
            << "Both Elements and Conditions are present which is not allowed!\n"
            << "Name of ModelPart: " << mrModelPartOrigin.Name()
            << "\nNumber of Elements: " << num_elements
            << "; Number of Condition: " << num_conditions << std::endl;

        KRATOS_ERROR_IF(num_elements_global+num_conditions_global == 0)
            << "No Elements and Conditions are present which is not allowed!\n"
            << "Name of ModelPart: " << mrModelPartOrigin.Name() << std::endl;

        mpInterfaceObjectsOrigin->resize(num_elements+num_conditions); // one of them has to be zero!!!

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(num_elements); ++i) {
            auto it_elem = elements_begin + i;
            (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceGeometryObject>((*it_elem)->pGetGeometry().get());
        }
        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(num_conditions); ++i) {
            auto it_cond = conditions_begin + i;
            (*mpInterfaceObjectsOrigin)[i] = Kratos::make_unique<InterfaceGeometryObject>((*it_cond)->pGetGeometry().get());
        }
    }
    else {
        KRATOS_ERROR << "Type of interface object construction not implemented" << std::endl;
    }

    // Making sure that the data-structure was correctly initialized
    int num_interface_objects = mpInterfaceObjectsOrigin->size(); // int bcs of MPI
    num_interface_objects = mrModelPartOrigin.GetCommunicator().GetDataCommunicator().SumAll(num_interface_objects);

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

    KRATOS_ERROR_IF(mSearchRadius < 0.0) << "Search-Radius has to be larger than 0.0!"
        << std::endl;

    if (num_interface_obj_bin > 0) { // this partition has a bin structure
        InterfaceObjectConfigure::ResultContainerType neighbor_results(num_interface_obj_bin);
        std::vector<double> neighbor_distances(num_interface_obj_bin);
        auto interface_obj(Kratos::make_shared<InterfaceObject>(array_1d<double, 3>(0.0)));

        for (auto& r_interface_infos_rank : mMapperInterfaceInfosContainer) { // loop the ranks
            // #pragma omp parallel for // TODO this requires to make some things thread-local!
            // it makes more sense to omp this loop even though it is not the outermost one ...
            for (IndexType i=0; i<r_interface_infos_rank.size(); ++i) {
                auto& r_interface_info = r_interface_infos_rank[i];

                interface_obj->Coordinates() = r_interface_info->Coordinates();
                double search_radius = mSearchRadius; // reset search radius // TODO check this

                // reset the containers
                auto results_itr = neighbor_results.begin();
                auto distance_itr = neighbor_distances.begin();

                const SizeType number_of_results = mpLocalBinStructure->SearchObjectsInRadius(
                    interface_obj, search_radius, results_itr,
                    distance_itr, num_interface_obj_bin);

                for (IndexType j=0; j<number_of_results; ++j) {
                    r_interface_info->ProcessSearchResult(*(neighbor_results[j]), neighbor_distances[j]);
                }

                // If the search did not result in a "valid" result (e.g. the projection fails)
                // we try to compute an approximation
                if (!r_interface_info->GetLocalSearchWasSuccessful()) {
                    for (IndexType j=0; j<number_of_results; ++j) {
                        r_interface_info->ProcessSearchResultForApproximation(
                            *(neighbor_results[j]), neighbor_distances[j]);
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

    for (const auto& local_sys : mrMapperLocalSystems) {
        if (!local_sys->HasInterfaceInfo()) {
            all_neighbors_found = 0;
            break;
        }
    }

    // This is necessary bcs not all partitions would start a new search iteration!
    all_neighbors_found = rComm.GetDataCommunicator().MinAll(all_neighbors_found);

    return all_neighbors_found > 0;
}

}  // namespace Kratos.
