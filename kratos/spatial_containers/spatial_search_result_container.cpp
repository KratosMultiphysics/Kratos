//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <numeric>
#include <functional>

// External includes

// Project includes
#include "includes/data_communicator.h"
#include "includes/node.h"
#include "includes/geometrical_object.h"
#include "spatial_containers/spatial_search_result_container.h"

namespace Kratos
{

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::SpatialSearchResultContainer()
{
    /// Some error in case asynchronous
    if constexpr (TSpatialSearchCommunication == SpatialSearchCommunication::ASYNCHRONOUS) {
        KRATOS_ERROR << "Asynchronous communication is not yet implemented for SpatialSearchResultContainer" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::AddResult(SpatialSearchResultType& rResult)
{
    // Check if the object has been found
    if (rResult.IsObjectFound()) {
        // Adding to the local results
        mLocalResults.emplace_back(rResult);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::AddResult(
    TObjectType* pResult,
    const int Rank
    )
{
    // Check if the object has been found (not nullptr)
    if (pResult != nullptr) {
        // Adding to the local results
        mLocalResults.emplace_back(pResult, Rank);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::AddResult(
    TObjectType* pResult,
    const double Distance,
    const int Rank
    )
{
    // Check if the object has been found (not nullptr)
    if (pResult != nullptr) {
        // Adding to the local results
        mLocalResults.emplace_back(pResult, Rank);
        (mLocalResults.back()).SetDistance(Distance);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::Clear()
{
    // Clear pointer
    mpGlobalPointerCommunicator = nullptr;

    // Clear local results
    mLocalResults.clear();

    // Clear global results
    mGlobalResults.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GenerateGlobalPointerCommunicator(const DataCommunicator& rDataCommunicator)
{
    // Generate the communicator
    mpGlobalPointerCommunicator = Kratos::make_shared<GlobalPointerCommunicatorType>(rDataCommunicator, mGlobalResults.ptr_begin(), mGlobalResults.ptr_end());
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::SynchronizeAll(const DataCommunicator& rDataCommunicator)
{
    // TODO: Try to avoid use simple calls of Send/Recv and try to use more efficient methods

    // Synchronize local results to global results
    if(rDataCommunicator.IsDistributed()) { // MPI code
        // MPI info
        const int world_size = rDataCommunicator.Size();
        const int rank = rDataCommunicator.Rank();
        const int root_rank = 0;

        // 1. All ranks determine their local result size.
        const int local_result_size = mLocalResults.size();

        // 2. Gather the result counts from every rank onto all ranks.
        // This replaces the separate SumAll and AllGather calls with a single AllGather.
        std::vector<int> result_counts(world_size);
        rDataCommunicator.AllGather({local_result_size}, result_counts);

        // Calculate the total size from the gathered counts.
        const int global_result_size = std::accumulate(result_counts.begin(), result_counts.end(), 0);
        mGlobalResults.reserve(global_result_size);

        // 3. The logic now splits between the root rank (which gathers) and worker ranks (which send).
        std::vector<GlobalPointerResultType> global_gp;
        if (rank == root_rank) {
            // --- ROOT PROCESS (GATHER & BROADCAST) ---
            global_gp.reserve(global_result_size);

            // Add its own local results first.
            for (auto& r_value : mLocalResults) {
                global_gp.emplace_back(&r_value, rank);
            }

            // Receive results from all other ranks that have data.
            for (int i_rank = 0; i_rank < world_size; ++i_rank) {
                if (i_rank != root_rank && result_counts[i_rank] > 0) {
                    std::vector<GlobalPointerResultType> received_gps;
                    rDataCommunicator.Recv(received_gps, i_rank);
                    // Move received elements for efficiency.
                    global_gp.insert(global_gp.end(), std::make_move_iterator(received_gps.begin()), std::make_move_iterator(received_gps.end()));
                }
            }
        } else {
            // --- WORKER PROCESSES (SEND & RECEIVE) ---

            // Send local results to the root rank, if any exist.
            if (local_result_size > 0) {
                std::vector<GlobalPointerResultType> local_gp;
                local_gp.reserve(local_result_size);
                for (auto& r_value : mLocalResults) {
                    local_gp.emplace_back(&r_value, rank);
                }
                rDataCommunicator.Send(local_gp, root_rank);
            }
        }

        // Receive the final, complete vector from the root/Broadcast the complete results to all worker ranks.
        rDataCommunicator.Broadcast(global_gp, root_rank);

        // Finally, populate the member variable.
        mGlobalResults.clear();
        for (auto& r_gp : global_gp) {
            mGlobalResults.push_back(std::move(r_gp));
        }
    } else { // Serial code
        // Fill global vector
        for (auto& r_value : mLocalResults) {
            mGlobalResults.push_back(&r_value);
        }
    }

    // Generate the communicator
    GenerateGlobalPointerCommunicator(rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::RemoveResultsFromRanksList(
    const std::vector<int>& rRanks,
    const DataCommunicator& rDataCommunicator
    )
{
    // Get current rank
    const int rank = rDataCommunicator.Rank();

    // Remove current rank if required
    auto it_find = std::find(rRanks.begin(), rRanks.end(), rank);
    if (it_find != rRanks.end()) {
        mLocalResults.clear();
    }

    // Prepare a map of current indexes for quick lookup
    std::vector<int> ranks;
    GetResultRank(ranks, rDataCommunicator);
    std::vector<IndexType> indexes_to_remove;
    indexes_to_remove.reserve(ranks.size());
    for (IndexType i = 0; i < ranks.size(); ++i) {
        auto it_find = std::find(rRanks.begin(), rRanks.end(), ranks[i]);
        if (it_find != rRanks.end()) {
            indexes_to_remove.push_back(i);
        }
    }
    indexes_to_remove.shrink_to_fit();

    // Determine and remove indexes from global results
    std::sort(indexes_to_remove.begin(), indexes_to_remove.end(), std::greater<IndexType>());
    for (const auto index : indexes_to_remove) {
        mGlobalResults.erase(mGlobalResults.begin() + index);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::RemoveResultsFromIndexesList(const std::vector<IndexType>& rIndexes)
{
    // Prepare for local remove
    std::unordered_map<IndexType, IndexType> local_indexes_map;
    const auto it_begin_local = mLocalResults.begin();
    for (IndexType i = 0; i < mLocalResults.size(); ++i) {
        local_indexes_map.insert({(it_begin_local + i)->Get()->Id(), i});
    }

    // Remove from local results (indexes may be repeated if for example a node is in more than one partition, using therefore a set)
    std::unordered_set<IndexType> local_indexes_to_remove_set;
    for (const auto index : rIndexes) {
        const auto it_find = local_indexes_map.find(index);
        if (it_find != local_indexes_map.end()) {
            local_indexes_to_remove_set.insert(it_find->second);
        }
    }
    std::vector<IndexType> local_indexes_to_remove(local_indexes_to_remove_set.begin(), local_indexes_to_remove_set.end());
    std::sort(local_indexes_to_remove.begin(), local_indexes_to_remove.end(), std::greater<IndexType>());
    for (const auto index : local_indexes_to_remove) {
        mLocalResults.erase(mLocalResults.begin() + index);
    }

    // Prepare a map of current indexes for quick lookup
    std::vector<IndexType> current_indexes;
    this->GetResultIndices(current_indexes);
    std::unordered_map<IndexType, IndexType> current_indexes_map;
    for (IndexType i = 0; i < current_indexes.size(); ++i) {
        current_indexes_map[current_indexes[i]] = i;
    }

    // Determine and remove indexes from global results
    std::vector<IndexType> indexes_to_remove;
    indexes_to_remove.reserve(rIndexes.size());
    for (auto i_remove : rIndexes) {
        indexes_to_remove.push_back(current_indexes_map[i_remove]);
    }
    std::sort(indexes_to_remove.begin(), indexes_to_remove.end(), std::greater<IndexType>());
    for (const auto index : indexes_to_remove) {
        mGlobalResults.erase(mGlobalResults.begin() + index);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::string SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::Info() const
{
    std::stringstream buffer;
    buffer << "SpatialSearchResultContainer" ;
    return buffer.str();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "SpatialSearchResultContainer" << "\n";
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::PrintData(std::ostream& rOStream) const
{
    rOStream << "SpatialSearchResultContainer data summary: " << "\n";
    rOStream << "\tNumber of local results: " << mLocalResults.size() << "\n";
    rOStream << "\tNumber of global results: " << mGlobalResults.size() << "\n";
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultRank(
    std::vector<int>& rResults,
    const DataCommunicator& rDataCommunicator
    )
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    if (rResults.size() != number_of_gp) {
        rResults.resize(number_of_gp);
    }

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> int {
        return rGP->Get().GetRank();
    });

    // Get the rank
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        rResults[i] = rDataCommunicator.MaxAll(proxy.Get(r_gp));
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultIndices(std::vector<IndexType>& rResults)
{
    // Define the indices vector
    const std::size_t number_of_gp = mGlobalResults.size();
    if (rResults.size() != number_of_gp) {
        rResults.resize(number_of_gp);
    }

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> std::size_t {
        auto p_object = rGP->Get();
        return p_object->Id();
    });

    // Get the indices
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        rResults[i] = proxy.Get(r_gp);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::save(Serializer& rSerializer) const
{
    rSerializer.save("LocalResults", mLocalResults);
    rSerializer.save("GlobalResults", mGlobalResults);
    //rSerializer.save("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::load(Serializer& rSerializer)
{
    rSerializer.load("LocalResults", mLocalResults);
    rSerializer.load("GlobalResults", mGlobalResults);
    //rSerializer.load("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
template class SpatialSearchResultContainer<Node, SpatialSearchCommunication::SYNCHRONOUS>;
template class SpatialSearchResultContainer<GeometricalObject, SpatialSearchCommunication::SYNCHRONOUS>;
template class SpatialSearchResultContainer<Element, SpatialSearchCommunication::SYNCHRONOUS>;
template class SpatialSearchResultContainer<Condition, SpatialSearchCommunication::SYNCHRONOUS>;

}  // namespace Kratos