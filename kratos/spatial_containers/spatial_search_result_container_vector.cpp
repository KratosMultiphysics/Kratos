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
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "spatial_containers/spatial_search_result_container_vector.h"

namespace Kratos
{
template <class TObjectType>
SpatialSearchResultContainerVector<TObjectType>::~SpatialSearchResultContainerVector()
{
    // Make sure to delete the pointers stored in the container
    block_for_each(mPointResults, [](auto p_result) {
        delete p_result;
    });
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::size_t SpatialSearchResultContainerVector<TObjectType>::NumberOfSearchResults() const
{
    // Some check
    KRATOS_DEBUG_ERROR_IF_NOT(mMapIdsStored.size() == mPointResults.size()) << "Allocation non consistent as mMapIdsStored size: " << mMapIdsStored.size() << " and mPointResults size: " << mPointResults.size() << std::endl;
    return mMapIdsStored.size();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
typename SpatialSearchResultContainerVector<TObjectType>::SpatialSearchResultContainerReferenceType SpatialSearchResultContainerVector<TObjectType>::InitializeResult(
    const IndexType Index,
    const DataCommunicator& rDataCommunicator
    )
{
    // Some check
    KRATOS_DEBUG_ERROR_IF_NOT(mMapIdsStored.size() == mPointResults.size()) << "Allocation non consistent as mMapIdsStored size: " << mMapIdsStored.size() << " and mPointResults size: " << mPointResults.size() << std::endl;

    // If doesn't exists, create it
    const auto it_find = mMapIdsStored.find(Index);
    if (it_find == mMapIdsStored.end()) {
        // Resize vector
        const std::size_t current_size = mPointResults.size();
        mPointResults.resize(current_size + 1);

        // Add index to map
        mMapIdsStored.insert({Index, current_size});

        // Create the result
        mPointResults[current_size] = new SpatialSearchResultContainer<TObjectType>(rDataCommunicator);
        return *mPointResults[current_size];
    } else {
        return *mPointResults[it_find->second];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::InitializeResults(
    const std::vector<IndexType>& rIndexes,
    const std::vector<const DataCommunicator*>& rDataCommunicators
    )
{
    // Some check
    KRATOS_DEBUG_ERROR_IF_NOT(mMapIdsStored.size() == mPointResults.size()) << "Allocation non consistent as mMapIdsStored size: " << mMapIdsStored.size() << " and mPointResults size: " << mPointResults.size() << std::endl;

    // Define counter
    std::size_t counter = 0;
    const std::size_t current_size = mPointResults.size();
    std::vector<IndexType> values_to_initialize;
    for (const auto index : rIndexes) {
        const auto it_find = mMapIdsStored.find(index);
        if (it_find == mMapIdsStored.end()) {
            // Add index to map
            mMapIdsStored.insert({index, current_size + counter});
            values_to_initialize.push_back(current_size + counter);

            // Update counters
            ++counter;
        }
    }
    // Resize vector
    mPointResults.resize(current_size + counter);

    // Create the results
    block_for_each(values_to_initialize, [this, &rDataCommunicators](const auto Index) {
        mPointResults[Index] = new SpatialSearchResultContainer<TObjectType>(*rDataCommunicators[Index]);
    });
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
bool SpatialSearchResultContainerVector<TObjectType>::HasResult(const IndexType Index) const
{
    // Some check
    KRATOS_DEBUG_ERROR_IF_NOT(mMapIdsStored.size() == mPointResults.size()) << "Allocation non consistent as mMapIdsStored size: " << mMapIdsStored.size() << " and mPointResults size: " << mPointResults.size() << std::endl;

    // Check if index is initialized
    const auto it_find = mMapIdsStored.find(Index);
    return (it_find != mMapIdsStored.end());
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::Clear()
{
    mPointResults.clear();
    mMapIdsStored.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::SynchronizeAll(const DataCommunicator& rDataCommunicator)
{
    // Synchronize local results to global results
    if(rDataCommunicator.IsDistributed()) { // MPI code
        // For the moment just sync manually
        for (auto p_result : mPointResults) {
            p_result->SynchronizeAll();
        }
        // TODO:: FIX MPI CODE. Requires coloring and maybe asynchronous communication
        // // First lets define the sizes of the local results to send and receive
        // std::vector<std::size_t> active_results;
        // active_results.reserve(mPointResults.size());
        // std::vector<std::size_t> result_global_size;
        // std::vector<GlobalPointerResultType> global_gp;

        // // MPI variables
        // const int world_size = rDataCommunicator.Size();
        // const int rank = rDataCommunicator.Rank();

        // // Define the total sizes
        // std::size_t total_local_size = 0;
        // std::size_t total_global_size = 0;

        // // Iterate over all the results
        // for (std::size_t i = 0; i < mPointResults.size(); ++i) {
        //     auto p_result = mPointResults[i];
        //     active_results.push_back(i);
        //     const int local_result_size = p_result->GetLocalResults().size();
        //     total_local_size += local_result_size;
        //     const auto& r_sub_data_communicator = p_result->GetDataCommunicator();
        //     const int global_result_size = r_sub_data_communicator.SumAll(local_result_size);
        //     result_global_size.push_back(global_result_size);
        //     total_global_size += global_result_size;
        //     p_result->GetGlobalResults().reserve(global_result_size);
        // }

        // // Gather sizes to send/receive
        // std::vector<int> send_buffer(1, total_local_size);
        // std::vector<int> recv_buffer(world_size);
        // rDataCommunicator.AllGather(send_buffer, recv_buffer);

        // // Transfer the data along all partitions
        // if (rank == 0) { // In rank 0
        //     // Prepare
        //     std::vector<GlobalPointerResultType> global_gp;
        //     global_gp.reserve(total_global_size);

        //     // Fill global vector with local result
        //     for (std::size_t i : active_results) {
        //         auto p_result = mPointResults[i];
        //         for (auto& r_value : p_result->GetLocalResults()) {
        //             global_gp.push_back(GlobalPointerResultType(&r_value, rank));
        //         }
        //     }

        //     // Call the lambda to generate the result vector of partitions with results
        //     std::vector<int> result_vector = SpatialSearchResultContainer<TObjectType>::GenerateGreaterThanZeroIndexes(recv_buffer);

        //     // Iterate over the ranks
        //     for (int rank_to_recv : result_vector) {
        //         std::vector<GlobalPointerResultType> recv_gps;
        //         rDataCommunicator.Recv(recv_gps, rank_to_recv);
        //         for (auto& r_value : recv_gps) {
        //             global_gp.push_back(r_value);
        //         }
        //     }

        //     // Send now to all ranks
        //     for (int i_rank = 1; i_rank < world_size; ++i_rank) {
        //         rDataCommunicator.Send(global_gp, i_rank);
        //     }
        // } else {
        //     // Sending local results if any
        //     if (total_local_size > 0) {
        //         std::vector<GlobalPointerResultType> local_gp;
        //         local_gp.reserve(total_local_size);
        //         for (std::size_t i : active_results) {
        //             auto p_result = mPointResults[i];
        //             for (auto& r_value : p_result->GetLocalResults()) {
        //                 local_gp.push_back(GlobalPointerResultType(&r_value, rank));
        //             }
        //         }
        //         rDataCommunicator.Send(local_gp, 0);
        //     }

        //     // Receiving synced result
        //     std::vector<GlobalPointerResultType> global_gp;
        //     rDataCommunicator.Recv(global_gp, 0);
        // }

        // // Copying values to global results
        // CopyingValuesToGlobalResultsVector(global_gp, active_results, result_global_size);
    } else { // Serial code
        // Iterate over all the results
        block_for_each(mPointResults, [](auto p_result) {
            // Fill global vector
            auto& r_global_results = p_result->GetGlobalResults();
            for (auto& r_value : p_result->GetLocalResults()) {
                r_global_results.push_back(&r_value);
            }
        });

        // Generate global pointer communicator
        block_for_each(mPointResults, [](auto p_result) {
            p_result->GenerateGlobalPointerCommunicator();
        });
    }

    // // Generate global pointer communicator
    // block_for_each(mPointResults, [](auto p_result) {
    //     p_result->GenerateGlobalPointerCommunicator();
    // });
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::string SpatialSearchResultContainerVector<TObjectType>::Info() const
{
    std::stringstream buffer;
    buffer << "SpatialSearchResultContainerVector" ;
    return buffer.str();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "SpatialSearchResultContainerVector" << "\n";
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::PrintData(std::ostream& rOStream) const
{
    // Print results
    rOStream << "SpatialSearchResultContainerVector data summary: " << "\n";
    for (auto& r_pair : mMapIdsStored) {
        auto p_result = mPointResults[r_pair.second];
        rOStream << "Point " << r_pair.first << ": ";
        p_result->PrintData(rOStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::CopyingValuesToGlobalResultsVector(
    const std::vector<GlobalPointerResultType>& rGlobalResults,
    const std::vector<std::size_t>& rActiveResults,
    const std::vector<std::size_t>& rResultGlobalSize
    )
{
    // Transfer to global pointer
    std::size_t counter = 1;
    std::size_t index_solution = 0;
    std::size_t global_results_number_current_result = rResultGlobalSize[index_solution]; // I can do rResultGlobalSize[0], for consistency
    auto p_result = mPointResults[rActiveResults[index_solution]]; // I can do mPointResults[rActiveResults[0]], for consistency
    const std::size_t global_results_size = rGlobalResults.size();
    for (std::size_t i_gp = 0; i_gp < global_results_size; ++i_gp) {
        auto& r_gp = rGlobalResults[i_gp];

        // Add to the result
        p_result->GetGlobalResults().push_back(r_gp);

        // Check if jumping to next vector
        if (counter == global_results_number_current_result && i_gp < global_results_size - 1) {
            counter = 0;
            ++index_solution;
            global_results_number_current_result = rResultGlobalSize[index_solution];
            p_result = mPointResults[rActiveResults[index_solution]];
        }

        // Update counter
        ++counter;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::save(Serializer& rSerializer) const
{
    rSerializer.save("PointResults", mPointResults);
    rSerializer.save("MapIdsStored", mMapIdsStored);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::load(Serializer& rSerializer)
{
    rSerializer.save("PointResults", mPointResults);
    rSerializer.save("MapIdsStored", mMapIdsStored);
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
template class SpatialSearchResultContainerVector<Node>;
template class SpatialSearchResultContainerVector<GeometricalObject>;
template class SpatialSearchResultContainerVector<Element>;
template class SpatialSearchResultContainerVector<Condition>;

}  // namespace Kratos