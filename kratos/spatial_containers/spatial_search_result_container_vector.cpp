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
    return block_for_each<SumReduction<IndexType>>(mPointResults, [](auto p_result) {
        if (p_result != nullptr) {
            return 1;
        }
        return 0;
    });
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
typename SpatialSearchResultContainerVector<TObjectType>::SpatialSearchResultContainerReferenceType SpatialSearchResultContainerVector<TObjectType>::InitializeResult(const IndexType Index)
{
    // If doesn't exists, create it
    if (!HasResult(Index)) {
        // Resize vector
        mPointResults.resize(Index + 1);

        // Create the result
        mPointResults[Index] = new SpatialSearchResultContainer<TObjectType>();
    }
    return *mPointResults[Index];
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::InitializeResults(const std::vector<IndexType>& rIndexes)
{
    // Get the max index
    const auto it_max_index = std::max_element(rIndexes.begin(), rIndexes.end());
    const IndexType max_index = *it_max_index;

    // Resize vector
    if (max_index >= mPointResults.size()) {
        mPointResults.resize(max_index + 1);
    }

    // Create the results
    block_for_each(rIndexes, [this](const IndexType Index) {
        if (!HasResult(Index)) {
            mPointResults[Index] = new SpatialSearchResultContainer<TObjectType>();
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
bool SpatialSearchResultContainerVector<TObjectType>::HasResult(const IndexType Index) const
{
    // Check size
    if (Index >= mPointResults.size()) {
        return false;
    } else {
        return mPointResults[Index] != nullptr;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::Clear()
{
    mPointResults.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::SynchronizeAll(const DataCommunicator& rDataCommunicator)
{
    // Synchronize local results to global results
    if(rDataCommunicator.IsDistributed()) { // MPI code
        // MPI variables
        const int world_size = rDataCommunicator.Size();
        const int rank = rDataCommunicator.Rank();

        // First lets define the sizes of the local results to send and receive
        std::vector<std::size_t> active_results;
        std::vector<std::size_t> result_global_size;

        // Define the total sizes
        std::size_t total_local_size = 0;
        std::size_t total_global_size = 0;

        // Iterate over all the results
        for (std::size_t i = 0; i < mPointResults.size(); ++i) {
            auto p_result = mPointResults[i];
            if (p_result != nullptr) {
                active_results.push_back(i);
                const int local_result_size = p_result->GetLocalResults().size();
                total_local_size += local_result_size;
                const int global_result_size = rDataCommunicator.SumAll(local_result_size);
                result_global_size.push_back(global_result_size);
                total_global_size += global_result_size;
                p_result->GetGlobalResults().reserve(global_result_size);
            }
        }

        // Gather sizes to send/receive
        std::vector<int> send_buffer(1, total_local_size);
        std::vector<int> recv_buffer(world_size);
        rDataCommunicator.AllGather(send_buffer, recv_buffer);

        // Transfer the data along all partitions
        if (rank == 0) { // In rank 0
            // Prepare
            std::vector<GlobalPointerResultType> global_gp;
            global_gp.reserve(total_global_size);

            // Fill global vector with local result
            for (std::size_t i : active_results) {
                auto p_result = mPointResults[i];
                for (auto& r_value : p_result->GetLocalResults()) {
                    global_gp.push_back(GlobalPointerResultType(&r_value, rank));
                }
            }

            // Call the lambda to generate the result vector of partitions with results
            std::vector<int> resultVector = SpatialSearchResultContainer<TObjectType>::GenerateGreaterThanZeroIndexes(recv_buffer);

            // Iterate over the ranks
            for (int rank_to_recv : resultVector) {
                std::vector<GlobalPointerResultType> recv_gps;
                rDataCommunicator.Recv(recv_gps, rank_to_recv);
                for (auto& r_value : recv_gps) {
                    global_gp.push_back(r_value);
                }
            }

            // Send now to all ranks
            for (int i_rank = 1; i_rank < world_size; ++i_rank) {
                rDataCommunicator.Send(global_gp, i_rank);
            }

            // Copying values to global results
            CopyingValuesToGlobalResultsVector(global_gp, active_results, result_global_size);
        } else {
            // Sending local results if any
            if (total_local_size > 0) {
                std::vector<GlobalPointerResultType> local_gp;
                local_gp.reserve(total_local_size);
                for (std::size_t i : active_results) {
                    auto p_result = mPointResults[i];
                    for (auto& r_value : p_result->GetLocalResults()) {
                        local_gp.push_back(GlobalPointerResultType(&r_value, rank));
                    }
                }
                rDataCommunicator.Send(local_gp, 0);
            }

            // Receiving synced result
            std::vector<GlobalPointerResultType> global_gp;
            rDataCommunicator.Recv(global_gp, 0);

            // Copying values to global results
            CopyingValuesToGlobalResultsVector(global_gp, active_results, result_global_size);
        }
    } else { // Serial code
        // Iterate over all the results
        block_for_each(mPointResults, [](auto p_result) {
            if (p_result != nullptr) {
                // Fill global vector
                auto& r_global_results = p_result->GetGlobalResults();
                for (auto& r_value : p_result->GetLocalResults()) {
                    r_global_results.push_back(&r_value);
                }
            }
        });
    }

    // Generate global pointer communicator
    block_for_each(mPointResults, [&rDataCommunicator](auto p_result) {
        if (p_result != nullptr) {
            p_result->GenerateGlobalPointerCommunicator(rDataCommunicator);
        }
    });
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
    for (auto p_result : mPointResults) {
        p_result->PrintData(rOStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::save(Serializer& rSerializer) const
{
    rSerializer.save("PointResults", mPointResults);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::load(Serializer& rSerializer)
{
    rSerializer.load("PointResults", mPointResults);
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
template class SpatialSearchResultContainerVector<Node>;
template class SpatialSearchResultContainerVector<GeometricalObject>;
template class SpatialSearchResultContainerVector<Element>;
template class SpatialSearchResultContainerVector<Condition>;

}  // namespace Kratos