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
template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::~SpatialSearchResultContainerVector()
{
    // Make sure to delete the pointers stored in the container
    block_for_each(mPointResults, [](auto p_result) {
        delete p_result;
    });
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::size_t SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::NumberOfSearchResults() const
{
    return mPointResults.size();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
typename SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::SpatialSearchResultContainerReferenceType SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::InitializeResult()
{
    // Resize vector
    const std::size_t current_size = mPointResults.size();
    mPointResults.resize(current_size + 1);

    // Create the result
    mPointResults[current_size] = new SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>();
    return *mPointResults[current_size];
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::InitializeResults(const std::size_t NumberOfResults)
{
    // Define counter
    std::size_t counter = 0;
    const std::size_t current_size = mPointResults.size();
    const std::size_t new_size = current_size + NumberOfResults;
    std::vector<IndexType> values_to_initialize(NumberOfResults, current_size);
    for (auto& r_index : values_to_initialize) {
        r_index += counter;
        ++counter;
    }
    // Resize vector
    mPointResults.resize(new_size);

    // Create the results
    block_for_each(values_to_initialize, [this](const auto Index) {
        mPointResults[Index] = new SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>();
    });
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
bool SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::HasResult(const IndexType Index) const
{
    // Check size
    return Index < mPointResults.size();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::Clear()
{
    mPointResults.clear();
    mpGlobalPointerCommunicator = nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GenerateGlobalPointerCommunicator(const DataCommunicator& rDataCommunicator)
{
    // Prepare the global results
    GlobalResultsVector all_global_results;
    const std::size_t size_vector = block_for_each<SumReduction<std::size_t>>(mPointResults, [&](auto p_result){
        return p_result->NumberOfGlobalResults();
    });
    all_global_results.reserve(size_vector);
    const std::size_t number_of_global_solutions = mPointResults.size();
    for(std::size_t i=0; i<number_of_global_solutions; ++i) {
        auto& r_global_results = mPointResults[i]->GetGlobalResults();
        const std::size_t number_of_gp = r_global_results.size();
        for(std::size_t i=0; i<number_of_gp; ++i) {
            auto& r_gp = r_global_results(i);
            all_global_results.push_back(r_gp);
        }
    }

    // Generate the communicator
    GenerateGlobalPointerCommunicator(all_global_results, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GenerateGlobalPointerCommunicator(
    const GlobalResultsVector& rAllGlobalResults,
    const DataCommunicator& rDataCommunicator
    )
{
    // Generate the communicator
    mpGlobalPointerCommunicator = Kratos::make_shared<GlobalPointerCommunicatorType>(rDataCommunicator, rAllGlobalResults.ptr_begin(), rAllGlobalResults.ptr_end());
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::SynchronizeAll(const DataCommunicator& rDataCommunicator)
{
    // Prepare the global results
    GlobalResultsVector all_global_results;

    // Synchronize local results to global results
    if(rDataCommunicator.IsDistributed()) { // MPI code
        // Get MPI info once before the loop
        const int world_size = rDataCommunicator.Size();
        const int rank = rDataCommunicator.Rank();
        const int root_rank = 0;

        // Declare reusable buffers outside the loop to avoid repeated memory allocations
        std::vector<int> result_counts(world_size);
        std::vector<GlobalPointerResultType> global_gp;

        // Iterate over all the results
        for (auto p_result : mPointResults) {
            auto& r_local_results = p_result->GetLocalResults();
            auto& r_global_results = p_result->GetGlobalResults();

            // 1. All ranks determine their local result size for the current iteration.
            const int local_result_size = r_local_results.size();

            // 2. Gather result counts from every rank.
            // This is more efficient than calling SumAll and then AllGather.
            rDataCommunicator.AllGather({local_result_size}, result_counts);

            // 3. Calculate total size locally from the gathered counts.
            const int global_result_size = std::accumulate(result_counts.begin(), result_counts.end(), 0);
            r_global_results.reserve(global_result_size);

            // 4. The "gather-to-root, then broadcast" pattern begins.
            global_gp.clear(); // Clear the temporary buffer for this iteration
            if (rank == root_rank) {
                // --- ROOT PROCESS (GATHER & PREPARE BROADCAST) ---
                global_gp.reserve(global_result_size);

                // Add its own local results.
                for (auto& r_value : r_local_results) {
                    global_gp.emplace_back(&r_value, rank);
                }

                // Receive results from all other ranks that have data.
                for (int i_rank = 0; i_rank < world_size; ++i_rank) {
                    if (i_rank != root_rank && result_counts[i_rank] > 0) {
                        std::vector<GlobalPointerResultType> received_gps;
                        rDataCommunicator.Recv(received_gps, i_rank);
                        global_gp.insert(global_gp.end(), std::make_move_iterator(received_gps.begin()), std::make_move_iterator(received_gps.end()));
                    }
                }
            } else {
                // --- WORKER PROCESSES (SEND TO ROOT) ---
                if (local_result_size > 0) {
                    std::vector<GlobalPointerResultType> local_gp;
                    local_gp.reserve(local_result_size);
                    for (auto& r_value : r_local_results) {
                        local_gp.emplace_back(&r_value, rank);
                    }
                    rDataCommunicator.Send(local_gp, root_rank);
                }
            }

            // 5. Broadcast the complete data from the root to ALL ranks (including itself).
            // This single collective call is more efficient and robust than a loop of Sends/Recvs.
            rDataCommunicator.Broadcast(global_gp, root_rank);

            // 6. Populate the final Kratos container from the temporary buffer.
            r_global_results.clear();
            const std::size_t number_of_gp = global_gp.size();
            r_global_results.reserve(number_of_gp);
            const std::size_t current_size = all_global_results.size();
            all_global_results.reserve(current_size + number_of_gp);
            for (auto& r_gp : global_gp) {
                r_global_results.push_back(std::move(r_gp));
                all_global_results.push_back(std::move(r_gp));
            }

            // 7. Mark the result as synchronized
            p_result->SetIsSynchronized(true);
        }
    } else { // Serial code
        // Iterate over all the results
        for (auto p_result : mPointResults) {
            // Fill global vector
            auto& r_global_results = p_result->GetGlobalResults();
            const std::size_t number_of_gp = p_result->NumberOfLocalResults();
            r_global_results.reserve(number_of_gp);
            const std::size_t current_size = all_global_results.size();
            all_global_results.reserve(current_size + number_of_gp);
            for (auto& r_value : p_result->GetLocalResults()) {
                r_global_results.push_back(&r_value);
                all_global_results.push_back(&r_value);
            }

            // Mark the result as synchronized
            p_result->SetIsSynchronized(true);
        }
    }

    // Generate the communicator
    GenerateGlobalPointerCommunicator(all_global_results, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GetDistances(std::vector<std::vector<double>>& rResults)
{
    // Define the coordinates vector
    const std::size_t number_of_global_solutions = mPointResults.size();
    if (rResults.size() != number_of_global_solutions) {
        rResults.resize(number_of_global_solutions);
    }

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> double {
        return rGP->GetDistance();
    });

    // Get the distances
    for(std::size_t i = 0; i < number_of_global_solutions; ++i) {
        auto& r_global_results = mPointResults[i]->GetGlobalResults();
        const std::size_t number_of_gp = r_global_results.size();
        auto& r_distances = rResults[i];
        r_distances.resize(number_of_gp);
        for(std::size_t j = 0; j < number_of_gp; ++j) {
            auto& r_gp = r_global_results(j);
            r_distances[j] = proxy.Get(r_gp);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GetResultIsLocal(
    std::vector<std::vector<bool>>& rResults,
    const DataCommunicator& rDataCommunicator
    )
{
    // Define the coordinates vector
    const std::size_t number_of_global_solutions = mPointResults.size();
    if (rResults.size() != number_of_global_solutions) {
        rResults.resize(number_of_global_solutions);
    }

    // Get the is local
    const int rank = rDataCommunicator.Rank();
    for(std::size_t i = 0; i < number_of_global_solutions; ++i) {
        auto& r_global_results = mPointResults[i]->GetGlobalResults();
        const std::size_t number_of_gp = r_global_results.size();
        auto& r_is_local = rResults[i];
        r_is_local.resize(number_of_gp);
        for(std::size_t j = 0; j < number_of_gp; ++j) {
            auto& r_gp = r_global_results(j);
            const int retrieved_rank = r_gp.GetRank();
            r_is_local[j] = (rank == retrieved_rank);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GetResultRank(std::vector<std::vector<int>>& rResults)
{
    // Define the coordinates vector
    const std::size_t number_of_global_solutions = mPointResults.size();
    if (rResults.size() != number_of_global_solutions) {
        rResults.resize(number_of_global_solutions);
    }

    // Get the ranks
    for(std::size_t i = 0; i < number_of_global_solutions; ++i) {
        auto& r_global_results = mPointResults[i]->GetGlobalResults();
        const std::size_t number_of_gp = r_global_results.size();
        auto& r_ranks = rResults[i];
        r_ranks.resize(number_of_gp);
        for(std::size_t j = 0; j < number_of_gp; ++j) {
            auto& r_gp = r_global_results(j);
            r_ranks[j] = r_gp.GetRank();
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GetResultIsActive(
    std::vector<std::vector<bool>>& rResults,
    const DataCommunicator& rDataCommunicator
    )
{
    // Define the coordinates vector
    const std::size_t number_of_global_solutions = mPointResults.size();
    if (rResults.size() != number_of_global_solutions) {
        rResults.resize(number_of_global_solutions);
    }

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> bool {
        auto p_object = rGP->Get();
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value || std::is_same<TObjectType, Node>::value) {
            return p_object->IsActive();
        } else {
            KRATOS_ERROR << "Not implemented yet. Not possible to compute is active for point." << std::endl;
            return false;
        }
    });

    // Get the is inside
    for(std::size_t i = 0; i < number_of_global_solutions; ++i) {
        auto& r_global_results = mPointResults[i]->GetGlobalResults();
        const std::size_t number_of_gp = r_global_results.size();
        auto& r_is_active = rResults[i];
        r_is_active.resize(number_of_gp);
        for(std::size_t j = 0; j < number_of_gp; ++j) {
            auto& r_gp = r_global_results(j);
            r_is_active[j] = rDataCommunicator.MaxAll(proxy.Get(r_gp));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GetResultIndices(std::vector<std::vector<IndexType>>& rResults)
{
    // Define the indices vector
    const std::size_t number_of_global_solutions = mPointResults.size();
    if (rResults.size() != number_of_global_solutions) {
        rResults.resize(number_of_global_solutions);
    }

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> std::size_t {
        auto p_object = rGP->Get();
        return p_object->Id();
    });

    // Get the indices
    for(std::size_t i = 0; i < number_of_global_solutions; ++i) {
        auto& r_global_results = mPointResults[i]->GetGlobalResults();
        const std::size_t number_of_gp = r_global_results.size();
        auto& r_indices = rResults[i];
        r_indices.resize(number_of_gp);
        for(std::size_t j = 0; j < number_of_gp; ++j) {
            auto& r_gp = r_global_results(j);
            r_indices[j] = proxy.Get(r_gp);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GetResultNodeIndices(std::vector<std::vector<std::vector<IndexType>>>& rResults)
{
    // Define the coordinates vector
    const std::size_t number_of_global_solutions = mPointResults.size();
    if (rResults.size() != number_of_global_solutions) {
        rResults.resize(number_of_global_solutions);
    }

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> std::vector<IndexType> {
        auto p_object = rGP->Get();
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {
            auto& r_geometry = p_object->GetGeometry();
            std::vector<IndexType> gp_indices(r_geometry.size());
            for (unsigned int i = 0; i < r_geometry.size(); ++i) {
                gp_indices[i] = r_geometry[i].Id();
            }
            return gp_indices;
        } else if constexpr (std::is_same<TObjectType, Node>::value) {
            std::vector<IndexType> gp_indices(1, p_object->Id());
            return gp_indices;
        } else {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
            std::vector<IndexType> gp_indices;
            return gp_indices;
        }
    });

    // Get the indices
    for(std::size_t i=0; i<number_of_global_solutions; ++i) {
        auto& r_global_results = mPointResults[i]->GetGlobalResults();
        const std::size_t number_of_gp = r_global_results.size();
        auto& r_indices = rResults[i];
        r_indices.resize(number_of_gp);
        for(std::size_t j=0; j<number_of_gp; ++j) {
            auto& r_gp = r_global_results(j);
            r_indices[j] = proxy.Get(r_gp);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GetResultPartitionIndices(std::vector<std::vector<std::vector<int>>>& rResults)
{
    // Define the coordinates vector
    const std::size_t number_of_global_solutions = mPointResults.size();
    if (rResults.size() != number_of_global_solutions) {
        rResults.resize(number_of_global_solutions);
    }

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> std::vector<int> {
        auto p_object = rGP->Get();
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {
            auto& r_geometry = p_object->GetGeometry();
            std::vector<int> gp_indices(r_geometry.size());
            for (unsigned int i = 0; i < r_geometry.size(); ++i) {
                gp_indices[i] = r_geometry[i].FastGetSolutionStepValue(PARTITION_INDEX);
            }
            return gp_indices;
        } else if constexpr (std::is_same<TObjectType, Node>::value) {
            std::vector<int> gp_indices(1, p_object->FastGetSolutionStepValue(PARTITION_INDEX));
            return gp_indices;
        } else {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
            std::vector<int> gp_indices;
            return gp_indices;
        }
    });

    // Get the indices
    for(std::size_t i=0; i<number_of_global_solutions; ++i) {
        auto& r_global_results = mPointResults[i]->GetGlobalResults();
        const std::size_t number_of_gp = r_global_results.size();
        auto& r_indices = rResults[i];
        r_indices.resize(number_of_gp);
        for(std::size_t j=0; j<number_of_gp; ++j) {
            auto& r_gp = r_global_results(j);
            r_indices[j] = proxy.Get(r_gp);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::GetResultCoordinates(std::vector<std::vector<std::vector<array_1d<double, 3>>>>& rResults)
{
    // Define the coordinates vector
    const std::size_t number_of_global_solutions = mPointResults.size();
    if (rResults.size() != number_of_global_solutions) {
        rResults.resize(number_of_global_solutions);
    }

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> std::vector<array_1d<double, 3>> {
        auto p_object = rGP->Get();
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {
            auto& r_geometry = p_object->GetGeometry();
            std::vector<array_1d<double, 3>> coordinates(r_geometry.size());
            for (unsigned int i = 0; i < r_geometry.size(); ++i) {
                coordinates[i] = r_geometry[i].Coordinates();
            }
            return coordinates;
        } else if constexpr (std::is_same<TObjectType, Node>::value) {
            std::vector<array_1d<double, 3>> coordinates(1, p_object->Coordinates());
            return coordinates;
        } else {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
            std::vector<array_1d<double, 3>> coordinates;
            return coordinates;
        }
    });

    // Get the coordinates
    for(std::size_t i=0; i<number_of_global_solutions; ++i) {
        auto& r_global_results = mPointResults[i]->GetGlobalResults();
        const std::size_t number_of_gp = r_global_results.size();
        auto& r_coordinates = rResults[i];
        r_coordinates.resize(number_of_gp);
        for(std::size_t j=0; j<number_of_gp; ++j) {
            auto& r_gp = r_global_results(j);
            r_coordinates[j] = proxy.Get(r_gp);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::string SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::Info() const
{
    std::stringstream buffer;
    buffer << "SpatialSearchResultContainerVector" ;
    return buffer.str();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "SpatialSearchResultContainerVector" << "\n";
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::PrintData(std::ostream& rOStream) const
{
    // Print results
    rOStream << "SpatialSearchResultContainerVector data summary: " << "\n";
    std::size_t counter = 0;
    for (auto p_result : mPointResults) {
        rOStream << "Point " << counter << ": ";
        p_result->PrintData(rOStream);
        ++counter;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::CopyingValuesToGlobalResultsVector(
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

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::save(Serializer& rSerializer) const
{
    //rSerializer.save("PointResults", mPointResults);
    //rSerializer.save("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>::load(Serializer& rSerializer)
{
    //rSerializer.load("PointResults", mPointResults);
    //rSerializer.load("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
template class SpatialSearchResultContainerVector<Node, SpatialSearchCommunication::SYNCHRONOUS>;
template class SpatialSearchResultContainerVector<GeometricalObject, SpatialSearchCommunication::SYNCHRONOUS>;
template class SpatialSearchResultContainerVector<Element, SpatialSearchCommunication::SYNCHRONOUS>;
template class SpatialSearchResultContainerVector<Condition, SpatialSearchCommunication::SYNCHRONOUS>;


}  // namespace Kratos