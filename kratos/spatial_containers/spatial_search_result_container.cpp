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
SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::SpatialSearchResultContainer(const DataCommunicator& rDataCommunicator)
    : mrDataCommunicator(rDataCommunicator)
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
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GenerateGlobalPointerCommunicator()
{
    // Generate the communicator
    mpGlobalPointerCommunicator = Kratos::make_shared<GlobalPointerCommunicatorType>(mrDataCommunicator, mGlobalResults.ptr_begin(), mGlobalResults.ptr_end());
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::Barrier()
{
    // Only in MPI code and in SYNCHRONOUS_HETEROGENEOUS
    if constexpr (TSpatialSearchCommunication == SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS) {
        if(mrDataCommunicator.IsDistributed()) {
            mrDataCommunicator.Barrier();
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::SynchronizeAll()
{
    // TODO: Try to avoid use simple calls of Send/Recv and try to use more efficient methods

    // Synchronize local results to global results
    if(mrDataCommunicator.IsDistributed()) { // MPI code
        const int local_result_size = mLocalResults.size();
        const int global_result_size = mrDataCommunicator.SumAll(local_result_size);
        mGlobalResults.reserve(global_result_size);
        const int world_size = mrDataCommunicator.Size();
        const int rank = mrDataCommunicator.Rank();
        std::vector<int> send_buffer(1, local_result_size);
        std::vector<int> recv_buffer(world_size);
        mrDataCommunicator.AllGather(send_buffer, recv_buffer);
        std::vector<int> send_rank(1, rank);
        std::vector<int> ranks(world_size);
        mrDataCommunicator.AllGather(send_rank, ranks);

        // Retrieve first rank
        const int first_rank = ranks[0];

        // In first rank
        if (rank == first_rank) {
            // Prepare
            std::vector<GlobalPointerResultType> global_gp;
            global_gp.reserve(global_result_size);

            // Fill global vector with local result
            for (auto& r_value : mLocalResults) {
                global_gp.push_back(GlobalPointerResultType(&r_value, rank));
            }

            // Call the lambda to generate the result vector of partitions with results
            std::vector<int> result_vector = GenerateGreaterThanZeroIndexes(recv_buffer);

            // Iterate over the ranks
            for (int rank_to_recv : result_vector) {
                std::vector<GlobalPointerResultType> recv_gps;
                mrDataCommunicator.Recv(recv_gps, rank_to_recv);
                for (auto& r_value : recv_gps) {
                    global_gp.push_back(r_value);
                }
            }

            // Send now to all ranks
            for (int i_rank = 1; i_rank < world_size; ++i_rank) {
                mrDataCommunicator.Send(global_gp, ranks[i_rank]);
            }

            // Transfer to global pointer
            for (auto& r_gp : global_gp) {
                mGlobalResults.push_back(r_gp);
            }
        } else {
            // Sending local results if any
            if (local_result_size > 0) {
                std::vector<GlobalPointerResultType> local_gp;
                local_gp.reserve(local_result_size);
                for (auto& r_value : mLocalResults) {
                    local_gp.push_back(GlobalPointerResultType(&r_value, rank));
                }
                mrDataCommunicator.Send(local_gp, first_rank);
            }

            // Receiving synced result
            std::vector<GlobalPointerResultType> global_gp;
            mrDataCommunicator.Recv(global_gp, first_rank);

            // Transfer to global pointer
            for (auto& r_gp : global_gp) {
                mGlobalResults.push_back(r_gp);
            }
        }
    } else { // Serial code
        // Fill global vector
        for (auto& r_value : mLocalResults) {
            mGlobalResults.push_back(&r_value);
        }
    }

    // Generate the communicator
    GenerateGlobalPointerCommunicator();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<double> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetDistances()
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<double> distances(number_of_gp);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> double {
        return rGP->GetDistance();
    });

    // Get the distances
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        distances[i] = proxy.Get(r_gp);
    }

    return distances;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<bool> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultIsLocal()
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<bool> is_local(number_of_gp, false);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> int {
        return rGP->Get().GetRank();
    });

    // Get the is local
    const int rank = mrDataCommunicator.Rank();
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        const int retrieved_rank = mrDataCommunicator.MaxAll(proxy.Get(r_gp));
        is_local[i] = (rank == retrieved_rank);
    }

    return is_local;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<int> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultRank()
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<int> ranks(number_of_gp, 0);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> int {
        return rGP->Get().GetRank();
    });

    // Get the rank
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        ranks[i] = mrDataCommunicator.MaxAll(proxy.Get(r_gp));
    }

    return ranks;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<bool> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultIsActive()
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<bool> is_active(number_of_gp, false);

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

    // Get the is active
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        is_active[i] = mrDataCommunicator.MaxAll(proxy.Get(r_gp));
    }

    return is_active;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<bool> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultIsInside(
    const array_1d<double, 3>& rPoint,
    const double Tolerance
    )
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<bool> is_inside(number_of_gp, false);

    // Call Apply to get the proxy
    auto proxy = this->Apply([&rPoint, &Tolerance](GlobalPointerResultType& rGP) -> bool {
        auto p_object = rGP->Get();
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {
            auto& r_geometry = p_object->GetGeometry();
            Point::CoordinatesArrayType aux_coords;
            return r_geometry.IsInside(rPoint, aux_coords, Tolerance);
        } else if constexpr (std::is_same<TObjectType, Node>::value) {
            KRATOS_ERROR << "Nodes do not provide is inside. Not possible to compute is inside for point: " << rPoint[0]<< "\t" << rPoint[1] << "\t" << rPoint[2] << " with tolerance " << Tolerance << std::endl;
            return false;
        } else {
            KRATOS_ERROR << "Not implemented yet. Not possible to compute is inside for point: " << rPoint[0]<< "\t" << rPoint[1] << "\t" << rPoint[2] << " with tolerance " << Tolerance << std::endl;
            return false;
        }
    });

    // Get the is inside
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        is_inside[i] = mrDataCommunicator.MaxAll(proxy.Get(r_gp));
    }

    return is_inside;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<Vector> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultShapeFunctions(const array_1d<double, 3>& rPoint)
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<Vector> shape_functions(number_of_gp);

    // Call Apply to get the proxy
    auto proxy = this->Apply([&rPoint](GlobalPointerResultType& rGP) -> Vector {
        auto p_object = rGP->Get();
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {
            auto& r_geometry = p_object->GetGeometry();
            Vector N(r_geometry.size());
            array_1d<double, 3> local_coordinates;
            r_geometry.PointLocalCoordinates(local_coordinates, rPoint);
            r_geometry.ShapeFunctionsValues(N, local_coordinates);
            return N;
        } else if constexpr (std::is_same<TObjectType, Node>::value) {
            KRATOS_ERROR << "Nodes do not provide shape functions. Not possible to compute shape functions for point: " << rPoint[0]<< "\t" << rPoint[1] << "\t" << rPoint[2] << std::endl;
            Vector N;
            return N;
        } else {
            KRATOS_ERROR << "Not implemented yet. Not possible to compute shape functions for point: " << rPoint[0]<< "\t" << rPoint[1] << "\t" << rPoint[2] << std::endl;
            Vector N;
            return N;
        }
    });

    // Get the shape functions
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        shape_functions[i] = proxy.Get(r_gp);
    }

    return shape_functions;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<IndexType> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultIndices()
{
    // Define the indices vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<IndexType> indices(number_of_gp);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> std::size_t {
        auto p_object = rGP->Get();
        return p_object->Id();
    });

    // Get the indices
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        indices[i] = proxy.Get(r_gp);
    }

    return indices;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<std::vector<IndexType>> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultNodeIndices()
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<std::vector<IndexType>> indices(number_of_gp);

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
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        indices[i] = proxy.Get(r_gp);
    }

    return indices;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<std::vector<int>> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultPartitionIndices()
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<std::vector<int>> indices(number_of_gp);

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
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        indices[i] = proxy.Get(r_gp);
    }

    return indices;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<std::vector<array_1d<double, 3>>> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GetResultCoordinates()
{
    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<std::vector<array_1d<double, 3>>> coordinates(number_of_gp);

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
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        coordinates[i] = proxy.Get(r_gp);
    }

    return coordinates;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::RemoveResultsFromRanksList(const std::vector<int>& rRanks)
{
    // Get current rank
    const int rank = mrDataCommunicator.Rank();

    // Remove current rank if required
    auto it_find = std::find(rRanks.begin(), rRanks.end(), rank);
    if (it_find != rRanks.end()) {
        mLocalResults.clear();
    }

    // Prepare a map of current indexes for quick lookup
    const auto ranks = this->GetResultRank();
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
    const auto current_indexes = this->GetResultIndices();
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
std::vector<int> SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::GenerateGreaterThanZeroIndexes(const std::vector<int>& rInputVector)
{
    std::vector<int> indexes;
    for (int i = 1; i < static_cast<int>(rInputVector.size()); ++i) {
        if (rInputVector[i] > 0) {
            indexes.push_back(i);
        }
    }
    return indexes;
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
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::save(Serializer& rSerializer) const
{
    //rSerializer.save("DataCommunicator", mrDataCommunicator); // TODO: DataCommunicator does not define Serializer
    rSerializer.save("LocalResults", mLocalResults);
    rSerializer.save("GlobalResults", mGlobalResults);
    //rSerializer.save("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::load(Serializer& rSerializer)
{
    //rSerializer.load("DataCommunicator", mrDataCommunicator); // TODO: DataCommunicator does not define Serializer
    rSerializer.load("LocalResults", mLocalResults);
    rSerializer.load("GlobalResults", mGlobalResults);
    //rSerializer.load("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
// SYNCHRONOUS_HOMOGENEOUS
template class SpatialSearchResultContainer<Node, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SpatialSearchResultContainer<GeometricalObject, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SpatialSearchResultContainer<Element, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SpatialSearchResultContainer<Condition, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;

// SYNCHRONOUS_HETEROGENEOUS
template class SpatialSearchResultContainer<Node, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SpatialSearchResultContainer<GeometricalObject, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SpatialSearchResultContainer<Element, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SpatialSearchResultContainer<Condition, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;

}  // namespace Kratos