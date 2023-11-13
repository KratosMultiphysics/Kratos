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

template <class TObjectType>
SpatialSearchResultContainer<TObjectType>::SpatialSearchResultContainer()
{
    // TODO: Add something if required
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::AddResult(SpatialSearchResultType& rResult)
{
    // Check if the object has been found
    if (rResult.IsObjectFound()) {
        // Adding to the local results
        mLocalResults.emplace_back(rResult);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::AddResult(
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

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::AddResult(
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

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::Clear()
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

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::SynchronizeAll(const DataCommunicator& rDataCommunicator)
{
    // Synchronize local results to global results
    if(rDataCommunicator.IsDistributed()) { // MPI code
        const int local_result_size = mLocalResults.size();
        const int global_result_size = rDataCommunicator.SumAll(local_result_size);
        mGlobalResults.reserve(global_result_size);
        const int world_size = rDataCommunicator.Size();
        const int rank = rDataCommunicator.Rank();
        std::vector<int> send_buffer(1, local_result_size);
        std::vector<int> recv_buffer(world_size);
        rDataCommunicator.AllGather(send_buffer, recv_buffer);

        // In rank 0
        if (rank == 0) {
            // Prepare
            std::vector<GlobalPointerResultType> global_gp;
            global_gp.reserve(global_result_size);

            // Fill global vector with local result
            for (auto& r_value : mLocalResults) {
                global_gp.push_back(GlobalPointerResultType(&r_value, rank));
            }

            // Create a lambda to generate the vector of indexes greater than zero
            auto GenerateIndexesLambda = [](const std::vector<int>& rInputVector) {
                std::vector<int> indexes;
                for (int i = 1; i < static_cast<int>(rInputVector.size()); ++i) {
                    if (rInputVector[i] > 0) {
                        indexes.push_back(i);
                    }
                }
                return indexes;
            };

            // Call the lambda to generate the result vector of partitions with results
            std::vector<int> resultVector = GenerateIndexesLambda(recv_buffer);

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
                rDataCommunicator.Send(local_gp, 0);
            }

            // Receiving synced result
            std::vector<GlobalPointerResultType> global_gp;
            rDataCommunicator.Recv(global_gp, 0);

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
    mpGlobalPointerCommunicator = Kratos::make_shared<GlobalPointerCommunicatorType>(rDataCommunicator, mGlobalResults.ptr_begin(), mGlobalResults.ptr_end());
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::vector<double> SpatialSearchResultContainer<TObjectType>::GetDistances()
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<double> distances(number_of_gp);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> double {
        return rGP->GetDistance();
    });

    // Get the coordinates
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        distances[i] = proxy.Get(r_gp);
    }

    return distances;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::vector<bool> SpatialSearchResultContainer<TObjectType>::GetResultIsLocal()
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<bool> is_local(number_of_gp, false);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> bool {
        return rGP->Get().GetRank();
    });

    // Get the is inside
    const auto& r_data_comm = mpGlobalPointerCommunicator->GetDataCommunicator();
    const int rank = r_data_comm.Rank();
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        is_local[i] = (rank == proxy.Get(r_gp));
    }

    return is_local;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::vector<bool> SpatialSearchResultContainer<TObjectType>::GetResultIsActive()
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

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

    // Get the is inside
    const auto& r_data_comm = mpGlobalPointerCommunicator->GetDataCommunicator();
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        is_active[i] = r_data_comm.MaxAll(proxy.Get(r_gp));
    }

    return is_active;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::vector<bool> SpatialSearchResultContainer<TObjectType>::GetResultIsInside(
    const array_1d<double, 3>& rPoint,
    const double Tolerance
    )
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

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
    const auto& r_data_comm = mpGlobalPointerCommunicator->GetDataCommunicator();
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        is_inside[i] = r_data_comm.MaxAll(proxy.Get(r_gp));
    }

    return is_inside;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::vector<Vector> SpatialSearchResultContainer<TObjectType>::GetResultShapeFunctions(const array_1d<double, 3>& rPoint)
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

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

template <class TObjectType>
std::vector<std::size_t> SpatialSearchResultContainer<TObjectType>::GetResultIndices()
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

    // Define the indices vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<std::size_t> indices(number_of_gp);

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

template <class TObjectType>
std::vector<std::vector<std::size_t>> SpatialSearchResultContainer<TObjectType>::GetResultNodeIndices()
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<std::vector<std::size_t>> indices(number_of_gp);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> std::vector<std::size_t> {
        auto p_object = rGP->Get();
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {
            auto& r_geometry = p_object->GetGeometry();
            std::vector<std::size_t> gp_indices(r_geometry.size());
            for (unsigned int i = 0; i < r_geometry.size(); ++i) {
                gp_indices[i] = r_geometry[i].Id();
            }
            return gp_indices;
        } else if constexpr (std::is_same<TObjectType, Node>::value) {
            std::vector<std::size_t> gp_indices(1, p_object->Id());
            return gp_indices;
        } else {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
            std::vector<std::size_t> gp_indices;
            return gp_indices;
        }
    });

    // Get the coordinates
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        indices[i] = proxy.Get(r_gp);
    }

    return indices;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::vector<std::vector<std::size_t>> SpatialSearchResultContainer<TObjectType>::GetResultPartitionIndices()
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalResults.size();
    std::vector<std::vector<std::size_t>> indices(number_of_gp);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointerResultType& rGP) -> std::vector<std::size_t> {
        auto p_object = rGP->Get();
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {
            auto& r_geometry = p_object->GetGeometry();
            std::vector<std::size_t> gp_indices(r_geometry.size());
            for (unsigned int i = 0; i < r_geometry.size(); ++i) {
                gp_indices[i] = r_geometry[i].FastGetSolutionStepValue(PARTITION_INDEX);
            }
            return gp_indices;
        } else if constexpr (std::is_same<TObjectType, Node>::value) {
            std::vector<std::size_t> gp_indices(1, p_object->FastGetSolutionStepValue(PARTITION_INDEX));
            return gp_indices;
        } else {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
            std::vector<std::size_t> gp_indices;
            return gp_indices;
        }
    });

    // Get the coordinates
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalResults(i);
        indices[i] = proxy.Get(r_gp);
    }

    return indices;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::vector<std::vector<array_1d<double, 3>>> SpatialSearchResultContainer<TObjectType>::GetResultCoordinates()
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

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

template <class TObjectType>
std::string SpatialSearchResultContainer<TObjectType>::Info() const
{
    std::stringstream buffer;
    buffer << "SpatialSearchResultContainer" ;
    return buffer.str();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "SpatialSearchResultContainer" << "\n";
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::PrintData(std::ostream& rOStream) const
{
    rOStream << "SpatialSearchResultContainer data summary: " << "\n";
    rOStream << "\tNumber of local results: " << mLocalResults.size() << "\n";
    rOStream << "\tNumber of global results: " << mGlobalResults.size() << "\n";
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::save(Serializer& rSerializer) const
{
    rSerializer.save("LocalResults", mLocalResults);
    rSerializer.save("GlobalResults", mGlobalResults);
    //rSerializer.save("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::load(Serializer& rSerializer)
{
    rSerializer.load("LocalResults", mLocalResults);
    rSerializer.load("GlobalResults", mGlobalResults);
    //rSerializer.load("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
template class SpatialSearchResultContainer<Node>;
template class SpatialSearchResultContainer<GeometricalObject>;
template class SpatialSearchResultContainer<Element>;
template class SpatialSearchResultContainer<Condition>;

}  // namespace Kratos