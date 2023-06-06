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
#include "utilities/global_pointer_utilities.h"
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
void SpatialSearchResultContainer<TObjectType>::AddResult(SpatialSearchResult<TObjectType>& rResult)
{
    // Check if the object has been found
    if (rResult.IsObjectFound()) {
        // Push_back in local pointers
        TObjectType* p_local_result = rResult.Get().get();
        mLocalPointers.push_back(p_local_result);

        // Add distances
        const IndexType id = p_local_result->Id();
        mLocalDistances.insert({id, rResult.GetDistance()});
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::AddResult(TObjectType* pResult)
{
    // Check if the object has been found (not nullptr)
    if (pResult != nullptr) {
        mLocalPointers.push_back(pResult);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::Clear()
{
    // Clear pointer
    mpGlobalPointerCommunicator = nullptr;

    // Clear local pointers
    mLocalPointers.clear();

    // Clear local distances
    mLocalDistances.clear();

    // Clear global pointers
    mGlobalPointers.clear();

    // Clear global distances
    mGlobalDistances.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::SynchronizeAll(const DataCommunicator& rDataCommunicator)
{
    // Synchronize local pointers to global pointers
    mGlobalPointers = GlobalPointerUtilities::GlobalRetrieveGlobalPointers(mLocalPointers, rDataCommunicator);

    // Synchronize local distances to global distances
    std::vector<double> local_values;
    for (const auto& r_pair : mLocalDistances) {
        local_values.push_back(r_pair.second);
    }

    // MPI information
    const int world_size = rDataCommunicator.Size();

    // Generate vectors with sizes for AllGatherv
    std::vector<int> recv_sizes(world_size);
    std::vector<int> send_points_per_partition(1, mLocalPointers.size());
    rDataCommunicator.AllGather(send_points_per_partition, recv_sizes);
    std::vector<int> recv_offsets(world_size, 0);
    for (int i_rank = 1; i_rank < world_size; ++i_rank) {
        recv_offsets[i_rank] = recv_offsets[i_rank - 1] + recv_sizes[i_rank - 1];
    }

    // Invoque AllGatherv
    mGlobalDistances.resize(std::accumulate(recv_sizes.begin(), recv_sizes.end(), 0));
    rDataCommunicator.AllGatherv(local_values, mGlobalDistances, recv_sizes, recv_offsets);

    // Generate the communicator
    mpGlobalPointerCommunicator = Kratos::make_shared<GlobalPointerCommunicator<TObjectType>>(rDataCommunicator, mGlobalPointers.ptr_begin(), mGlobalPointers.ptr_end());
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::vector<Vector> SpatialSearchResultContainer<TObjectType>::GetResultShapeFunctions(const array_1d<double, 3>& rPoint)
{
    // Check if the communicator has been created
    KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

    // Define the coordinates vector
    const std::size_t number_of_gp = mGlobalPointers.size();
    std::vector<Vector> shape_functions(number_of_gp);

    // Call Apply to get the proxy
    auto proxy = this->Apply([&rPoint](GlobalPointer<TObjectType>& rGP) -> Vector {
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {   
            auto& r_geometry = rGP->GetGeometry();
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
        auto& r_gp = mGlobalPointers(i);
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
    const std::size_t number_of_gp = mGlobalPointers.size();
    std::vector<std::size_t> indices(number_of_gp);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointer<TObjectType>& rGP) -> std::size_t {
        return rGP->Id();
    });

    // Get the indices
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalPointers(i);
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
    const std::size_t number_of_gp = mGlobalPointers.size();
    std::vector<std::vector<array_1d<double, 3>>> coordinates(number_of_gp);

    // Call Apply to get the proxy
    auto proxy = this->Apply([](GlobalPointer<TObjectType>& rGP) -> std::vector<array_1d<double, 3>> {
        if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {   
            auto& r_geometry = rGP->GetGeometry();
            std::vector<array_1d<double, 3>> coordinates(r_geometry.size());
            for (unsigned int i = 0; i < r_geometry.size(); ++i) {
                coordinates[i] = r_geometry[i].Coordinates();
            }
            return coordinates;
        } else if constexpr (std::is_same<TObjectType, Node>::value) {   
            std::vector<array_1d<double, 3>> coordinates(1, rGP->Coordinates());
            return coordinates;
        } else {   
            KRATOS_ERROR << "Not implemented yet" << std::endl;
            std::vector<array_1d<double, 3>> coordinates;
            return coordinates;
        }
    });

    // Get the coordinates
    for(std::size_t i=0; i<number_of_gp; ++i) {
        auto& r_gp = mGlobalPointers(i);
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
    rOStream << "SpatialSearchResultContainer";
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::PrintData(std::ostream& rOStream) const
{
    rOStream << "SpatialSearchResultContainer data summary: " << "\n";
    rOStream << "\tNumber of local pointers: " << mLocalPointers.size() << "\n";
    rOStream << "\tNumber of global pointers: " << mGlobalPointers.size() << "\n";
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::save(Serializer& rSerializer) const
{ 
    rSerializer.save("LocalPointers", mLocalPointers);
    rSerializer.save("GlobalPointers", mGlobalPointers);
    rSerializer.save("LocalDistances", mLocalDistances);
    rSerializer.save("GlobalDistances", mGlobalDistances);
    //rSerializer.save("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::load(Serializer& rSerializer)
{
    rSerializer.load("LocalPointers", mLocalPointers);
    rSerializer.load("GlobalPointers", mGlobalPointers);
    rSerializer.load("LocalDistances", mLocalDistances);
    rSerializer.load("GlobalDistances", mGlobalDistances);
    //rSerializer.load("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
template class SpatialSearchResultContainer<Node>;
template class SpatialSearchResultContainer<GeometricalObject>;

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
SpatialSearchResultContainer<TObjectType>& SpatialSearchResultContainerMap<TObjectType>::InitializeResult(const array_1d<double, 3>& rCoordinates)
{
    const HashType hash = Hash(rCoordinates);
    // If doesn't exists, create it
    if (!HasResult(rCoordinates)) {
        mPointResults.insert({hash, SpatialSearchResultContainer<TObjectType>()});
    }
    return mPointResults[hash];
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
bool SpatialSearchResultContainerMap<TObjectType>::HasResult(const array_1d<double, 3>& rCoordinates) const
{
    const HashType hash = Hash(rCoordinates);
    if (mPointResults.find(hash) != mPointResults.end()) {
        return true;
    }
    return false;       
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerMap<TObjectType>::Clear()
{
    mPointResults.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerMap<TObjectType>::SynchronizeAll(const DataCommunicator& rDataCommunicator)
{
    // Synchronize all the results
    for (auto& r_point_result : mPointResults) {
        r_point_result.second.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
typename SpatialSearchResultContainerMap<TObjectType>::HashType SpatialSearchResultContainerMap<TObjectType>::Hash(const array_1d<double, 3>& rCoordinates) const
{
    std::hash<array_1d<double, 3>> hasher;
    return hasher(rCoordinates);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::string SpatialSearchResultContainerMap<TObjectType>::Info() const
{
    std::stringstream buffer;
    buffer << "SpatialSearchResultContainerMap" ;
    return buffer.str();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerMap<TObjectType>::PrintInfo(std::ostream& rOStream) const 
{
    rOStream << "SpatialSearchResultContainerMap";
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerMap<TObjectType>::PrintData(std::ostream& rOStream) const
{
    // Print results
    rOStream << "SpatialSearchResultContainerMap data summary: " << "\n";
    for (auto it = mPointResults.begin(); it != mPointResults.end(); ++it) {
        rOStream << "Hash " << it->first << ":\n";
        it->second.PrintData(rOStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerMap<TObjectType>::save(Serializer& rSerializer) const
{ 
    rSerializer.save("PointResults", mPointResults);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerMap<TObjectType>::load(Serializer& rSerializer)
{
    rSerializer.load("PointResults", mPointResults);
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
template class SpatialSearchResultContainerMap<Node>;
template class SpatialSearchResultContainerMap<GeometricalObject>;

}  // namespace Kratos