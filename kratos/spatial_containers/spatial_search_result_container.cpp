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
void SpatialSearchResultContainer<TObjectType>::AddResult(SpatialSearchResultType& rResult)
{
    // Check if the object has been found
    if (rResult.IsObjectFound()) {
        // Adding to the local results
        auto p_result = Kratos::make_shared<SpatialSearchResultType>(rResult);
        mLocalResults.push_back(p_result);
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
        auto p_result = Kratos::make_shared<SpatialSearchResultType>(pResult, Rank);
        mLocalResults.push_back(p_result);
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
        auto p_result = Kratos::make_shared<SpatialSearchResultType>(pResult, Rank);
        p_result->SetDistance(Distance);
        mLocalResults.push_back(p_result);
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
    mGlobalResults = GlobalPointerUtilities::GlobalRetrieveGlobalPointers(mLocalResults, rDataCommunicator);

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

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
SpatialSearchResultContainer<TObjectType>& SpatialSearchResultContainerMap<TObjectType>::InitializeResult(const IndexType Index)
{
    // If doesn't exists, create it
    if (!HasResult(Index)) {
        mPointResults.insert({Index, SpatialSearchResultContainer<TObjectType>()});
    }
    return mPointResults[Index];
}

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
bool SpatialSearchResultContainerMap<TObjectType>::HasResult(const IndexType Index) const
{
    if (mPointResults.find(Index) != mPointResults.end()) {
        return true;
    }
    return false;       
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
    return reinterpret_cast<std::size_t>(&rCoordinates);
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
    rOStream << "SpatialSearchResultContainerMap" << "\n";
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