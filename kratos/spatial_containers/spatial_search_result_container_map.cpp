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
#include "spatial_containers/spatial_search_result_container_map.h"

namespace Kratos
{

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