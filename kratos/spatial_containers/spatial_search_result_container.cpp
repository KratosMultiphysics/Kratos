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
#include <functional>

// External includes

// Project includes
#include "includes/data_communicator.h"
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
    // Push_back in local pointers
    TObjectType* p_local_result = rResult.Get().get();
    mLocalPointers.push_back(Kratos::intrusive_ptr<TObjectType>(p_local_result));

    // Add distances
    const IndexType id = p_local_result->Id();
    mLocalDistances.insert({id, rResult.GetDistance()});
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

    // Clear global pointers
    mGlobalPointers.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::SynchronizeAll(const DataCommunicator& rDataCommunicator)
{
    // Synchronize local pointers to global pointers
    mGlobalPointers = GlobalPointerUtilities::GlobalRetrieveGlobalPointers(mLocalPointers, rDataCommunicator);

    // Generate the communicator
    mpGlobalPointerCommunicator = Kratos::make_shared<GlobalPointerCommunicator<TObjectType>>(rDataCommunicator, mGlobalPointers.ptr_begin(), mGlobalPointers.ptr_end());
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
    //rSerializer.load("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

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

template class SpatialSearchResultContainerMap<GeometricalObject>;

}  // namespace Kratos