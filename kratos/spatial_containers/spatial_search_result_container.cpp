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

// External includes

// Project includes
#include "includes/data_communicator.h"
#include "includes/geometrical_object.h"
#include "utilities/global_pointer_utilities.h"
#include "spatial_containers/spatial_search_result_container.h"

namespace Kratos
{

template <class TObjectType>
SpatialSearchResultContainer<TObjectType>::SpatialSearchResultContainer(const DataCommunicator& rDataCommunicator)
    : mrDataCommunicator(rDataCommunicator)
{
    // TODO: Add something if required
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::AddResult(SpatialSearchResult<TObjectType>& rResult)
{
    // Push_back in local pointers
    mLocalPointers.push_back(Kratos::intrusive_ptr<TObjectType>(rResult.Get().get()));
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
void SpatialSearchResultContainer<TObjectType>::SynchronizeAll()
{
    // Synchronize local pointers to global pointers
    mGlobalPointers = GlobalPointerUtilities::GlobalRetrieveGlobalPointers(mLocalPointers, mrDataCommunicator);

    // Generate the communicator
    mpGlobalPointerCommunicator = Kratos::make_shared<GlobalPointerCommunicator<TObjectType>>(mrDataCommunicator, mGlobalPointers.ptr_begin(), mGlobalPointers.ptr_end());
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
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::save(Serializer& rSerializer) const
{ 
    rSerializer.save("LocalPointers", mLocalPointers);
    rSerializer.save("GlobalPointers", mGlobalPointers);
    //rSerializer.save("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainer<TObjectType>::load(Serializer& rSerializer)
{
    rSerializer.load("LocalPointers", mLocalPointers);
    rSerializer.load("GlobalPointers", mGlobalPointers);
    //rSerializer.load("GlobalPointerCommunicator", mpGlobalPointerCommunicator); // Not necessary, is created and filled during use
}

/***********************************************************************************/
/***********************************************************************************/

template class SpatialSearchResultContainer<GeometricalObject>;

}  // namespace Kratos