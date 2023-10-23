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
std::size_t SpatialSearchResultContainerVector<TObjectType>::NumberOfSearchResults() const
{
    return block_for_each<SumReduction<IndexType>>(mPointResults, [&](auto p_result) {
        if (p_result != nullptr) {
            return 1;
        }
        return 0;
    });
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
SpatialSearchResultContainer<TObjectType>& SpatialSearchResultContainerVector<TObjectType>::InitializeResult(const IndexType Index)
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
    // Synchronize all the results
    for (auto p_result : mPointResults) {
        p_result->SynchronizeAll(rDataCommunicator);
    }
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

}  // namespace Kratos