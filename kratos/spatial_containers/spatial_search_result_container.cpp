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
SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::SpatialSearchResultContainer()
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

        // Data is not synchronized anymore
        mIsSynchronized = false;
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

        // Data is not synchronized anymore
        mIsSynchronized = false;
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

        // Data is not synchronized anymore
        mIsSynchronized = false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::Clear()
{
    // Clear synchronization flag
    mIsSynchronized = false;

    // Clear local results
    mLocalResults.clear();

    // Clear global results
    mGlobalResults.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::RemoveResultsFromRanksList(
    const std::vector<int>& rRanks,
    const std::vector<int>& rAllRanks,
    const DataCommunicator& rDataCommunicator
    )
{
    // Get current rank
    const int rank = rDataCommunicator.Rank();

    // Remove current rank if required
    auto it_find = std::find(rRanks.begin(), rRanks.end(), rank);
    if (it_find != rRanks.end()) {
        mLocalResults.clear();
    }

    // Prepare a map of current indexes for quick lookup
    std::vector<IndexType> indexes_to_remove;
    indexes_to_remove.reserve(rAllRanks.size());
    for (IndexType i = 0; i < rAllRanks.size(); ++i) {
        auto it_find = std::find(rRanks.begin(), rRanks.end(), rAllRanks[i]);
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
    rOStream << "\tLocal index: " << mLocalIndex << "\n";
    rOStream << "\tGlobal index: " << mGlobalIndex << "\n";
    rOStream << "\tIs data synchronized: " << (mIsSynchronized ? "Yes" : "No") << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::save(Serializer& rSerializer) const
{
    rSerializer.save("LocalResults", mLocalResults);
    rSerializer.save("GlobalResults", mGlobalResults);
    rSerializer.save("LocalIndex", mLocalIndex);
    rSerializer.save("GlobalIndex", mGlobalIndex);
    rSerializer.save("Synchronized", mIsSynchronized);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
void SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>::load(Serializer& rSerializer)
{
    rSerializer.load("LocalResults", mLocalResults);
    rSerializer.load("GlobalResults", mGlobalResults);
    rSerializer.load("LocalIndex", mLocalIndex);
    rSerializer.load("GlobalIndex", mGlobalIndex);
    rSerializer.load("Synchronized", mIsSynchronized);
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
template class SpatialSearchResultContainer<Node, SpatialSearchCommunication::SYNCHRONOUS>;
template class SpatialSearchResultContainer<GeometricalObject, SpatialSearchCommunication::SYNCHRONOUS>;
template class SpatialSearchResultContainer<Element, SpatialSearchCommunication::SYNCHRONOUS>;
template class SpatialSearchResultContainer<Condition, SpatialSearchCommunication::SYNCHRONOUS>;

}  // namespace Kratos