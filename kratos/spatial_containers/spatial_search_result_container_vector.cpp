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
SpatialSearchResultContainerVector<TObjectType>::~SpatialSearchResultContainerVector()
{
    // Make sure to delete the pointers stored in the container
    block_for_each(mPointResults, [](auto p_result) {
        delete p_result;
    });
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
std::size_t SpatialSearchResultContainerVector<TObjectType>::NumberOfSearchResults() const
{
    return mPointResults.size();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
typename SpatialSearchResultContainerVector<TObjectType>::SpatialSearchResultContainerReferenceType SpatialSearchResultContainerVector<TObjectType>::InitializeResult(const DataCommunicator& rDataCommunicator)
{
    // Resize vector
    const std::size_t current_size = mPointResults.size();
    mPointResults.resize(current_size + 1);

    // Create the result
    mPointResults[current_size] = new SpatialSearchResultContainer<TObjectType>(rDataCommunicator);
    return *mPointResults[current_size];
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::InitializeResults(const std::vector<const DataCommunicator*>& rDataCommunicators)
{
    // Define counter
    std::size_t counter = 0;
    const std::size_t current_size = mPointResults.size();
    const std::size_t to_be_added = rDataCommunicators.size();
    const std::size_t new_size = current_size + to_be_added;
    std::vector<IndexType> values_to_initialize(to_be_added, current_size);
    for (auto& r_index : values_to_initialize) {
        r_index += counter;
        ++counter;
    }
    // Resize vector
    mPointResults.resize(new_size);

    // Create the results
    block_for_each(values_to_initialize, [this, &rDataCommunicators](const auto Index) {
        mPointResults[Index] = new SpatialSearchResultContainer<TObjectType>(*rDataCommunicators[Index]);
    });
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
        return true;
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
    // Synchronize local results to global results
    if(rDataCommunicator.IsDistributed()) { // MPI code
        // For the moment just sync manually
        for (auto p_result : mPointResults) {
            p_result->SynchronizeAll();
        }
        // TODO: FIX MPI CODE. Requires coloring and maybe asynchronous communication
    } else { // Serial code
        // Iterate over all the results
        block_for_each(mPointResults, [](auto p_result) {
            // Fill global vector
            auto& r_global_results = p_result->GetGlobalResults();
            for (auto& r_value : p_result->GetLocalResults()) {
                r_global_results.push_back(&r_value);
            }
        });

        // Generate global pointer communicator
        block_for_each(mPointResults, [](auto p_result) {
            p_result->GenerateGlobalPointerCommunicator();
        });
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
    std::size_t counter = 0;
    for (auto p_result : mPointResults) {
        rOStream << "Point " << counter << ": ";
        p_result->PrintData(rOStream);
        ++counter;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TObjectType>
void SpatialSearchResultContainerVector<TObjectType>::CopyingValuesToGlobalResultsVector(
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
    rSerializer.save("PointResults", mPointResults);
}

/***********************************************************************************/
/***********************************************************************************/

/// Template instantiation
template class SpatialSearchResultContainerVector<Node>;
template class SpatialSearchResultContainerVector<GeometricalObject>;
template class SpatialSearchResultContainerVector<Element>;
template class SpatialSearchResultContainerVector<Condition>;

}  // namespace Kratos