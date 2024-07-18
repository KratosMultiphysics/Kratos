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
#include "includes/key_hash.h"
#include "utilities/search_utilities.h"
#include "spatial_containers/search_wrapper.h"

namespace Kratos
{

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
SearchWrapper<TSearchObject, TSpatialSearchCommunication>::~SearchWrapper()
{
    // Cleanup subdatacommunicators leftovers
    if constexpr (TSpatialSearchCommunication == SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS) {
        for (auto& r_name : mSubDataCommunicatorNames) {
            ParallelEnvironment::UnregisterDataCommunicator(r_name);
        }
    }
}

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
BoundingBox<Point> SearchWrapper<TSearchObject, TSpatialSearchCommunication>::GetGlobalBoundingBox() const
{
    // Generate BB
    BoundingBox<Point> bb;
    auto& r_max = bb.GetMaxPoint();
    auto& r_min = bb.GetMinPoint();

    // We get the global bounding box
    array_1d<double, 3> local_max, local_min;
    if (mpSearchObject) {
        const auto& r_local_bb = mpSearchObject->GetBoundingBox();
        noalias(local_max) = r_local_bb.GetMaxPoint().Coordinates();
        noalias(local_min) = r_local_bb.GetMinPoint().Coordinates();
    } else {
        noalias(local_max) = ZeroVector(3);
        noalias(local_min) = ZeroVector(3);
    }

    // Getting max values
    r_max.Coordinates() = mrDataCommunicator.MaxAll(local_max);

    // Getting min values
    r_min.Coordinates() = mrDataCommunicator.MinAll(local_min);

    return bb;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
int SearchWrapper<TSearchObject, TSpatialSearchCommunication>::GetRank() const
{
    return mrDataCommunicator.Rank();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
int SearchWrapper<TSearchObject, TSpatialSearchCommunication>::GetWorldSize() const
{
    return mrDataCommunicator.Size();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void SearchWrapper<TSearchObject, TSpatialSearchCommunication>::InitializeGlobalBoundingBoxes()
{
    // We get the world size
    const int world_size = GetWorldSize();

    // Set up the global bounding boxes
    if (static_cast<int>(mGlobalBoundingBoxes.size()) != 6*world_size) {
        mGlobalBoundingBoxes.resize(6*world_size);
    }

    // Set up the local bounding boxes
    std::vector<double> local_bounding_box(6);
    array_1d<double, 3> local_max, local_min;
    if (mpSearchObject) {
        const auto& r_local_bb = mpSearchObject->GetBoundingBox();
        noalias(local_max) = r_local_bb.GetMaxPoint().Coordinates();
        noalias(local_min) = r_local_bb.GetMinPoint().Coordinates();
    } else {
        noalias(local_max) = ZeroVector(3);
        noalias(local_min) = ZeroVector(3);
    }

    // Assign local BB
    for (int i = 0; i < 3; ++i) {
        local_bounding_box[2 * i] = local_max[i];
        local_bounding_box[2 * i + 1] = local_min[i];
    }

    // Gather all bounding boxes
    mrDataCommunicator.AllGather(local_bounding_box, mGlobalBoundingBoxes);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<int> SearchWrapper<TSearchObject, TSpatialSearchCommunication>::RansksPointIsInsideBoundingBox(const array_1d<double, 3>& rCoords)
{
    std::vector<int> ranks;
    const int world_size = GetWorldSize();
    std::array<double, 6> local_bb;
    const auto it_begin = mGlobalBoundingBoxes.begin();
    for (int i = 0; i < world_size; ++i) {
        auto vec_it = it_begin + 6 * i;
        for (unsigned int j = 0; j < 6; ++j, ++vec_it) {
            local_bb[j] = *vec_it;
        }
        if (SearchUtilities::PointIsInsideBoundingBox(local_bb, rCoords)) {
            ranks.push_back(i);
        }
    }

    return ranks;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<int> SearchWrapper<TSearchObject, TSpatialSearchCommunication>::RansksPointIsInsideBoundingBoxWithTolerance(
    const array_1d<double, 3>& rCoords,
    const double Tolerance
    )
{
    std::vector<int> ranks;
    const int world_size = GetWorldSize();
    std::array<double, 6> local_bb;
    std::vector<double> bb_tolerance(mGlobalBoundingBoxes.size());
    SearchUtilities::ComputeBoundingBoxesWithToleranceCheckingNullBB(mGlobalBoundingBoxes, Tolerance, bb_tolerance);
    const auto it_begin = bb_tolerance.begin();
    for (int i = 0; i < world_size; ++i) {
        auto vec_it = it_begin + 6 * i;
        for (unsigned int j = 0; j < 6; ++j, ++vec_it) {
            local_bb[j] = *vec_it;
        }
        if (SearchUtilities::PointIsInsideBoundingBox(local_bb, rCoords)) {
            ranks.push_back(i);
        }
    }

    return ranks;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void SearchWrapper<TSearchObject, TSpatialSearchCommunication>::LocalSearchInRadius(
    const PointType& rPoint,
    const double Radius,
    std::vector<ResultType>& rResults,
    const int Rank,
    const int AllocationSize
    )
{
    // If we are using GeometricalObjectBins we can use the optimized search
    if constexpr (IsGeometricalObjectBins) {
        mpSearchObject->SearchInRadius(rPoint, Radius, rResults);
        for(auto& r_result : rResults) {
            r_result.Get().SetRank(Rank);
        }
    } else { // Using trees
        // We search if geometrical objects are provided
        if (mpPointVector) {
            PointVector results(AllocationSize);
            DistanceVector results_distances(AllocationSize);
            const std::size_t number_of_results = mpSearchObject->SearchInRadius(rPoint, Radius, results.begin(), results_distances.begin(), AllocationSize);
            if (number_of_results > 0) {
                // Set the results
                rResults.reserve(number_of_results);
                for (std::size_t i = 0; i < number_of_results; ++i) {
                    auto p_point = results[i];
                    const double distance = results_distances[i];
                    rResults.emplace_back((p_point->pGetObject()).get(), Rank);
                    rResults[i].SetDistance(distance);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void SearchWrapper<TSearchObject, TSpatialSearchCommunication>::LocalSearchNearestInRadius(
    const PointType& rPoint,
    const double Radius,
    ResultType& rResult,
    const int Rank,
    const int AllocationSize
    )
{
    // If we are using GeometricalObjectBins we can use the optimized search
    if constexpr (IsGeometricalObjectBins) {
        rResult = mpSearchObject->SearchNearestInRadius(rPoint, Radius);
        rResult.Get().SetRank(Rank);
    } else { // Using trees
        // We search if geometrical objects are provided
        if (mpPointVector) {
            PointVector results(AllocationSize);
            DistanceVector results_distances(AllocationSize);
            const std::size_t number_of_results = mpSearchObject->SearchInRadius(rPoint, Radius, results.begin(), results_distances.begin(), AllocationSize);
            if (number_of_results > 0) {
                // Resize the results
                results_distances.resize(number_of_results);

                // Find the iterator pointing to the smallest value
                auto it_min = std::min_element(results_distances.begin(), results_distances.end());

                // Calculate the index
                IndexType index = std::distance(results_distances.begin(), it_min);

                // Set the result
                rResult = ResultType((results[index]->pGetObject()).get(), Rank);
                rResult.SetDistance(results_distances[index]);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void SearchWrapper<TSearchObject, TSpatialSearchCommunication>::LocalSearchNearest(
    const PointType& rPoint,
    ResultType& rResult,
    const int Rank
    )
{
    // If we are using GeometricalObjectBins we can use the optimized search
    if constexpr (IsGeometricalObjectBins) {
        rResult = mpSearchObject->SearchNearest(rPoint);
        rResult.Get().SetRank(Rank);
    } else { // Using trees
        // We search if geometrical objects are provided
        if (mpPointVector) {
            // Search nearest
            double distance = 0.0;
            auto p_point = mpSearchObject->SearchNearestPoint(rPoint, distance);

            // Set the result
            rResult = ResultType((p_point->pGetObject()).get(), Rank);
            rResult.SetDistance(distance);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void SearchWrapper<TSearchObject, TSpatialSearchCommunication>::LocalSearchIsInside(
    const PointType& rPoint,
    ResultType& rResult,
    const int Rank
    )
{
    // If we are using GeometricalObjectBins we can use the optimized search
    if constexpr (IsGeometricalObjectBins) {
        rResult = mpSearchObject->SearchIsInside(rPoint);
        rResult.Get().SetRank(Rank);
    } else { // Using trees
        KRATOS_ERROR << "SearchIsInside not compatible with Search trees" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void SearchWrapper<TSearchObject, TSpatialSearchCommunication>::KeepOnlyClosestResult(ResultContainerVectorType& rResults)
{
    auto distance_lambda = [](ResultContainerType& rResult) -> std::vector<double> {
        return rResult.GetDistances();
    };
    KeepOnlyGivenLambdaResult(rResults, distance_lambda);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void SearchWrapper<TSearchObject, TSpatialSearchCommunication>::KeepOnlyLowestRankResult(ResultContainerVectorType& rResults)
{
    auto rank_lambda = [](ResultContainerType& rResult) -> std::vector<int> {
        return rResult.GetResultRank();
    };
    KeepOnlyGivenLambdaResult(rResults, rank_lambda);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
std::string SearchWrapper<TSearchObject, TSpatialSearchCommunication>::GenerateNameFromRanks(
    const std::string& rBaseName,
    const std::vector<int>& rRanks
    )
{
    std::stringstream ss;
    ss << rBaseName;
    for (const int rank : rRanks) {
        ss << rank; // Convert int to string and append to stringstream
    }
    return ss.str();
}

/***********************************************************************************/
/***********************************************************************************/

/* TODO: Move it to key_hash.h, in more generic way */
struct VectorHash {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t seed = 0;
        for (int i : v) {
            HashCombine(seed, i);
        }
        return seed;
    }
};

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void SearchWrapper<TSearchObject, TSpatialSearchCommunication>::PrepareResultsInProperRanks(
    ResultContainerVectorType& rResults,
    const DistributedSearchInformation& rSearchInfo
    )
{
    // If considering global data communicator
    if constexpr (ConsiderGlobalDataCommunicator) {
        // Prepare the data communicators
        std::vector<const DataCommunicator*> data_communicators(rSearchInfo.TotalNumberOfPoints, &mrDataCommunicator);
        // Initialize results
        rResults.InitializeResults(data_communicators);
    } else { // If not considering global data communicator
        // Get the ranks and prepare the data communicators
        const auto& r_ranks = rSearchInfo.Ranks;
        std::vector<const DataCommunicator*> data_communicators(r_ranks.size(), &mrDataCommunicator);

        // The base sub data communicator name
        const std::string base_name = "SubCommunicator_";
        std::unordered_map<std::vector<int>, const DataCommunicator*, VectorHash> data_communicators_database; // NOTE: WE use this to avoid the creating of strings concatenating integers and the search of std::string that is expensive
        for (std::size_t i = 0; i < r_ranks.size(); ++i) {
            const auto& r_current_ranks = r_ranks[i];
            auto it_find = data_communicators_database.find(r_current_ranks);
            // Found
            if (it_find != data_communicators_database.end()) {
                data_communicators[i] = it_find->second;
            } else { // Not found
                // Generate the name
                const std::string name = GenerateNameFromRanks(base_name, r_current_ranks);
                const DataCommunicator& r_sub_communicator = mrDataCommunicator.GetSubDataCommunicator(r_current_ranks, name);
                if constexpr (TSpatialSearchCommunication == SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS) {
                    mSubDataCommunicatorNames.push_back(name);
                }
                data_communicators[i] = &r_sub_communicator;
                data_communicators_database.insert({r_current_ranks, &r_sub_communicator});
            }
        }

        // Initialize results
        rResults.InitializeResults(data_communicators);
    }

    // Set some values
    const auto& r_local_indices = rSearchInfo.LocalIndices;
    const auto& r_global_indices = rSearchInfo.GlobalIndices;
    auto& r_results_vector = rResults.GetContainer();
    if constexpr (ConsiderGlobalDataCommunicator) {
        const auto& r_global_position = rSearchInfo.GlobalPosition;
        IndexPartition<IndexType>(r_local_indices.size()).for_each([&r_results_vector, &r_local_indices, &r_global_indices, &r_global_position](const IndexType Index) {
            auto& r_point_result = *(r_results_vector[r_global_position[Index]]);
            r_point_result.SetLocalIndex(r_local_indices[Index]);
            r_point_result.SetGlobalIndex(r_global_indices[Index]);
        });
    } else {
        IndexPartition<IndexType>(r_local_indices.size()).for_each([&r_results_vector, &r_local_indices, &r_global_indices](const IndexType Index) {
            auto& r_point_result = *(r_results_vector[Index]);
            r_point_result.SetLocalIndex(r_local_indices[Index]);
            r_point_result.SetGlobalIndex(r_global_indices[Index]);
        });
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
const Parameters SearchWrapper<TSearchObject, TSpatialSearchCommunication>::GetDefaultParameters() const
{
    return Parameters(R"({
        "allocation_size"   : 1000,
        "bucket_size"       : 10
    })");
}

/***********************************************************************************/
/***********************************************************************************/

// GeometricalObjectsBins
// SYNCHRONOUS_HOMOGENEOUS
template class SearchWrapper<GeometricalObjectsBins, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
// SYNCHRONOUS_HETEROGENEOUS
template class SearchWrapper<GeometricalObjectsBins, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;

// KDTree
// SYNCHRONOUS_HOMOGENEOUS
template class SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
// SYNCHRONOUS_HETEROGENEOUS
template class SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;

// OCTree
// SYNCHRONOUS_HOMOGENEOUS
template class SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
// SYNCHRONOUS_HETEROGENEOUS
template class SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;

// StaticBinsTree
// SYNCHRONOUS_HOMOGENEOUS
template class SearchWrapper<Tree<Bins<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SearchWrapper<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SearchWrapper<Tree<Bins<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
// SYNCHRONOUS_HETEROGENEOUS
template class SearchWrapper<Tree<Bins<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SearchWrapper<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SearchWrapper<Tree<Bins<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;

// DynamicBins
// SYNCHRONOUS_HOMOGENEOUS
template class SearchWrapper<BinsDynamic<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SearchWrapper<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
template class SearchWrapper<BinsDynamic<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>;
// SYNCHRONOUS_HETEROGENEOUS
template class SearchWrapper<BinsDynamic<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SearchWrapper<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;
template class SearchWrapper<BinsDynamic<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS_HETEROGENEOUS>;

}  // namespace Kratos.