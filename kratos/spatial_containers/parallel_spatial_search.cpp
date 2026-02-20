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
#include "spatial_containers/parallel_spatial_search.h"

namespace Kratos
{

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::~ParallelSpatialSearch()
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
BoundingBox<Point> ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::GetGlobalBoundingBox() const
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
int ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::GetRank() const
{
    return mrDataCommunicator.Rank();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
int ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::GetWorldSize() const
{
    return mrDataCommunicator.Size();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::InitializeGlobalBoundingBoxes()
{
    // We get the world size
    const int world_size = GetWorldSize();

    // Set up the global bounding boxes
    std::vector<double> global_bounding_boxes(6*world_size);
    if (mGlobalBoundingBoxes.size() != static_cast<std::size_t>(world_size)) {
        mGlobalBoundingBoxes.resize(world_size);
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
    mrDataCommunicator.AllGather(local_bounding_box, global_bounding_boxes);

    // Translate to BoundingBox type
    for (int i = 0; i < world_size; ++i) {
        auto& r_bb = mGlobalBoundingBoxes[i];
        auto& r_min_point = r_bb.GetMinPoint();
        auto& r_max_point = r_bb.GetMaxPoint();
        r_max_point[0] = global_bounding_boxes[6 * i + 0];
        r_min_point[0] = global_bounding_boxes[6 * i + 1];
        r_max_point[1] = global_bounding_boxes[6 * i + 2];
        r_min_point[1] = global_bounding_boxes[6 * i + 3];
        r_max_point[2] = global_bounding_boxes[6 * i + 4];
        r_min_point[2] = global_bounding_boxes[6 * i + 5];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<int> ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::RansksPointIsInsideBoundingBox(const array_1d<double, 3>& rCoords)
{
    std::vector<int> ranks;
    const int world_size = GetWorldSize();
    for (int i = 0; i < world_size; ++i) {
        if (SearchUtilities::PointIsInsideBoundingBox(mGlobalBoundingBoxes[i], rCoords)) {
            ranks.push_back(i);
        }
    }

    return ranks;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
std::vector<int> ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::RansksPointIsInsideBoundingBoxWithTolerance(
    const array_1d<double, 3>& rCoords,
    const double Tolerance
    )
{
    std::vector<int> ranks;
    const int world_size = GetWorldSize();
    std::vector<BoundingBox<Point>> bb_tolerance(mGlobalBoundingBoxes.size());
    SearchUtilities::ComputeBoundingBoxesWithToleranceCheckingNullBB(mGlobalBoundingBoxes, Tolerance, bb_tolerance);
    for (int i = 0; i < world_size; ++i) {
        if (SearchUtilities::PointIsInsideBoundingBox(bb_tolerance[i], rCoords)) {
            ranks.push_back(i);
        }
    }

    return ranks;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::LocalSearchInRadius(
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
void ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::LocalSearchNearestInRadius(
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
void ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::LocalSearchNearest(
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
void ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::LocalSearchIsInside(
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
void ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::KeepOnlyClosestResult(ResultContainerVectorType& rResults)
{
    // MPI only
    if (mrDataCommunicator.IsDistributed()) {
        // Getting current rank
        const int rank = mrDataCommunicator.Rank();

        // The values
        std::vector<std::vector<double>> all_distances;
        rResults.GetDistances(all_distances);

        // The ranks
        std::vector<std::vector<int>> all_ranks;
        rResults.GetResultRank(all_ranks);

        // Retrieve the solution
        auto& r_results_vector = rResults.GetContainer();
        for (std::size_t i = 0; i < r_results_vector.size(); ++i) {
            auto& p_partial_result = r_results_vector[i];
            auto& r_partial_result = *p_partial_result;
            // Then must have at least one solution, but just filter if at least 2
            const std::size_t number_of_global_results = r_partial_result.NumberOfGlobalResults();
            if (number_of_global_results > 1) {
                // The values
                const auto& r_values = all_distances[i];

                // The indexes
                std::vector<int> ranks = all_ranks[i];
                std::vector<int>& r_ranks = all_ranks[i];

                // Find the index of the minimum value
                auto it_min_distance = std::min_element(r_values.begin(), r_values.end());

                // Check if the values vector is not empty
                if (it_min_distance != r_values.end()) {
                    // Calculate the position
                    const IndexType pos = std::distance(r_values.begin(), it_min_distance);
                    if (rank == ranks[pos]) {
                        KRATOS_ERROR_IF(r_partial_result.NumberOfLocalResults() > 1) << "The rank criteria to filter results assumes that one rank only holds one local result. This is not true for " << r_partial_result.GetGlobalIndex() << " in rank " << rank << std::endl;
                    }

                    // Remove the index from the ranks vector
                    ranks.erase(ranks.begin() + pos);

                    // Remove all results but the closest one
                    r_partial_result.RemoveResultsFromRanksList(ranks, r_ranks, mrDataCommunicator);
                } else {
                    KRATOS_ERROR << "Distances vector is empty." << std::endl;
                }
            }
        }

        // Checking that is properly cleaned
    #ifdef KRATOS_DEBUG
        for (auto& p_partial_result : r_results_vector) {
            auto& r_partial_result = *p_partial_result;
            // Check that the number of results is 0 or 1
            KRATOS_ERROR_IF(r_partial_result.NumberOfGlobalResults() > 1) << "Cleaning has not been done properly. Number of results: " << r_partial_result.NumberOfGlobalResults() << std::endl;
            // Check that is not empty locally
            if (r_partial_result.NumberOfGlobalResults() == 1) {
                const unsigned int number_of_local_results = r_partial_result.NumberOfLocalResults();
                KRATOS_ERROR_IF(mrDataCommunicator.SumAll(number_of_local_results) == 0) << "Local results also removed in result " << r_partial_result.GetGlobalIndex() << std::endl;
            }
        }
    #endif
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
std::string ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::GenerateNameFromRanks(
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

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::PrepareResultsInProperRanks(
    ResultContainerVectorType& rResults,
    const DistributedSearchInformation& rSearchInfo
    )
{
    // Initialize results
    rResults.InitializeResults(rSearchInfo.TotalNumberOfPoints);

    // Set some values
    const auto& r_local_indices = rSearchInfo.LocalIndices;
    const auto& r_global_indices = rSearchInfo.GlobalIndices;
    auto& r_results_vector = rResults.GetContainer();
    const auto& r_global_position = rSearchInfo.GlobalPosition;
    IndexPartition<IndexType>(r_local_indices.size()).for_each([&r_results_vector, &r_local_indices, &r_global_indices, &r_global_position](const IndexType Index) {
        auto& r_point_result = *(r_results_vector[r_global_position[Index]]);
        r_point_result.SetLocalIndex(r_local_indices[Index]);
        r_point_result.SetGlobalIndex(r_global_indices[Index]);
    });
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
const Parameters ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::GetDefaultParameters() const
{
    return Parameters(R"({
        "allocation_size" : 1000,
        "bucket_size"     : 10
    })");
}

/***********************************************************************************/
/***********************************************************************************/

// GeometricalObjectsBins
template class ParallelSpatialSearch<GeometricalObjectsBins, SpatialSearchCommunication::SYNCHRONOUS>;

// KDTree
template class ParallelSpatialSearch<Tree<KDTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>;
template class ParallelSpatialSearch<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>;
template class ParallelSpatialSearch<Tree<KDTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>;

// OCTree
template class ParallelSpatialSearch<Tree<OCTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>;
template class ParallelSpatialSearch<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>;
template class ParallelSpatialSearch<Tree<OCTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>;

// StaticBinsTree
template class ParallelSpatialSearch<Tree<Bins<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS>;
template class ParallelSpatialSearch<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS>;
template class ParallelSpatialSearch<Tree<Bins<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS>;

// DynamicBins
template class ParallelSpatialSearch<BinsDynamic<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS>;
template class ParallelSpatialSearch<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS>;
template class ParallelSpatialSearch<BinsDynamic<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS>;

}  // namespace Kratos.