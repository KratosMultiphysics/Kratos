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
#include "utilities/search_utilities.h"
#include "spatial_containers/search_wrapper.h"

namespace Kratos
{

template<class TSearchObject>
BoundingBox<Point> SearchWrapper<TSearchObject>::GetBoundingBox() const
{
    // Generate BB
    BoundingBox<Point> bb;

    // We get the global bounding box
    if (mpSearchObject) {
        auto& r_max = bb.GetMaxPoint();
        auto& r_min = bb.GetMinPoint();

        const auto& r_local_bb = mpSearchObject->GetBoundingBox();
        const auto& r_local_max = r_local_bb.GetMaxPoint();
        const auto& r_local_min = r_local_bb.GetMinPoint();

        // Getting max values
        r_max[0] = mrDataCommunicator.MaxAll(r_local_max[0]);
        r_max[1] = mrDataCommunicator.MaxAll(r_local_max[1]);
        r_max[2] = mrDataCommunicator.MaxAll(r_local_max[2]);

        // Getting min values
        r_min[0] = mrDataCommunicator.MinAll(r_local_min[0]);
        r_min[1] = mrDataCommunicator.MinAll(r_local_min[1]);
        r_min[2] = mrDataCommunicator.MinAll(r_local_min[2]);
    }

    return bb;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
int SearchWrapper<TSearchObject>::GetRank() const
{
    return mrDataCommunicator.Rank();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
int SearchWrapper<TSearchObject>::GetWorldSize() const
{
    return mrDataCommunicator.Size();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
void SearchWrapper<TSearchObject>::InitializeGlobalBoundingBoxes()
{
    // Just executed in MPI
    if (mrDataCommunicator.IsDistributed() && mpSearchObject) {
        // We get the world size
        const int world_size = GetWorldSize();

        // Set up the global bounding boxes
        if (static_cast<int>(mGlobalBoundingBoxes.size()) != 6*world_size) {
            mGlobalBoundingBoxes.resize(6*world_size);
        }

        // Set up the local bounding boxes
        std::vector<double> local_bounding_box(6);
        const auto& r_local_bb = mpSearchObject->GetBoundingBox();
        const auto& r_max = r_local_bb.GetMaxPoint();
        const auto& r_min = r_local_bb.GetMinPoint();
        for (int i = 0; i < 3; ++i) {
            local_bounding_box[2 * i] = r_max[i];
            local_bounding_box[2 * i + 1] = r_min[i];
        }

        // Gather all bounding boxes
        mrDataCommunicator.AllGather(local_bounding_box, mGlobalBoundingBoxes);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
std::vector<int> SearchWrapper<TSearchObject>::RansksPointIsInsideBoundingBox(const array_1d<double, 3>& rCoords)
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

template<class TSearchObject>
std::vector<int> SearchWrapper<TSearchObject>::RansksPointIsInsideBoundingBoxWithTolerance(
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

template<class TSearchObject>
void SearchWrapper<TSearchObject>::LocalSearchInRadius(
    const PointType& rPoint,
    const double Radius,
    std::vector<ResultType>& rResults,
    const int AllocationSize
    )
{
    // If we are using GeometricalObjectBins we can use the optimized search
    if constexpr (IsGeometricalObjectBins) {
        mpSearchObject->SearchInRadius(rPoint, Radius, rResults);
    } else { // Using trees
        // We search if geometrical objects are provided
        if (mpPointVector) {
            PointVector results(AllocationSize);
            DistanceVector results_distances(AllocationSize);
            const std::size_t number_of_results = mpSearchObject->SearchInRadius(rPoint, Radius, results.begin(), results_distances.begin(), AllocationSize);
            if (number_of_results > 0) {
                // Get the rank
                const int rank = GetRank();

                // Set the results
                rResults.reserve(number_of_results);
                for (std::size_t i = 0; i < number_of_results; ++i) {
                    auto p_point = results[i];
                    const double distance = results_distances[i];
                    rResults.emplace_back((p_point->pGetObject()).get(), rank);
                    rResults[i].SetDistance(distance);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
void SearchWrapper<TSearchObject>::LocalSearchNearestInRadius(
    const PointType& rPoint,
    const double Radius,
    ResultType& rResult,
    const int AllocationSize
    )
{
    // If we are using GeometricalObjectBins we can use the optimized search
    if constexpr (IsGeometricalObjectBins) {
        rResult = mpSearchObject->SearchNearestInRadius(rPoint, Radius);
    } else { // Using trees
        // We search if geometrical objects are provided
        if (mpPointVector) {
            PointVector results(AllocationSize);
            DistanceVector results_distances(AllocationSize);
            const std::size_t number_of_results = mpSearchObject->SearchInRadius(rPoint, Radius, results.begin(), results_distances.begin(), AllocationSize);
            if (number_of_results > 0) {
                // Resize the results
                results_distances.resize(number_of_results);

                // Get the rank
                const int rank = GetRank();

                // Find the iterator pointing to the smallest value
                auto it_min = std::min_element(results_distances.begin(), results_distances.end());

                // Calculate the index
                IndexType index = std::distance(results_distances.begin(), it_min);

                // Set the result
                rResult = ResultType((results[index]->pGetObject()).get(), rank);
                rResult.SetDistance(results_distances[index]);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
void SearchWrapper<TSearchObject>::LocalSearchNearest(
    const PointType& rPoint,
    ResultType& rResult
    )
{
    // If we are using GeometricalObjectBins we can use the optimized search
    if constexpr (IsGeometricalObjectBins) {
        rResult = mpSearchObject->SearchNearest(rPoint);
    } else { // Using trees
        // We search if geometrical objects are provided
        if (mpPointVector) {
            // Get the rank
            const int rank = GetRank();

            // Search nearest
            double distance = 0.0;
            auto p_point = mpSearchObject->SearchNearestPoint(rPoint, distance);

            // Set the result
            rResult = ResultType((p_point->pGetObject()).get(), rank);
            rResult.SetDistance(distance);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
void SearchWrapper<TSearchObject>::LocalSearchIsInside(
    const PointType& rPoint,
    ResultType& rResult
    )
{
    // If we are using GeometricalObjectBins we can use the optimized search
    if constexpr (IsGeometricalObjectBins) {
        rResult = mpSearchObject->SearchIsInside(rPoint);
    } else { // Using trees
        KRATOS_ERROR << "SearchIsInside not compatible with Search trees" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
const Parameters SearchWrapper<TSearchObject>::GetDefaultParameters() const 
{
    const Parameters default_parameters = Parameters(R"(
    {
        "allocation_size" : 1000,
        "bucket_size"     : 4
    })" );
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

// GeometricalObjectsBins
template class SearchWrapper<GeometricalObjectsBins>;

// KDTree
template class SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>>;
template class SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>>;
template class SearchWrapper<Tree<KDTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>>;

// OCTree
template class SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>>;
template class SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>>;
template class SearchWrapper<Tree<OCTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>>;

// StaticBinsTree
template class SearchWrapper<Tree<Bins<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>;
template class SearchWrapper<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>;
template class SearchWrapper<Tree<Bins<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>;

// DynamicBins
template class SearchWrapper<BinsDynamic<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>;
template class SearchWrapper<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>;
template class SearchWrapper<BinsDynamic<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>;

}  // namespace Kratos.