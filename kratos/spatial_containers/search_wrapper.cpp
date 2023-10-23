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

template<class TSearchObject, class TObjectType>
BoundingBox<Point> SearchWrapper<TSearchObject, TObjectType>::GetBoundingBox() const
{
    // Generate BB
    BoundingBox<Point> bb;
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

    return bb;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::SearchInRadius(
    const Point& rPoint,
    const double Radius,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
#ifdef KRATOS_USING_MPI
    // Call distributed version
    if (mrDataCommunicator.IsDistributed()) {
        DistributedSearchInRadius(rPoint, Radius, rResults, SyncronizeResults);
    } else { // Call serial version
        SerialSearchInRadius(rPoint, Radius, rResults, SyncronizeResults);
    }
#else // Calling just serial method
    SerialSearchInRadius(rPoint, Radius, rResults, SyncronizeResults);
#endif
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::SearchNearestInRadius(
    const Point& rPoint,
    const double Radius,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
#ifdef KRATOS_USING_MPI
    // Call distributed version
    if (mrDataCommunicator.IsDistributed()) {
        DistributedSearchNearestInRadius(rPoint, Radius, rResults, SyncronizeResults);
    } else { // Call serial version
        SerialSearchNearestInRadius(rPoint, Radius, rResults, SyncronizeResults);
    }
#else // Calling just serial method
    SerialSearchNearestInRadius(rPoint, Radius, rResults, SyncronizeResults);
#endif
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::SearchNearest(
    const Point& rPoint,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
#ifdef KRATOS_USING_MPI
    // Call distributed version
    if (mrDataCommunicator.IsDistributed()) {
        DistributedSearchNearest(rPoint, rResults, SyncronizeResults);
    } else { // Call serial version
        SerialSearchNearest(rPoint, rResults, SyncronizeResults);
    }
#else // Calling just serial method
    SerialSearchNearest(rPoint, rResults, SyncronizeResults);
#endif
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::SearchIsInside(
    const Point& rPoint,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
#ifdef KRATOS_USING_MPI
    // Call distributed version
    if (mrDataCommunicator.IsDistributed()) {
        DistributedSearchIsInside(rPoint, rResults, SyncronizeResults);
    } else { // Call serial version
        SerialSearchIsInside(rPoint, rResults, SyncronizeResults);
    }
#else // Calling just serial method
    SerialSearchIsInside(rPoint, rResults, SyncronizeResults);
#endif
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::SerialSearchInRadius(
    const Point& rPoint,
    const double Radius,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
    // Search
    std::vector<ResultType> results;
    mpSearchObject->SearchInRadius(rPoint, Radius, results);
    for (auto& r_result : results) {
        rResults.AddResult(r_result);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::SerialSearchNearestInRadius(
    const Point& rPoint,
    const double Radius,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
    // Search
    auto result = mpSearchObject->SearchNearestInRadius(rPoint, Radius);
    rResults.AddResult(result);

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::SerialSearchNearest(
    const Point& rPoint,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
    // Search
    auto result = mpSearchObject->SearchNearest(rPoint);
    rResults.AddResult(result);

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::SerialSearchIsInside(
    const Point& rPoint,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
    // Search
    auto result = mpSearchObject->SearchIsInside(rPoint);
    rResults.AddResult(result);

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/
#ifdef KRATOS_USING_MPI

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::DistributedSearchInRadius(
    const Point& rPoint,
    const double Radius,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
    // Check if the point is inside the set
    if (SearchUtilities::PointIsInsideBoundingBox(mpSearchObject->GetBoundingBox(), rPoint, Radius)) {
        // Call local search
        SerialSearchInRadius(rPoint, Radius, rResults, false);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::DistributedSearchNearestInRadius(
    const Point& rPoint,
    const double Radius,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
    // Result to return
    ResultType local_result;

    // Get the rank
    const int current_rank = GetRank();

    // Check if the point is inside the set
    if (SearchUtilities::PointIsInsideBoundingBox(mpSearchObject->GetBoundingBox(), rPoint, Radius)) {
        // Call local search
        local_result = mpSearchObject->SearchNearestInRadius(rPoint, Radius);
    }

    /* Now sync results between partitions */

    // Get the distance
    const double local_distance = (local_result.IsObjectFound() && local_result.IsDistanceCalculated()) ? local_result.GetDistance() : std::numeric_limits<double>::max();

    // Find the minimum value and the rank that holds it
    const auto global_min = mrDataCommunicator.MinLocAll(local_distance);

    // Get the solution from the computed_rank
    if (global_min.second == current_rank) {
        // Add the local search
        rResults.AddResult(local_result);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::DistributedSearchNearest(
    const Point& rPoint,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
    ResultType current_result;
    const auto bb = GetBoundingBox();
    const array_1d<double, 3> box_size = bb.GetMaxPoint() - bb.GetMinPoint();
    const double max_radius= *std::max_element(box_size.begin(), box_size.end());
    this->SearchNearestInRadius(rPoint, max_radius, rResults, SyncronizeResults);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::DistributedSearchIsInside(
    const Point& rPoint,
    ResultContainerType& rResults,
    const bool SyncronizeResults
    )
{
    // Result to return
    ResultType local_result;

    // Get the rank
    const int current_rank = GetRank();

    // Check if the point is inside the set
    int computed_rank = std::numeric_limits<int>::max();
    if (SearchUtilities::PointIsInsideBoundingBox(mpSearchObject->GetBoundingBox(), rPoint)) {
        // Call local search
        local_result = mpSearchObject->SearchIsInside(rPoint);

        // Set current rank
        computed_rank = current_rank;
    }

    // Now sync results between partitions
    //KRATOS_WARNING("SearchWrapper.SearchIsInside") << "This assumes that first one of the viable results is the closest one. The algorithm  sets distance to 0, so it is ambiguous which is the closest one." << std::endl;

    // Min rank of all ranks
    computed_rank = mrDataCommunicator.MinAll(computed_rank);

    // Get the solution from the computed_rank
    if (computed_rank == current_rank) {
        // Add the local search
        rResults.AddResult(local_result);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

#endif

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
int SearchWrapper<TSearchObject, TObjectType>::GetRank() const
{
    return mrDataCommunicator.Rank();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
int SearchWrapper<TSearchObject, TObjectType>::GetWorldSize() const
{
    return mrDataCommunicator.Size();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject, class TObjectType>
void SearchWrapper<TSearchObject, TObjectType>::InitializeGlobalBoundingBoxes()
{
    // Just executed in MPI
    if (mrDataCommunicator.IsDistributed()) {
        // We get the world size
        const int world_size = GetWorldSize();

        // Set up the global bounding boxes
        if (static_cast<int>(mGlobalBoundingBoxes.size()) != 6*world_size) {
            mGlobalBoundingBoxes.resize(6*world_size);
        }

        // Set up the local bounding boxes
        std::vector<double> local_bounding_box(6);
        const auto& r_bb = mpSearchObject->GetBoundingBox();
        const auto& r_max = r_bb.GetMaxPoint();
        const auto& r_min = r_bb.GetMinPoint();
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

template<class TSearchObject, class TObjectType>
std::vector<int> SearchWrapper<TSearchObject, TObjectType>::RansksPointIsInsideBoundingBox(const array_1d<double, 3>& rCoords)
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

template<class TSearchObject, class TObjectType>
std::vector<int> SearchWrapper<TSearchObject, TObjectType>::RansksPointIsInsideBoundingBoxWithTolerance(
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

template class SearchWrapper<GeometricalObjectsBins, GeometricalObject>;

}  // namespace Kratos.