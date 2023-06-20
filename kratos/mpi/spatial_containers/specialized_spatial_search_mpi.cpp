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
#include "utilities/parallel_utilities.h"
#include "mpi/spatial_containers/specialized_spatial_search_mpi.h"
#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos
{

template<SpatialContainer TSearchBackend>
typename SpecializedSpatialSearchMPI<TSearchBackend>::ElementSpatialSearchResultContainerMapType SpecializedSpatialSearchMPI<TSearchBackend>::SearchElementsInRadiusExclusive (
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    // The points vector
    const auto it_elem_begin = rInputElements.begin();
    std::vector<Point> points(rInputElements.size());
    IndexPartition<std::size_t>(rInputElements.size()).for_each([&](std::size_t i) {
        points[i] = Point((it_elem_begin + i)->GetGeometry().Center());
    });

    return SearchElementsOverPointsInRadius(rStructureElements, points.begin(), points.end(), rRadius, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
typename SpecializedSpatialSearchMPI<TSearchBackend>::NodeSpatialSearchResultContainerMapType SpecializedSpatialSearchMPI<TSearchBackend>::SearchNodesInRadiusExclusive (
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return SearchNodesOverPointsInRadius(rStructureNodes, rInputNodes.begin(), rInputNodes.end(), rRadius, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
typename SpecializedSpatialSearchMPI<TSearchBackend>::ConditionSpatialSearchResultContainerMapType SpecializedSpatialSearchMPI<TSearchBackend>::SearchConditionsInRadiusExclusive (
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    // The points vector
    const auto it_cond_begin = rInputConditions.begin();
    std::vector<Point> points(rInputConditions.size());
    IndexPartition<std::size_t>(rInputConditions.size()).for_each([&](std::size_t i) {
        points[i] = Point((it_cond_begin + i)->GetGeometry().Center());
    });

    return SearchConditionsOverPointsInRadius(rStructureConditions, points.begin(), points.end(), rRadius, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchNodesOverPointInRadius (
    const NodesContainerType& rStructureNodes,
    const array_1d<double,3>& rPoint,
    const double Radius,
    NodeSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{    
    // Initialize the BB is required
    if (!mBoundingBoxesInitialized) {
        InitializeLocalBoundingBox(rStructureNodes);
    }

    // Check if the point is inside the set
    if (SearchUtilities::PointIsInsideBoundingBox(mLocalBoundingBox, rPoint, Radius)) {
        // Call local search
        BaseType::SearchNodesOverPointInRadius(rStructureNodes, rPoint, Radius, rResults, rDataCommunicator, false);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchNodesOverPointNearestPoint (
    const NodesContainerType& rStructureNodes,
    const array_1d<double,3>& rPoint,
    NodeSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{    
    // Initialize the BB is required
    if (!mBoundingBoxesInitialized) {
        InitializeLocalBoundingBox(rStructureNodes);
    }

    // Compute max radius
    const array_1d<double, 3> box_size = mLocalBoundingBox.GetMaxPoint() - mLocalBoundingBox.GetMinPoint();
    const double max_radius= *std::max_element(box_size.begin(), box_size.end());

    // Get the rank
    const int current_rank = rDataCommunicator.Rank();

    // Check if the point is inside the set
    NodeSpatialSearchResultContainerType local_result;
    if (SearchUtilities::PointIsInsideBoundingBox(mLocalBoundingBox, rPoint, max_radius)) {
        // Call local search
        BaseType::SearchNodesOverPointNearestPoint(rStructureNodes, rPoint, local_result, rDataCommunicator, false);
    }

    /* Now sync results between partitions */

    // Get the distance
    const double local_distance = local_result.NumberOfLocalResults() > 0 ? local_result.GetLocalDistances().begin()->second : std::numeric_limits<double>::max();

    // Find the minimum value and the rank that holds it
    struct {
        double value;
        int rank;
    } local_min, global_min;

    local_min.value = local_distance;
    local_min.rank = current_rank;
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPIDataCommunicator::GetMPICommunicator(rDataCommunicator));

    // Get the solution from the computed_rank
    if (global_min.rank == current_rank) {
        // Add the local search
        rResults.AddResult(*(local_result.GetLocalPointers().begin().base()));
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchElementsOverPointInRadius (
    const ElementsContainerType& rStructureElements,
    const array_1d<double,3>& rPoint,
    const double Radius,
    ElementSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{    
    // Initialize the BB is required
    if (!mBoundingBoxesInitialized) {
        InitializeLocalBoundingBox(rStructureElements.begin(), rStructureElements.end());
    }

    // Check if the point is inside the set
    if (SearchUtilities::PointIsInsideBoundingBox(mLocalBoundingBox, rPoint, Radius)) {
        // Call local search
        BaseType::SearchElementsOverPointInRadius(rStructureElements, rPoint, Radius, rResults, rDataCommunicator, false);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchElementsOverPointNearestPoint (
    const ElementsContainerType& rStructureElements,
    const array_1d<double,3>& rPoint,
    ElementSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    // Initialize the BB is required
    if (!mBoundingBoxesInitialized) {
        InitializeLocalBoundingBox(rStructureElements.begin(), rStructureElements.end());
    }

    // Compute max radius
    const array_1d<double, 3> box_size = mLocalBoundingBox.GetMaxPoint() - mLocalBoundingBox.GetMinPoint();
    const double max_radius= *std::max_element(box_size.begin(), box_size.end());

    // Get the rank
    const int current_rank = rDataCommunicator.Rank();

    // Check if the point is inside the set
    ElementSpatialSearchResultContainerType local_result;
    if (SearchUtilities::PointIsInsideBoundingBox(mLocalBoundingBox, rPoint, max_radius)) {
        // Call local search
        BaseType::SearchElementsOverPointNearestPoint(rStructureElements, rPoint, local_result, rDataCommunicator, false);
    }

    /* Now sync results between partitions */

    // Get the distance
    const double local_distance = local_result.NumberOfLocalResults() > 0 ? local_result.GetLocalDistances().begin()->second : std::numeric_limits<double>::max();

    // Find the minimum value and the rank that holds it
    struct {
        double value;
        int rank;
    } local_min, global_min;

    local_min.value = local_distance;
    local_min.rank = current_rank;
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPIDataCommunicator::GetMPICommunicator(rDataCommunicator));

    // Get the solution from the computed_rank
    if (global_min.rank == current_rank) {
        // Add the local search
        rResults.AddResult(*(local_result.GetLocalPointers().begin().base()));
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchConditionsOverPointInRadius (
    const ConditionsContainerType& rStructureConditions,
    const array_1d<double,3>& rPoint,
    const double Radius,
    ConditionSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    // Initialize the BB is required
    if (!mBoundingBoxesInitialized) {
        InitializeLocalBoundingBox(rStructureConditions.begin(), rStructureConditions.end());
    }

    // Check if the point is inside the set
    if (SearchUtilities::PointIsInsideBoundingBox(mLocalBoundingBox, rPoint, Radius)) {
        // Call local search
        BaseType::SearchConditionsOverPointInRadius(rStructureConditions, rPoint, Radius, rResults, rDataCommunicator, false);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchConditionsOverPointNearestPoint (
    const ConditionsContainerType& rStructureConditions,
    const array_1d<double,3>& rPoint,
    ConditionSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    // Initialize the BB is required
    if (!mBoundingBoxesInitialized) {
        InitializeLocalBoundingBox(rStructureConditions.begin(), rStructureConditions.end());
    }

    // Compute max radius
    const array_1d<double, 3> box_size = mLocalBoundingBox.GetMaxPoint() - mLocalBoundingBox.GetMinPoint();
    const double max_radius= *std::max_element(box_size.begin(), box_size.end());

    // Get the rank
    const int current_rank = rDataCommunicator.Rank();

    // Check if the point is inside the set
    ConditionSpatialSearchResultContainerType local_result;
    if (SearchUtilities::PointIsInsideBoundingBox(mLocalBoundingBox, rPoint, max_radius)) {
        // Call local search
        BaseType::SearchConditionsOverPointNearestPoint(rStructureConditions, rPoint, local_result, rDataCommunicator, false);
    }

    /* Now sync results between partitions */

    // Get the distance
    const double local_distance = local_result.NumberOfLocalResults() > 0 ? local_result.GetLocalDistances().begin()->second : std::numeric_limits<double>::max();

    // Find the minimum value and the rank that holds it
    struct {
        double value;
        int rank;
    } local_min, global_min;

    local_min.value = local_distance;
    local_min.rank = current_rank;
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPIDataCommunicator::GetMPICommunicator(rDataCommunicator));

    // Get the solution from the computed_rank
    if (global_min.rank == current_rank) {
        // Add the local search
        rResults.AddResult(*(local_result.GetLocalPointers().begin().base()));
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::InitializeLocalBoundingBox(const NodesContainerType& rStructureNodes)
{
    const std::size_t number_of_nodes = rStructureNodes.size();
    if (number_of_nodes > 0) {
        mLocalBoundingBox.Set(rStructureNodes.begin(), rStructureNodes.end());
        mLocalBoundingBox.Extend(Tolerance);
    }
    mBoundingBoxesInitialized = true;
}

/***********************************************************************************/
/***********************************************************************************/

template class SpecializedSpatialSearchMPI<SpatialContainer::KDTree>;
template class SpecializedSpatialSearchMPI<SpatialContainer::Octree>;
template class SpecializedSpatialSearchMPI<SpatialContainer::BinsStatic>;
template class SpecializedSpatialSearchMPI<SpatialContainer::BinsDynamic>;

} // namespace Kratos