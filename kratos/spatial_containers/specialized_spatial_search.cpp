//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |  (   | |  (   |\__ `
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
#include "spatial_containers/specialized_spatial_search.h"
#include "spatial_containers/spatial_containers.h"

namespace Kratos
{

template class PointObject<Node>;
template class PointObject<GeometricalObject>;
template class PointObject<Condition>;
template class PointObject<Element>;

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchElementsInRadiusExclusive(
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Element>;
    using PointVector = std::vector<PointType::Pointer>;

    // Defining the PointVector
    PointVector points = PrepareSearch(rStructureElements, rInputElements, rResults, rResultsDistance);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        ParallelSearchInRadius(rInputElements, rRadius, kd_tree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), mBucketSize);

        // Performing search
        ParallelSearchInRadius(rInputElements, rRadius, octree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        ParallelSearchInRadius(rInputElements, rRadius, static_bins_tree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        ParallelSearchInRadius(rInputElements, rRadius, dynamic_bins, rResults, rResultsDistance);
    } else {
        static_assert(TSearchBackend == SpatialContainer::KDTree || TSearchBackend == SpatialContainer::Octree || TSearchBackend == SpatialContainer::BinsStatic || TSearchBackend == SpatialContainer::BinsDynamic, "Unknown search backend");
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchElementsInRadiusExclusive(
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchElementsInRadiusExclusive(rStructureElements, rInputElements, rRadius, rResults, distances);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
typename SpecializedSpatialSearch<TSearchBackend>::ElementSpatialSearchResultContainerMapType SpecializedSpatialSearch<TSearchBackend>::SearchElementsInRadiusExclusive (
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return BaseType::SearchElementsInRadiusExclusive(rStructureElements, rInputElements, rRadius, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchNodesInRadiusExclusive(
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Node>;
    using PointVector = std::vector<PointType::Pointer>;

    // Defining the PointVector
    PointVector points = PrepareSearch(rStructureNodes, rInputNodes, rResults, rResultsDistance);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        ParallelSearchInRadius(rInputNodes, rRadius, kd_tree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), mBucketSize);

        // Performing search
        ParallelSearchInRadius(rInputNodes, rRadius, octree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        ParallelSearchInRadius(rInputNodes, rRadius, static_bins_tree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        ParallelSearchInRadius(rInputNodes, rRadius, dynamic_bins, rResults, rResultsDistance);
    } else {
        static_assert(TSearchBackend == SpatialContainer::KDTree || TSearchBackend == SpatialContainer::Octree || TSearchBackend == SpatialContainer::BinsStatic || TSearchBackend == SpatialContainer::BinsDynamic, "Unknown search backend");
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchNodesInRadiusExclusive(
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, rRadius, rResults, distances);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
typename SpecializedSpatialSearch<TSearchBackend>::NodeSpatialSearchResultContainerMapType SpecializedSpatialSearch<TSearchBackend>::SearchNodesInRadiusExclusive (
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return BaseType::SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, rRadius, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchConditionsInRadiusExclusive(
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Condition>;
    using PointVector = std::vector<PointType::Pointer>;

    // Defining the PointVector
    PointVector points = PrepareSearch(rStructureConditions, rInputConditions, rResults, rResultsDistance);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        ParallelSearchInRadius(rInputConditions, rRadius, kd_tree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), mBucketSize);

        // Performing search
        ParallelSearchInRadius(rInputConditions, rRadius, octree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        ParallelSearchInRadius(rInputConditions, rRadius, static_bins_tree, rResults, rResultsDistance);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        ParallelSearchInRadius(rInputConditions, rRadius, dynamic_bins, rResults, rResultsDistance);
    } else {
        static_assert(TSearchBackend == SpatialContainer::KDTree || TSearchBackend == SpatialContainer::Octree || TSearchBackend == SpatialContainer::BinsStatic || TSearchBackend == SpatialContainer::BinsDynamic, "Unknown search backend");
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchConditionsInRadiusExclusive(
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchConditionsInRadiusExclusive(rStructureConditions, rInputConditions, rRadius, rResults, distances);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
typename SpecializedSpatialSearch<TSearchBackend>::ConditionSpatialSearchResultContainerMapType SpecializedSpatialSearch<TSearchBackend>::SearchConditionsInRadiusExclusive (
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return BaseType::SearchConditionsInRadiusExclusive(rStructureConditions, rInputConditions, rRadius, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchNodesOverPointInRadius (
    const NodesContainerType& rStructureNodes,
    const array_1d<double,3>& rPoint,
    const double Radius,
    NodeSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Node>;
    using PointVector = std::vector<PointType::Pointer>;

    // Defining the PointVector
    PointVector points = PrepareSearchPoints(rStructureNodes);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchInRadius(rStructureNodes, rPoint, Radius, kd_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchInRadius(rStructureNodes, rPoint, Radius, octree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchInRadius(rStructureNodes, rPoint, Radius, static_bins_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        SearchInRadius(rStructureNodes, rPoint, Radius, dynamic_bins, rResults, rDataCommunicator);
    } else {
        static_assert(TSearchBackend == SpatialContainer::KDTree || TSearchBackend == SpatialContainer::Octree || TSearchBackend == SpatialContainer::BinsStatic || TSearchBackend == SpatialContainer::BinsDynamic, "Unknown search backend");
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchNodesOverPointNearestPoint (
    const NodesContainerType& rStructureNodes,
    const array_1d<double,3>& rPoint,
    NodeSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Node>;
    using PointVector = std::vector<PointType::Pointer>;

    // Defining the PointVector
    PointVector points = PrepareSearchPoints(rStructureNodes);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchNearestPoint(rStructureNodes, rPoint, kd_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchNearestPoint(rStructureNodes, rPoint, octree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchNearestPoint(rStructureNodes, rPoint, static_bins_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        SearchNearestPoint(rStructureNodes, rPoint, dynamic_bins, rResults, rDataCommunicator);
    } else {
        static_assert(TSearchBackend == SpatialContainer::KDTree || TSearchBackend == SpatialContainer::Octree || TSearchBackend == SpatialContainer::BinsStatic || TSearchBackend == SpatialContainer::BinsDynamic, "Unknown search backend");
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchElementsOverPointInRadius (
    const ElementsContainerType& rStructureElements,
    const array_1d<double,3>& rPoint,
    const double Radius,
    ElementSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Element>;
    using PointVector = std::vector<PointType::Pointer>;

    // Defining the PointVector
    PointVector points = PrepareSearchPoints(rStructureElements);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchInRadius(rStructureElements, rPoint, Radius, kd_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchInRadius(rStructureElements, rPoint, Radius, octree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchInRadius(rStructureElements, rPoint, Radius, static_bins_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        SearchInRadius(rStructureElements, rPoint, Radius, dynamic_bins, rResults, rDataCommunicator);
    } else {
        static_assert(TSearchBackend == SpatialContainer::KDTree || TSearchBackend == SpatialContainer::Octree || TSearchBackend == SpatialContainer::BinsStatic || TSearchBackend == SpatialContainer::BinsDynamic, "Unknown search backend");
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchElementsOverPointNearestPoint (
    const ElementsContainerType& rStructureElements,
    const array_1d<double,3>& rPoint,
    ElementSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Element>;
    using PointVector = std::vector<PointType::Pointer>;

    // Defining the PointVector
    PointVector points = PrepareSearchPoints(rStructureElements);
    
    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchNearestPoint(rStructureElements, rPoint, kd_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchNearestPoint(rStructureElements, rPoint, octree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchNearestPoint(rStructureElements, rPoint, static_bins_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        SearchNearestPoint(rStructureElements, rPoint, dynamic_bins, rResults, rDataCommunicator);
    } else {
        static_assert(TSearchBackend == SpatialContainer::KDTree || TSearchBackend == SpatialContainer::Octree || TSearchBackend == SpatialContainer::BinsStatic || TSearchBackend == SpatialContainer::BinsDynamic, "Unknown search backend");
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchConditionsOverPointInRadius (
    const ConditionsContainerType& rStructureConditions,
    const array_1d<double,3>& rPoint,
    const double Radius,
    ConditionSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Condition>;
    using PointVector = std::vector<PointType::Pointer>;

    // Defining the PointVector
    PointVector points = PrepareSearchPoints(rStructureConditions);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchInRadius(rStructureConditions, rPoint, Radius, kd_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchInRadius(rStructureConditions, rPoint, Radius, octree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchInRadius(rStructureConditions, rPoint, Radius, static_bins_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        SearchInRadius(rStructureConditions, rPoint, Radius, dynamic_bins, rResults, rDataCommunicator);
    } else {
        static_assert(TSearchBackend == SpatialContainer::KDTree || TSearchBackend == SpatialContainer::Octree || TSearchBackend == SpatialContainer::BinsStatic || TSearchBackend == SpatialContainer::BinsDynamic, "Unknown search backend");
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(rDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchConditionsOverPointNearestPoint (
    const ConditionsContainerType& rStructureConditions,
    const array_1d<double,3>& rPoint,
    ConditionSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Condition>;
    using PointVector = std::vector<PointType::Pointer>;

    // Defining the PointVector
    PointVector points = PrepareSearchPoints(rStructureConditions);

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchNearestPoint(rStructureConditions, rPoint, kd_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchNearestPoint(rStructureConditions, rPoint, octree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), mBucketSize);

        // Performing search
        SearchNearestPoint(rStructureConditions, rPoint, static_bins_tree, rResults, rDataCommunicator);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        SearchNearestPoint(rStructureConditions, rPoint, dynamic_bins, rResults, rDataCommunicator);
    } else {
        static_assert(TSearchBackend == SpatialContainer::KDTree || TSearchBackend == SpatialContainer::Octree || TSearchBackend == SpatialContainer::BinsStatic || TSearchBackend == SpatialContainer::BinsDynamic, "Unknown search backend");
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
Parameters SpecializedSpatialSearch<TSearchBackend>::GetDefaultParameters() const
{
    Parameters default_parameters = Parameters(R"(
    {
        "allocation_size" : 1000,
        "bucket_size"     : 4
    })" );

    return default_parameters;
}

template class SpecializedSpatialSearch<SpatialContainer::KDTree>;
template class SpecializedSpatialSearch<SpatialContainer::Octree>;
template class SpecializedSpatialSearch<SpatialContainer::BinsStatic>;
template class SpecializedSpatialSearch<SpatialContainer::BinsDynamic>;
} // namespace Kratos.