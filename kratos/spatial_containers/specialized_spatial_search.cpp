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
#include "utilities/search_utilities.h"
#include "spatial_containers/specialized_spatial_search.h"
#include "spatial_containers/spatial_containers.h"

namespace Kratos
{

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
    PointVector points = SearchUtilities::PrepareSearch(rStructureElements, rInputElements, rResults, rResultsDistance);

    // Retrieving parameters
    const int allocation_size = mParameters["allocation_size"].GetInt();

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        // Retrieving parameters
        const int bucket_size = mParameters["bucket_size"].GetInt();

        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        SearchUtilities::ParallelSearch(rInputElements, rRadius, kd_tree, rResults, rResultsDistance, allocation_size);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        // Retrieving parameters
        const int bucket_size = mParameters["bucket_size"].GetInt();

        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), bucket_size);

        // Performing search
        SearchUtilities::ParallelSearch(rInputElements, rRadius, octree, rResults, rResultsDistance, allocation_size);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        // Retrieving parameters
        const int bucket_size = mParameters["bucket_size"].GetInt();

        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        SearchUtilities::ParallelSearch(rInputElements, rRadius, static_bins_tree, rResults, rResultsDistance, allocation_size);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        SearchUtilities::ParallelSearch(rInputElements, rRadius, dynamic_bins, rResults, rResultsDistance, allocation_size);
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
    PointVector points = SearchUtilities::PrepareSearch(rStructureNodes, rInputNodes, rResults, rResultsDistance);

    // Retrieving parameters
    const int allocation_size = mParameters["allocation_size"].GetInt();

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        // Retrieving parameters
        const int bucket_size = mParameters["bucket_size"].GetInt();

        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        SearchUtilities::ParallelSearch(rInputNodes, rRadius, kd_tree, rResults, rResultsDistance, allocation_size);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        // Retrieving parameters
        const int bucket_size = mParameters["bucket_size"].GetInt();

        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), bucket_size);

        // Performing search
        SearchUtilities::ParallelSearch(rInputNodes, rRadius, octree, rResults, rResultsDistance, allocation_size);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        // Retrieving parameters
        const int bucket_size = mParameters["bucket_size"].GetInt();

        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        SearchUtilities::ParallelSearch(rInputNodes, rRadius, static_bins_tree, rResults, rResultsDistance, allocation_size);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        SearchUtilities::ParallelSearch(rInputNodes, rRadius, dynamic_bins, rResults, rResultsDistance, allocation_size);
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
    PointVector points = SearchUtilities::PrepareSearch(rStructureConditions, rInputConditions, rResults, rResultsDistance);

    // Retrieving parameters
    const int allocation_size = mParameters["allocation_size"].GetInt();

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        // Retrieving parameters
        const int bucket_size = mParameters["bucket_size"].GetInt();

        /// KDtree definitions
        using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        SearchUtilities::ParallelSearch(rInputConditions, rRadius, kd_tree, rResults, rResultsDistance, allocation_size);
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        // Retrieving parameters
        const int bucket_size = mParameters["bucket_size"].GetInt();

        /// Octree definitions
        using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

        // Creating the tree
        Octree octree(points.begin(), points.end(), bucket_size);

        // Performing search
        SearchUtilities::ParallelSearch(rInputConditions, rRadius, octree, rResults, rResultsDistance, allocation_size);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        // Retrieving parameters
        const int bucket_size = mParameters["bucket_size"].GetInt();

        /// StaticBins definitions
        using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

        // Creating the tree
        StaticBinsTree static_bins_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        SearchUtilities::ParallelSearch(rInputConditions, rRadius, static_bins_tree, rResults, rResultsDistance, allocation_size);
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

        // Creating the bins
        DynamicBins dynamic_bins(points.begin(), points.end());

        // Performing search
        SearchUtilities::ParallelSearch(rInputConditions, rRadius, dynamic_bins, rResults, rResultsDistance, allocation_size);
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