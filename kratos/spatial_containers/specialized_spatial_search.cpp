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
#include "utilities/parallel_utilities.h"
#include "spatial_containers/specialized_spatial_search.h"
#include "spatial_containers/spatial_containers.h"

namespace Kratos
{

template<>
void PointObject<Node<3>>::UpdatePoint()
{
    noalias(this->Coordinates()) = mpObject->Coordinates();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void PointObject<Condition>::UpdatePoint()
{
    noalias(this->Coordinates()) = mpObject->GetGeometry().Center().Coordinates();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void PointObject<Element>::UpdatePoint()
{
    noalias(this->Coordinates()) = mpObject->GetGeometry().Center().Coordinates();
}

template class PointObject<Node<3>>;
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
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;

    // Retrieving parameters
    const int allocation_size = mParameters["allocation_size"].GetInt();
    const int bucket_size = mParameters["bucket_size"].GetInt();

    // Defining the PointVector
    PointVector points;
    const std::size_t structure_size = rStructureElements.size();
    points.reserve(structure_size);
    const auto it_begin = rStructureElements.begin();
    for (std::size_t i = 0; i < structure_size; ++i) {
        auto it_elem = it_begin + i;
        points.push_back(PointTypePointer(new PointType(*(it_elem.base()))));
    }

    // Resizing the results
    const std::size_t input_size = rInputElements.size();
    if (rResults.size() != input_size) {
        rResults.resize(input_size);
    }
    if (rResultsDistance.size() != input_size) {
        rResultsDistance.resize(input_size);
    }

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using KDTree = Tree<KDTreePartition<BucketType>>;

        // Creating the tree
        KDTree kd_tree(points.begin(), points.end(), bucket_size);

        // Performing search
        IndexPartition<std::size_t>(input_size).for_each([&](std::size_t i) {
            auto it_elem = rInputElements.begin() + i;
            PointType aux_point(*(it_elem.base()));
            PointVector results;
            DistanceVector results_distances;
            const std::size_t number_of_results = kd_tree.SearchInRadius(aux_point, rRadius[i], results.begin(), results_distances.begin(), allocation_size);
            if (number_of_results > 0) {
                auto& r_results = rResults[i];
                auto& r_results_distance = rResultsDistance[i];
                r_results.reserve(number_of_results);
                r_results_distance.reserve(number_of_results);
                for (std::size_t j = 0; j < number_of_results; ++j) {
                    PointType::Pointer p_point = results[j];
                    r_results.push_back(p_point->pGetObject());
                    r_results_distance.push_back(results_distances[j]);
                }
            }
        });
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using Octree = Tree<OCTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBins = Bins<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsStaticObjects) {
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamicObjects) {
    } else {
        KRATOS_ERROR << "Unknown search backend" << std::endl;
    }

}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchElementsInRadiusInclusive(
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Element>;
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;

    // Retrieving parameters
    const int allocation_size = mParameters["allocation_size"].GetInt();
    const int bucket_size = mParameters["bucket_size"].GetInt();

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using KDTree = Tree<KDTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using Octree = Tree<OCTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBins = Bins<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsStaticObjects) {
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamicObjects) {
    } else {
        KRATOS_ERROR << "Unknown search backend" << std::endl;
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
void SpecializedSpatialSearch<TSearchBackend>::SearchElementsInRadiusInclusive(
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchElementsInRadiusInclusive(rStructureElements, rInputElements, rRadius, rResults, distances);
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
    using PointType = PointObject<Node<3>>;
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;

    // Retrieving parameters
    const int allocation_size = mParameters["allocation_size"].GetInt();
    const int bucket_size = mParameters["bucket_size"].GetInt();

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using KDTree = Tree<KDTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using Octree = Tree<OCTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBins = Bins<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsStaticObjects) {
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamicObjects) {
    } else {
        KRATOS_ERROR << "Unknown search backend" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchNodesInRadiusInclusive(
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Node<3>>;
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;

    // Retrieving parameters
    const int allocation_size = mParameters["allocation_size"].GetInt();
    const int bucket_size = mParameters["bucket_size"].GetInt();

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using KDTree = Tree<KDTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using Octree = Tree<OCTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBins = Bins<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsStaticObjects) {
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamicObjects) {
    } else {
        KRATOS_ERROR << "Unknown search backend" << std::endl;
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
void SpecializedSpatialSearch<TSearchBackend>::SearchNodesInRadiusInclusive(
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchNodesInRadiusInclusive(rStructureNodes, rInputNodes, rRadius, rResults, distances);
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
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;

    // Retrieving parameters
    const int allocation_size = mParameters["allocation_size"].GetInt();
    const int bucket_size = mParameters["bucket_size"].GetInt();

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using KDTree = Tree<KDTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using Octree = Tree<OCTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBins = Bins<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsStaticObjects) {
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamicObjects) {
    } else {
        KRATOS_ERROR << "Unknown search backend" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearch<TSearchBackend>::SearchConditionsInRadiusInclusive(
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    // Defining the point type for the search
    using PointType = PointObject<Condition>;
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;

    // Retrieving parameters
    const int allocation_size = mParameters["allocation_size"].GetInt();
    const int bucket_size = mParameters["bucket_size"].GetInt();

    // Defining the search structure
    if constexpr (TSearchBackend == SpatialContainer::KDTree) {
        /// KDtree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using KDTree = Tree<KDTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::Octree) {
        /// Octree definitions
        using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using Octree = Tree<OCTreePartition<BucketType>>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsStatic) {
        /// StaticBins definitions
        using StaticBins = Bins<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamic) {
        /// BinsDynamic definitions
        using DynamicBins = BinsDynamic<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsStaticObjects) {
    // } else if constexpr (TSearchBackend == SpatialContainer::BinsDynamicObjects) {
    } else {
        KRATOS_ERROR << "Unknown search backend" << std::endl;
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
void SpecializedSpatialSearch<TSearchBackend>::SearchConditionsInRadiusInclusive(
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchConditionsInRadiusInclusive(rStructureConditions, rInputConditions, rRadius, rResults, distances);
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
// template class SpecializedSpatialSearch<SpatialContainer::BinsStaticObjects>;
// template class SpecializedSpatialSearch<SpatialContainer::BinsDynamicObjects>;
} // namespace Kratos.