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
#include "utilities/parallel_utilities.h"
#include "mpi/spatial_containers/specialized_spatial_search_mpi.h"

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
    const DataCommunicator& rDataCommunicator
    )
{
    
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchNodesOverPointNearestPoint (
    const NodesContainerType& rStructureNodes,
    const array_1d<double,3>& rPoint,
    NodeSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchElementsOverPointInRadius (
    const ElementsContainerType& rStructureElements,
    const array_1d<double,3>& rPoint,
    const double Radius,
    ElementSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchElementsOverPointNearestPoint (
    const ElementsContainerType& rStructureElements,
    const array_1d<double,3>& rPoint,
    ElementSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchConditionsOverPointInRadius (
    const ConditionsContainerType& rStructureConditions,
    const array_1d<double,3>& rPoint,
    const double Radius,
    ConditionSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
void SpecializedSpatialSearchMPI<TSearchBackend>::SearchConditionsOverPointNearestPoint (
    const ConditionsContainerType& rStructureConditions,
    const array_1d<double,3>& rPoint,
    ConditionSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template class SpecializedSpatialSearchMPI<SpatialContainer::KDTree>;
template class SpecializedSpatialSearchMPI<SpatialContainer::Octree>;
template class SpecializedSpatialSearchMPI<SpatialContainer::BinsStatic>;
template class SpecializedSpatialSearchMPI<SpatialContainer::BinsDynamic>;

} // namespace Kratos