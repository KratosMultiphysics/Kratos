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
    std::vector<Point> points(rInputElements.size());
    for (std::size_t i = 0; i < rInputElements.size(); ++i) {
        points[i] = Point((rInputElements.begin() + i)->GetGeometry().Center());
    }

    return SearchElementsOverPointsInRadius(rStructureElements, points.begin(), points.end(), rRadius, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
typename SpecializedSpatialSearchMPI<TSearchBackend>::NodeSpatialSearchResultContainerMapType SpecializedSpatialSearchMPI<TSearchBackend>::SearchElementsInRadiusInclusive (
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    NodeSpatialSearchResultContainerMapType results;

    return results;
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
typename SpecializedSpatialSearchMPI<TSearchBackend>::NodeSpatialSearchResultContainerMapType SpecializedSpatialSearchMPI<TSearchBackend>::SearchNodesInRadiusInclusive (
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    NodeSpatialSearchResultContainerMapType results;

    return results;
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
    std::vector<Point> points(rInputConditions.size());
    for (std::size_t i = 0; i < rInputConditions.size(); ++i) {
        points[i] = Point((rInputConditions.begin() + i)->GetGeometry().Center());
    }

    return SearchConditionsOverPointsInRadius(rStructureConditions, points.begin(), points.end(), rRadius, rDataCommunicator);
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearchBackend>
typename SpecializedSpatialSearchMPI<TSearchBackend>::NodeSpatialSearchResultContainerMapType SpecializedSpatialSearchMPI<TSearchBackend>::SearchConditionsInRadiusInclusive (
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    NodeSpatialSearchResultContainerMapType results;

    return results;
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