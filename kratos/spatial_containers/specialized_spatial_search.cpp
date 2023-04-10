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

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchElementsInRadiusExclusive(
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchElementsInRadiusInclusive(
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchElementsInRadiusExclusive(
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

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchElementsInRadiusInclusive(
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

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchNodesInRadiusExclusive(
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchNodesInRadiusInclusive(
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchNodesInRadiusExclusive(
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

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchNodesInRadiusInclusive(
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

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchConditionsInRadiusExclusive(
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchConditionsInRadiusInclusive(
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchConditionsInRadiusExclusive(
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

template<SpatialContainer TSearhcBackend>
void SpecializedSpatialSearch<TSearhcBackend>::SearchConditionsInRadiusInclusive(
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    VectorDistanceType distances;
    SearchConditionsInRadiusInclusive(rStructureConditions, rInputConditions, rRadius, rResults, distances);
}

template class SpecializedSpatialSearch<SpatialContainer::KDTree>;
template class SpecializedSpatialSearch<SpatialContainer::Octree>;
template class SpecializedSpatialSearch<SpatialContainer::BinsStatic>;
template class SpecializedSpatialSearch<SpatialContainer::BinsDynamic>;
template class SpecializedSpatialSearch<SpatialContainer::BinsStaticObjects>;
template class SpecializedSpatialSearch<SpatialContainer::BinsDynamicObjects>;
} // namespace Kratos.