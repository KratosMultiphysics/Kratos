//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//

// System includes

// External includes

// Project includes
#include "spatial_containers/spatial_search.h"

namespace Kratos
{

void SpatialSearch::SearchElementsInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    ModelPart& rModelPart,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        InputElements,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    ElementsContainerType const& StructureElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusExclusive(StructureElements,
                                        StructureElements,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    ElementsContainerType const& StructureElements,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    ModelPart& rModelPart,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        InputElements,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    ElementsContainerType const& StructureElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusInclusive(StructureElements,
                                        StructureElements,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    ElementsContainerType const& StructureElements,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults
    )
{
    this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    ModelPart& rModelPart,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults
    )
{
    this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        InputElements,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    ElementsContainerType const& StructureElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults
    )
{
    this->SearchElementsInRadiusExclusive(StructureElements,
                                        StructureElements,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    ElementsContainerType const& StructureElements,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    ModelPart& rModelPart,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        InputElements,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    ElementsContainerType const& StructureElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchElementsInRadiusInclusive(StructureElements,
                                        StructureElements,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    ElementsContainerType const& StructureElements,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    ModelPart& rModelPart,
    NodesContainerType const& InputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        InputNodes,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    NodesContainerType const& StructureNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusExclusive(StructureNodes,
                                        StructureNodes,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    NodesContainerType const& StructureNodes,
    NodesContainerType const& InputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    ModelPart& rModelPart,
    NodesContainerType const& InputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        InputNodes,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    NodesContainerType const& StructureNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusInclusive(StructureNodes,
                                        StructureNodes,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    NodesContainerType const& StructureNodes,
    NodesContainerType const& InputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    ModelPart& rModelPart,
    NodesContainerType const& InputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        InputNodes,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    NodesContainerType const& StructureNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusExclusive(StructureNodes,
                                        StructureNodes,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    NodesContainerType const& StructureNodes,
    NodesContainerType const& InputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    ModelPart& rModelPart,
    NodesContainerType const& InputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        InputNodes,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    NodesContainerType const& StructureNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusInclusive(StructureNodes,
                                        StructureNodes,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    NodesContainerType const& StructureNodes,
    NodesContainerType const& InputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    ModelPart& rModelPart,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        InputConditions,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    ConditionsContainerType const& StructureConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusExclusive(StructureConditions,
                                        StructureConditions,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    ConditionsContainerType const& StructureConditions,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    ModelPart& rModelPart,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        InputConditions,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    ConditionsContainerType const& StructureConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusInclusive(StructureConditions,
                                        StructureConditions,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    ConditionsContainerType const& StructureConditions,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults
    )
{
    this->SearchConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    ModelPart& rModelPart,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults
    )
{
    this->SearchConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        InputConditions,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    ConditionsContainerType const& StructureConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults
    )
{
    this->SearchConditionsInRadiusExclusive(StructureConditions,
                                        StructureConditions,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    ConditionsContainerType const& StructureConditions,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    ModelPart& rModelPart,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        InputConditions,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    ConditionsContainerType const& StructureConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchConditionsInRadiusInclusive(StructureConditions,
                                        StructureConditions,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    ConditionsContainerType const& StructureConditions,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverElementsInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsOverElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                        rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverElementsInRadiusExclusive (
    ModelPart& rModelPart,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsOverElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                        InputConditions,
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverElementsInRadiusExclusive (
    ElementsContainerType const& StructureElements,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverElementsInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsOverElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                        rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverElementsInRadiusInclusive (
    ModelPart& rModelPart,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsOverElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                        InputConditions,
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverElementsInRadiusInclusive (
    ElementsContainerType const& StructureElements,
    ConditionsContainerType const& InputConditions,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverConditionsInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsOverConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                        rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverConditionsInRadiusExclusive (
    ModelPart& rModelPart,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsOverConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                        InputElements,
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverConditionsInRadiusExclusive (
    ConditionsContainerType const& StructureElements,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverConditionsInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsOverConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                        rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverConditionsInRadiusInclusive (
    ModelPart& rModelPart,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsOverConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                        InputElements,
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverConditionsInRadiusInclusive (
    ConditionsContainerType const& StructureElements,
    ElementsContainerType const& InputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

}

