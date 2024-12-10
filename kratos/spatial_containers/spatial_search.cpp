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
    std::cout << "In SearchElementsInRadiusExclusive" << std::endl;
    this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rRadius,rResults,rResultsDistance);
    std::cout << "End SearchElementsInRadiusExclusive" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    ModelPart& rModelPart,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rInputElements,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    const ElementsContainerType& rStructureElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusExclusive(rStructureElements,
                                        rStructureElements,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
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
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rInputElements,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    const ElementsContainerType& rStructureElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsInRadiusInclusive(rStructureElements,
                                        rStructureElements,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
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
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults
    )
{
    this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rInputElements,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    const ElementsContainerType& rStructureElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults
    )
{
    this->SearchElementsInRadiusExclusive(rStructureElements,
                                        rStructureElements,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusExclusive (
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
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
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                        rInputElements,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    const ElementsContainerType& rStructureElements,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchElementsInRadiusInclusive(rStructureElements,
                                        rStructureElements,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsInRadiusInclusive (
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
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
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                       rInputNodes,
                                       rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    const NodesContainerType& rStructureNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusExclusive(rStructureNodes,
                                        rStructureNodes,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
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
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rInputNodes,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    const NodesContainerType& rStructureNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchNodesInRadiusInclusive(rStructureNodes,
                                        rStructureNodes,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
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
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rInputNodes,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    const NodesContainerType& rStructureNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusExclusive(rStructureNodes,
                                        rStructureNodes,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusExclusive (
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
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
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                        rInputNodes,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    const NodesContainerType& rStructureNodes,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchNodesInRadiusInclusive(rStructureNodes,
                                        rStructureNodes,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesInRadiusInclusive (
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
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
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rInputConditions,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    const ConditionsContainerType& rStructureConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusExclusive(rStructureConditions,
                                        rStructureConditions,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
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
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rInputConditions,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    const ConditionsContainerType& rStructureConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsInRadiusInclusive(rStructureConditions,
                                        rStructureConditions,
                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
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
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults
    )
{
    this->SearchConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rInputConditions,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    const ConditionsContainerType& rStructureConditions,
    const RadiusArrayType& rRadius,
    VectorResultConditionsContainerType& rResults
    )
{
    this->SearchConditionsInRadiusExclusive(rStructureConditions,
                                        rStructureConditions,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusExclusive (
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
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
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                        rInputConditions,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    const ConditionsContainerType& rStructureConditions,
    const RadiusArrayType& rRadius,
    VectorResultNodesContainerType& rResults
    )
{
    this->SearchConditionsInRadiusInclusive(rStructureConditions,
                                        rStructureConditions,
                                        rRadius,rResults);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsInRadiusInclusive (
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
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
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsOverElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                        rInputConditions,
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverElementsInRadiusExclusive (
    const ElementsContainerType& rStructureElements,
    const ConditionsContainerType& rInputConditions,
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
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchConditionsOverElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                        rInputConditions,
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverElementsInRadiusInclusive (
    const ElementsContainerType& rStructureElements,
    const ConditionsContainerType& rInputConditions,
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
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsOverConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                        rInputElements,
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverConditionsInRadiusExclusive (
    const ConditionsContainerType& rStructureConditions,
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
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    VectorResultElementsContainerType& rResults,
    VectorDistanceType& rResultsDistance
    )
{
    this->SearchElementsOverConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                        rInputElements,
                                                        rRadius,rResults,rResultsDistance);
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverConditionsInRadiusInclusive (
    const ConditionsContainerType& rStructureConditions,
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

void SpatialSearch::InitializeBoundingBox(const NodesContainerType& rStructureNodes)
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::InitializeBoundingBox(const ElementsContainerType& rStructureElements)
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::InitializeBoundingBox(const ConditionsContainerType& rStructureConditions)
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

}