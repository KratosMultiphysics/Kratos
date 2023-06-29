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

SpatialSearch::ElementSpatialSearchResultContainerMapType SpatialSearch::SearchElementsInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), rModelPart.GetCommunicator().LocalMesh().Elements(), rRadius, rDataCommunicator);
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

SpatialSearch::ElementSpatialSearchResultContainerMapType SpatialSearch::SearchElementsInRadiusExclusive (
    ModelPart& rModelPart,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), rInputElements, rRadius, rDataCommunicator);
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

SpatialSearch::ElementSpatialSearchResultContainerMapType SpatialSearch::SearchElementsInRadiusExclusive (
    const ElementsContainerType& rStructureElements,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchElementsInRadiusExclusive(rStructureElements, rStructureElements, rRadius, rDataCommunicator);
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

SpatialSearch::ElementSpatialSearchResultContainerMapType SpatialSearch::SearchElementsInRadiusExclusive (
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    // First search in the structure
    VectorResultElementsContainerType results;
    VectorDistanceType distance;
    this->SearchElementsInRadiusExclusive(rStructureElements, rInputElements, rRadius, results, distance);
    
    // Now pass this results to the container
    ElementSpatialSearchResultContainerMapType result;
    const SizeType number_of_results = results.size();
    for (IndexType i = 0; i < number_of_results; ++i) {
        // Partial results
        auto& r_partial_results = results[i];
        auto& r_partial_distances = distance[i];
        
        // Getting id of the element
        const IndexType id = (rInputElements.begin() + i)->Id();
        
        // Adding partial results
        auto& r_result_i = result.InitializeResult(id);
        IndexType j = 0;
        for (auto& r_partial_result : r_partial_results) {
            r_result_i.AddResult(r_partial_result.get(), r_partial_distances[j]);
            ++j;
        }
    }

    // Synchronize all the results
    result.SynchronizeAll(rDataCommunicator);
    return result;
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchElementsInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), rModelPart.GetCommunicator().LocalMesh().Elements(), rRadius,rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchElementsInRadiusInclusive (
    ModelPart& rModelPart,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), rInputElements, rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchElementsInRadiusInclusive (
    const ElementsContainerType& rStructureElements,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchElementsInRadiusInclusive(rStructureElements, rStructureElements, rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchElementsInRadiusInclusive (
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    // First search in the structure
    VectorResultNodesContainerType results;
    VectorDistanceType distance;
    this->SearchElementsInRadiusInclusive(rStructureElements, rInputElements, rRadius, results, distance);
    
    // Now pass this results to the container
    NodeSpatialSearchResultContainerMapType result;
    const SizeType number_of_results = results.size();
    for (IndexType i = 0; i < number_of_results; ++i) {
        // Partial results
        auto& r_partial_results = results[i];
        auto& r_partial_distances = distance[i];
        
        // Getting id of the element
        const IndexType id = (rInputElements.begin() + i)->Id();
        
        // Adding partial results
        auto& r_result_i = result.InitializeResult(id);
        IndexType j = 0;
        for (auto& r_partial_result : r_partial_results) {
            r_result_i.AddResult(r_partial_result.get(), r_partial_distances[j]);
            ++j;
        }
    }

    // Synchronize all the results
    result.SynchronizeAll(rDataCommunicator);
    return result;
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchNodesInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(), rModelPart.GetCommunicator().LocalMesh().Nodes(), rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchNodesInRadiusExclusive (
    ModelPart& rModelPart,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(), rInputNodes, rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchNodesInRadiusExclusive (
    const NodesContainerType& rStructureNodes,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchNodesInRadiusExclusive(rStructureNodes, rStructureNodes, rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchNodesInRadiusExclusive (
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    // First search in the structure
    VectorResultNodesContainerType results;
    VectorDistanceType distance;
    this->SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, rRadius, results, distance);
    
    // Now pass this results to the container
    NodeSpatialSearchResultContainerMapType result;
    const SizeType number_of_results = results.size();
    for (IndexType i = 0; i < number_of_results; ++i) {
        // Partial results
        auto& r_partial_results = results[i];
        auto& r_partial_distances = distance[i];
        
        // Getting id of the element
        const IndexType id = (rInputNodes.begin() + i)->Id();
        
        // Adding partial results
        auto& r_result_i = result.InitializeResult(id);
        IndexType j = 0;
        for (auto& r_partial_result : r_partial_results) {
            r_result_i.AddResult(r_partial_result.get(), r_partial_distances[j]);
            ++j;
        }
    }

    // Synchronize all the results
    result.SynchronizeAll(rDataCommunicator);
    return result;
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchNodesInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(), rModelPart.GetCommunicator().LocalMesh().Nodes(), rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchNodesInRadiusInclusive (
    ModelPart& rModelPart,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(), rInputNodes, rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchNodesInRadiusInclusive (
    const NodesContainerType& rStructureNodes,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchNodesInRadiusInclusive(rStructureNodes, rStructureNodes, rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchNodesInRadiusInclusive (
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    // First search in the structure
    VectorResultNodesContainerType results;
    VectorDistanceType distance;
    this->SearchNodesInRadiusInclusive(rStructureNodes, rInputNodes, rRadius, results, distance);
    
    // Now pass this results to the container
    NodeSpatialSearchResultContainerMapType result;
    const SizeType number_of_results = results.size();
    for (IndexType i = 0; i < number_of_results; ++i) {
        // Partial results
        auto& r_partial_results = results[i];
        auto& r_partial_distances = distance[i];
        
        // Getting id of the element
        const IndexType id = (rInputNodes.begin() + i)->Id();
        
        // Adding partial results
        auto& r_result_i = result.InitializeResult(id);
        IndexType j = 0;
        for (auto& r_partial_result : r_partial_results) {
            r_result_i.AddResult(r_partial_result.get(), r_partial_distances[j]);
            ++j;
        }
    }

    // Synchronize all the results
    result.SynchronizeAll(rDataCommunicator);
    return result;
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

SpatialSearch::ConditionSpatialSearchResultContainerMapType SpatialSearch::SearchConditionsInRadiusExclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(), rModelPart.GetCommunicator().LocalMesh().Conditions(), rRadius, rDataCommunicator);
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

SpatialSearch::ConditionSpatialSearchResultContainerMapType SpatialSearch::SearchConditionsInRadiusExclusive (
    ModelPart& rModelPart,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(), rInputConditions, rRadius, rDataCommunicator);
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

SpatialSearch::ConditionSpatialSearchResultContainerMapType SpatialSearch::SearchConditionsInRadiusExclusive (
    const ConditionsContainerType& rStructureConditions,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchConditionsInRadiusExclusive(rStructureConditions, rStructureConditions, rRadius, rDataCommunicator);
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

SpatialSearch::ConditionSpatialSearchResultContainerMapType SpatialSearch::SearchConditionsInRadiusExclusive (
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    // First search in the structure
    VectorResultConditionsContainerType results;
    VectorDistanceType distance;
    this->SearchConditionsInRadiusExclusive(rStructureConditions, rInputConditions, rRadius, results, distance);
    
    // Now pass this results to the container
    ConditionSpatialSearchResultContainerMapType result;
    const SizeType number_of_results = results.size();
    for (IndexType i = 0; i < number_of_results; ++i) {
        // Partial results
        auto& r_partial_results = results[i];
        auto& r_partial_distances = distance[i];
        
        // Getting id of the element
        const IndexType id = (rInputConditions.begin() + i)->Id();
        
        // Adding partial results
        auto& r_result_i = result.InitializeResult(id);
        IndexType j = 0;
        for (auto& r_partial_result : r_partial_results) {
            r_result_i.AddResult(r_partial_result.get(), r_partial_distances[j]);
            ++j;
        }
    }

    // Synchronize all the results
    result.SynchronizeAll(rDataCommunicator);
    return result;
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchConditionsInRadiusInclusive (
    ModelPart& rModelPart,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(), rModelPart.GetCommunicator().LocalMesh().Conditions(), rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchConditionsInRadiusInclusive (
    ModelPart& rModelPart,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(), rInputConditions, rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchConditionsInRadiusInclusive (
    const ConditionsContainerType& rStructureConditions,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    return this->SearchConditionsInRadiusInclusive(rStructureConditions, rStructureConditions, rRadius, rDataCommunicator);
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

SpatialSearch::NodeSpatialSearchResultContainerMapType SpatialSearch::SearchConditionsInRadiusInclusive (
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions,
    const RadiusArrayType& rRadius,
    const DataCommunicator& rDataCommunicator
    )
{
    // First search in the structure
    VectorResultNodesContainerType results;
    VectorDistanceType distance;
    this->SearchConditionsInRadiusInclusive(rStructureConditions, rInputConditions, rRadius, results, distance);
    
    // Now pass this results to the container
    NodeSpatialSearchResultContainerMapType result;
    const SizeType number_of_results = results.size();
    for (IndexType i = 0; i < number_of_results; ++i) {
        // Partial results
        auto& r_partial_results = results[i];
        auto& r_partial_distances = distance[i];
        
        // Getting id of the element
        const IndexType id = (rInputConditions.begin() + i)->Id();
        
        // Adding partial results
        auto& r_result_i = result.InitializeResult(id);
        IndexType j = 0;
        for (auto& r_partial_result : r_partial_results) {
            r_result_i.AddResult(r_partial_result.get(), r_partial_distances[j]);
            ++j;
        }
    }

    // Synchronize all the results
    result.SynchronizeAll(rDataCommunicator);
    return result;
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

void SpatialSearch::SearchNodesOverPointInRadius (
    const NodesContainerType& rStructureNodes,
    const array_1d<double,3>& rPoint,
    const double Radius,
    NodeSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchNodesOverPointNearestPoint (
    const NodesContainerType& rStructureNodes,
    const array_1d<double,3>& rPoint,
    NodeSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverPointInRadius (
    const ElementsContainerType& rStructureElements,
    const array_1d<double,3>& rPoint,
    const double Radius,
    ElementSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchElementsOverPointNearestPoint (
    const ElementsContainerType& rStructureElements,
    const array_1d<double,3>& rPoint,
    ElementSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverPointInRadius (
    const ConditionsContainerType& rStructureConditions,
    const array_1d<double,3>& rPoint,
    const double Radius,
    ConditionSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
    )
{
    KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SpatialSearch::SearchConditionsOverPointNearestPoint (
    const ConditionsContainerType& rStructureConditions,
    const array_1d<double,3>& rPoint,
    ConditionSpatialSearchResultContainerType& rResults,
    const DataCommunicator& rDataCommunicator,
    const bool SyncronizeResults
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
