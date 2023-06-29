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

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "spatial_containers/spatial_search_result_container.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
* @class SpatialSearch
* @ingroup KratosCore
* @brief This class is used to search for elements, conditions and nodes in a given model part
* @author Carlos Roig
*/
class KRATOS_API(KRATOS_CORE) SpatialSearch
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialSearch
    KRATOS_CLASS_POINTER_DEFINITION(SpatialSearch);

    /// The size type definition
    using SizeType = std::size_t;

    /// The index type definition
    using IndexType = std::size_t;

    /// Defining dimension
    static constexpr std::size_t Dimension = 3;

    /// Maximum level of the tree
    static constexpr std::size_t MAX_LEVEL = 16;

    /// Minimum level of the tree
    static constexpr std::size_t MIN_LEVEL = 2;

    /// The point type
    using PointType = Point;

    /// The node type
    using NodeType = Node;

    /// Nodes classes
    using NodesContainerType = ModelPart::NodesContainerType;
    using ResultNodesContainerType = NodesContainerType::ContainerType;
    using VectorResultNodesContainerType = std::vector<ResultNodesContainerType>;
    using NodeSpatialSearchResultContainerType = SpatialSearchResultContainer<Node>;
    using NodeSpatialSearchResultContainerMapType = SpatialSearchResultContainerMap<Node>;

    /// Elements classes
    using ElementsContainerType = ModelPart::ElementsContainerType;
    using ResultElementsContainerType = ElementsContainerType::ContainerType;
    using VectorResultElementsContainerType = std::vector<ResultElementsContainerType>;
    using ElementSpatialSearchResultContainerType = SpatialSearchResultContainer<GeometricalObject>;
    using ElementSpatialSearchResultContainerMapType = SpatialSearchResultContainerMap<GeometricalObject>;

    /// Conditions classes
    using ConditionsContainerType = ModelPart::ConditionsContainerType;
    using ResultConditionsContainerType = ConditionsContainerType::ContainerType;
    using VectorResultConditionsContainerType = std::vector<ResultConditionsContainerType>;
    using ConditionSpatialSearchResultContainerType = SpatialSearchResultContainer<GeometricalObject>;
    using ConditionSpatialSearchResultContainerMapType = SpatialSearchResultContainerMap<GeometricalObject>;

    /// Input/output types
    using RadiusArrayType = std::vector<double>;
    using DistanceType = std::vector<double>;
    using VectorDistanceType = std::vector<DistanceType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpatialSearch(){}

    /// Destructor.
    virtual ~SpatialSearch(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    //************************************************************************
    // Elemental Exclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual ElementSpatialSearchResultContainerMapType SearchElementsInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ModelPart& rModelPart,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual ElementSpatialSearchResultContainerMapType SearchElementsInRadiusExclusive (
        ModelPart& rModelPart,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" excluding itself
    * @param rStructureElements  List of elements modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusExclusive (
        const ElementsContainerType& rStructureElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" excluding itself
    * @param rStructureElements  List of elements modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual ElementSpatialSearchResultContainerMapType SearchElementsInRadiusExclusive (
        const ElementsContainerType& rStructureElements,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusExclusive (
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual ElementSpatialSearchResultContainerMapType SearchElementsInRadiusExclusive (
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    //************************************************************************
    // Elemental Inclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchElementsInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ModelPart& rModelPart,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchElementsInRadiusInclusive (
        ModelPart& rModelPart,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" including itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusInclusive (
        const ElementsContainerType& rStructureElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" including itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchElementsInRadiusInclusive (
        const ElementsContainerType& rStructureElements,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusInclusive (
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchElementsInRadiusInclusive (
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    //************************************************************************
    // Elemental Exclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ModelPart& rModelPart,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" excluding itself
    * @param rStructureElements  List of nodes against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusExclusive (
        const ElementsContainerType& rStructureElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusExclusive (
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        );

    //************************************************************************
    // Elemental Inclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ModelPart& rModelPart,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" including itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusInclusive (
        const ElementsContainerType& rStructureElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusInclusive (
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    //************************************************************************
    // Nodal Exclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every node in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchNodesInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusExclusive (
        ModelPart& rModelPart,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchNodesInRadiusExclusive (
        ModelPart& rModelPart,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every node in "rStructureNodes" excluding itself
    * @param rStructureNodes     List of nodes against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusExclusive (
        const NodesContainerType& rStructureNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "rStructureNodes" excluding itself
    * @param rStructureNodes     List of nodes against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchNodesInRadiusExclusive (
        const NodesContainerType& rStructureNodes,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" excluding itself
    * @param rModelPart          List of nodes against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusExclusive (
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" excluding itself
    * @param rModelPart          List of nodes against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchNodesInRadiusExclusive (
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    //************************************************************************
    // Nodal Inclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every node in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchNodesInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusInclusive (
        ModelPart& rModelPart,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchNodesInRadiusInclusive (
        ModelPart& rModelPart,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every node in "rStructureNodes" including itself
    * @param rStructureNodes     List of nodes against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusInclusive (
        const NodesContainerType& rStructureNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "rStructureNodes" including itself
    * @param rStructureNodes     List of nodes against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchNodesInRadiusInclusive (
        const NodesContainerType& rStructureNodes,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" including itself
    * @param rStructureNodes     List of nodes against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusInclusive (
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" including itself
    * @param rStructureNodes     List of nodes against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchNodesInRadiusInclusive (
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    //************************************************************************
    // Nodal Exclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every node in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusExclusive (
        ModelPart& rModelPart,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "rStructureNodes" excluding itself
    * @param rStructureNodes     List of nodes against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusExclusive (
        const NodesContainerType& rStructureNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" excluding itself
    * @param rModelPart          List of nodes against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusExclusive (
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    //************************************************************************
    // Nodal Inclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every node in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusInclusive (
        ModelPart& rModelPart,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "rStructureNodes" including itself
    * @param rStructureNodes     List of nodes against which the neighbours are searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusInclusive (
        const NodesContainerType& rStructureNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" including itself
    * @param rStructureNodes     List of nodes against which the neighbours are searched
    * @param rInputNodes         List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusInclusive (
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    //************************************************************************
    // Conditional Exclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every Condition in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    * @param rResultDistance     Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every Condition
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual ConditionSpatialSearchResultContainerMapType SearchConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    * @param rResultDistance     Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every Condition
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual ConditionSpatialSearchResultContainerMapType SearchConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every Condition in "rStructureConditions" excluding itself
    * @param rStructureConditions  List of conditions modelpart against which the neighbours are searched
    * @param rRadius               List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    * @param rResultDistance       Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        const ConditionsContainerType& rStructureConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "rStructureConditions" excluding itself
    * @param rStructureConditions  List of conditions modelpart against which the neighbours are searched
    * @param rRadius               List of search radius for every Condition
    * @param rDataCommunicator     The data communicator
    * @return                      The results maps
    */
    virtual ConditionSpatialSearchResultContainerMapType SearchConditionsInRadiusExclusive (
        const ConditionsContainerType& rStructureConditions,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    * @param rStructureConditions  List of conditions against which the neighbours are searched
    * @param rInputConditions      List of conditions to be searched
    * @param rRadius               List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    * @param rResultDistance       Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    * @param rStructureConditions  List of conditions against which the neighbours are searched
    * @param rInputConditions      List of conditions to be searched
    * @param rRadius               List of search radius for every Condition
    * @param rDataCommunicator     The data communicator
    * @return                      The results maps
    */
    virtual ConditionSpatialSearchResultContainerMapType SearchConditionsInRadiusExclusive (
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    //************************************************************************
    // Conditional Inclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every Condition in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    * @param rResultDistance     Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every Condition
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    * @param rResultDistance     Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every Condition
    * @param rDataCommunicator   The data communicator
    * @return                    The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every Condition in "rStructureConditions" including itself
    * @param rStructureConditions  List of conditions against which the neighbours are searched
    * @param rRadius               List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    * @param rResultDistance       Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        const ConditionsContainerType& rStructureConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "rStructureConditions" including itself
    * @param rStructureConditions  List of conditions against which the neighbours are searched
    * @param rRadius               List of search radius for every Condition
    * @param rDataCommunicator     The data communicator
    * @return                      The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchConditionsInRadiusInclusive (
        const ConditionsContainerType& rStructureConditions,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" including itself
    * @param rStructureConditions  List of conditions against which the neighbours are searched
    * @param rInputConditions      List of conditions to be searched
    * @param rRadius               List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    * @param rResultDistance       Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" including itself
    * @param rStructureConditions  List of conditions against which the neighbours are searched
    * @param rInputConditions      List of conditions to be searched
    * @param rRadius               List of search radius for every Condition
    * @param rDataCommunicator     The data communicator
    * @return                      The results maps
    */
    virtual NodeSpatialSearchResultContainerMapType SearchConditionsInRadiusInclusive (
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        );

    //************************************************************************
    // Conditional Exclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every Condition in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "rStructureConditions" excluding itself
    * @param rStructureConditions  List of nodes against which the neighbours are searched
    * @param rRadius               List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        const ConditionsContainerType& rStructureConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    * @param rStructureConditions  List of conditions against which the neighbours are searched
    * @param rInputConditions      List of conditions to be searched
    * @param rRadius               List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        );

    //************************************************************************
    // Conditional Inclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every Condition in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "rStructureConditions" including itself
    * @param rStructureConditions  List of conditions against which the neighbours are searched
    * @param rRadius               List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        const ConditionsContainerType& rStructureConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" including itself
    * @param rStructureConditions  List of conditions against which the neighbours are searched
    * @param rInputConditions      List of conditions to be searched
    * @param rRadius               List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    //************************************************************************
    // Elemental vs Condition Exclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusExclusive (
        ModelPart& rModelPart,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusExclusive (
        const ElementsContainerType& rStructureElements,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    //************************************************************************
    // Elemental vs Condition Inclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusInclusive (
        ModelPart& rModelPart,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rStructureElements  List of elements against which the neighbours are searched
    * @param rInputConditions    List of conditions to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusInclusive (
        const ElementsContainerType& rStructureElements,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    //************************************************************************
    // Condition vs Elemental Exclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rStructureConditions List of conditions against which the neighbours are searched
    * @param InputElements        List of elements to be searched
    * @param rRadius              List of search radius for every element
    * @param rResults             Array of results for each element
    * @param rResultDistance      Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusExclusive (
        const ConditionsContainerType& rStructureConditions,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    //************************************************************************
    // Condition vs Elemental Inclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param rInputElements      List of elements to be searched
    * @param rRadius             List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rStructureConditions List of conditions against which the neighbours are searched
    * @param InputElements        List of elements to be searched
    * @param rRadius              List of search radius for every element
    * @param rResults             Array of results for each element
    * @param rResultDistance      Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusInclusive (
        const ConditionsContainerType& rStructureConditions,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    //************************************************************************
    // Point vs Entities (these are new interfaces and already use the new containers)
    //************************************************************************

    /**
     * @brief Search neighbours nodes for one point in a given radius
     * @param rStructureNodes      List of nodes to be searched
     * @param rPoint               Point to be searched
     * @param Radius               Radius of the search
     * @param rResults             Results of the search
     * @param rDataCommunicator    The data communicator
     * @param SyncronizeResults    If true, the results are synchronized
     */
    virtual void SearchNodesOverPointInRadius (
        const NodesContainerType& rStructureNodes,
        const array_1d<double,3>& rPoint,
        const double Radius,
        NodeSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        );

    /**
     * @brief Search neighbours nodes for several points in a given radius
     * @param rStructureNodes      List of nodes to be searched
     * @param itPointBegin         Iterator to the first point
     * @param itPointEnd           Iterator to the last point
     * @param rRadius              List of search radius for every node
     * @param rDataCommunicator    The data communicator
     */
    template<typename TPointIteratorType>    
    NodeSpatialSearchResultContainerMapType SearchNodesOverPointsInRadius (
        const NodesContainerType& rStructureNodes,
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        )
    {
        NodeSpatialSearchResultContainerMapType results;
        int i = 0;
        for (auto it_point = itPointBegin; it_point != itPointEnd; ++it_point) {
            const Point& r_point = *it_point;
            auto& r_partial_result = results.InitializeResult(r_point);
            SearchNodesOverPointInRadius(rStructureNodes, r_point, rRadius[i], r_partial_result, rDataCommunicator);
            ++i;
        }
        return results;
    }

    /**
     * @brief Search nearest neighbour node for one point
     * @param rStructureNodes      List of nodes to be searched
     * @param rPoint               Point to be searched
     * @param rResults             Results of the search
     * @param rDataCommunicator    The data communicator
     * @param SyncronizeResults    If true, the results are synchronized
     */
    virtual void SearchNodesOverPointNearestPoint (
        const NodesContainerType& rStructureNodes,
        const array_1d<double,3>& rPoint,
        NodeSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        );

    /**
     * @brief Search nearest neighbour node for several points
     * @param rStructureNodes      List of nodes to be searched
     * @param itPointBegin         Iterator to the first point
     * @param itPointEnd           Iterator to the last point
     * @param rDataCommunicator    The data communicator
     */
    template<typename TPointIteratorType>    
    NodeSpatialSearchResultContainerMapType SearchNodesOverPointsNearestPoint (
        const NodesContainerType& rStructureNodes,
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const DataCommunicator& rDataCommunicator
        )
    {
        NodeSpatialSearchResultContainerMapType results;
        for (auto it_point = itPointBegin; it_point != itPointEnd; ++it_point) {
            const Point& r_point = *it_point;
            auto& r_partial_result = results.InitializeResult(r_point);
            SearchNodesOverPointNearestPoint(rStructureNodes, r_point, r_partial_result, rDataCommunicator);
        }
        return results;
    }

    /**
     * @brief Search neighbours elements for one point in a given radius
     * @param rStructureElements   List of elements to be searched
     * @param rPoint               Point to be searched
     * @param Radius               Radius of the search
     * @param rResults             Results of the search
     * @param rDataCommunicator    The data communicator
     * @param SyncronizeResults    If true, the results are synchronized
     */
    virtual void SearchElementsOverPointInRadius (
        const ElementsContainerType& rStructureElements,
        const array_1d<double,3>& rPoint,
        const double Radius,
        ElementSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        );

    /**
     * @brief Search neighbours elements for several points in a given radius
     * @param rStructureElements   List of elements to be searched
     * @param itPointBegin         Iterator to the first point
     * @param itPointEnd           Iterator to the last point
     * @param rRadius              List of search radius for every element
     * @param rDataCommunicator    The data communicator
     */
    template<typename TPointIteratorType>    
    ElementSpatialSearchResultContainerMapType SearchElementsOverPointsInRadius (
        const ElementsContainerType& rStructureElements,
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        )
    {
        ElementSpatialSearchResultContainerMapType results;
        int i = 0;
        for (auto it_point = itPointBegin; it_point != itPointEnd; ++it_point) {
            const Point& r_point = *it_point;
            auto& r_partial_result = results.InitializeResult(r_point);
            SearchElementsOverPointInRadius(rStructureElements, r_point, rRadius[i], r_partial_result, rDataCommunicator);
            ++i;
        }
        return results;
    }

    /**
     * @brief Search nearest neighbour element for one point
     * @param rStructureElements   List of elements to be searched
     * @param rPoint               Point to be searched
     * @param rResults             Results of the search
     * @param rDataCommunicator    The data communicator
     * @param SyncronizeResults    If true, the results are synchronized
     */
    virtual void SearchElementsOverPointNearestPoint (
        const ElementsContainerType& rStructureElements,
        const array_1d<double,3>& rPoint,
        ElementSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        );

    /**
     * @brief Search nearest neighbour element for several points
     * @param rStructureElements   List of elements to be searched
     * @param itPointBegin         Iterator to the first point
     * @param itPointEnd           Iterator to the last point
     * @param rDataCommunicator    The data communicator
     */
    template<typename TPointIteratorType>    
    ElementSpatialSearchResultContainerMapType SearchElementsOverPointsNearestPoint (
        const ElementsContainerType& rStructureElements,
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const DataCommunicator& rDataCommunicator
        )
    {
        ElementSpatialSearchResultContainerMapType results;
        for (auto it_point = itPointBegin; it_point != itPointEnd; ++it_point) {
            const Point& r_point = *it_point;
            auto& r_partial_result = results.InitializeResult(r_point);
            SearchElementsOverPointNearestPoint(rStructureElements, r_point, r_partial_result, rDataCommunicator);
        }
        return results;
    }

    /**
     * @brief Search neighbours conditions for one point in a given radius
     * @param rStructureConditions List of conditions to be searched
     * @param rPoint               Point to be searched
     * @param Radius               Radius of the search
     * @param rResults             Results of the search
     * @param rDataCommunicator    The data communicator
     * @param SyncronizeResults    If true, the results are synchronized
     */
    virtual void SearchConditionsOverPointInRadius (
        const ConditionsContainerType& rStructureConditions,
        const array_1d<double,3>& rPoint,
        const double Radius,
        ConditionSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        );

    /**
     * @brief Search neighbours conditions for several points in a given radius
     * @param rStructureConditions List of conditions to be searched
     * @param itPointBegin         Iterator to the first point
     * @param itPointEnd           Iterator to the last point
     * @param rRadius              List of search radius for every condition
     * @param rDataCommunicator    The data communicator
     */
    template<typename TPointIteratorType>    
    ConditionSpatialSearchResultContainerMapType SearchConditionsOverPointsInRadius (
        const ConditionsContainerType& rStructureConditions,
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        )
    {
        ConditionSpatialSearchResultContainerMapType results;
        int i = 0;
        for (auto it_point = itPointBegin; it_point != itPointEnd; ++it_point) {
            const Point& r_point = *it_point;
            auto& r_partial_result = results.InitializeResult(r_point);
            SearchConditionsOverPointInRadius(rStructureConditions, r_point, rRadius[i], r_partial_result, rDataCommunicator);
            ++i;
        }
        return results;
    }

    /**
     * @brief Search nearest neighbour condition for one point
     * @param rStructureConditions List of conditions to be searched
     * @param rPoint               Point to be searched
     * @param rResults             Results of the search
     * @param rDataCommunicator    The data communicator
     * @param SyncronizeResults    If true, the results are synchronized
     */
    virtual void SearchConditionsOverPointNearestPoint (
        const ConditionsContainerType& rStructureConditions,
        const array_1d<double,3>& rPoint,
        ConditionSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        );

    /**
     * @brief Search nearest neighbour condition for several points
     * @param rStructureConditions List of conditions to be searched
     * @param itPointBegin         Iterator to the first point
     * @param itPointEnd           Iterator to the last point
     * @param rDataCommunicator    The data communicator
     */
    template<typename TPointIteratorType>    
    ConditionSpatialSearchResultContainerMapType SearchConditionsOverPointsNearestPoint (
        const ConditionsContainerType& rStructureConditions,
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const DataCommunicator& rDataCommunicator
        )
    {
        ConditionSpatialSearchResultContainerMapType results;
        for (auto it_point = itPointBegin; it_point != itPointEnd; ++it_point) {
            const Point& r_point = *it_point;
            auto& r_partial_result = results.InitializeResult(r_point);
            SearchConditionsOverPointNearestPoint(rStructureConditions, r_point, r_partial_result, rDataCommunicator);
        }
        return results;
    }

    //************************************************************************
    // Bounding box methods
    //************************************************************************

   /**
     * @brief This method allows to initialize the local bounding box (for nodes)
     * @param rStructureNodes The container of nodes
     */
    virtual void InitializeBoundingBox(const NodesContainerType& rStructureNodes);

    /**
     * @param rStructureElements The container of elements
     * @brief This method allows to initialize the local bounding box (for elements)
     */
    virtual void InitializeBoundingBox(const ElementsContainerType& rStructureElements);

    /**
     * @brief This method allows to initialize the local bounding box (for conditions)
     * @param rStructureConditions The container of conditions
     */
    virtual void InitializeBoundingBox(const ConditionsContainerType& rStructureConditions);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "SpatialSearch" ;

        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SpatialSearch";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SpatialSearch& operator=(SpatialSearch const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    SpatialSearch(SpatialSearch const& rOther)
    {
        *this = rOther;
    }

    ///@}
}; // Class SpatialSearch

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, 
                const SpatialSearch& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}addtogroup block

}  // namespace Kratos.