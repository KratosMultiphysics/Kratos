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

    /// Defining dimension
    static constexpr std::size_t Dimension = 3;

    /// Maximum level of the tree
    static constexpr std::size_t MAX_LEVEL = 16;

    /// Minimum level of the tree
    static constexpr std::size_t MIN_LEVEL = 2;

    /// Common Defines
    using PointType = Point;

    using ElementsContainerType = ModelPart::ElementsContainerType;
    using ElementType = ModelPart::ElementType;
    using ElementPointerType = ModelPart::ElementType::Pointer;
    using ResultElementsContainerType = ElementsContainerType::ContainerType;
    using VectorResultElementsContainerType = std::vector<ResultElementsContainerType>;

    using NodesContainerType = ModelPart::NodesContainerType;
    using NodeType = ModelPart::NodeType;
    using NodePointerType = ModelPart::NodeType::Pointer;
    using ResultNodesContainerType = NodesContainerType::ContainerType;
    using VectorResultNodesContainerType = std::vector<ResultNodesContainerType>;

    using ConditionsContainerType = ModelPart::ConditionsContainerType;
    using ConditionType = ModelPart::ConditionType;
    using ConditionPointerType = ModelPart::ConditionType::Pointer;
    using ResultConditionsContainerType = ConditionsContainerType::ContainerType;
    using VectorResultConditionsContainerType = std::vector<ResultConditionsContainerType>;

    using RadiusArrayType = std::vector<double>;
    using DistanceType = std::vector<double>;
    using VectorDistanceType = std::vector<DistanceType>;

    using ResultIteratorType = ElementsContainerType::ContainerType::iterator;

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
    * @param Radius              List of search radius for every element
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
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ModelPart& rModelPart,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" excluding itself
    * @param StructureElements   List of elements modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ElementsContainerType const& StructureElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ElementsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    //************************************************************************
    // Elemental Inclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every element
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
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ModelPart& rModelPart,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" including itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ElementsContainerType const& StructureElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ElementsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    //************************************************************************
    // Elemental Exclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ModelPart& rModelPart,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ModelPart& rModelPart,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" excluding itself
    * @param StructureElements   List of nodes against which the neighbours are searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ElementsContainerType const& StructureElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusExclusive (
        ElementsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        );

    //************************************************************************
    // Elemental Inclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every element
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
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ModelPart& rModelPart,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "StructureElements" including itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ElementsContainerType const& StructureElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    */
    virtual void SearchElementsInRadiusInclusive (
        ElementsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    //************************************************************************
    // Nodal Exclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every node in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every node
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
    * @brief Search neighbours for every node in "InputNodes" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputNodes          List of nodes to be searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusExclusive (
        ModelPart& rModelPart,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "StructureNodes" excluding itself
    * @param StructureNodes      Lis of nodes against which the neighbours are searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusExclusive (
        NodesContainerType const& StructureNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" excluding itself
    * @param rModelPart          List of nodes against which the neighbours are searched
    * @param InputNodes          List of nodes to be searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusExclusive (
        NodesContainerType const& StructureNodes,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    //************************************************************************
    // Nodal Inclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every node in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every node
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
    * @brief Search neighbours for every node in "InputNodes" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputNodes          List of nodes to be searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusInclusive (
        ModelPart& rModelPart,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "StructureNodes" including itself
    * @param StructureNodes      List of nodes against which the neighbours are searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusInclusive (
        NodesContainerType const& StructureNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" including itself
    * @param StructureNodes      List of nodes against which the neighbours are searched
    * @param InputNodes          List of nodes to be searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    virtual void SearchNodesInRadiusInclusive (
        NodesContainerType const& StructureNodes,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );


    //************************************************************************
    // Nodal Exclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every node in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every node
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
    * @param InputNodes          List of nodes to be searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusExclusive (
        ModelPart& rModelPart,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "StructureNodes" excluding itself
    * @param StructureNodes      Lis of nodes against which the neighbours are searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusExclusive (
        NodesContainerType const& StructureNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" excluding itself
    * @param rModelPart          List of nodes against which the neighbours are searched
    * @param InputNodes          List of nodes to be searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusExclusive (
        NodesContainerType const& StructureNodes,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    //************************************************************************
    // Nodal Inclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every node in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every node
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
    * @param InputNodes          List of nodes to be searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusInclusive (
        ModelPart& rModelPart,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "StructureNodes" including itself
    * @param StructureNodes      List of nodes against which the neighbours are searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusInclusive (
        NodesContainerType const& StructureNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every node in "InputNodes" including itself
    * @param StructureNodes      List of nodes against which the neighbours are searched
    * @param InputNodes          List of nodes to be searched
    * @param Radius              List of search radius for every node
    * @param rResults            Array of results for each node
    */
    virtual void SearchNodesInRadiusInclusive (
        NodesContainerType const& StructureNodes,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    //************************************************************************
    // Conditional Exclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every Condition in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every Condition
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
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputConditions     List of Conditions to be searched
    * @param Radius              List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    * @param rResultDistance     Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "StructureConditions" excluding itself
    * @param StructureConditions   List of Conditions modelpart against which the neighbours are searched
    * @param Radius                List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    * @param rResultDistance       Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ConditionsContainerType const& StructureConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    * @param StructureConditions   List of Conditions against which the neighbours are searched
    * @param InputConditions       List of Conditions to be searched
    * @param Radius                List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    * @param rResultDistance       Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ConditionsContainerType const& StructureConditions,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    //************************************************************************
    // Conditional Inclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every Condition in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every Condition
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
    * @brief Search neighbours for every Condition in "InputConditions" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputConditions     List of Conditions to be searched
    * @param Radius              List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    * @param rResultDistance     Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "StructureConditions" including itself
    * @param StructureConditions   List of Conditions against which the neighbours are searched
    * @param Radius                List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    * @param rResultDistance       Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ConditionsContainerType const& StructureConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" including itself
    * @param StructureConditions   List of Conditions against which the neighbours are searched
    * @param InputConditions       List of Conditions to be searched
    * @param Radius                List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    * @param rResultDistance       Array of distances for each result of each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ConditionsContainerType const& StructureConditions,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    //************************************************************************
    // Conditional Exclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every Condition in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every Condition
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
    * @param InputConditions     List of Conditions to be searched
    * @param Radius              List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "StructureConditions" excluding itself
    * @param StructureConditions   List of nodes against which the neighbours are searched
    * @param Radius                List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ConditionsContainerType const& StructureConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    * @param StructureConditions   List of Conditions against which the neighbours are searched
    * @param InputConditions       List of Conditions to be searched
    * @param Radius                List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusExclusive (
        ConditionsContainerType const& StructureConditions,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        );

    //************************************************************************
    // Conditional Inclusive search without distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every Condition in "rModelpart" including itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every Condition
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
    * @param InputConditions     List of Conditions to be searched
    * @param Radius              List of search radius for every Condition
    * @param rResults            Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "StructureConditions" including itself
    * @param StructureConditions   List of Conditions against which the neighbours are searched
    * @param Radius                List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ConditionsContainerType const& StructureConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    /**
    * @brief Search neighbours for every Condition in "InputConditions" including itself
    * @param StructureConditions   List of Conditions against which the neighbours are searched
    * @param InputConditions       List of Conditions to be searched
    * @param Radius                List of search radius for every Condition
    * @param rResults              Array of results for each Condition
    */
    virtual void SearchConditionsInRadiusInclusive (
        ConditionsContainerType const& StructureConditions,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        );

    //************************************************************************
    // Elemental vs Condition Exclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every element in "rModelpart" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param Radius              List of search radius for every element
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
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputConditions     List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusExclusive (
        ModelPart& rModelPart,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputConditions     List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusExclusive (
        ElementsContainerType const& StructureElements,
        ConditionsContainerType const& InputConditions,
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
    * @param Radius              List of search radius for every element
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
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputConditions     List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusInclusive (
        ModelPart& rModelPart,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputConditions     List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchConditionsOverElementsInRadiusInclusive (
        ElementsContainerType const& StructureElements,
        ConditionsContainerType const& InputConditions,
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
    * @param Radius              List of search radius for every element
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
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusExclusive (
        ModelPart& rModelPart,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusExclusive (
        ConditionsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
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
    * @param Radius              List of search radius for every element
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
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param rModelPart          Input modelpart against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusInclusive (
        ModelPart& rModelPart,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

    /**
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    virtual void SearchElementsOverConditionsInRadiusInclusive (
        ConditionsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        );

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


