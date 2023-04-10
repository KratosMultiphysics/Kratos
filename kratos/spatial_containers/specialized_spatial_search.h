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

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/point.h"
#include "includes/kratos_parameters.h"
#include "spatial_containers/spatial_search.h"

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

/**
 * @brief This enum defines the different spatial containers available
 * @details The different spatial containers available are:
 * - KDTree: This is a tree based on the k-d tree algorithm
 * - Octree: This is a tree based on the octree algorithm
 * - BinsStatic: This is a bin based on the static bins algorithm
 * - BinsDynamic: This is a bin based on the dynamic bins algorithm
 * - BinsStaticObjects: This is a bin based on the static bins algorithm, but with objects
 * - BinsDynamicObjects: This is a bin based on the dynamic bins algorithm, but with objects
 */
enum class SpatialContainer
{
    KDTree,
    Octree,
    BinsStatic,
    BinsDynamic,
    BinsStaticObjects,
    BinsDynamicObjects
};

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class PointObject
 * @ingroup KratosCore
 * @brief Custom Point container to be used by the search
 * @details It stores the pointer of a certain object
 * @author Vicente Mataix Ferrandiz
 */
template<class TObject>
class PointObject
    : public Point
{
public:

    ///@name Type Definitions
    ///@{

    /// Base class definition
    typedef Point BaseType;

    /// Counted pointer of PointObject
    KRATOS_CLASS_POINTER_DEFINITION( PointObject );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointObject():
        BaseType()
    {}

    PointObject(const double X, const double Y, const double Z)
        : BaseType(X, Y, Z)
    {

    }

    PointObject(typename TObject::Pointer pObject):
        mpObject(pObject)
    {
        UpdatePoint();
    }

    ///Copy constructor  (not really required)
    PointObject(const PointObject& rRHS):
        BaseType(rRHS),
        mpObject(rRHS.mpObject)
    {
    }

    /// Destructor.
    ~PointObject() override= default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the geometry associated to the point
     * @return mrGeometry The reference to the geometry associated to the point
     */
    typename TObject::Pointer pGetObject()
    {
        return mpObject;
    }

    /**
     * @brief This function updates the database, using as base for the coordinates the condition center
     */
    void UpdatePoint();

private:
    ///@}
    ///@name Member Variables
    ///@{

    typename TObject::Pointer mpObject = nullptr;

    ///@}

}; // Class PointObject

/**
* @class SpecializedSpatialSearch
* @ingroup KratosCore
* @brief This class is used to search for elements, conditions and nodes in a given model part
* @details In order to perform the search it uses as backend some of the the spatial containers defined `spatial_containers` folder
* @tparam TSearchBackend The spatial container to be used as backend
* @author Vicente Mataix Ferrandiz
*/
template<SpatialContainer TSearchBackend>
class KRATOS_API(KRATOS_CORE) SpecializedSpatialSearch
    : public SpatialSearch
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpecializedSpatialSearch
    KRATOS_CLASS_POINTER_DEFINITION(SpecializedSpatialSearch);

    /// The base type
    using BaseType = SpatialSearch;

    /// Common Defines
    using PointType = BaseType::PointType;

    using ElementsContainerType = BaseType::ElementsContainerType;
    using ElementType = BaseType::ElementType;
    using ElementPointerType = BaseType::ElementPointerType;
    using ResultElementsContainerType = BaseType::ResultElementsContainerType;
    using VectorResultElementsContainerType = BaseType::VectorResultElementsContainerType;

    using NodesContainerType = BaseType::NodesContainerType;
    using NodeType = BaseType::NodeType;
    using NodePointerType = BaseType::NodePointerType;
    using ResultNodesContainerType = BaseType::ResultNodesContainerType;
    using VectorResultNodesContainerType = BaseType::VectorResultNodesContainerType;

    using ConditionsContainerType = BaseType::ConditionsContainerType;
    using ConditionType = BaseType::ConditionType;
    using ConditionPointerType = BaseType::ConditionPointerType;
    using ResultConditionsContainerType = BaseType::ResultConditionsContainerType;
    using VectorResultConditionsContainerType = BaseType::VectorResultConditionsContainerType;

    using RadiusArrayType = BaseType::RadiusArrayType;
    using DistanceType = BaseType::DistanceType;
    using VectorDistanceType = BaseType::VectorDistanceType;

    using ResultIteratorType = BaseType::ResultIteratorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SpecializedSpatialSearch()
    {
        // We asssign default parameters
        mParameters = GetDefaultParameters();
    }

    /**
     * @brief Constructor with parameters
     * @param ThisParameters The parameters to be considered
     */
    SpecializedSpatialSearch(Parameters ThisParameters)
        : mParameters(ThisParameters)
    {
        // We asssign default parameters
        const Parameters default_parameters = GetDefaultParameters();

        // We update the parameters
        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    }

    /// Destructor.
    ~SpecializedSpatialSearch() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rStructureElements   List of elements against which the neighbours are searched
    * @param rInputElements       List of elements to be searched
    * @param rRadius              List of search radius for every element
    * @param rResults             Array of results for each element
    * @param rResultDistance      Array of distances for each result of each element
    */
    void SearchElementsInRadiusExclusive(
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override;

    /**
    * @brief Search neighbours for every element in "rInputElements" including itself
    * @param rStructureElements   List of elements against which the neighbours are searched
    * @param rInputElements       List of elements to be searched
    * @param rRadius              List of search radius for every element
    * @param rResults             Array of results for each element
    * @param rResultDistance      Array of distances for each result of each element
    */
    void SearchElementsInRadiusInclusive(
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override;

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    * @param rStructureElements   List of elements against which the neighbours are searched
    * @param rInputElements       List of elements to be searched
    * @param rRadius              List of search radius for every element
    * @param rResults             Array of results for each element
    */
    void SearchElementsInRadiusExclusive(
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        ) override;

    /**
    * @brief Search neighbours for every element in "InputElements" including itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param rRadius              List of search radius for every element
    * @param rResults            Array of results for each element
    */
    void SearchElementsInRadiusInclusive(
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        ) override;

    /**
    * @brief Search neighbours for every node in "rInputNodes" excluding itself
    * @param rModelPart          List of nodes against which the neighbours are searched
    * @param rInputNodes          List of nodes to be searched
    * @param rRadius             List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    void SearchNodesInRadiusExclusive(
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override;

    /**
    * @brief Search neighbours for every node in "rInputNodes" including itself
    * @param rStructureNodes      List of nodes against which the neighbours are searched
    * @param rInputNodes          List of nodes to be searched
    * @param rRadius              List of search radius for every node
    * @param rResults            Array of results for each node
    * @param rResultDistance     Array of distances for each result of each node
    */
    void SearchNodesInRadiusInclusive(
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override;

    /**
    * @brief Search neighbours for every node in "rInputNodes" excluding itself
    * @param rModelPart          List of nodes against which the neighbours are searched
    * @param rInputNodes          List of nodes to be searched
    * @param rRadius              List of search radius for every node
    * @param rResults            Array of results for each node
    */
    void SearchNodesInRadiusExclusive(
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        ) override;

    /**
    * @brief Search neighbours for every node in "rInputNodes" including itself
    * @param rStructureNodes      List of nodes against which the neighbours are searched
    * @param rInputNodes          List of nodes to be searched
    * @param rRadius              List of search radius for every node
    * @param rResults             Array of results for each node
    */
    void SearchNodesInRadiusInclusive(
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        ) override;

    /**
    * @brief Search neighbours for every Condition in "rInputConditions" excluding itself
    * @param rStructureConditions   List of Conditions against which the neighbours are searched
    * @param rInputConditions       List of Conditions to be searched
    * @param rRadius                List of search radius for every Condition
    * @param rResults               Array of results for each Condition
    * @param rResultDistance        Array of distances for each result of each Condition
    */
    void SearchConditionsInRadiusExclusive(
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override;

    /**
    * @brief Search neighbours for every Condition in "rInputConditions" including itself
    * @param rStructureConditions   List of Conditions against which the neighbours are searched
    * @param rInputConditions       List of Conditions to be searched
    * @param rRadius                List of search radius for every Condition
    * @param rResults               Array of results for each Condition
    * @param rResultDistance        Array of distances for each result of each Condition
    */
    void SearchConditionsInRadiusInclusive(
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override;

    /**
    * @brief Search neighbours for every Condition in "rInputConditions" excluding itself
    * @param rStructureConditions   List of Conditions against which the neighbours are searched
    * @param rInputConditions       List of Conditions to be searched
    * @param rRadius                List of search radius for every Condition
    * @param rResults               Array of results for each Condition
    */
    void SearchConditionsInRadiusExclusive(
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        ) override;

    /**
    * @brief Search neighbours for every Condition in "rInputConditions" including itself
    * @param rStructureConditions   List of Conditions against which the neighbours are searched
    * @param rInputConditions       List of Conditions to be searched
    * @param rRadius                List of search radius for every Condition
    * @param rResults               Array of results for each Condition
    */
    void SearchConditionsInRadiusInclusive(
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        ) override;

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SpecializedSpatialSearch" ;

        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SpecializedSpatialSearch";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {

    }

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

    Parameters mParameters; /// The configuration parameters

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    Parameters GetDefaultParameters() const;

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
    SpecializedSpatialSearch& operator=(SpecializedSpatialSearch const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    SpecializedSpatialSearch(SpecializedSpatialSearch const& rOther)
    {
        *this = rOther;
    }

    ///@}
}; // Class SpecializedSpatialSearch

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template<SpatialContainer TSearchBackend>
inline std::ostream& operator << (std::ostream& rOStream, 
                const SpecializedSpatialSearch<TSearchBackend>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}addtogroup block

}  // namespace Kratos.