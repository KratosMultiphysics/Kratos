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

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
* @class SpecializedSpatialSearch
* @ingroup KratosCore
* @brief This class is used to search for elements, conditions and nodes in a given model part
* @details In order to perform the search it uses as backend some of the the spatial containers defined `spatial_containers` folder
* @author Vicente Mataix Ferrandiz
*/
template<class TSearhcBackend>
class KRATOS_API(KRATOS_CORE) SpecializedSpatialSearch
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpecializedSpatialSearch
    KRATOS_CLASS_POINTER_DEFINITION(SpecializedSpatialSearch);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpecializedSpatialSearch()
    {
    }

    /// Destructor.
    virtual ~SpecializedSpatialSearch(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Search neighbours for every element in "Inputelements" excluding itself
    * @param StructureElements   List of elements against which the neighbours are searched
    * @param InputElements       List of elements to be searched
    * @param Radius              List of search radius for every element
    * @param rResults            Array of results for each element
    * @param rResultDistance     Array of distances for each result of each element
    */
    void SearchElementsInRadiusExclusive (
        ElementsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchElementsInRadiusInclusive (
        ElementsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchElementsInRadiusExclusive (
        ElementsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchElementsInRadiusInclusive (
        ElementsContainerType const& StructureElements,
        ElementsContainerType const& InputElements,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchNodesInRadiusExclusive (
        NodesContainerType const& StructureNodes,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchNodesInRadiusInclusive (
        NodesContainerType const& StructureNodes,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchNodesInRadiusExclusive (
        NodesContainerType const& StructureNodes,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchNodesInRadiusInclusive (
        NodesContainerType const& StructureNodes,
        NodesContainerType const& InputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchConditionsInRadiusExclusive (
        ConditionsContainerType const& StructureConditions,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchConditionsInRadiusInclusive (
        ConditionsContainerType const& StructureConditions,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchConditionsInRadiusExclusive (
        ConditionsContainerType const& StructureConditions,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

    void SearchConditionsInRadiusInclusive (
        ConditionsContainerType const& StructureConditions,
        ConditionsContainerType const& InputConditions,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        ) override
    {
        KRATOS_ERROR << "Direct call of an abstract method" << std::endl;
    }

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
        buffer << "SpecializedSpatialSearch" ;

        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SpecializedSpatialSearch";}

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

    TSearhcBackend mSearchStructure; /// The search structure

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
inline std::ostream& operator << (std::ostream& rOStream, 
                const SpecializedSpatialSearch& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}addtogroup block

}  // namespace Kratos.