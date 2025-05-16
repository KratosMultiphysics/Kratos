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
#include "spatial_containers/specialized_spatial_search.h"

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
* @class SpecializedSpatialSearchFactory
* @ingroup KratosCore
* @brief Factory for the specialized spatial search
* @author Vicente Mataix Ferrandiz
*/
class SpecializedSpatialSearchFactory
    : public SpatialSearch
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SpecializedSpatialSearchFactory);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpecializedSpatialSearchFactory()
    {
        Parameters default_parameters = GetDefaultParameters();
        Parameters search_parameters(default_parameters["search_parameters"].WriteJsonString());
        mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearch<SpatialContainer::KDTree>(search_parameters));
    }

    /// Constructor with parameters
    SpecializedSpatialSearchFactory(Parameters ThisParameters)
    {
        Parameters default_parameters = GetDefaultParameters();
        ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
        const std::string& r_container_type = ThisParameters["container_type"].GetString();
        Parameters search_parameters(ThisParameters["search_parameters"].WriteJsonString());
        if (r_container_type == "KDTree" || r_container_type == "kd_tree") {
            mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearch<SpatialContainer::KDTree>(search_parameters));
        } else if (r_container_type == "Octree" || r_container_type == "octree") {
            mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearch<SpatialContainer::Octree>(search_parameters));
        } else if (r_container_type == "BinsStatic" || r_container_type == "bins_static") {
            mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearch<SpatialContainer::BinsStatic>(search_parameters));
        } else if (r_container_type == "BinsDynamic" || r_container_type == "bins_dynamic") {
            mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearch<SpatialContainer::BinsDynamic>(search_parameters));
        } else {
            KRATOS_ERROR << "Unknown container type: " << r_container_type << std::endl;
        }
    }

    /// Destructor.
    ~SpecializedSpatialSearchFactory() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    */
    void SearchElementsInRadiusExclusive(
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override
    {
        mpSpatialSearch->SearchElementsInRadiusExclusive(rStructureElements, rInputElements, rRadius, rResults, rResultsDistance);
    }

    /**
    * @brief Search neighbours for every element in "rInputElements" excluding itself
    */
    void SearchElementsInRadiusExclusive(
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        VectorResultElementsContainerType& rResults
        ) override
    {
        mpSpatialSearch->SearchElementsInRadiusExclusive(rStructureElements, rInputElements, rRadius, rResults);
    }

    /**
    * @brief Search neighbours for every node in "rInputNodes" excluding itself
    */
    void SearchNodesInRadiusExclusive(
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override
    {
        mpSpatialSearch->SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, rRadius, rResults, rResultsDistance);
    }

    /**
    * @brief Search neighbours for every node in "rInputNodes" excluding itself
    */
    void SearchNodesInRadiusExclusive(
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        VectorResultNodesContainerType& rResults
        ) override
    {
        mpSpatialSearch->SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, rRadius, rResults);
    }

    /**
    * @brief Search neighbours for every Condition in "rInputConditions" excluding itself
    */
    void SearchConditionsInRadiusExclusive(
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults,
        VectorDistanceType& rResultsDistance
        ) override
    {
        mpSpatialSearch->SearchConditionsInRadiusExclusive(rStructureConditions, rInputConditions, rRadius, rResults, rResultsDistance);
    }

    /**
    * @brief Search neighbours for every Condition in "rInputConditions" excluding itself
    */
    void SearchConditionsInRadiusExclusive(
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        VectorResultConditionsContainerType& rResults
        ) override
    {
        mpSpatialSearch->SearchConditionsInRadiusExclusive(rStructureConditions, rInputConditions, rRadius, rResults);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SpecializedSpatialSearchFactory" ;

        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SpecializedSpatialSearchFactory";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {

    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    SpatialSearch::Pointer mpSpatialSearch = nullptr; /// The spatial search object

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method returns the default parameters
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const
    {
        Parameters default_parameters = Parameters(R"(
        {   "container_type"  : "KDTree",
            "search_parameters" : {
                "bucket_size"     : 4,
                "allocation_size" : 1000
            }
        })" );

        return default_parameters;
    }

    ///@}
}; // Class SpecializedSpatialSearchFactory

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const SpecializedSpatialSearchFactory& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.