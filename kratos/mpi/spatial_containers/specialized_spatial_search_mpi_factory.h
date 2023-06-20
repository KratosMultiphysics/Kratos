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
#include "mpi/spatial_containers/specialized_spatial_search_mpi.h"

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
* @class SpecializedSpatialSearchMPIFactory
* @ingroup KratosCore
* @brief Factory for the specialized spatial search
* @author Vicente Mataix Ferrandiz
*/
class SpecializedSpatialSearchMPIFactory 
    : public SpatialSearch
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SpecializedSpatialSearchMPIFactory);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpecializedSpatialSearchMPIFactory()
    {
        Parameters default_parameters = GetDefaultParameters();
        Parameters search_parameters(default_parameters["search_parameters"].WriteJsonString());
        mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearchMPI<SpatialContainer::KDTree>(search_parameters));
    }

    /// Constructor with parameters
    SpecializedSpatialSearchMPIFactory(Parameters ThisParameters)
    {
        Parameters default_parameters = GetDefaultParameters();
        ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
        const std::string& r_container_type = ThisParameters["container_type"].GetString();
        Parameters search_parameters(ThisParameters["search_parameters"].WriteJsonString());
        if (r_container_type == "KDTree" || r_container_type == "kd_tree") {
            mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearchMPI<SpatialContainer::KDTree>(search_parameters));
        } else if (r_container_type == "Octree" || r_container_type == "octree") {
            mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearchMPI<SpatialContainer::Octree>(search_parameters));
        } else if (r_container_type == "BinsStatic" || r_container_type == "bins_static") {
            mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearchMPI<SpatialContainer::BinsStatic>(search_parameters));
        } else if (r_container_type == "BinsDynamic" || r_container_type == "bins_dynamic") {
            mpSpatialSearch = SpatialSearch::Pointer(new SpecializedSpatialSearchMPI<SpatialContainer::BinsDynamic>(search_parameters));
        } else {
            KRATOS_ERROR << "Unknown container type: " << r_container_type << std::endl;
        }
    }

    /// Destructor.
    ~SpecializedSpatialSearchMPIFactory() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Search neighbours nodes for one point in a given radius
     */
    void SearchNodesOverPointInRadius (
        const NodesContainerType& rStructureNodes,
        const array_1d<double,3>& rPoint,
        const double Radius,
        NodeSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        ) override
    {
        mpSpatialSearch->SearchNodesOverPointInRadius(rStructureNodes, rPoint, Radius, rResults, rDataCommunicator, SyncronizeResults);
    }

    /**
     * @brief Search nearest neighbour node for one point
     */
    void SearchNodesOverPointNearestPoint (
        const NodesContainerType& rStructureNodes,
        const array_1d<double,3>& rPoint,
        NodeSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        ) override
    {
        mpSpatialSearch->SearchNodesOverPointNearestPoint(rStructureNodes, rPoint, rResults, rDataCommunicator, SyncronizeResults);
    }

    /**
     * @brief Search neighbours elements for one point in a given radius
     */
    void SearchElementsOverPointInRadius (
        const ElementsContainerType& rStructureElements,
        const array_1d<double,3>& rPoint,
        const double Radius,
        ElementSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        ) override
    {
        mpSpatialSearch->SearchElementsOverPointInRadius(rStructureElements, rPoint, Radius, rResults, rDataCommunicator, SyncronizeResults);
    }

    /**
     * @brief Search nearest neighbour element for one point
     */
    void SearchElementsOverPointNearestPoint (
        const ElementsContainerType& rStructureElements,
        const array_1d<double,3>& rPoint,
        ElementSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        ) override
    {
        mpSpatialSearch->SearchElementsOverPointNearestPoint(rStructureElements, rPoint, rResults, rDataCommunicator, SyncronizeResults);
    }

    /**
     * @brief Search neighbours conditions for one point in a given radius
     */
    void SearchConditionsOverPointInRadius (
        const ConditionsContainerType& rStructureConditions,
        const array_1d<double,3>& rPoint,
        const double Radius,
        ConditionSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        ) override
    {
        mpSpatialSearch->SearchConditionsOverPointInRadius(rStructureConditions, rPoint, Radius, rResults, rDataCommunicator, SyncronizeResults);
    }

    /**
     * @brief Search nearest neighbour condition for one point
     */
    void SearchConditionsOverPointNearestPoint (
        const ConditionsContainerType& rStructureConditions,
        const array_1d<double,3>& rPoint,
        ConditionSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator,
        const bool SyncronizeResults = true
        ) override
    {
        mpSpatialSearch->SearchConditionsOverPointNearestPoint(rStructureConditions, rPoint, rResults, rDataCommunicator, SyncronizeResults);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SpecializedSpatialSearchMPIFactory" ;

        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SpecializedSpatialSearchMPIFactory";
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
}; // Class SpecializedSpatialSearchMPIFactory

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const SpecializedSpatialSearchMPIFactory& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@}addtogroup block

}  // namespace Kratos.