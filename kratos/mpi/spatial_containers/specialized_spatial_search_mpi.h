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
#include "mpi/utilities/mpi_search_utilities.h"
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
* @class SpecializedSpatialSearchMPI
* @ingroup KratosCore
* @brief This class is used to search for elements, conditions and nodes in a given model part (MPI version)
* @details In order to perform the search it uses as backend some of the the spatial containers defined `spatial_containers` folder
* @tparam TSearchBackend The spatial container to be used as backend
* @author Vicente Mataix Ferrandiz
*/
template<SpatialContainer TSearchBackend>
class KRATOS_API(KRATOS_CORE) SpecializedSpatialSearchMPI
    : public SpecializedSpatialSearch<TSearchBackend>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpecializedSpatialSearchMPI
    KRATOS_CLASS_POINTER_DEFINITION(SpecializedSpatialSearchMPI);

    /// The base type
    using BaseType = SpecializedSpatialSearch<TSearchBackend>;

    /// Nodes classes
    using NodesContainerType = typename BaseType::NodesContainerType;
    using ResultNodesContainerType = typename BaseType::ResultNodesContainerType;
    using VectorResultNodesContainerType = typename BaseType::VectorResultNodesContainerType;
    using NodeSpatialSearchResultContainerType = typename BaseType::NodeSpatialSearchResultContainerType;
    using NodeSpatialSearchResultContainerMapType = typename SpatialSearch::NodeSpatialSearchResultContainerMapType;

    /// Elements classes
    using ElementsContainerType = typename BaseType::ElementsContainerType;
    using ResultElementsContainerType = typename BaseType::ResultElementsContainerType;
    using VectorResultElementsContainerType = typename BaseType::VectorResultElementsContainerType;
    using ElementSpatialSearchResultContainerType = typename BaseType::ElementSpatialSearchResultContainerType;
    using ElementSpatialSearchResultContainerMapType = typename SpatialSearch::ElementSpatialSearchResultContainerMapType;

    /// Conditions classes
    using ConditionsContainerType = typename BaseType::ConditionsContainerType;
    using ResultConditionsContainerType = typename BaseType::ResultConditionsContainerType;
    using VectorResultConditionsContainerType = typename BaseType::VectorResultConditionsContainerType;
    using ConditionSpatialSearchResultContainerType = typename BaseType::ConditionSpatialSearchResultContainerType;
    using ConditionSpatialSearchResultContainerMapType = typename SpatialSearch::ConditionSpatialSearchResultContainerMapType;

    /// Input/output types
    using RadiusArrayType = typename BaseType::RadiusArrayType;
    using DistanceType = typename BaseType::DistanceType;
    using VectorDistanceType = typename BaseType::VectorDistanceType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SpecializedSpatialSearchMPI()
        : BaseType()
    {
    }

    /**
     * @brief Constructor with parameters
     * @param ThisParameters The parameters to be considered
     */
    SpecializedSpatialSearchMPI(Parameters ThisParameters)
        : BaseType(ThisParameters)
    {
    }

    /// Destructor.
    ~SpecializedSpatialSearchMPI() override = default;

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
     * @brief Search neighbours for every element in "Inputelements" excluding itself
     */
    ElementSpatialSearchResultContainerMapType SearchElementsInRadiusExclusive (
        const ElementsContainerType& rStructureElements,
        const ElementsContainerType& rInputElements,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        ) override;

    //************************************************************************
    // Nodal Exclusive search with distance calculation
    //************************************************************************

    /**
     * @brief Search neighbours for every node in "InputNodes" excluding itself
     */
    NodeSpatialSearchResultContainerMapType SearchNodesInRadiusExclusive (
        const NodesContainerType& rStructureNodes,
        const NodesContainerType& rInputNodes,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        ) override;

    //************************************************************************
    // Conditional Exclusive search with distance calculation
    //************************************************************************

    /**
    * @brief Search neighbours for every Condition in "InputConditions" excluding itself
    */
    ConditionSpatialSearchResultContainerMapType SearchConditionsInRadiusExclusive (
        const ConditionsContainerType& rStructureConditions,
        const ConditionsContainerType& rInputConditions,
        const RadiusArrayType& rRadius,
        const DataCommunicator& rDataCommunicator
        ) override;

    //************************************************************************
    // Point vs Entities (these are new interfaces and already use the new containers)
    //************************************************************************

    /**
     * @brief Search neighbours nodes for one point in a given radius
     */
    void SearchNodesOverPointInRadius (
        const NodesContainerType& rStructureNodes,
        const array_1d<double,3>& rPoint,
        const double Radius,
        NodeSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator
        ) override;

    /**
     * @brief Search neighbours nodes for several points in a given radius
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
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        const auto all_points_distances = MPISearchUtilities::MPISynchronousPointSynchronizationWithDistances(itPointBegin, itPointEnd, all_points_coordinates, rRadius, rDataCommunicator);

        // Perform the corresponding searchs
        NodeSpatialSearchResultContainerMapType results;
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = results.InitializeResult(point);
            SearchNodesOverPointInRadius(rStructureNodes, point, all_points_distances[i_node], r_partial_result, rDataCommunicator);
        }
        return results;
    }

    /**
     * @brief Search nearest neighbour node for one point
     */
    void SearchNodesOverPointNearestPoint (
        const NodesContainerType& rStructureNodes,
        const array_1d<double,3>& rPoint,
        NodeSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator
        ) override;

    /**
     * @brief Search nearest neighbour node for several points
     */
    template<typename TPointIteratorType>    
    NodeSpatialSearchResultContainerMapType SearchNodesOverPointsNearestPoint (
        const NodesContainerType& rStructureNodes,
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const DataCommunicator& rDataCommunicator
        )
    {
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        MPISearchUtilities::MPISynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, rDataCommunicator);

        // Perform the corresponding searchs
        NodeSpatialSearchResultContainerMapType results;
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = results.InitializeResult(point);
            SearchNodesOverPointNearestPoint(rStructureNodes, point, r_partial_result, rDataCommunicator);
        }
        return results;
    }

    /**
     * @brief Search neighbours elements for one point in a given radius
     */
    void SearchElementsOverPointInRadius (
        const ElementsContainerType& rStructureElements,
        const array_1d<double,3>& rPoint,
        const double Radius,
        ElementSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator
        ) override;

    /**
     * @brief Search neighbours elements for several points in a given radius
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
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        const auto all_points_distances = MPISearchUtilities::MPISynchronousPointSynchronizationWithDistances(itPointBegin, itPointEnd, all_points_coordinates, rRadius, rDataCommunicator);

        // Perform the corresponding searchs
        ElementSpatialSearchResultContainerMapType results;
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = results.InitializeResult(point);
            SearchElementsOverPointInRadius(rStructureElements, point, all_points_distances[i_node], r_partial_result, rDataCommunicator);
        }
        return results;
    }

    /**
     * @brief Search nearest neighbour element for one point
     */
    void SearchElementsOverPointNearestPoint (
        const ElementsContainerType& rStructureElements,
        const array_1d<double,3>& rPoint,
        ElementSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator
        ) override;

    /**
     * @brief Search nearest neighbour element for several points
     */
    template<typename TPointIteratorType>    
    ElementSpatialSearchResultContainerMapType SearchElementsOverPointsNearestPoint (
        const ElementsContainerType& rStructureElements,
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const DataCommunicator& rDataCommunicator
        )
    {
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        MPISearchUtilities::MPISynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, rDataCommunicator);

        // Perform the corresponding searchs
        ElementSpatialSearchResultContainerMapType results;
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = results.InitializeResult(point);
            SearchElementsOverPointNearestPoint(rStructureElements, point, r_partial_result, rDataCommunicator);
        }
        return results;
    }

    /**
     * @brief Search neighbours conditions for one point in a given radius
     */
    void SearchConditionsOverPointInRadius (
        const ConditionsContainerType& rStructureConditions,
        const array_1d<double,3>& rPoint,
        const double Radius,
        ConditionSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator
        ) override;

    /**
     * @brief Search neighbours conditions for several points in a given radius
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
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        const auto all_points_distances = MPISearchUtilities::MPISynchronousPointSynchronizationWithDistances(itPointBegin, itPointEnd, all_points_coordinates, rRadius, rDataCommunicator);

        // Perform the corresponding searchs
        ConditionSpatialSearchResultContainerMapType results;
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = results.InitializeResult(point);
            SearchConditionsOverPointInRadius(rStructureConditions, point, all_points_distances[i_node], r_partial_result, rDataCommunicator);
        }
        return results;
    }

    /**
     * @brief Search nearest neighbour condition for one point
     */
    void SearchConditionsOverPointNearestPoint (
        const ConditionsContainerType& rStructureConditions,
        const array_1d<double,3>& rPoint,
        ConditionSpatialSearchResultContainerType& rResults,
        const DataCommunicator& rDataCommunicator
        ) override;

    /**
     * @brief Search nearest neighbour condition for several points
     */
    template<typename TPointIteratorType>    
    ConditionSpatialSearchResultContainerMapType SearchConditionsOverPointsNearestPoint (
        const ConditionsContainerType& rStructureConditions,
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const DataCommunicator& rDataCommunicator
        )
    {        
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        MPISearchUtilities::MPISynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, rDataCommunicator);

        // Perform the corresponding searchs
        ConditionSpatialSearchResultContainerMapType results;
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = results.InitializeResult(point);
            SearchConditionsOverPointNearestPoint(rStructureConditions, point, r_partial_result, rDataCommunicator);
        }
        return results;
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SpecializedSpatialSearchMPI" ;

        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SpecializedSpatialSearchMPI";
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

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method returns the default parameters
     * @return The default parameters
     */
    Kratos::Parameters GetDefaultParameters() const
    {
        return BaseType::GetDefaultParameters();
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SpecializedSpatialSearchMPI& operator=(SpecializedSpatialSearchMPI const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    SpecializedSpatialSearchMPI(SpecializedSpatialSearchMPI const& rOther)
    {
        *this = rOther;
    }

    ///@}
}; // Class SpecializedSpatialSearchMPI

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template<SpatialContainer TSearchBackend>
inline std::ostream& operator << (std::ostream& rOStream,
                const SpecializedSpatialSearchMPI<TSearchBackend>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@}addtogroup block

}  // namespace Kratos.