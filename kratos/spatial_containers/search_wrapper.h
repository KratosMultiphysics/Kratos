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
#include <numeric>

// External includes

// Project includes
#include "includes/geometrical_object.h"
#ifdef KRATOS_USING_MPI
#include "mpi/utilities/mpi_search_utilities.h"
#endif
#include "spatial_containers/geometrical_objects_bins.h"
#include "spatial_containers/spatial_search_result_container.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class SearchWrapper
 * @ingroup KratosCore
 * @brief This is a search wrapper ready for MPI searches
 * @details Must be adapted and specialized for every search object
 * @author Vicente Mataix Ferrandiz
 * @tparam TSearchObject The seach object considered
 * @tparam TObjectType The objec type considered
 */
template<class TSearchObject, class TObjectType = GeometricalObject>
class KRATOS_API(KRATOS_CORE) SearchWrapper
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SearchWrapper
    KRATOS_CLASS_POINTER_DEFINITION(SearchWrapper);

    /// The buffer type definition
    using BufferTypeDouble = std::vector<std::vector<double>>;
    using BufferTypeChar = std::vector<std::vector<char>>;

    /// The type of geometrical object to be stored in the bins
    using CellType = GeometricalObjectsBins::CellType;
    using ResultType = GeometricalObjectsBins::ResultType;

    // Search containers
    using ResultContainerType = SpatialSearchResultContainer<GeometricalObject>;
    using ResultContainerVectorType = SpatialSearchResultContainerVector<GeometricalObject>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    SearchWrapper() = delete;

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param GeometricalObjectsBegin The begin iterator of the geometries to be stored
     * @param GeometricalObjectsEnd The end iterator of the geometries to be stored
     * @param rDataCommunicator The data communicator
     * @tparam TIteratorType The type of the iterator
     */
    template<typename TIteratorType>
    SearchWrapper(
        TIteratorType GeometricalObjectsBegin,
        TIteratorType GeometricalObjectsEnd,
        const DataCommunicator& rDataCommunicator
        ) : mrDataCommunicator(rDataCommunicator)
    {
        // Create base search object
        if constexpr (std::is_same_v<TSearchObject, GeometricalObjectsBins>) {
            mSearchObject = Kratos::make_shared<TSearchObject>(GeometricalObjectsBegin, GeometricalObjectsEnd);
        }

        // Set up the global bounding boxes
        InitializeGlobalBoundingBoxes();
    }

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param rGeometricalObjectsVector The geometries to be stored
     * @param rDataCommunicator The data communicator
     * @tparam TContainer The container type
     */
    template<typename TContainer>
    SearchWrapper(
        TContainer& rGeometricalObjectsVector,
        const DataCommunicator& rDataCommunicator
        ) : SearchWrapper(rGeometricalObjectsVector.begin(), rGeometricalObjectsVector.end(), rDataCommunicator)
    {
    }

    /// Destructor.
    ~SearchWrapper() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Getting the bins bounding box
     * @return The bounding box of the bins
     */
    BoundingBox<Point> GetBoundingBox() const;

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it.
     * @details The result contains the object and also its distance to the point.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
     */
    void SearchInRadius(
        const Point& rPoint,
        const double Radius,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it.
     * @details The result contains the object and also its distance to the point.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    void SearchInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true
        )
    {
        // Clear current solution
        if (ClearSolution) {
            rResults.Clear();
        }

        // TODO: Check if TPointIteratorType is node interator and check the PARTITION_INDEX

        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        MPISearchUtilities::MPISynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, mrDataCommunicator);

        // Perform the corresponding searches
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = rResults.InitializeResult(i_node);
            this->SearchInRadius(point, Radius, r_partial_result);
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius.
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
     */
    void SearchNearestInRadius(
        const Point& rPoint,
        const double Radius,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius.
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    void SearchNearestInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true
        )
    {
        // Clear current solution
        if (ClearSolution) {
            rResults.Clear();
        }

        // TODO: Check if TPointIteratorType is node interator and check the PARTITION_INDEX

        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        MPISearchUtilities::MPISynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, mrDataCommunicator);

        // Perform the corresponding searches
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            // Perform local search
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = rResults.InitializeResult(i_node);
            this->SearchNearestInRadius(point, Radius, r_partial_result);
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it.
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
    */
    void SearchNearest(
        const Point& rPoint,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and finds the nearest object to it.
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param rResults The results of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    void SearchNearest(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true
        )
    {
        // Clear current solution
        if (ClearSolution) {
            rResults.Clear();
        }
        
        // TODO: Check if TPointIteratorType is node interator and check the PARTITION_INDEX

        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        MPISearchUtilities::MPISynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, mrDataCommunicator);

        // Perform the corresponding searches
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            // Perform local search
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = rResults.InitializeResult(i_node);
            this->SearchNearest(point, r_partial_result);
        }
    }

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain.
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param rPoint The point to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
     */
    void SearchIsInside(
        const Point& rPoint,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain (iterative version).
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param rResults The results of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    void SearchIsInside(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true
        )
    {
        // Clear current solution
        if (ClearSolution) {
            rResults.Clear();
        }
        
        // TODO: Check if TPointIteratorType is node interator and check the PARTITION_INDEX

        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        MPISearchUtilities::MPISynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, mrDataCommunicator);

        // Perform the corresponding searches
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            // Perform local search
            const Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            auto& r_partial_result = rResults.InitializeResult(i_node);
            this->SearchIsInside(point, r_partial_result);
        }
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
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "SearchWrapper" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SearchWrapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Local Search for Rank " << GetRank() + 1 << "/" << GetWorldSize() << "\n";
        BaseType::PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    TSearchObject::Pointer mSearchObject;       /// The pointer to the base search considered

    std::vector<double> mGlobalBoundingBoxes;   /// All the global BB, data is xmax, xmin,  ymax, ymin,  zmax, zmin

    const DataCommunicator& mrDataCommunicator; /// The data communicator

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Returns the current rank
     * @return The current rank
     */
    int GetRank() const;

    /**
     * @brief Returns the world size
     * @return The world size
     */
    int GetWorldSize() const;

    /**
     * @brief Initializes the global bounding boxes
     */
    void InitializeGlobalBoundingBoxes();

    /**
     * @brief This method checks if a point is inside any bounding box of the global bounding boxes
     * @param rCoords The coordinates of the point
     * @return The index of the ranks inside the bounding box
     */
    std::vector<int> RansksPointIsInsideBoundingBox(const array_1d<double, 3>& rCoords);

    /**
     * @brief This method checks if a point is inside any bounding box of the global bounding boxes considering a certain tolerance
     * @param rCoords The coordinates of the point
     * @param Tolerance The tolerance
     * @return The index of the ranks inside the bounding box
     */
    std::vector<int> RansksPointIsInsideBoundingBoxWithTolerance(
        const array_1d<double, 3>& rCoords,
        const double Tolerance
        );

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it (serial version).
     * @details The result contains the object and also its distance to the point.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
     */
    void SerialSearchInRadius(
        const Point& rPoint,
        const double Radius,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius (serial version).
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
     */
    void SerialSearchNearestInRadius(
        const Point& rPoint,
        const double Radius,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and finds the nearest object to it (serial version).
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
    */
    void SerialSearchNearest(
        const Point& rPoint,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain (serial version).
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param rPoint The point to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
     */
    void SerialSearchIsInside(
        const Point& rPoint,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );


    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it (MPI version).
     * @details The result contains the object and also its distance to the point.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
     */
    void DistributedSearchInRadius(
        const Point& rPoint,
        const double Radius,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius (MPI version).
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
     */
    void DistributedSearchNearestInRadius(
        const Point& rPoint,
        const double Radius,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and finds the nearest object to it (MPI version).
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
    */
    void DistributedSearchNearest(
        const Point& rPoint,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain (MPI version).
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param rPoint The point to be checked
     * @param rResults The results of the search
     * @param SyncronizeResults If the results should be synchronized or not
     */
    void DistributedSearchIsInside(
        const Point& rPoint,
        ResultTypeContainer& rResults,
        const bool SyncronizeResults = true
        );

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator deleted.
    SearchWrapper& operator=(SearchWrapper const& rOther) = delete;

    /// Copy constructor deleted.
    SearchWrapper(SearchWrapper const& rOther) = delete;

    ///@}

}; // Class SearchWrapper

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@} addtogroup block

}  // namespace Kratos.