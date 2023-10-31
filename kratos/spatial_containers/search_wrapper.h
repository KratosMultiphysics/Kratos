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
#include "includes/kratos_parameters.h"
#include "utilities/search_utilities.h"
#include "spatial_containers/geometrical_objects_bins.h"
#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/spatial_search_result_container.h"
#include "spatial_containers/spatial_search_result_container_vector.h"

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
 */
template<class TSearchObject>
class KRATOS_API(KRATOS_CORE) SearchWrapper
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SearchWrapper
    KRATOS_CLASS_POINTER_DEFINITION(SearchWrapper);

    /// The type of geometrical object to be stored in the bins
    using ObjectType = typename TSearchObject::ObjectType;

    /// The result type definition
    using ResultType = SpatialSearchResult<ObjectType>;

    /// Search containers
    using ResultContainerType = SpatialSearchResultContainer<ObjectType>;
    using ResultContainerVectorType = SpatialSearchResultContainerVector<ObjectType>;

    /// Defining the point type for the search
    using PointType = typename TSearchObject::PointType;
    using PointVector = std::vector<typename PointType::Pointer>;

    /// Trees types definitions
    using DistanceVector = std::vector<double>;

    /// Some constexpr flags
    static constexpr bool IsGeometricalObjectBins = std::is_same_v<TSearchObject, GeometricalObjectsBins>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    SearchWrapper() = delete;

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param rGeometricalObjectsVector The geometries to be stored
     * @param rDataCommunicator The data communicator
     * @param Settings The settings of the search
     * @tparam TContainer The container type
     */
    template<typename TContainer>
    SearchWrapper(
        TContainer& rGeometricalObjectsVector,
        const DataCommunicator& rDataCommunicator,
        Parameters Settings = Parameters(R"({})")
        ) : mrDataCommunicator(rDataCommunicator),
            mSettings(Settings)
    {
        // Validate and assign defaults
        mSettings.ValidateAndAssignDefaults(GetDefaultParameters());

        // Create base search object
        if constexpr (IsGeometricalObjectBins) {
            mpSearchObject = Kratos::make_shared<TSearchObject>(rGeometricalObjectsVector.begin(), rGeometricalObjectsVector.end());
        } else { // Otherwise we assume it will be some kind of spatial search
            // Defining the PointVector
            mpPointVector = Kratos::make_unique<PointVector>(SearchUtilities::PreparePointsSearch(rGeometricalObjectsVector));
            
            // Retrieving parameters
            const int bucket_size = mSettings["bucket_size"].GetInt();

            // Create the search
            mpSearchObject = Kratos::make_shared<TSearchObject>(mpPointVector->begin(), mpPointVector->end(), bucket_size);
        }

    #ifdef KRATOS_USING_MPI
        // Set up the global bounding boxes
        InitializeGlobalBoundingBoxes();
    #endif
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
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<class TPointIteratorType>
    void SearchInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true
        )
    {
#ifdef KRATOS_USING_MPI
    // Call distributed version
    if (mrDataCommunicator.IsDistributed()) {
        DistributedSearchInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution);
    } else { // Call serial version
        SerialSearchInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution);
    }
#else // Calling just serial method
    SerialSearchInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution);
#endif
    }

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
    template<class TPointIteratorType>
    void SearchNearestInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true
        )
    {
#ifdef KRATOS_USING_MPI
    // Call distributed version
    if (mrDataCommunicator.IsDistributed()) {
        DistributedSearchNearestInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution);
    } else { // Call serial version
        SerialSearchNearestInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution);
    }
#else // Calling just serial method
    SerialSearchNearestInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution);
#endif
    }

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
    template<class TPointIteratorType>
    void SearchNearest(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true
        )
    {
#ifdef KRATOS_USING_MPI
    // Call distributed version
    if (mrDataCommunicator.IsDistributed()) {
        DistributedSearchNearest(itPointBegin, itPointEnd, rResults, ClearSolution);
    } else { // Call serial version
        SerialSearchNearest(itPointBegin, itPointEnd, rResults, ClearSolution);
    }
#else // Calling just serial method
    SerialSearchNearest(itPointBegin, itPointEnd, rResults, ClearSolution);
#endif
    }

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
    template<class TPointIteratorType>
    void SearchIsInside(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true
        )
    {
#ifdef KRATOS_USING_MPI
    // Call distributed version
    if (mrDataCommunicator.IsDistributed()) {
        DistributedSearchIsInside(itPointBegin, itPointEnd, rResults, ClearSolution);
    } else { // Call serial version
        SerialSearchIsInside(itPointBegin, itPointEnd, rResults, ClearSolution);
    }
#else // Calling just serial method
    SerialSearchIsInside(itPointBegin, itPointEnd, rResults, ClearSolution);
#endif
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
        mpSearchObject->PrintData(rOStream);
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

    typename TSearchObject::Pointer mpSearchObject;            /// The pointer to the base search considered
    Kratos::unique_ptr<PointVector> mpPointVector =  nullptr;  /// The point vector considered in the search trees
    std::vector<double> mGlobalBoundingBoxes;                  /// All the global BB, data is xmax, xmin,  ymax, ymin,  zmax, zmin
    const DataCommunicator& mrDataCommunicator;                /// The data communicator
    Parameters mSettings;                                      /// The settings considered

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
     * @brief This method takes a point and finds all of the objects in the given radius to it (iterative version) (serial version).
     * @details The result contains the object and also its distance to the point.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<class TPointIteratorType>
    void SerialSearchInRadius(
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

        // Retrieving parameters
        const int allocation_size = mSettings["allocation_size"].GetInt();

        // Adding the results to the container
        std::size_t counter = 0, id = 0;
        // TODO: Use ParallelUtilities
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++) {
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value || std::is_same<TPointIteratorType, ModelPart::NodeConstantIterator>::value) {
                id = it_point->Id();
            } else {
                id = counter;
            }
            auto& r_partial_result = rResults.InitializeResult(id);

            // Search
            std::vector<ResultType> results;
            LocalSearchInRadius(*it_point, Radius, results, allocation_size);
            for (auto& r_result : results) {
                r_partial_result.AddResult(r_result);
            }

            // Synchronize
            r_partial_result.SynchronizeAll(mrDataCommunicator);

            // Update counter
            ++counter;
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius (iterative version) (serial version).
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param Radius The radius to be checked
     * @param rResults The result of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<class TPointIteratorType>
    void SerialSearchNearestInRadius(
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

        // Retrieving parameters
        const int allocation_size = mSettings["allocation_size"].GetInt();

        // Adding the results to the container
        std::size_t counter = 0, id = 0;
        // TODO: Use ParallelUtilities
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++) {
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value || std::is_same<TPointIteratorType, ModelPart::NodeConstantIterator>::value) {
                id = it_point->Id();
            } else {
                id = counter;
            }
            auto& r_point_result = rResults.InitializeResult(id);
            ResultType result;
            LocalSearchNearestInRadius(*it_point, Radius, result, allocation_size);
            r_point_result.AddResult(result);
            
            // Synchronize
            r_point_result.SynchronizeAll(mrDataCommunicator);

            // Update counter
            ++counter;
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it (iterative version)  (serial version).
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param rResults The result of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<class TPointIteratorType>
    void SerialSearchNearest(
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

        // Adding the results to the container
        std::size_t counter = 0, id = 0;
        // TODO: Use ParallelUtilities
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++) {
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value || std::is_same<TPointIteratorType, ModelPart::NodeConstantIterator>::value) {
                id = it_point->Id();
            } else {
                id = counter;
            }
            auto& r_point_result = rResults.InitializeResult(id);
            ResultType result;
            LocalSearchNearest(*it_point, result);
            r_point_result.AddResult(result);
            
            // Synchronize
            r_point_result.SynchronizeAll(mrDataCommunicator);
            
            // Update counter
            ++counter;
        }
    }

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain (iterative version)  (serial version).
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param rResults The result of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<class TPointIteratorType>
    void SerialSearchIsInside(
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

        // Adding the results to the container
        std::size_t counter = 0, id = 0;
        // TODO: Use ParallelUtilities
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++) {
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value || std::is_same<TPointIteratorType, ModelPart::NodeConstantIterator>::value) {
                id = it_point->Id();
            } else {
                id = counter;
            }
            auto& r_point_result = rResults.InitializeResult(id);
            ResultType result;
            LocalSearchIsInside(*it_point, result);
            r_point_result.AddResult(result);

            // Synchronize
            r_point_result.SynchronizeAll(mrDataCommunicator);

            // Update counter
            ++counter;
        }
    }

#ifdef KRATOS_USING_MPI

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it (MPI version).
     * @details The result contains the object and also its distance to the point.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<class TPointIteratorType>
    void DistributedSearchInRadius(
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

        // Retrieving parameters
        const int allocation_size = mSettings["allocation_size"].GetInt();

        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        std::vector<IndexType> all_points_ids;
        SearchUtilities::SynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, all_points_ids, mrDataCommunicator);

        // The local bounding box
        const auto& r_local_bb = mpSearchObject->GetBoundingBox(); 

        // Perform the corresponding searches
        const std::size_t total_number_of_points = all_points_coordinates.size()/3;
        // TODO: Use ParallelUtilities
        for (std::size_t i_point = 0; i_point < total_number_of_points; ++i_point) {
            const Point point(all_points_coordinates[i_point * 3 + 0], all_points_coordinates[i_point * 3 + 1], all_points_coordinates[i_point * 3 + 2]);
            auto& r_partial_result = rResults.InitializeResult(all_points_ids[i_point]);

            // Check if the point is inside the set
            if (SearchUtilities::PointIsInsideBoundingBox(r_local_bb, point, Radius)) {
                // Search
                std::vector<ResultType> results;
                LocalSearchInRadius(point, Radius, results, allocation_size);
                for (auto& r_result : results) {
                    r_partial_result.AddResult(r_result);
                }
            }

            // Synchronize
            r_partial_result.SynchronizeAll(mrDataCommunicator);
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius (MPI version).
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
    template<class TPointIteratorType>
    void DistributedSearchNearestInRadius(
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

        // Retrieving parameters
        const int allocation_size = mSettings["allocation_size"].GetInt();

        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        std::vector<IndexType> all_points_ids;
        SearchUtilities::SynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, all_points_ids, mrDataCommunicator);

        // Get the rank
        const int current_rank = GetRank();

        // The local bounding box
        const auto& r_local_bb = mpSearchObject->GetBoundingBox(); 

        // Perform the corresponding searches
        const std::size_t total_number_of_points = all_points_coordinates.size()/3;
        // TODO: Use ParallelUtilities
        for (std::size_t i_point = 0; i_point < total_number_of_points; ++i_point) {
            // Perform local search
            const Point point(all_points_coordinates[i_point * 3 + 0], all_points_coordinates[i_point * 3 + 1], all_points_coordinates[i_point * 3 + 2]);
            auto& r_partial_result = rResults.InitializeResult(all_points_ids[i_point]);
            
            // Result of search
            ResultType local_result;

            // Check if the point is inside the set
            if (SearchUtilities::PointIsInsideBoundingBox(r_local_bb, point, Radius)) {
                // Call local search
                LocalSearchNearestInRadius(point, Radius, local_result, allocation_size);
            }

            /* Now sync results between partitions */

            // Get the distance
            const double local_distance = (local_result.IsObjectFound() && local_result.IsDistanceCalculated()) ? local_result.GetDistance() : std::numeric_limits<double>::max();

            // Find the minimum value and the rank that holds it
            const auto global_min = mrDataCommunicator.MinLocAll(local_distance);

            // Get the solution from the computed_rank
            if (global_min.second == current_rank) {
                // Add the local search
                r_partial_result.AddResult(local_result);
            }

            // Synchronize
            r_partial_result.SynchronizeAll(mrDataCommunicator);
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it (MPI version).
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param rResults The results of the search
     * @param ClearSolution Clear the current solution
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<class TPointIteratorType>
    void DistributedSearchNearest(
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

        // Retrieving parameters
        const int allocation_size = mSettings["allocation_size"].GetInt();
        
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        std::vector<IndexType> all_points_ids;
        SearchUtilities::SynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, all_points_ids, mrDataCommunicator);

        // Get the rank
        const int current_rank = GetRank();

        // Get the maximum radius
        const auto bb = GetBoundingBox();
        const array_1d<double, 3> box_size = bb.GetMaxPoint() - bb.GetMinPoint();
        const double max_radius= *std::max_element(box_size.begin(), box_size.end());

        // The local bounding box
        const auto& r_local_bb = mpSearchObject->GetBoundingBox(); 

        // Perform the corresponding searches
        const std::size_t total_number_of_points = all_points_coordinates.size()/3;
        // TODO: Use ParallelUtilities
        for (std::size_t i_point = 0; i_point < total_number_of_points; ++i_point) {
            // Perform local search
            const Point point(all_points_coordinates[i_point * 3 + 0], all_points_coordinates[i_point * 3 + 1], all_points_coordinates[i_point * 3 + 2]);
            auto& r_partial_result = rResults.InitializeResult(all_points_ids[i_point]);

            // Check if the point is inside the set
            ResultType local_result;
            if (SearchUtilities::PointIsInsideBoundingBox(r_local_bb, point, max_radius)) {
                // Call local search
                LocalSearchNearestInRadius(point, max_radius, local_result, allocation_size);
            }

            /* Now sync results between partitions */

            // Get the distance
            const double local_distance = (local_result.IsObjectFound() && local_result.IsDistanceCalculated()) ? local_result.GetDistance() : std::numeric_limits<double>::max();

            // Find the minimum value and the rank that holds it
            const auto global_min = mrDataCommunicator.MinLocAll(local_distance);

            // Get the solution from the computed_rank
            if (global_min.second == current_rank) {
                // Add the local search
                r_partial_result.AddResult(local_result);
            }

            // Synchronize
            r_partial_result.SynchronizeAll(mrDataCommunicator);
        }
    }

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain (iterative version) (MPI version).
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
    template<class TPointIteratorType>
    void DistributedSearchIsInside(
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
        
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        std::vector<IndexType> all_points_ids;
        SearchUtilities::SynchronousPointSynchronization(itPointBegin, itPointEnd, all_points_coordinates, all_points_ids, mrDataCommunicator);

        // Get the rank
        const int current_rank = GetRank();

        // The local bounding box
        const auto& r_local_bb = mpSearchObject->GetBoundingBox(); 

        // Perform the corresponding searches
        const std::size_t total_number_of_points = all_points_coordinates.size()/3;
        // TODO: Use ParallelUtilities
        for (std::size_t i_point = 0; i_point < total_number_of_points; ++i_point) {
            // Perform local search
            const Point point(all_points_coordinates[i_point * 3 + 0], all_points_coordinates[i_point * 3 + 1], all_points_coordinates[i_point * 3 + 2]);
            auto& r_partial_result = rResults.InitializeResult(all_points_ids[i_point]);

            // Check if the point is inside the set
            int computed_rank = std::numeric_limits<int>::max();
            ResultType local_result;
            if (SearchUtilities::PointIsInsideBoundingBox(r_local_bb, point)) {
                // Call local search
                LocalSearchIsInside(point, local_result);

                // Set current rank
                computed_rank = current_rank;
            }

            // Now sync results between partitions
            //KRATOS_WARNING("SearchWrapper.SearchIsInside") << "This assumes that first one of the viable results is the closest one. The algorithm  sets distance to 0, so it is ambiguous which is the closest one." << std::endl;

            // Min rank of all ranks
            computed_rank = mrDataCommunicator.MinAll(computed_rank);

            // Get the solution from the computed_rank
            if (computed_rank == current_rank) {
                // Add the local search
                r_partial_result.AddResult(local_result);
            }

            // Synchronize
            r_partial_result.SynchronizeAll(mrDataCommunicator);
        }
    }

#endif

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it.
     * @details The result contains the object and also its distance to the point.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @param AllocationSize The allocation size of the results
     */
    void LocalSearchInRadius(
        const PointType& rPoint,
        const double Radius,
        std::vector<ResultType>& rResults,
        const int AllocationSize = 1000
        );

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius.
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResult The result of the search
     * @param AllocationSize The allocation size of the results
     */
    void LocalSearchNearestInRadius(
        const PointType& rPoint,
        const double Radius, 
        ResultType& rResult,
        const int AllocationSize = 1000
        );

    /**
     * @brief This method takes a point and finds the nearest object to it.
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param rResult The result of the search
    */
    void LocalSearchNearest(
        const PointType& rPoint, 
        ResultType& rResult
        );

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain.
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param rPoint The point to be checked
     * @param rResult The result of the search
     */
    void LocalSearchIsInside(
        const PointType& rPoint, 
        ResultType& rResult
        );

    /**
     * @brief This method returns the default parameters of the search
     * @return The default parameters
     */
    const Parameters GetDefaultParameters() const;

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