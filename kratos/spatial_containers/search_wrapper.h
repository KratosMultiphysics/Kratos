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
#include "utilities/timer.h"
#include "utilities/search_utilities.h"
#include "utilities/parallel_utilities.h"
#include "spatial_containers/geometrical_objects_bins.h"
#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/spatial_search_result_container.h"
#include "spatial_containers/spatial_search_result_container_vector.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SearchWrapper
 * @ingroup KratosCore
 * @brief This is a search wrapper ready for MPI searches
 * @details Must be adapted and specialized for every search object
 * @author Vicente Mataix Ferrandiz
 * @tparam TSearchObject The seach object considered
 * @tparam TSpatialSearchCommunication The communication type to be used
 */
template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication = SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS>
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

    /// If considering the global data communicator
    static constexpr bool ConsiderGlobalDataCommunicator = TSpatialSearchCommunication == SpatialSearchCommunication::SYNCHRONOUS_HOMOGENEOUS;

    /// Search containers
    using ResultContainerType = SpatialSearchResultContainer<ObjectType, TSpatialSearchCommunication>;
    using ResultContainerVectorType = SpatialSearchResultContainerVector<ObjectType, TSpatialSearchCommunication>;

    /// Defining the point type for the search
    using PointType = typename TSearchObject::PointType;
    using PointVector = std::vector<typename PointType::Pointer>;

    /// Trees types definitions
    using DistanceVector = std::vector<double>;

    /// Define Zero tolerance as epsilon
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Some constexpr flags
    static constexpr bool IsGeometricalObjectBins = std::is_same_v<TSearchObject, GeometricalObjectsBins>;
    static constexpr bool IsDynamicBins = std::is_same_v<TSearchObject, BinsDynamic<3ul, PointType, PointVector>>;

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
        // Checking we are using a global communicator
        KRATOS_ERROR_IF(mrDataCommunicator.IsNullOnThisRank()) << "The data communicator is null on this rank. Try to use a global communicator" << std::endl;

        // Validate and assign defaults
        mSettings.ValidateAndAssignDefaults(GetDefaultParameters());

        // Create base search object
        if constexpr (IsGeometricalObjectBins) {
            mpSearchObject = Kratos::make_shared<TSearchObject>(rGeometricalObjectsVector.begin(), rGeometricalObjectsVector.end());
        } else { // Otherwise we assume it will be some kind of spatial search
            // Prepare the point vector if geometrical objects are provided
            if (rGeometricalObjectsVector.size() > 0) {
                // Defining the PointVector
                const int rank = mrDataCommunicator.IsDistributed() ? mrDataCommunicator.Rank() : -1;
                const auto preprocessed_points = SearchUtilities::PreparePointsSearch(rGeometricalObjectsVector, rank);
                // Check that is greater tha zero (may happen that is empty due to non-local nodes container)
                if (preprocessed_points.size() > 0) {
                    mpPointVector = Kratos::make_unique<PointVector>(preprocessed_points);

                    // Create the search object
                    if constexpr (!IsDynamicBins) {
                        const int bucket_size = mSettings["bucket_size"].GetInt();
                        mpSearchObject = Kratos::make_shared<TSearchObject>(mpPointVector->begin(), mpPointVector->end(), bucket_size);
                    } else {
                        mpSearchObject = Kratos::make_shared<TSearchObject>(mpPointVector->begin(), mpPointVector->end());
                    }
                }
            }
        }

        // Set up the global bounding boxes
        InitializeGlobalBoundingBoxes();
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
     * @brief Getting the bins global bounding box
     * @return The global bounding box of the bins
     */
    BoundingBox<Point> GetGlobalBoundingBox() const;

    /**
     * @brief This method takes several points and finds all of the objects in the given radius to it.
     * @details The result contains the object and also its distance to the point.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param Radius The radius to be checked.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @tparam TPointIteratorType The type of the point iterator.
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
        // Clear current solution
        if (ClearSolution) {
            rResults.Clear();
        }

        // Retrieving parameters
        const int allocation_size = mSettings["allocation_size"].GetInt();

        // The local bounding box
        const auto& r_local_bb = mpSearchObject ? mpSearchObject->GetBoundingBox() : BoundingBox<PointType>();

        // Prepare MPI search
        Timer::Start("SearchWrapper::SynchronousPointSynchronizationWithBoundingBox");
        DistributedSearchInformation search_info;
        SearchUtilities::SynchronousPointSynchronizationWithBoundingBox(itPointBegin, itPointEnd, search_info, r_local_bb, Radius, mrDataCommunicator, true);
        Timer::Stop("SearchWrapper::SynchronousPointSynchronizationWithBoundingBox");

        // Initialize results
        Timer::Start("SearchWrapper::PrepareResultsInProperRanks");
        PrepareResultsInProperRanks(rResults, search_info);
        Timer::Stop("SearchWrapper::PrepareResultsInProperRanks");

        // Perform the corresponding searches
        Timer::Start("SearchWrapper::Search");
        const std::size_t total_number_of_points = search_info.LocalIndices.size();
        struct TLS {
            Point point;
            IndexType global_position;
        };
        const bool is_distributed = mrDataCommunicator.IsDistributed();
        IndexPartition<IndexType>(total_number_of_points).for_each(TLS(),[this, &is_distributed, &search_info, &rResults, &Radius, &allocation_size](std::size_t i_point, TLS& rTLS) {
            rTLS.point[0] = search_info.PointCoordinates[i_point * 3 + 0];
            rTLS.point[1] = search_info.PointCoordinates[i_point * 3 + 1];
            rTLS.point[2] = search_info.PointCoordinates[i_point * 3 + 2];
            if constexpr (ConsiderGlobalDataCommunicator) {
                rTLS.global_position = search_info.GlobalPosition[i_point];
            } else {
                rTLS.global_position = i_point;
            }
            auto& r_point_result = rResults[rTLS.global_position];

            // Search
            std::vector<ResultType> results;
            const int rank = is_distributed ? r_point_result.GetDataCommunicator().Rank() : -1;
            LocalSearchInRadius(rTLS.point, Radius, results, rank, allocation_size);
            for (auto& r_result : results) {
                r_point_result.AddResult(r_result);
            }
        });
        Timer::Stop("SearchWrapper::Search");

        // Synchronize
        Timer::Start("SearchWrapper::SynchronizeAll");
        rResults.SynchronizeAll(mrDataCommunicator);
        Timer::Stop("SearchWrapper::SynchronizeAll");
    }

    /**
     * @brief This method takes several points and finds the nearest object to it in a given radius.
     * @details If there are more than one object in the same minimum distance only one is returned.
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param Radius The radius to be checked.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @tparam TPointIteratorType The type of the point iterator.
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
        // Clear current solution
        if (ClearSolution) {
            rResults.Clear();
        }

        // Retrieving parameters
        const int allocation_size = mSettings["allocation_size"].GetInt();

        // The local bounding box
        const auto& r_local_bb = mpSearchObject ? mpSearchObject->GetBoundingBox() : BoundingBox<PointType>();

        // Prepare MPI search
        Timer::Start("SearchWrapper::SynchronousPointSynchronizationWithBoundingBox");
        DistributedSearchInformation search_info;
        SearchUtilities::SynchronousPointSynchronizationWithBoundingBox(itPointBegin, itPointEnd, search_info, r_local_bb, Radius, mrDataCommunicator, true);
        Timer::Stop("SearchWrapper::SynchronousPointSynchronizationWithBoundingBox");

        // Initialize results
        Timer::Start("SearchWrapper::PrepareResultsInProperRanks");
        PrepareResultsInProperRanks(rResults, search_info);
        Timer::Stop("SearchWrapper::PrepareResultsInProperRanks");

        // Perform the corresponding searches
        Timer::Start("SearchWrapper::Search");
        const std::size_t total_number_of_points = search_info.LocalIndices.size();
        struct TLS {
            Point point;
            IndexType global_position;
        };
        const bool is_distributed = mrDataCommunicator.IsDistributed();
        IndexPartition<IndexType>(total_number_of_points).for_each(TLS(),[this, &is_distributed, &search_info, &rResults, &Radius, &allocation_size](std::size_t i_point, TLS& rTLS) {
            rTLS.point[0] = search_info.PointCoordinates[i_point * 3 + 0];
            rTLS.point[1] = search_info.PointCoordinates[i_point * 3 + 1];
            rTLS.point[2] = search_info.PointCoordinates[i_point * 3 + 2];
            if constexpr (ConsiderGlobalDataCommunicator) {
                rTLS.global_position = search_info.GlobalPosition[i_point];
            } else {
                rTLS.global_position = i_point;
            }
            auto& r_point_result = rResults[rTLS.global_position];

            // Result of search
            ResultType local_result;
            const int rank = is_distributed ? r_point_result.GetDataCommunicator().Rank() : -1;
            LocalSearchNearestInRadius(rTLS.point, Radius, local_result, rank, allocation_size);
            if (local_result.GetIsObjectFound()) {
                r_point_result.AddResult(local_result);
            }
        });
        Timer::Stop("SearchWrapper::Search");

        // Synchronize
        Timer::Start("SearchWrapper::SynchronizeAll");
        rResults.SynchronizeAll(mrDataCommunicator);
        Timer::Stop("SearchWrapper::SynchronizeAll");

        // Remove the non closest results
        Timer::Start("SearchWrapper::KeepOnlyClosestResult");
        KeepOnlyClosestResult(rResults);
        Timer::Stop("SearchWrapper::KeepOnlyClosestResult");
    }

    /**
     * @brief This method takes several points and finds the nearest object to it.
     * @details If there are more than one object in the same minimum distance only one is returned.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
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
        // Clear current solution
        if (ClearSolution) {
            rResults.Clear();
        }

        // Get the maximum radius
        const auto global_bb = GetGlobalBoundingBox();
        const array_1d<double, 3> box_size = global_bb.GetMaxPoint() - global_bb.GetMinPoint();
        const double max_radius= *std::max_element(box_size.begin(), box_size.end());

        // The local bounding box
        const auto& r_local_bb = mpSearchObject ? mpSearchObject->GetBoundingBox() : BoundingBox<PointType>();

        // Prepare MPI search
        Timer::Start("SearchWrapper::SynchronousPointSynchronizationWithBoundingBox");
        DistributedSearchInformation search_info;
        SearchUtilities::SynchronousPointSynchronizationWithBoundingBox(itPointBegin, itPointEnd, search_info, r_local_bb, max_radius, mrDataCommunicator, true);
        Timer::Stop("SearchWrapper::SynchronousPointSynchronizationWithBoundingBox");

        // Initialize results
        Timer::Start("SearchWrapper::PrepareResultsInProperRanks");
        PrepareResultsInProperRanks(rResults, search_info);
        Timer::Stop("SearchWrapper::PrepareResultsInProperRanks");

        // Perform the corresponding searches
        Timer::Start("SearchWrapper::Search");
        const std::size_t total_number_of_points = search_info.LocalIndices.size();
        struct TLS {
            Point point;
            IndexType global_position;
        };
        const bool is_distributed = mrDataCommunicator.IsDistributed();
        IndexPartition<IndexType>(total_number_of_points).for_each(TLS(),[this, &is_distributed, &search_info, &rResults](std::size_t i_point, TLS& rTLS) {
            rTLS.point[0] = search_info.PointCoordinates[i_point * 3 + 0];
            rTLS.point[1] = search_info.PointCoordinates[i_point * 3 + 1];
            rTLS.point[2] = search_info.PointCoordinates[i_point * 3 + 2];
            if constexpr (ConsiderGlobalDataCommunicator) {
                rTLS.global_position = search_info.GlobalPosition[i_point];
            } else {
                rTLS.global_position = i_point;
            }
            auto& r_point_result = rResults[rTLS.global_position];

            // Result of search
            ResultType local_result;
            const int rank = is_distributed ? r_point_result.GetDataCommunicator().Rank() : -1;
            LocalSearchNearest(rTLS.point, local_result, rank);
            if (local_result.GetIsObjectFound()) {
                r_point_result.AddResult(local_result);
            }
        });
        Timer::Stop("SearchWrapper::Search");

        // Synchronize
        Timer::Start("SearchWrapper::SynchronizeAll");
        rResults.SynchronizeAll(mrDataCommunicator);
        Timer::Stop("SearchWrapper::SynchronizeAll");

        // Remove the non closest results
        Timer::Start("SearchWrapper::KeepOnlyClosestResult");
        KeepOnlyClosestResult(rResults);
        Timer::Stop("SearchWrapper::KeepOnlyClosestResult");
    }

    /**
     * @brief This method takes several points and search if it's inside an geometrical object of the domain (iterative version).
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @tparam TPointIteratorType The type of the point iterator.
     */
    template<class TPointIteratorType>
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

        // The local bounding box
        const auto& r_local_bb = mpSearchObject ? mpSearchObject->GetBoundingBox() : BoundingBox<PointType>();

        // Prepare MPI search
        Timer::Start("SearchWrapper::SynchronousPointSynchronizationWithBoundingBox");
        DistributedSearchInformation search_info;
        SearchUtilities::SynchronousPointSynchronizationWithBoundingBox(itPointBegin, itPointEnd, search_info, r_local_bb, 0.0, mrDataCommunicator, true);
        Timer::Stop("SearchWrapper::SynchronousPointSynchronizationWithBoundingBox");

        // Initialize results
        Timer::Start("SearchWrapper::PrepareResultsInProperRanks");
        PrepareResultsInProperRanks(rResults, search_info);
        Timer::Stop("SearchWrapper::PrepareResultsInProperRanks");

        // Perform the corresponding searches
        Timer::Start("SearchWrapper::Search");
        const std::size_t total_number_of_points = search_info.LocalIndices.size();
        struct TLS {
            Point point;
            IndexType global_position;
        };
        const bool is_distributed = mrDataCommunicator.IsDistributed();
        IndexPartition<IndexType>(total_number_of_points).for_each(TLS(),[this, &is_distributed, &search_info, &rResults](std::size_t i_point, TLS& rTLS) {
            rTLS.point[0] = search_info.PointCoordinates[i_point * 3 + 0];
            rTLS.point[1] = search_info.PointCoordinates[i_point * 3 + 1];
            rTLS.point[2] = search_info.PointCoordinates[i_point * 3 + 2];
            if constexpr (ConsiderGlobalDataCommunicator) {
                rTLS.global_position = search_info.GlobalPosition[i_point];
            } else {
                rTLS.global_position = i_point;
            }
            auto& r_point_result = rResults[rTLS.global_position];

            // Result of search
            ResultType local_result;
            const int rank = is_distributed ? r_point_result.GetDataCommunicator().Rank() : -1;
            LocalSearchIsInside(rTLS.point, local_result, rank);
            if (local_result.GetIsObjectFound()) {
                r_point_result.AddResult(local_result);
            }
        });
        Timer::Stop("SearchWrapper::Search");

        // Synchronize
        Timer::Start("SearchWrapper::SynchronizeAll");
        rResults.SynchronizeAll(mrDataCommunicator);
        Timer::Stop("SearchWrapper::SynchronizeAll");

        // Remove the non lowest rank results
        Timer::Start("SearchWrapper::KeepOnlyLowestRankResult");
        KeepOnlyLowestRankResult(rResults);
        Timer::Stop("SearchWrapper::KeepOnlyLowestRankResult");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Retrieves the static instance of SpatialSearchCommunication.
     * @details This function returns the static instance of SpatialSearchCommunication used for spatial searching.
     * @return The static instance of SpatialSearchCommunication.
     */
    static constexpr SpatialSearchCommunication GetSpatialSearchCommunication()
    {
        return TSpatialSearchCommunication;
    }

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

    typename TSearchObject::Pointer mpSearchObject = nullptr;   /// The pointer to the base search considered
    Kratos::unique_ptr<PointVector> mpPointVector = nullptr;    /// The point vector considered in the search trees
    std::vector<double> mGlobalBoundingBoxes;                   /// All the global BB, data is xmax, xmin,  ymax, ymin,  zmax, zmin
    const DataCommunicator& mrDataCommunicator;                 /// The data communicator
    Parameters mSettings;                                       /// The settings considered

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
     * @brief This method takes several points and finds all of the objects in the given radius to it.
     * @details The result contains the object and also its distance to the point.
     * @param rPoint The point to be checked.
     * @param Radius The radius to be checked.
     * @param rResults The results of the search.
     * @param Rank The rank of the search.
     * @param AllocationSize The allocation size of the results.
     */
    void LocalSearchInRadius(
        const PointType& rPoint,
        const double Radius,
        std::vector<ResultType>& rResults,
        const int Rank,
        const int AllocationSize = 1000
        );

    /**
     * @brief This method takes several points and finds the nearest object to it in a given radius.
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked.
     * @param Radius The radius to be checked.
     * @param rResult The result of the search.
     * @param Rank The rank of the search.
     * @param AllocationSize The allocation size of the results.
     */
    void LocalSearchNearestInRadius(
        const PointType& rPoint,
        const double Radius,
        ResultType& rResult,
        const int Rank,
        const int AllocationSize = 1000
        );

    /**
     * @brief This method takes several points and finds the nearest object to it.
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked.
     * @param rResult The result of the search.
     * @param Rank The rank of the search.
    */
    void LocalSearchNearest(
        const PointType& rPoint,
        ResultType& rResult,
        const int Rank
        );

    /**
     * @brief This method takes several points and search if it's inside an geometrical object of the domain.
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked.
     * @param rResult The result of the search.
     * @param Rank The rank of the search.
     */
    void LocalSearchIsInside(
        const PointType& rPoint,
        ResultType& rResult,
        const int Rank
        );

    /**
     * @brief This method removes tle solutions according to a given lambda
     * @param rResults Results to be simplified
     * @param rLambda Lambda to be considered as a criteria
     */
    template<class TLambda>
    void KeepOnlyGivenLambdaResult(
        ResultContainerVectorType& rResults,
        const TLambda& rLambda
        )
    {
        // MPI only
        if (mrDataCommunicator.IsDistributed()) {
            // Getting current rank
            const int rank = mrDataCommunicator.Rank();

            // Retrieve the solution
            auto& r_results_vector = rResults.GetContainer();
            for (auto& p_partial_result : r_results_vector) {
                auto& r_partial_result = *p_partial_result;
                // Then must have at least one solution, but just filter if at least 2
                const std::size_t number_of_global_results = r_partial_result.NumberOfGlobalResults();
                if (number_of_global_results > 1) {
                    // The values
                    const auto values = rLambda(r_partial_result);

                    // The indexes
                    std::vector<int> ranks = r_partial_result.GetResultRank();

                    // Find the index of the minimum value
                    auto it_min_distance = std::min_element(values.begin(), values.end());

                    // Check if the values vector is not empty
                    if (it_min_distance != values.end()) {
                        // Calculate the position
                        const IndexType pos = std::distance(values.begin(), it_min_distance);
                        if (rank == ranks[pos]) {
                            KRATOS_ERROR_IF(r_partial_result.NumberOfLocalResults() > 1) << "The rank criteria to filter results assumes that one rank only holds one local result. This is not true for " << r_partial_result.GetGlobalIndex() << " in rank " << rank << std::endl;
                        }

                        // Remove the index from the ranks vector
                        ranks.erase(ranks.begin() + pos);

                        // Remove all results but the closest one
                        r_partial_result.RemoveResultsFromRanksList(ranks);
                    } else {
                        KRATOS_ERROR << "Distances vector is empty." << std::endl;
                    }
                }
            }

            // Checking that is properly cleaned
            for (auto& p_partial_result : r_results_vector) {
                auto& r_partial_result = *p_partial_result;
                // Check that the number of results is 0 or 1
                KRATOS_ERROR_IF(r_partial_result.NumberOfGlobalResults() > 1) << "Cleaning has not been done properly. Number of results: " << r_partial_result.NumberOfGlobalResults() << std::endl;
                // Check that is not empty locally
                if (r_partial_result.NumberOfGlobalResults() == 1) {
                    r_partial_result.Barrier();
                    const unsigned int number_of_local_results = r_partial_result.NumberOfLocalResults();
                    const auto& r_sub_data_communicator = r_partial_result.GetDataCommunicator();
                    KRATOS_ERROR_IF(r_sub_data_communicator.SumAll(number_of_local_results) == 0) << "Local results also removed in result " << r_partial_result.GetGlobalIndex() << std::endl;
                }
            }
        }
    }

    /**
     * @brief This method removes the solutions if there are not the closest
     * @param rResults Results to be simplified
     */
    void KeepOnlyClosestResult(ResultContainerVectorType& rResults);

    /**
     * @brief This method removes the solutions not in the lowest rank
     * @param rResults Results to be simplified
     */
    void KeepOnlyLowestRankResult(ResultContainerVectorType& rResults);

    /**
     * @brief Generates a string by appending integer ranks to a base name.
     * @details This function takes a base name and a vector of integer ranks and concatenates  them to generate a single string. The integer ranks are converted to strings and appended to the base name in the order they appear in the vector.
     * @param rBaseName The base name to which the integer ranks will be appended.
     * @param rRanks    A vector of integer ranks to be appended to the base name.
     * @return A string containing the base name followed by the integer ranks.
     */
    std::string GenerateNameFromRanks(
        const std::string& rBaseName,
        const std::vector<int>& rRanks
        );

    /**
     * @brief Prepare the search results in proper ranks.
     * @details This function prepares the search results in proper ranks based on the provided search information.
     * @param rResults The vector containing the search results to be prepared.
     * @param rSearchInfo The distributed search information to determine the proper ranks.
     */
    void PrepareResultsInProperRanks(
        ResultContainerVectorType& rResults,
        const DistributedSearchInformation& rSearchInfo
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