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
                mpPointVector = Kratos::make_unique<PointVector>(SearchUtilities::PreparePointsSearch(rGeometricalObjectsVector));

                // Create the search object
                if constexpr (!IsDynamicBins) {
                    const int bucket_size = mSettings["bucket_size"].GetInt();
                    mpSearchObject = Kratos::make_shared<TSearchObject>(mpPointVector->begin(), mpPointVector->end(), bucket_size);
                } else {
                    mpSearchObject = Kratos::make_shared<TSearchObject>(mpPointVector->begin(), mpPointVector->end());
                }
            }
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
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param Radius The radius to be checked.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @param ConsiderGlobalDataCommunicator If a sub DataCommunicator or a global DataCommunicator is considered.
     * @tparam TPointIteratorType The type of the point iterator.
     */
    template<class TPointIteratorType>
    void SearchInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
        )
    {
        // Call distributed version
        if (mrDataCommunicator.IsDistributed()) {
            DistributedSearchInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution, ConsiderGlobalDataCommunicator);
        } else { // Call serial version
            SerialSearchInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution);
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius.
     * @details If there are more than one object in the same minimum distance only one is returned.
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param Radius The radius to be checked.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @param ConsiderGlobalDataCommunicator If a sub DataCommunicator or a global DataCommunicator is considered.
     * @tparam TPointIteratorType The type of the point iterator.
     */
    template<class TPointIteratorType>
    void SearchNearestInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
        )
    {
        // Call distributed version
        if (mrDataCommunicator.IsDistributed()) {
            DistributedSearchNearestInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution, ConsiderGlobalDataCommunicator);
        } else { // Call serial version
            SerialSearchNearestInRadius(itPointBegin, itPointEnd, Radius, rResults, ClearSolution);
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it.
     * @details If there are more than one object in the same minimum distance only one is returned.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @param ConsiderGlobalDataCommunicator If a sub DataCommunicator or a global DataCommunicator is considered.
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<class TPointIteratorType>
    void SearchNearest(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
        )
    {
        // Call distributed version
        if (mrDataCommunicator.IsDistributed()) {
            DistributedSearchNearest(itPointBegin, itPointEnd, rResults, ClearSolution, ConsiderGlobalDataCommunicator);
        } else { // Call serial version
            SerialSearchNearest(itPointBegin, itPointEnd, rResults, ClearSolution);
        }
    }

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain (iterative version).
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @param ConsiderGlobalDataCommunicator If a sub DataCommunicator or a global DataCommunicator is considered.
     * @tparam TPointIteratorType The type of the point iterator.
     */
    template<class TPointIteratorType>
    void SearchIsInside(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
        )
    {
        // Call distributed version
        if (mrDataCommunicator.IsDistributed()) {
            DistributedSearchIsInside(itPointBegin, itPointEnd, rResults, ClearSolution, ConsiderGlobalDataCommunicator);
        } else { // Call serial version
            SerialSearchIsInside(itPointBegin, itPointEnd, rResults, ClearSolution);
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

    typename TSearchObject::Pointer mpSearchObject =  nullptr;  /// The pointer to the base search considered
    Kratos::unique_ptr<PointVector> mpPointVector =  nullptr;   /// The point vector considered in the search trees
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

        // Retrieving the number of points
        const std::size_t number_of_points = itPointEnd - itPointBegin;

        // Initialize results
        Timer::Start("SearchWrapper::InitializeResults");
        {
            const std::vector<const DataCommunicator*> data_communicators(number_of_points, &mrDataCommunicator);
            rResults.InitializeResults(data_communicators);
        }
        Timer::Stop("SearchWrapper::InitializeResults");

        // Adding the results to the container
        Timer::Start("SearchWrapper::Search");
        IndexPartition<IndexType>(number_of_points).for_each([this, &itPointBegin, &rResults, &Radius, &allocation_size](std::size_t index) {
            auto it_point = itPointBegin + index;
            auto& r_point_result = rResults[index];

            // Set some values
            r_point_result.SetLocalIndex(index);
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value || std::is_same<TPointIteratorType, ModelPart::NodeConstantIterator>::value) {
                r_point_result.SetGlobalIndex(it_point->Id());
            } else {
                r_point_result.SetGlobalIndex(index);
            }

            // Search
            std::vector<ResultType> results;
            LocalSearchInRadius(*it_point, Radius, results, 0, allocation_size);
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

        // Retrieving the number of points
        const std::size_t number_of_points = itPointEnd - itPointBegin;

        // Initialize results
        Timer::Start("SearchWrapper::InitializeResults");
        {
            const std::vector<const DataCommunicator*> data_communicators(number_of_points, &mrDataCommunicator);
            rResults.InitializeResults(data_communicators);
        }
        Timer::Stop("SearchWrapper::InitializeResults");

        // Adding the results to the container
        Timer::Start("SearchWrapper::Search");
        IndexPartition<IndexType>(number_of_points).for_each([this, &itPointBegin, &rResults, &Radius, &allocation_size](std::size_t index) {
            auto it_point = itPointBegin + index;
            auto& r_point_result = rResults[index];

            // Set some values
            r_point_result.SetLocalIndex(index);
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value || std::is_same<TPointIteratorType, ModelPart::NodeConstantIterator>::value) {
                r_point_result.SetGlobalIndex(it_point->Id());
            } else {
                r_point_result.SetGlobalIndex(index);
            }

            // Search
            ResultType result;
            LocalSearchNearestInRadius(*it_point, Radius, result, 0, allocation_size);
            if (result.GetIsObjectFound()) {
                r_point_result.AddResult(result);
            }
        });
        Timer::Stop("SearchWrapper::Search");

        // Synchronize
        Timer::Start("SearchWrapper::SynchronizeAll");
        rResults.SynchronizeAll(mrDataCommunicator);
        Timer::Stop("SearchWrapper::SynchronizeAll");
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

        // Retrieving the number of points
        const std::size_t number_of_points = itPointEnd - itPointBegin;

        // Initialize results
        Timer::Start("SearchWrapper::InitializeResults");
        {
            const std::vector<const DataCommunicator*> data_communicators(number_of_points, &mrDataCommunicator);
            rResults.InitializeResults(data_communicators);
        }
        Timer::Stop("SearchWrapper::InitializeResults");

        // Adding the results to the container
        Timer::Start("SearchWrapper::Search");
        IndexPartition<IndexType>(number_of_points).for_each([this, &itPointBegin, &rResults](std::size_t index) {
            auto it_point = itPointBegin + index;
            auto& r_point_result = rResults[index];

            // Set some values
            r_point_result.SetLocalIndex(index);
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value || std::is_same<TPointIteratorType, ModelPart::NodeConstantIterator>::value) {
                r_point_result.SetGlobalIndex(it_point->Id());
            } else {
                r_point_result.SetGlobalIndex(index);
            }

            // Search
            ResultType result;
            LocalSearchNearest(*it_point, result, 0);
            if (result.GetIsObjectFound()) {
                r_point_result.AddResult(result);
            }
        });
        Timer::Stop("SearchWrapper::Search");

        // Synchronize
        Timer::Start("SearchWrapper::SynchronizeAll");
        rResults.SynchronizeAll(mrDataCommunicator);
        Timer::Stop("SearchWrapper::SynchronizeAll");
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

        // Retrieving the number of points
        const std::size_t number_of_points = itPointEnd - itPointBegin;

        // Initialize results
        Timer::Start("SearchWrapper::InitializeResults");
        {
            const std::vector<const DataCommunicator*> data_communicators(number_of_points, &mrDataCommunicator);
            rResults.InitializeResults(data_communicators);
        }
        Timer::Stop("SearchWrapper::InitializeResults");

        // Adding the results to the container
        Timer::Start("SearchWrapper::Search");
        IndexPartition<IndexType>(number_of_points).for_each([this, &itPointBegin, &rResults](std::size_t index) {
            auto it_point = itPointBegin + index;
            auto& r_point_result = rResults[index];

            // Set some values
            r_point_result.SetLocalIndex(index);
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value || std::is_same<TPointIteratorType, ModelPart::NodeConstantIterator>::value) {
                r_point_result.SetGlobalIndex(it_point->Id());
            } else {
                r_point_result.SetGlobalIndex(index);
            }

            // Search
            ResultType result;
            LocalSearchIsInside(*it_point, result, 0);
            if (result.GetIsObjectFound()) {
                r_point_result.AddResult(result);
            }
        });
        Timer::Stop("SearchWrapper::Search");

        // Synchronize
        Timer::Start("SearchWrapper::SynchronizeAll");
        rResults.SynchronizeAll(mrDataCommunicator);
        Timer::Stop("SearchWrapper::SynchronizeAll");
    }

#ifdef KRATOS_USING_MPI

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it (MPI version).
     * @details The result contains the object and also its distance to the point.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param Radius The radius to be checked.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @param ConsiderGlobalDataCommunicator If a sub DataCommunicator or a global DataCommunicator is considered.
     * @tparam TPointIteratorType The type of the point iterator.
     */
    template<class TPointIteratorType>
    void DistributedSearchInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
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
        PrepareResultsInProperRanks(rResults, search_info, ConsiderGlobalDataCommunicator);
        Timer::Stop("SearchWrapper::PrepareResultsInProperRanks");

        // Perform the corresponding searches
        Timer::Start("SearchWrapper::Search");
        const std::size_t total_number_of_points = search_info.Ranks.size();
        struct TLS {
            Point point;
        };
        IndexPartition<IndexType>(total_number_of_points).for_each(TLS(),[this, &search_info, &rResults, &Radius, &allocation_size](std::size_t i_point, TLS& rTLS) {
            rTLS.point[0] = search_info.PointCoordinates[i_point * 3 + 0];
            rTLS.point[1] = search_info.PointCoordinates[i_point * 3 + 1];
            rTLS.point[2] = search_info.PointCoordinates[i_point * 3 + 2];
            auto& r_point_result = rResults[i_point];

            // Search
            std::vector<ResultType> results;
            LocalSearchInRadius(rTLS.point, Radius, results, r_point_result.GetDataCommunicator().Rank(), allocation_size);
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
     * @brief This method takes a point and finds the nearest object to it in a given radius (MPI version).
     * @details If there are more than one object in the same minimum distance only one is returned.
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param Radius The radius to be checked.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @param ConsiderGlobalDataCommunicator If a sub DataCommunicator or a global DataCommunicator is considered.
     * @tparam TPointIteratorType The type of the point iterator.
     */
    template<class TPointIteratorType>
    void DistributedSearchNearestInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
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
        PrepareResultsInProperRanks(rResults, search_info, ConsiderGlobalDataCommunicator);
        Timer::Stop("SearchWrapper::PrepareResultsInProperRanks");

        // Perform the corresponding searches
        Timer::Start("SearchWrapper::Search");
        const std::size_t total_number_of_points = search_info.Ranks.size();
        struct TLS {
            Point point;
        };
        IndexPartition<IndexType>(total_number_of_points).for_each(TLS(),[this, &search_info, &rResults, &Radius, &allocation_size](std::size_t i_point, TLS& rTLS) {
            rTLS.point[0] = search_info.PointCoordinates[i_point * 3 + 0];
            rTLS.point[1] = search_info.PointCoordinates[i_point * 3 + 1];
            rTLS.point[2] = search_info.PointCoordinates[i_point * 3 + 2];
            auto& r_point_result = rResults[i_point];

            // Result of search
            ResultType local_result;
            LocalSearchNearestInRadius(rTLS.point, Radius, local_result, r_point_result.GetDataCommunicator().Rank(), allocation_size);
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
     * @brief This method takes a point and finds the nearest object to it (MPI version).
     * @details If there are more than one object in the same minimum distance only one is returned.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @param ConsiderGlobalDataCommunicator If a sub DataCommunicator or a global DataCommunicator is considered.
     * @tparam TPointIteratorType The type of the point iterator.
     */
    template<class TPointIteratorType>
    void DistributedSearchNearest(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
        )
    {
        // Clear current solution
        if (ClearSolution) {
            rResults.Clear();
        }

        // Get the maximum radius
        const auto bb = GetBoundingBox();
        const array_1d<double, 3> box_size = bb.GetMaxPoint() - bb.GetMinPoint();
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
        PrepareResultsInProperRanks(rResults, search_info, ConsiderGlobalDataCommunicator);
        Timer::Stop("SearchWrapper::PrepareResultsInProperRanks");

        // Perform the corresponding searches
        Timer::Start("SearchWrapper::Search");
        const std::size_t total_number_of_points = search_info.Ranks.size();
        struct TLS {
            Point point;
        };
        IndexPartition<IndexType>(total_number_of_points).for_each(TLS(),[this, &search_info, &rResults](std::size_t i_point, TLS& rTLS) {
            rTLS.point[0] = search_info.PointCoordinates[i_point * 3 + 0];
            rTLS.point[1] = search_info.PointCoordinates[i_point * 3 + 1];
            rTLS.point[2] = search_info.PointCoordinates[i_point * 3 + 2];
            auto& r_point_result = rResults[i_point];

            // Result of search
            ResultType local_result;
            LocalSearchNearest(rTLS.point, local_result, r_point_result.GetDataCommunicator().Rank());
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
     * @brief This method takes a point and search if it's inside an geometrical object of the domain (iterative version) (MPI version).
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param itPointBegin The first point iterator.
     * @param itPointEnd The last point iterator.
     * @param rResults The results of the search.
     * @param ClearSolution Clear the current solution.
     * @param ConsiderGlobalDataCommunicator If a sub DataCommunicator or a global DataCommunicator is considered.
     * @tparam TPointIteratorType The type of the point iterator.
     */
    template<class TPointIteratorType>
    void DistributedSearchIsInside(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
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
        PrepareResultsInProperRanks(rResults, search_info, ConsiderGlobalDataCommunicator);
        Timer::Stop("SearchWrapper::PrepareResultsInProperRanks");

        // Perform the corresponding searches
        Timer::Start("SearchWrapper::Search");
        const std::size_t total_number_of_points = search_info.Ranks.size();
        struct TLS {
            Point point;
        };
        IndexPartition<IndexType>(total_number_of_points).for_each(TLS(),[this, &search_info, &rResults](std::size_t i_point, TLS& rTLS) {
            rTLS.point[0] = search_info.PointCoordinates[i_point * 3 + 0];
            rTLS.point[1] = search_info.PointCoordinates[i_point * 3 + 1];
            rTLS.point[2] = search_info.PointCoordinates[i_point * 3 + 2];
            auto& r_point_result = rResults[i_point];

            // Result of search
            ResultType local_result;
            LocalSearchIsInside(rTLS.point, local_result, r_point_result.GetDataCommunicator().Rank());
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

#else
    template<class TPointIteratorType>
    void DistributedSearchInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
        )
    {
        KRATOS_ERROR << "Running distributed method requires to compile with MPI support" << std::endl;
    }

    template<class TPointIteratorType>
    void DistributedSearchNearestInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
        )
    {
        KRATOS_ERROR << "Running distributed method requires to compile with MPI support" << std::endl;
    }

    template<class TPointIteratorType>
    void DistributedSearchNearest(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
        )
    {
        KRATOS_ERROR << "Running distributed method requires to compile with MPI support" << std::endl;
    }

    template<class TPointIteratorType>
    void DistributedSearchIsInside(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        ResultContainerVectorType& rResults,
        const bool ClearSolution = true,
        const bool ConsiderGlobalDataCommunicator = false
        )
    {
        KRATOS_ERROR << "Running distributed method requires to compile with MPI support" << std::endl;
    }

#endif

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it.
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
     * @brief This method takes a point and finds the nearest object to it in a given radius.
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
     * @brief This method takes a point and finds the nearest object to it.
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
     * @brief This method takes a point and search if it's inside an geometrical object of the domain.
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
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
        // Retrieve th solution
        auto& r_results_vector = rResults.GetContainer();
        for (auto& p_partial_result : r_results_vector) {
            auto& r_partial_result = *p_partial_result;
            // Then must have at least one solution, but just filter if at least 2
            const std::size_t number_of_global_results = r_partial_result.NumberOfGlobalResults();
            if (number_of_global_results > 1) {
                // The values
                const auto values = rLambda(r_partial_result);

                // The indexes
                std::vector<IndexType> indexes = r_partial_result.GetResultIndices();

                // Find the index of the minimum value
                auto it_min_distance = std::min_element(values.begin(), values.end());

                // Check if the values vector is not empty
                if (it_min_distance != values.end()) {
                    // Calculate the position
                    const IndexType pos = std::distance(values.begin(), it_min_distance);

                    // Retrieve the corresponding index from indexes vector
                    const IndexType index_to_remove = indexes[pos];

                    // Remove the index from the indexes vector
                    indexes.erase(std::remove(indexes.begin(), indexes.end(), index_to_remove), indexes.end());

                    // Remove all results but the closest one
                    r_partial_result.RemoveResultsFromIndexesList(indexes);
                } else {
                    KRATOS_ERROR << "Distances vector is empty." << std::endl;
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
     * @param ConsiderGlobalDataCommunicator If a sub DataCommunicator or a global DataCommunicator is considered.
     */
    void PrepareResultsInProperRanks(
        ResultContainerVectorType& rResults,
        const DistributedSearchInformation& rSearchInfo,
        const bool ConsiderGlobalDataCommunicator = false
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