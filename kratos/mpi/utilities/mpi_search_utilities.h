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
#include <vector>
#include <numeric>

// External includes

// Project includes
#include "includes/data_communicator.h"
#include "includes/variables.h"
#include "includes/model_part.h"

namespace Kratos
{

/**
 * @class MPISearchUtilities
 * @ingroup KratosCore
 * @brief This class provides the utilities for MPI search
 * @author Vicente Mataix Ferrandiz
 */
class MPISearchUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// The index definition
    using IndexType = std::size_t;

    /// Define zero tolerance as Epsilon
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief MPISynchronousPointSynchronization prepares synchronously the coordinates of the points for MPI search.
     * @param itPointBegin Iterator to the beginning of the points range
     * @param itPointEnd Iterator to the end of the points range
     * @param rAllPointsCoordinates vector where the computed coordinates will be stored
     * @param rAllPointsIds The ids of all the points (just a counter for points, and ids for nodes)
     * @param rDataCommunicator The data communicator
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static void MPISynchronousPointSynchronization(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        std::vector<IndexType>& rAllPointsIds,
        const DataCommunicator& rDataCommunicator
        )
    {
        // First check that the points are the same in all processes
        int number_of_points, total_number_of_points;
        CalculateNumberOfPoints(itPointBegin, itPointEnd, number_of_points, total_number_of_points, rDataCommunicator);

        KRATOS_DEBUG_ERROR_IF(number_of_points < 0) << "The number of points is negative" << std::endl;
        KRATOS_DEBUG_ERROR_IF(total_number_of_points < 0) << "The total number of points is negative" << std::endl;

        // We synchronize the points
        SynchronizePoints(itPointBegin, itPointEnd, rAllPointsCoordinates, rAllPointsIds, rDataCommunicator, number_of_points, total_number_of_points);
    }

    /**
     * @brief MPISynchronousPointSynchronizationWithRecvSizes prepares synchronously the coordinates of the points for MPI search including the recv sizes
     * @details With recv sizes
     * @param itPointBegin Iterator to the beginning of the points range
     * @param itPointEnd Iterator to the end of the points range
     * @param rAllPointsCoordinates vector where the computed coordinates will be stored
     * @param rAllPointsIds The ids of all the points (just a counter for points, and ids for nodes)
     * @param rDataCommunicator The data communicator
     * @return The resulting whole radius vector
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static std::vector<int> MPISynchronousPointSynchronizationWithRecvSizes(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        std::vector<IndexType>& rAllPointsIds,
        const DataCommunicator& rDataCommunicator
        )
    {
        // First check that the points are the same in all processes
        int number_of_points, total_number_of_points;
        CalculateNumberOfPoints(itPointBegin, itPointEnd, number_of_points, total_number_of_points, rDataCommunicator);

        KRATOS_DEBUG_ERROR_IF(number_of_points < 0) << "The number of points is negative" << std::endl;
        KRATOS_DEBUG_ERROR_IF(total_number_of_points < 0) << "The total number of points is negative" << std::endl;

        // We synchronize the points
        return SynchronizePoints(itPointBegin, itPointEnd, rAllPointsCoordinates, rAllPointsIds, rDataCommunicator, number_of_points, total_number_of_points);
    }

    /**
     * @brief MPISynchronousPointSynchronizationWithRadius prepares synchronously the coordinates of the points for MPI search including radius
     * @details With radius
     * @param itPointBegin Iterator to the beginning of the points range
     * @param itPointEnd Iterator to the end of the points range
     * @param rAllPointsCoordinates vector where the computed coordinates will be stored
     * @param rAllPointsIds The ids of all the points (just a counter for points, and ids for nodes)
     * @param rRadius The radius of the points
     * @param rDataCommunicator The data communicator
     * @return The resulting whole radius vector
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static std::vector<double> MPISynchronousPointSynchronizationWithRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        std::vector<IndexType>& rAllPointsIds,
        const std::vector<double>& rRadius,
        const DataCommunicator& rDataCommunicator
        )
    {
        // First check that the points are the same in all processes
        int number_of_points, total_number_of_points;
        CalculateNumberOfPoints(itPointBegin, itPointEnd, number_of_points, total_number_of_points, rDataCommunicator);

        KRATOS_DEBUG_ERROR_IF(number_of_points < 0) << "The number of points is negative" << std::endl;
        KRATOS_DEBUG_ERROR_IF(total_number_of_points < 0) << "The total number of points is negative" << std::endl;

        // We synchronize the points
        const auto recv_sizes = SynchronizePoints(itPointBegin, itPointEnd, rAllPointsCoordinates, rAllPointsIds, rDataCommunicator, number_of_points, total_number_of_points);

        // Get radius
        const auto radius = SynchronizeRadius(recv_sizes, rRadius, rDataCommunicator);
        return radius;
    }

    ///@}
private:
    ///@name Private Operations
    ///@{  

    /**
     * @details Synchronizes points between different processes. 
     * @details Synchonously
     * @param itPointBegin Iterator pointing to the beginning of the range of points
     * @param itPointEnd Iterator pointing to the end of the range of points
     * @param rAllPointsCoordinates Vector to store the synchronized points' coordinates
     * @param rAllPointsIds The ids of all the points (just a counter for points, and ids for nodes)
     * @param rDataCommunicator Object for data communication between processes
     * @param NumberOfPoints Local number of points to be synchronized
     * @param TotalNumberOfPoints Total number of points across all processes
     * @return A vector containing the sizes of data for each process
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static std::vector<int> SynchronizePoints(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        std::vector<IndexType>& rAllPointsIds,
        const DataCommunicator& rDataCommunicator,
        const int NumberOfPoints,
        const int TotalNumberOfPoints
        )
    {
        // Get the World Size and rank in MPI
        const int world_size = rDataCommunicator.Size();
        const int rank = rDataCommunicator.Rank();

        // Getting global number of points
        std::vector<int> points_per_partition(world_size);
        std::vector<int> send_points_per_partition(1, NumberOfPoints);
        rDataCommunicator.AllGather(send_points_per_partition, points_per_partition);
        std::size_t initial_id = 0;
        if constexpr (!std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value) {
            initial_id = std::accumulate(points_per_partition.begin(), points_per_partition.begin() + rank, 0);
        }

        // Initialize local points coordinates
        std::size_t counter = 0;
        array_1d<double, 3> coordinates;
        unsigned int i_coord;
        std::vector<double> send_points_coordinates(NumberOfPoints * 3);
        std::vector<IndexType> send_points_ids(NumberOfPoints);
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; ++it_point) {
            // In case of considering nodes
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value) {
                if (it_point->FastGetSolutionStepValue(PARTITION_INDEX) != rank) {
                    continue; // Skip if not local
                }
            }
            noalias(coordinates) = it_point->Coordinates();
            if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value) {
                send_points_ids[counter] = it_point->Id();
            } else {
                send_points_ids[counter] = initial_id + counter;
            }
            for (i_coord = 0; i_coord < 3; ++i_coord) {
                send_points_coordinates[3 * counter + i_coord] = coordinates[i_coord];
            }
            ++counter;
        }

        // Resizing
        rAllPointsCoordinates.resize(TotalNumberOfPoints * 3);
        rAllPointsIds.resize(TotalNumberOfPoints);

        /* Sync coordinates */

        // Generate vectors with sizes for AllGatherv
        std::vector<int> recv_sizes(world_size, 0);
        for (int i_rank = 0; i_rank < world_size; ++i_rank) {
            recv_sizes[i_rank] = 3 * points_per_partition[i_rank];
        }
        std::vector<int> recv_offsets(world_size, 0);
        for (int i_rank = 1; i_rank < world_size; ++i_rank) {
            recv_offsets[i_rank] = recv_offsets[i_rank - 1] + recv_sizes[i_rank - 1];
        }

        // Invoke AllGatherv
        rDataCommunicator.AllGatherv(send_points_coordinates, rAllPointsCoordinates, recv_sizes, recv_offsets);

        /* Sync Ids */

        // Generate vectors with sizes for AllGatherv
        for (int i_rank = 0; i_rank < world_size; ++i_rank) {
            recv_sizes[i_rank] = points_per_partition[i_rank];
        }
        for (int i_rank = 1; i_rank < world_size; ++i_rank) {
            recv_offsets[i_rank] = recv_offsets[i_rank - 1] + recv_sizes[i_rank - 1];
        }

        // Invoke AllGatherv
        rDataCommunicator.AllGatherv(send_points_ids, rAllPointsIds, recv_sizes, recv_offsets);

        // Return the recv sizes
        return recv_sizes;
    }

    /**
     * @brief Synchronizes the radius of all points in a distributed system.
     * @param rRecvSizes a vector containing the number of points to be received from each rank
     * @param rRadius a vector containing the radius of each point
     * @param rDataCommunicator the communication object used for data exchange
     * @return A vector containing the synchronized radius of all points
     */
    static std::vector<double> SynchronizeRadius(
        const std::vector<int>& rRecvSizes,
        const std::vector<double>& rRadius,
        const DataCommunicator& rDataCommunicator
        )
    {
        // First we calculate the total number of points to communicate
        const int total_number_of_points = std::accumulate(rRecvSizes.begin(), rRecvSizes.end(), 0);

        // Synchronize radius
        if (total_number_of_points == 0) { // If all points are the same
            return rRadius;
        } else {                           // If not
            // The resulting radius
            std::vector<double> all_points_radius(total_number_of_points);

            // MPI information
            const int world_size = rDataCommunicator.Size();

            // Generate vectors with sizes for AllGatherv
            std::vector<int> recv_offsets(world_size, 0);
            for (int i_rank = 1; i_rank < world_size; ++i_rank) {
                recv_offsets[i_rank] = recv_offsets[i_rank - 1] + rRecvSizes[i_rank - 1];
            }

            // Invoke AllGatherv
            rDataCommunicator.AllGatherv(rRadius, all_points_radius, rRecvSizes, recv_offsets);

            return all_points_radius;
        }
    }

    /**
     * @brief This method checks if all nodes are the same across all partitions
     * @param itPointBegin Iterator to the beginning of the points range
     * @param itPointEnd Iterator to the end of the points range
     * @param rNumberOfPoints Number of points in the range
     * @param rTotalNumberOfPoints Total number of points in all partitions
     * @return true if all points are the same in all partitions
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static void CalculateNumberOfPoints(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        int& rNumberOfPoints,
        int& rTotalNumberOfPoints,
        const DataCommunicator& rDataCommunicator
        )
    {
        // Get the World Size and rank in MPI
        const int rank = rDataCommunicator.Rank();

        // Getting local number of points
        rNumberOfPoints = std::distance(itPointBegin, itPointEnd);
        
        // In case of considering nodes
        if constexpr (std::is_same<TPointIteratorType, ModelPart::NodeIterator>::value) {
            for (auto it_node = itPointBegin; it_node < itPointEnd; ++it_node) {
                if (it_node->FastGetSolutionStepValue(PARTITION_INDEX) != rank) {
                    --rNumberOfPoints; // Remove is not local
                }
            }
        }

        // Sum all the nodes
        rTotalNumberOfPoints = rDataCommunicator.SumAll(rNumberOfPoints);
    }
    ///@}
};

}  // namespace Kratos