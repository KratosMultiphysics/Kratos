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
     * @param rDataCommunicator The data communicator
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static void MPISynchronousPointSynchronization(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        const DataCommunicator& rDataCommunicator
        )
    {
        // First check that the points are the same in all processes
        int number_of_points, total_number_of_points;
        const bool all_points_are_the_same = CheckAllPointsAreTheSameAndCalculateNumberOfPoints(itPointBegin, itPointEnd, number_of_points, total_number_of_points, rDataCommunicator);

        KRATOS_DEBUG_ERROR_IF(number_of_points < 0) << "The number of points is negative" << std::endl;
        KRATOS_DEBUG_ERROR_IF(total_number_of_points < 0) << "The total number of points is negative" << std::endl;

        // We synchronize the points
        SynchronizePoints(itPointBegin, itPointEnd, rAllPointsCoordinates, rDataCommunicator, all_points_are_the_same, number_of_points, total_number_of_points);
    }

    /**
     * @brief MPISynchronousPointSynchronizationWithRecvSizes prepares synchronously the coordinates of the points for MPI search including the recv sizes
     * @details With recv sizes
     * @param itPointBegin Iterator to the beginning of the points range
     * @param itPointEnd Iterator to the end of the points range
     * @param rAllPointsCoordinates vector where the computed coordinates will be stored
     * @param rDataCommunicator The data communicator
     * @return The resulting whole radius vector
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static std::vector<int> MPISynchronousPointSynchronizationWithRecvSizes(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        const DataCommunicator& rDataCommunicator
        )
    {
        // First check that the points are the same in all processes
        int number_of_points, total_number_of_points;
        const bool all_points_are_the_same = CheckAllPointsAreTheSameAndCalculateNumberOfPoints(itPointBegin, itPointEnd, number_of_points, total_number_of_points, rDataCommunicator);

        KRATOS_DEBUG_ERROR_IF(number_of_points < 0) << "The number of points is negative" << std::endl;
        KRATOS_DEBUG_ERROR_IF(total_number_of_points < 0) << "The total number of points is negative" << std::endl;

        // We synchronize the points
        SynchronizePoints(itPointBegin, itPointEnd, rAllPointsCoordinates, rDataCommunicator, all_points_are_the_same, number_of_points, total_number_of_points);

        // Get recv_sizes
        const auto recv_sizes = SynchronizeRecvSizes(number_of_points, all_points_are_the_same, rDataCommunicator);
        return recv_sizes;
    }

    /**
     * @brief MPISynchronousPointSynchronizationWithRadius prepares synchronously the coordinates of the points for MPI search including radius
     * @details With radius
     * @param itPointBegin Iterator to the beginning of the points range
     * @param itPointEnd Iterator to the end of the points range
     * @param rAllPointsCoordinates vector where the computed coordinates will be stored
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
        const std::vector<double>& rRadius,
        const DataCommunicator& rDataCommunicator
        )
    {
        // First check that the points are the same in all processes
        int number_of_points, total_number_of_points;
        const bool all_points_are_the_same = CheckAllPointsAreTheSameAndCalculateNumberOfPoints(itPointBegin, itPointEnd, number_of_points, total_number_of_points, rDataCommunicator);

        KRATOS_DEBUG_ERROR_IF(number_of_points < 0) << "The number of points is negative" << std::endl;
        KRATOS_DEBUG_ERROR_IF(total_number_of_points < 0) << "The total number of points is negative" << std::endl;

        // We synchronize the points
        SynchronizePoints(itPointBegin, itPointEnd, rAllPointsCoordinates, rDataCommunicator, all_points_are_the_same, number_of_points, total_number_of_points);

        // Get recv_sizes
        const auto recv_sizes = SynchronizeRecvSizes(number_of_points, all_points_are_the_same, rDataCommunicator);

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
     * @param rDataCommunicator Object for data communication between processes
     * @param AllPointsAreTheSame Flag indicating if all points are the same
     * @param NumberOfPoints Local number of points to be synchronized
     * @param TotalNumberOfPoints Total number of points across all processes
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static void SynchronizePoints(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        const DataCommunicator& rDataCommunicator,
        const bool AllPointsAreTheSame,
        const int NumberOfPoints,
        const int TotalNumberOfPoints
        )
    {
        // Initialize local points coordinates
        std::size_t counter = 0;
        array_1d<double, 3> coordinates;
        unsigned int i_coord;
        std::vector<double> send_points_coordinates(NumberOfPoints * 3);
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; ++it_point) {
            noalias(coordinates) = it_point->Coordinates();
            for (i_coord = 0; i_coord < 3; ++i_coord) {
                send_points_coordinates[3 * counter + i_coord] = coordinates[i_coord];
            }
            ++counter;
        }

        // If all points are the same
        if (AllPointsAreTheSame) {
            rAllPointsCoordinates = send_points_coordinates;
        } else { // If not
            // MPI information
            const int world_size = rDataCommunicator.Size();

            // Getting global number of points
            std::vector<int> points_per_partition(world_size);
            std::vector<int> send_points_per_partition(1, NumberOfPoints);
            rDataCommunicator.AllGather(send_points_per_partition, points_per_partition);

            // Getting global coordinates
            rAllPointsCoordinates.resize(TotalNumberOfPoints * 3);

            // Generate vectors with sizes for AllGatherv
            std::vector<int> recv_sizes(world_size, 0);
            for (int i_rank = 0; i_rank < world_size; ++i_rank) {
                recv_sizes[i_rank] = 3 * points_per_partition[i_rank];
            }
            std::vector<int> recv_offsets(world_size, 0);
            for (int i_rank = 1; i_rank < world_size; ++i_rank) {
                recv_offsets[i_rank] = recv_offsets[i_rank - 1] + recv_sizes[i_rank - 1];
            }

            // Invoque AllGatherv
            rDataCommunicator.AllGatherv(send_points_coordinates, rAllPointsCoordinates, recv_sizes, recv_offsets);
        }
    }

    /**
     * @brief Synchronizes the sizes of data among multiple processes using MPI.
     * @param NumberOfPoints The number of local points in the data
     * @param AllPointsAreTheSame Flag indicating if all points are the same
     * @param rDataCommunicator The data communicator for MPI communication
     * @return A vector containing the sizes of data for each process
     */
    static std::vector<int> SynchronizeRecvSizes(
        const int NumberOfPoints,
        const bool AllPointsAreTheSame,
        const DataCommunicator& rDataCommunicator
        )
    {
        // MPI information
        const int world_size = rDataCommunicator.Size();

        // Define recv_sizes
        std::vector<int> recv_sizes(world_size, 0);
        if (!AllPointsAreTheSame) { // If not all points are the same
            // Getting global number of points
            std::vector<int> points_per_partition(world_size);
            std::vector<int> send_points_per_partition(1, NumberOfPoints);
            rDataCommunicator.AllGather(send_points_per_partition, points_per_partition);

            // Generate vectors with sizes for AllGatherv
            for (int i_rank = 0; i_rank < world_size; ++i_rank) {
                recv_sizes[i_rank] = points_per_partition[i_rank];
            }
        } else if (world_size == 1) { // In case only one process is used
            recv_sizes[0] = NumberOfPoints;
        }

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

        // Synchonize radius
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

            // Invoque AllGatherv
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
    static bool CheckAllPointsAreTheSameAndCalculateNumberOfPoints(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        int& rNumberOfPoints,
        int& rTotalNumberOfPoints,
        const DataCommunicator& rDataCommunicator
        )
    {
        // Get the World Size in MPI
        const int world_size = rDataCommunicator.Size();

        // Getting local number of points
        rNumberOfPoints = std::distance(itPointBegin, itPointEnd);

        // First we check if all the partitions have the same number of points
        rTotalNumberOfPoints = rDataCommunicator.SumAll(rNumberOfPoints);
        if (rNumberOfPoints * world_size != rTotalNumberOfPoints) {
            return false;
        }

        // Now we check if all the points are the same (we are going to do a gross assumtion that the sum of coordinates in al partitions does not compensate between them)
        double x_sum, y_sum, z_sum;
        array_1d<double, 3> coordinates;
        bool all_points_are_the_same = true;
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++) {
            noalias(coordinates) = it_point->Coordinates();
            x_sum = rDataCommunicator.SumAll(coordinates[0]);
            if (std::abs(coordinates[0] - (x_sum/world_size)) > ZeroTolerance) all_points_are_the_same = false;
            all_points_are_the_same = rDataCommunicator.AndReduceAll(all_points_are_the_same);
            if (!all_points_are_the_same) break;
            y_sum = rDataCommunicator.SumAll(coordinates[1]);
            if (std::abs(coordinates[1] - (y_sum/world_size)) > ZeroTolerance) all_points_are_the_same = false;
            all_points_are_the_same = rDataCommunicator.AndReduceAll(all_points_are_the_same);
            if (!all_points_are_the_same) break;
            z_sum = rDataCommunicator.SumAll(coordinates[2]);
            if (std::abs(coordinates[2] - (z_sum/world_size)) > ZeroTolerance) all_points_are_the_same = false;
            all_points_are_the_same = rDataCommunicator.AndReduceAll(all_points_are_the_same);
            if (!all_points_are_the_same) break;
        }

        return all_points_are_the_same;
    }
    ///@}
};

}  // namespace Kratos