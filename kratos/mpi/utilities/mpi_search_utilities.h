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
     * @param AllPointsAreTheSame If all points are the same (-1 to be computed, 0 if not, 1 if yes)
     * @param NumberOfPoints The number of points in the range
     * @param TotalNumberOfPoints The total number of points in all the ranges
     * @return The number of points in the range
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static int MPISynchronousPointSynchronization(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        const DataCommunicator& rDataCommunicator,
        int AllPointsAreTheSame = -1,
        int NumberOfPoints = -1,
        int TotalNumberOfPoints = -1
        )
    {
        // First check that the points are the same in all processes
        if (AllPointsAreTheSame < 0) {
            AllPointsAreTheSame = static_cast<int>(CheckAllPointsAreTheSame(itPointBegin, itPointEnd, NumberOfPoints, TotalNumberOfPoints, rDataCommunicator));
        }

        KRATOS_DEBUG_ERROR_IF(NumberOfPoints < 0) << "The number of points is negative" << std::endl;
        KRATOS_DEBUG_ERROR_IF(TotalNumberOfPoints < 0) << "The total number of points is negative" << std::endl;

        // If all points are the same
        if (AllPointsAreTheSame) {
            rAllPointsCoordinates.resize(NumberOfPoints * 3);
            std::size_t counter = 0;
            array_1d<double, 3> coordinates;
            unsigned int i_coord;
            for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++) {
                noalias(coordinates) = it_point->Coordinates();
                for (i_coord = 0; i_coord < 3; ++i_coord) {
                    rAllPointsCoordinates[3 * counter + i_coord] = coordinates[i_coord];
                }
                ++counter;
            }
        } else { // If not
            // MPI information
            const int world_size = rDataCommunicator.Size();

            // Getting global number of points
            std::vector<int> points_per_partition(world_size);
            std::vector<int> send_points_per_partition(1, NumberOfPoints);
            rDataCommunicator.AllGather(send_points_per_partition, points_per_partition);

            // Getting global coordinates
            rAllPointsCoordinates.resize(TotalNumberOfPoints * 3);
            std::vector<double> send_points_coordinates(NumberOfPoints * 3);
            std::size_t counter = 0;
            array_1d<double, 3> coordinates;
            unsigned int i_coord;
            for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++) {
                noalias(coordinates) = it_point->Coordinates();
                for (i_coord = 0; i_coord < 3; ++i_coord) {
                    send_points_coordinates[3 * counter + i_coord] = coordinates[i_coord];
                }
                ++counter;
            }

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

        return NumberOfPoints;
    }

    /**
     * @brief MPISynchronousPointSynchronization prepares synchronously the coordinates of the points for MPI search.
     * @details With radius
     * @param itPointBegin Iterator to the beginning of the points range
     * @param itPointEnd Iterator to the end of the points range
     * @param rAllPointsCoordinates vector where the computed coordinates will be stored
     * @return  The resulting whole radius vector
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    static std::vector<double> MPISynchronousPointSynchronizationWithDistances(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        const std::vector<double>& rRadius,
        const DataCommunicator& rDataCommunicator
        )
    {
        // First check that the points are the same in all processes
        int number_of_points, total_number_of_points;
        const bool all_points_are_the_same = CheckAllPointsAreTheSame(itPointBegin, itPointEnd, number_of_points, total_number_of_points, rDataCommunicator);

        // Synchronize points
        MPISynchronousPointSynchronization(itPointBegin, itPointEnd, rAllPointsCoordinates, rDataCommunicator, static_cast<int>(all_points_are_the_same), number_of_points, total_number_of_points);

        // Synchonize radius
        if (all_points_are_the_same) { // If all points are the same
            return rRadius;
        } else {                       // If not
            // The resulting distances
            std::vector<double> all_points_distances(total_number_of_points);

            // MPI information
            const int world_size = rDataCommunicator.Size();

            // Getting global number of points
            std::vector<int> points_per_partition(world_size);
            std::vector<int> send_points_per_partition(1, number_of_points);
            rDataCommunicator.AllGather(send_points_per_partition, points_per_partition);

            // Generate vectors with sizes for AllGatherv
            std::vector<int> recv_sizes(world_size, 0);
            for (int i_rank = 0; i_rank < world_size; ++i_rank) {
                recv_sizes[i_rank] = points_per_partition[i_rank];
            }
            std::vector<int> recv_offsets(world_size, 0);
            for (int i_rank = 1; i_rank < world_size; ++i_rank) {
                recv_offsets[i_rank] = recv_offsets[i_rank - 1] + recv_sizes[i_rank - 1];
            }

            // Invoque AllGatherv
            rDataCommunicator.AllGatherv(rRadius, all_points_distances, recv_sizes, recv_offsets);

            return all_points_distances;
        }
    }

    ///@}
private:
    ///@name Private Operations
    ///@{  

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
    static bool CheckAllPointsAreTheSame(
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