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
#include "includes/data_communicator.h"
#include "spatial_containers/geometrical_objects_bins.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

class GeometricalObject; // forward declaration, to be included in the cpp. This is needed to reduce the compilation time. Can be done as we consider the GeometricalObject as a pointer

/**
 * @class GeometricalObjectsBinsMPI
 * @ingroup KratosCore
 * @brief A bins container for 3 dimensional GeometricalObject entities (MPI version)
 * @details This is the MPI version of the GeometricalObjectsBins, which is a container for geometrical objects. It is used to perform fast search of geometrical objects in a given space.
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(KRATOS_CORE) GeometricalObjectsBinsMPI
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalObjectsBinsMPI
    KRATOS_CLASS_POINTER_DEFINITION(GeometricalObjectsBinsMPI);

    /// The buffer type definition
    using BufferTypeDouble = std::vector<std::vector<double>>;
    using BufferTypeChar = std::vector<std::vector<char>>;

    /// The type of geometrical object to be stored in the bins
    using CellType = GeometricalObjectsBins::CellType;
    using ResultType = GeometricalObjectsBins::ResultType;

    /// Define zero tolerance as Epsilon
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    GeometricalObjectsBinsMPI() = delete;

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param GeometricalObjectsBegin The begin iterator of the geometries to be stored
     * @param GeometricalObjectsEnd The end iterator of the geometries to be stored
     * @param rDataCommunicator The data communicator
     * @tparam TIteratorType The type of the iterator
     */
    template<typename TIteratorType>
    GeometricalObjectsBinsMPI(
        TIteratorType GeometricalObjectsBegin,
        TIteratorType GeometricalObjectsEnd,
        const DataCommunicator& rDataCommunicator
        ) : mLocalGeometricalObjectsBins(GeometricalObjectsBegin, GeometricalObjectsEnd),
            mrDataCommunicator(rDataCommunicator)
    {
        // We get the world size
        const int world_size = GetWorldSize();

        // Set up the buffers
        mSendSizes.resize(world_size);
        mRecvSizes.resize(world_size);

        mSendBufferDouble.resize(world_size);
        mRecvBufferDouble.resize(world_size);

        mSendBufferChar.resize(world_size);
        mRecvBufferChar.resize(world_size);

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
    GeometricalObjectsBinsMPI(
        TContainer& rGeometricalObjectsVector,
        const DataCommunicator& rDataCommunicator
        ) : GeometricalObjectsBinsMPI(rGeometricalObjectsVector.begin(), rGeometricalObjectsVector.end(), rDataCommunicator)
    {
    }

    /// Destructor.
    ~GeometricalObjectsBinsMPI() = default;

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
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    void SearchInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        std::vector<std::vector<ResultType>>& rResults
        )
    {
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        std::array<int, 2> limits;
        const int number_of_points = PreparePointsForMPISearch(itPointBegin, itPointEnd, all_points_coordinates, limits);
        rResults.resize(number_of_points);

        // Performa the corresponding searchs
        const int lower_limit = limits[0];
        const int upper_limit = limits[1];
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            std::vector<ResultType> result;
            ImplSearchInRadius(point, Radius, result);
            // Added only in the corresponding partition
            if (i_node >= lower_limit && i_node < upper_limit) {
                rResults[i_node - lower_limit] = result;
            }
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius.
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param Radius The radius to be checked
     * @return ResultType The result of the search
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    std::vector<ResultType> SearchNearestInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius
        )
    {
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        std::array<int, 2> limits;
        const int number_of_points = PreparePointsForMPISearch(itPointBegin, itPointEnd, all_points_coordinates, limits);
        std::vector<ResultType> results(number_of_points);

        // Performa the corresponding searchs
        const int lower_limit = limits[0];
        const int upper_limit = limits[1];
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            const auto result = ImplSearchNearestInRadius(point, Radius);
            // Added only in the corresponding partition
            if (i_node >= lower_limit && i_node < upper_limit) {
                results[i_node - lower_limit] = result;
            }
        }
        return results;
    }

    /**
     * @brief This method takes a point and finds the nearest object to it.
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @return ResultType The result of the search
     */
    template<typename TPointIteratorType>
    std::vector<ResultType> SearchNearest(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd
        )
    {
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        std::array<int, 2> limits;
        const int number_of_points = PreparePointsForMPISearch(itPointBegin, itPointEnd, all_points_coordinates, limits);
        std::vector<ResultType> results(number_of_points);

        // Performa the corresponding searchs
        const int lower_limit = limits[0];
        const int upper_limit = limits[1];
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            const auto result = ImplSearchNearest(point);
            // Added only in the corresponding partition
            if (i_node >= lower_limit && i_node < upper_limit) {
                results[i_node - lower_limit] = result;
            }
        }
        return results;
    }

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain (iterative version).
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @return std::vector<ResultType> The result of the search
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    std::vector<ResultType> SearchIsInside(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd
        )
    {
        // Prepare MPI search
        std::vector<double> all_points_coordinates;
        std::array<int, 2> limits;
        const int number_of_points = PreparePointsForMPISearch(itPointBegin, itPointEnd, all_points_coordinates, limits);
        std::vector<ResultType> results(number_of_points);

        // Performa the corresponding searchs
        const int lower_limit = limits[0];
        const int upper_limit = limits[1];
        const int total_number_of_points = all_points_coordinates.size()/3;
        for (int i_node = 0; i_node < total_number_of_points; ++i_node) {
            Point point(all_points_coordinates[i_node * 3 + 0], all_points_coordinates[i_node * 3 + 1], all_points_coordinates[i_node * 3 + 2]);
            const auto result = ImplSearchIsInside(point);
            // Added only in the corresponding partition
            if (i_node >= lower_limit && i_node < upper_limit) {
                results[i_node - lower_limit] = result;
            }
        }
        return results;
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GeometricalObjectsBinsMPI" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GeometricalObjectsBinsMPI";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Local Geometrical Objects Bins for Rank " << GetRank() + 1 << "/" << GetWorldSize() << "\n";
        mLocalGeometricalObjectsBins.PrintData(rOStream);
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

    std::vector<double> mGlobalBoundingBoxes; /// All the global BB, data is xmax, xmin,  ymax, ymin,  zmax, zmin

    /// TODO: Replace with a vector and serialize it to all partitions
    GeometricalObjectsBins mLocalGeometricalObjectsBins; /// The local bins

    const DataCommunicator& mrDataCommunicator; /// The data communicator

    // TODO: Check what is necessary after final implementation on MPI is done
    std::vector<int> mSendSizes; /// The sizes of the send buffers
    std::vector<int> mRecvSizes; /// The sizes of the recv buffers

    BufferTypeDouble mSendBufferDouble; /// The send buffer (double)
    BufferTypeDouble mRecvBufferDouble; /// The recv buffer (double)

    BufferTypeChar mSendBufferChar; /// The send buffer (char)
    BufferTypeChar mRecvBufferChar; /// The recv buffer (char)

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief PreparePointsForMPISearch prepares the points for MPI search.
     * @param itPointBegin Iterator to the beginning of the points range
     * @param itPointEnd Iterator to the end of the points range
     * @param rAllPointsCoordinates vector where the computed coordinates will be stored
     * @param rLimits array with the lower and upper limits of the current partition
     * @return The number of points in the range
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    int PreparePointsForMPISearch(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        std::vector<double>& rAllPointsCoordinates,
        std::array<int, 2>& rLimits
        )
    {
        // First check that the points are the same in all processes
        int number_of_points, total_number_of_points;
        const bool all_points_are_the_same = CheckAllPointsAreTheSame(itPointBegin, itPointEnd, number_of_points, total_number_of_points);

        // If all points are the same
        if (all_points_are_the_same) {
            rAllPointsCoordinates.resize(number_of_points * 3);
            std::size_t counter = 0;
            array_1d<double, 3> coordinates;
            for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++) {
                noalias(coordinates) = it_point->Coordinates();
                rAllPointsCoordinates[3 * counter + 0] = coordinates[0];
                rAllPointsCoordinates[3 * counter + 1] = coordinates[1];
                rAllPointsCoordinates[3 * counter + 2] = coordinates[2];
                ++counter;
            }

            // Define limits
            rLimits[0] = 0;
            rLimits[1] = number_of_points;
        } else { // If not
            // MPI information
            const int rank = GetRank();
            const int world_size = GetWorldSize();

            // Getting global number of points
            std::vector<int> points_per_partition(world_size);
            std::vector<int> send_points_per_partition(1, number_of_points);
            mrDataCommunicator.AllGather(send_points_per_partition, points_per_partition);

            // Getting global coordinates
            rAllPointsCoordinates.resize(total_number_of_points * 3);
            std::vector<double> send_points_coordinates(number_of_points * 3);
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
            std::vector<int> recv_sizes(world_size);
            for (int i_rank = 1; i_rank < world_size; i_rank++) {
                recv_sizes[i_rank] = 3 * points_per_partition[i_rank];
            }
            int message_size = 0;
            std::vector<int> recv_offsets(world_size, 0);
            for (int i_rank = 1; i_rank < world_size; i_rank++) {
                recv_offsets[i_rank] = message_size;
                message_size += recv_sizes[i_rank];
            }

            // Invoque AllGatherv
            mrDataCommunicator.AllGatherv(send_points_coordinates, rAllPointsCoordinates, recv_sizes, recv_offsets);

            // Define limits
            const auto it_point_begin = points_per_partition.begin();
            rLimits[0] = std::reduce(it_point_begin, it_point_begin + rank + 1);
            rLimits[1] = std::reduce(it_point_begin, it_point_begin + rank + 2);
        }

        return total_number_of_points;
    }

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it.
     * @details The result contains the object and also its distance to the point.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     */
    void ImplSearchInRadius(
        const Point& rPoint,
        const double Radius,
        std::vector<ResultType>& rResults
        );

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius.
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @return ResultType The result of the search
     */
    ResultType ImplSearchNearestInRadius(
        const Point& rPoint,
        const double Radius
        );

    /**
     * @brief This method takes a point and finds the nearest object to it.
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @return ResultType The result of the search
    */
    ResultType ImplSearchNearest(const Point& rPoint);

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain.
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param rPoint The point to be checked
     * @return ResultType The result of the search
     */
    ResultType ImplSearchIsInside(const Point& rPoint);

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
     * @brief This method checks if all nodes are the same across all partitions
     * @param itPointBegin Iterator to the beginning of the points range
     * @param itPointEnd Iterator to the end of the points range
     * @param rNumberOfPoints Number of points in the range
     * @param rTotalNumberOfPoints Total number of points in all partitions
     * @return true if all points are the same in all partitions
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    bool CheckAllPointsAreTheSame(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        int& rNumberOfPoints,
        int& rTotalNumberOfPoints
        )
    {
        // Get the World Size in MPI
        const int world_size = GetWorldSize();

        // Getting local number of points
        rNumberOfPoints = std::distance(itPointBegin, itPointEnd);

        // First we check if all the partitions have the same number of points
        rTotalNumberOfPoints = mrDataCommunicator.SumAll(rNumberOfPoints);
        if (rNumberOfPoints * world_size != rTotalNumberOfPoints) {
            return false;
        }

        // Now we check if all the points are the same (we are going to do a gross assumtion that the sum of coordinates in al partitions does not compensate between them)
        double x_sum, y_sum, z_sum;
        array_1d<double, 3> coordinates;
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++) {
            noalias(coordinates) = it_point->Coordinates();
            x_sum = mrDataCommunicator.SumAll(coordinates[0]);
            if (std::abs(coordinates[0] - (x_sum/world_size)) > ZeroTolerance) return false;
            y_sum = mrDataCommunicator.SumAll(coordinates[1]);
            if (std::abs(coordinates[1] - (y_sum/world_size)) > ZeroTolerance) return false;
            z_sum = mrDataCommunicator.SumAll(coordinates[2]);
            if (std::abs(coordinates[2] - (z_sum/world_size)) > ZeroTolerance) return false;
        }

        return true;
    }

    // /**
    //  * @brief This method prepares the buffer for the result
    //  * @details Values are set in member variables
    //  * @param rLocalResult The local result
    //  */
    // void PrepareBufferResultType(const ResultType& rLocalResult);

    // /**
    //  * @brief This method deserializes the buffer for the result
    //  * @details Values are get from member variables
    //  * @param rGlobalResult The global result
    //  */
    // void DeserializeResultType(const ResultType& rGlobalResult);

    /**
     * @brief This method does the serialization/desiralization of the results and rerieves the result from a given partition
     * @param rLocalResult The local result
     * @param Rank The rank in MPI
     * @return The result from the given partition
     */
    ResultType GetResultFromGivenPartition(
        const ResultType& rLocalResult,
        const int Rank
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
    GeometricalObjectsBinsMPI& operator=(GeometricalObjectsBinsMPI const& rOther) = delete;

    /// Copy constructor deleted.
    GeometricalObjectsBinsMPI(GeometricalObjectsBinsMPI const& rOther) = delete;

    ///@}

}; // Class GeometricalObjectsBinsMPI

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@} addtogroup block

}  // namespace Kratos.