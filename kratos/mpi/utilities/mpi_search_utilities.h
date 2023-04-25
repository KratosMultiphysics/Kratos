//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Vicente Mataix Ferrandiz

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "containers/array_1d.h"

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
 * @class MPISearchData
 * @ingroup KratosMPI
 * @brief MPI search data
 * @details Original implementation from MappingUtilities
 * @author Philipp Bucher (moved by Vicente Mataix Ferrandiz)
 */
struct MPISearchData
{
    ///@name Type Definitions
    ///@{

    /// The buffer type for doubles
    using BufferTypeDouble = std::vector<std::vector<double>>;

    /// The buffer type for integers (char)
    using BufferTypeChar = std::vector<std::vector<char>>;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initializes the MPISearchData object by setting up the necessary buffers
     * and obtaining the rank and size of the MPI communicator.
     */
    void Initialize();

    ///@}
    ///@name Member Variables
    ///@{

    int CommRank;                      /// The rank of the current processor
    int CommSize;                      /// The size of the communicator

    std::vector<int> SendSizes;        /// The size of the send buffer
    std::vector<int> RecvSizes;        /// The size of the receive buffer

    BufferTypeDouble SendBufferDouble; /// The send buffer (double)
    BufferTypeDouble RecvBufferDouble; /// The receive buffer (double)

    BufferTypeChar SendBufferChar;     /// The send buffer (char)
    BufferTypeChar RecvBufferChar;     /// The receive buffer (char)

    ///@}
};

/**
 * @class MPISearchUtilities
 * @ingroup KratosMPI
 * @brief MPI utilities for searching geometrical objects
 * @details Original implementation from MappingUtilities 
 * @author Philipp Bucher (moved by Vicente Mataix Ferrandiz)
 */
class MPISearchUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// The Bounding Box type
    using BoundingBoxType = std::array<double, 6>;

    /// The index type definition
    using IndexType = std::size_t;

    /// The size type definition
    using SizeType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Check if a point is inside a bounding box
     * @param rBoundingBox The bounding box
     * @param rCoords The point
     * @return true if the point is inside the bounding box
     */
    static bool PointIsInsideBoundingBox(
        const BoundingBoxType& rBoundingBox,
        const array_1d<double, 3>& rCoords
        )
    {
        // The Bounding Box should have some tolerance already!
        if (rCoords[0] < rBoundingBox[0] && rCoords[0] > rBoundingBox[1])   // check x-direction
            if (rCoords[1] < rBoundingBox[2] && rCoords[1] > rBoundingBox[3])   // check y-direction
                if (rCoords[2] < rBoundingBox[4] && rCoords[2] > rBoundingBox[5])   // check z-direction
                    return true;
        return false;
    }

    /**
     * @brief Compute the bounding boxes of the given bounding boxes from a given tolerance
     * @param rBoundingBoxes The bounding boxes
     * @param Tolerance The tolerance
     * @param rBoundingBoxesWithTolerance The resulting bounding boxes with the applied tolerance
     */
    static void ComputeBoundingBoxesWithTolerance(
        const std::vector<double>& rBoundingBoxes,
        const double Tolerance,
        std::vector<double>& rBoundingBoxesWithTolerance
        );

    /**
     * @brief This method exchanges data asynchronously
     * @tparam TDataType The type of data to be exchanged
     * @param rSendBuffer The send buffer
     * @param rRecvBuffer The receive buffer
     * @param rSearchData The search data
     * @todo This method should be moved to the communicator
     */
    template<typename TDataType>
    static int ExchangeDataAsync(
        const std::vector<std::vector<TDataType>>& rSendBuffer,
        std::vector<std::vector<TDataType>>& rRecvBuffer,
        MPISearchData& rSearchData
        );


    ///@}
};

} // namespace Kratos