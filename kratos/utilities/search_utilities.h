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
#include "geometries/bounding_box.h"
#include "geometries/point.h"

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
 * @class SearchUtilities
 * @ingroup KratosMPI
 * @brief MPI utilities for searching geometrical objects
 * @details Original implementation from MappingUtilities 
 * @author Philipp Bucher (moved by Vicente Mataix Ferrandiz)
 */
class SearchUtilities
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
     * @details Bounding box class implementation
     * @param rBoundingBox The bounding box
     * @param rCoords The point
     * @return true if the point is inside the bounding box
     */
    static bool PointIsInsideBoundingBox(
        const BoundingBox<Point>& rBoundingBox,
        const array_1d<double, 3>& rCoords
        )
    {
        // Get the bounding box points
        const auto& r_max_point = rBoundingBox.GetMaxPoint();
        const auto& r_min_point = rBoundingBox.GetMinPoint();

        // The Bounding Box check
        if (rCoords[0] < r_max_point[0] && rCoords[0] > r_min_point[0])           // check x-direction
            if (rCoords[1] < r_max_point[1] && rCoords[1] > r_min_point[1])       // check y-direction
                if (rCoords[2] < r_max_point[2] && rCoords[2] > r_min_point[2])   // check z-direction
                    return true;
        return false;
    }

    /**
     * @brief Check if a point is inside a bounding box
     * @details Bounding box array of 6 doubles implementation
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
        if (rCoords[0] < rBoundingBox[0] && rCoords[0] > rBoundingBox[1])           // check x-direction
            if (rCoords[1] < rBoundingBox[2] && rCoords[1] > rBoundingBox[3])       // check y-direction
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
     * @brief This method checks if a point is inside any bounding box of the global bounding boxes considering a certain tolerance
     * @param rBoundingBox The bounding box
     * @param rCoords The coordinates of the point
     * @param Tolerance The tolerance
     * @return True if the point is inside the bounding box
     */
    static bool PointIsInsideBoundingBoxWithTolerance(
        const BoundingBox<Point>& rBoundingBox,
        const array_1d<double, 3>& rCoords,
        const double Tolerance
        )
    {
        // Get the bounding box points
        auto max_point = rBoundingBox.GetMaxPoint();
        auto min_point = rBoundingBox.GetMinPoint();
        
        // Apply Tolerances (only in non zero BB cases)
        const double epsilon = std::numeric_limits<double>::epsilon();
        if (norm_2(max_point) > epsilon && norm_2(min_point) > epsilon) {
            for (unsigned int i=0; i<3; ++i) {
                max_point[i] += Tolerance;
                min_point[i] -= Tolerance;
            }
        }

        // The Bounding Box check
        if (rCoords[0] < max_point[0] && rCoords[0] > min_point[0])           // check x-direction
            if (rCoords[1] < max_point[1] && rCoords[1] > min_point[1])       // check y-direction
                if (rCoords[2] < max_point[2] && rCoords[2] > min_point[2])   // check z-direction
                    return true;
        return false;
    }

    /**
     * @brief Compute the bounding boxes of the given bounding boxes from a given tolerance, additionally checking if the bounding boxes are initialized
     * @details This method is used when the bounding boxes are not initialized
     * @param rBoundingBoxes The bounding boxes
     * @param Tolerance The tolerance
     * @param rBoundingBoxesWithTolerance The resulting bounding boxes with the applied tolerance
     */
    static void ComputeBoundingBoxesWithToleranceCheckingNullBB(
        const std::vector<double>& rBoundingBoxes,
        const double Tolerance,
        std::vector<double>& rBoundingBoxesWithTolerance
        );

    ///@}
};

} // namespace Kratos