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
 * @class MPISearchUtilities
 * @ingroup KratosMPI
 * @brief MPI utilities for searching geometrical objects
 * @details Original implementation fron MappingUtilities 
 * @author Philipp Bucher
 */
class MPISearchUtilities
{
public:
    ///@name Type Definitions
    ///@{

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

    ///@}
};

} // namespace Kratos