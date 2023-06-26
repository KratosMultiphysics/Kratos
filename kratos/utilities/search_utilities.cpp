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

// System includes
#include <cmath>

// External includes

// Project includes
#include "utilities/search_utilities.h"

namespace Kratos
{

bool SearchUtilities::PointIsInsideBoundingBox(
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

/***********************************************************************************/
/***********************************************************************************/

bool SearchUtilities::PointIsInsideBoundingBox(
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

/***********************************************************************************/
/***********************************************************************************/

void SearchUtilities::ComputeBoundingBoxesWithTolerance(
    const std::vector<double>& rBoundingBoxes,
    const double Tolerance,
    std::vector<double>& rBoundingBoxesWithTolerance
    )
{
    const SizeType size_vec = rBoundingBoxes.size();

    KRATOS_ERROR_IF_NOT(std::fmod(size_vec, 6) == 0) << "Bounding Boxes size has to be a multiple of 6!" << std::endl;

    if (rBoundingBoxesWithTolerance.size() != size_vec) {
        rBoundingBoxesWithTolerance.resize(size_vec);
    }

    // Apply Tolerances
    for (IndexType i=0; i<size_vec; i+=2) {
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] + Tolerance;
    }

    for (IndexType i=1; i<size_vec; i+=2) {
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] - Tolerance;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SearchUtilities::ComputeBoundingBoxesWithToleranceCheckingNullBB(
    const std::vector<double>& rBoundingBoxes,
    const double Tolerance,
    std::vector<double>& rBoundingBoxesWithTolerance
    )
{
    const SizeType size_vec = rBoundingBoxes.size();

    KRATOS_ERROR_IF_NOT(std::fmod(size_vec, 6) == 0) << "Bounding Boxes size has to be a multiple of 6!" << std::endl;

    if (rBoundingBoxesWithTolerance.size() != size_vec) {
        rBoundingBoxesWithTolerance.resize(size_vec);
    }

    // Apply Tolerances if the BB is not null
    const double zero_tolerance = std::numeric_limits<double>::epsilon();
    if (std::abs(rBoundingBoxes[0] - rBoundingBoxes[1]) > zero_tolerance && 
        std::abs(rBoundingBoxes[2] - rBoundingBoxes[3]) > zero_tolerance &&
        std::abs(rBoundingBoxes[4] - rBoundingBoxes[5]) > zero_tolerance) {
        for (IndexType i=0; i<size_vec; i+=2) {
            rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] + Tolerance;
        }

        for (IndexType i=1; i<size_vec; i+=2) {
            rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] - Tolerance;
        }
    }
}

} // namespace Kratos