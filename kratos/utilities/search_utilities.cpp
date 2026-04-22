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

void SearchUtilities::ComputeBoundingBoxesWithTolerance(
    const std::vector<BoundingBox<Point>>& rBoundingBoxes,
    const double Tolerance,
    std::vector<BoundingBox<Point>>& rBoundingBoxesWithTolerance
    )
{
    const SizeType size_vec = rBoundingBoxes.size();

    if (rBoundingBoxesWithTolerance.size() != size_vec) {
        rBoundingBoxesWithTolerance.resize(size_vec);
    }

    // Apply Tolerances
    for (IndexType i = 0; i < size_vec; i++) {
        const auto& r_bounding_box = rBoundingBoxes[i];
        const auto& r_min_point = r_bounding_box.GetMinPoint();
        const auto& r_max_point = r_bounding_box.GetMaxPoint();
        auto& r_bb_with_tolerance = rBoundingBoxesWithTolerance[i];
        auto& r_min_point_with_tolerance = r_bb_with_tolerance.GetMinPoint();
        auto& r_max_point_with_tolerance = r_bb_with_tolerance.GetMaxPoint();
        for (unsigned int j = 0; j < 3; ++j) {
            r_min_point_with_tolerance[j] = r_min_point[j] - Tolerance;
            r_max_point_with_tolerance[j] = r_max_point[j] + Tolerance;
        }
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

/***********************************************************************************/
/***********************************************************************************/

void SearchUtilities::ComputeBoundingBoxesWithToleranceCheckingNullBB(
    const std::vector<BoundingBox<Point>>& rBoundingBoxes,
    const double Tolerance,
    std::vector<BoundingBox<Point>>& rBoundingBoxesWithTolerance
    )
{
    const SizeType size_vec = rBoundingBoxes.size();

    if (rBoundingBoxesWithTolerance.size() != size_vec) {
        rBoundingBoxesWithTolerance.resize(size_vec);
    }

    // Apply Tolerances if the BB is not null
    const double zero_tolerance = std::numeric_limits<double>::epsilon();
    const auto& r_first_bb = rBoundingBoxes[0];
    const auto& r_first_bb_min_point = r_first_bb.GetMinPoint();
    const auto& r_first_bb_max_point = r_first_bb.GetMaxPoint();
    if (std::abs(r_first_bb_max_point[0] - r_first_bb_min_point[0]) > zero_tolerance && 
        std::abs(r_first_bb_max_point[1] - r_first_bb_min_point[1]) > zero_tolerance &&
        std::abs(r_first_bb_max_point[2] - r_first_bb_min_point[2]) > zero_tolerance) {
        for (IndexType i = 0; i < size_vec; i++) {
            const auto& r_bounding_box = rBoundingBoxes[i];
            const auto& r_min_point = r_bounding_box.GetMinPoint();
            const auto& r_max_point = r_bounding_box.GetMaxPoint();
            auto& r_bb_with_tolerance = rBoundingBoxesWithTolerance[i];
            auto& r_min_point_with_tolerance = r_bb_with_tolerance.GetMinPoint();
            auto& r_max_point_with_tolerance = r_bb_with_tolerance.GetMaxPoint();
            for (unsigned int j = 0; j < 3; ++j) {
                r_min_point_with_tolerance[j] = r_min_point[j] - Tolerance;
                r_max_point_with_tolerance[j] = r_max_point[j] + Tolerance;
            }
        }
    }
}

} // namespace Kratos