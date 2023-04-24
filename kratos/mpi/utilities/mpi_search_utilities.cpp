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
#include "mpi.h"

// Project includes
#include "includes/exception.h"
#include "mpi/utilities/mpi_search_utilities.h"

namespace Kratos
{
void MPISearchUtilities::ComputeBoundingBoxesWithTolerance(
    const std::vector<double>& rBoundingBoxes,
    const double Tolerance,
    std::vector<double>& rBoundingBoxesWithTolerance
    )
{
    const SizeType size_vec = rBoundingBoxes.size();

    KRATOS_DEBUG_ERROR_IF_NOT(std::fmod(size_vec, 6) == 0) << "Bounding Boxes size has to be a multiple of 6!" << std::endl;

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

} // namespace Kratos