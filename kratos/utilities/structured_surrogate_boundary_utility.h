//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#pragma once

// System includes
#include <array>
#include <cstddef>
#include <vector>

// Project includes
#include "includes/define.h"

namespace Kratos {

/**
 * @class StructuredSurrogateBoundaryUtility
 * @brief Dimension-independent topology kernel for structured surrogate meshes.
 * @details The utility deliberately knows nothing about NURBS, finite elements,
 *          ModelParts or boundary-condition formulations. It classifies cells
 *          through a caller-provided predicate and extracts oriented faces
 *          between active and inactive cells. IGA and FEM modelers can therefore
 *          share the same topology without sharing their entity construction.
 */
template <std::size_t TDim>
class StructuredSurrogateBoundaryUtility
{
public:
    static_assert(TDim == 2 || TDim == 3, "Only 2D and 3D structured grids are supported.");

    using IndexArrayType = std::array<std::size_t, TDim>;

    struct BoundaryFace
    {
        IndexArrayType OwnerCell{};
        std::size_t Direction = 0;
        bool IsPositiveSide = false;
    };

    template<class TCellClassifier>
    static std::vector<bool> ClassifyCells(
        const IndexArrayType& rNumberOfCells,
        TCellClassifier&& rCellClassifier)
    {
        std::vector<bool> active_cells(NumberOfCells(rNumberOfCells), false);
        IndexArrayType index{};
        ClassifyRecursive<0>(active_cells, index, rNumberOfCells, rCellClassifier);
        return active_cells;
    }

    static std::vector<BoundaryFace> ExtractBoundaryFaces(
        const std::vector<bool>& rActiveCells,
        const IndexArrayType& rNumberOfCells,
        const bool IncludeGridBoundary = false)
    {
        KRATOS_ERROR_IF(rActiveCells.size() != NumberOfCells(rNumberOfCells))
            << "Active-cell mask size does not match the structured grid." << std::endl;

        std::vector<BoundaryFace> faces;
        IndexArrayType index{};
        ExtractRecursive<0>(faces, index, rActiveCells, rNumberOfCells, IncludeGridBoundary);
        return faces;
    }

    static std::size_t LinearIndex(
        const IndexArrayType& rIndex,
        const IndexArrayType& rNumberOfCells)
    {
        std::size_t stride = 1;
        std::size_t result = 0;
        for (std::size_t d = 0; d < TDim; ++d) {
            result += rIndex[d] * stride;
            stride *= rNumberOfCells[d];
        }
        return result;
    }

private:
    static std::size_t NumberOfCells(const IndexArrayType& rNumberOfCells)
    {
        std::size_t result = 1;
        for (const auto number : rNumberOfCells) {
            result *= number;
        }
        return result;
    }

    template<std::size_t TDirection, class TCellClassifier>
    static void ClassifyRecursive(
        std::vector<bool>& rActiveCells,
        IndexArrayType& rIndex,
        const IndexArrayType& rNumberOfCells,
        TCellClassifier& rCellClassifier)
    {
        if constexpr (TDirection == TDim) {
            rActiveCells[LinearIndex(rIndex, rNumberOfCells)] = rCellClassifier(rIndex);
        } else {
            for (rIndex[TDirection] = 0; rIndex[TDirection] < rNumberOfCells[TDirection]; ++rIndex[TDirection]) {
                ClassifyRecursive<TDirection + 1>(rActiveCells, rIndex, rNumberOfCells, rCellClassifier);
            }
        }
    }

    template<std::size_t TDirection>
    static void ExtractRecursive(
        std::vector<BoundaryFace>& rFaces,
        IndexArrayType& rIndex,
        const std::vector<bool>& rActiveCells,
        const IndexArrayType& rNumberOfCells,
        const bool IncludeGridBoundary)
    {
        if constexpr (TDirection == TDim) {
            if (!rActiveCells[LinearIndex(rIndex, rNumberOfCells)]) {
                return;
            }

            for (std::size_t direction = 0; direction < TDim; ++direction) {
                for (const bool positive_side : {false, true}) {
                    const bool on_grid_boundary = positive_side
                        ? rIndex[direction] + 1 == rNumberOfCells[direction]
                        : rIndex[direction] == 0;

                    bool neighbour_is_active = false;
                    if (!on_grid_boundary) {
                        IndexArrayType neighbour = rIndex;
                        if (positive_side) {
                            ++neighbour[direction];
                        } else {
                            --neighbour[direction];
                        }
                        neighbour_is_active = rActiveCells[LinearIndex(neighbour, rNumberOfCells)];
                    }

                    if (!neighbour_is_active && (!on_grid_boundary || IncludeGridBoundary)) {
                        rFaces.push_back(BoundaryFace{rIndex, direction, positive_side});
                    }
                }
            }
        } else {
            for (rIndex[TDirection] = 0; rIndex[TDirection] < rNumberOfCells[TDirection]; ++rIndex[TDirection]) {
                ExtractRecursive<TDirection + 1>(rFaces, rIndex, rActiveCells, rNumberOfCells, IncludeGridBoundary);
            }
        }
    }
};

} // namespace Kratos
