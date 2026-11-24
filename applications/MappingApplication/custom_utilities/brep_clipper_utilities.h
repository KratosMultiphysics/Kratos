//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti

#pragma once

// Project includes
#include "geometries/brep_surface.h"
#include "mapping_application.h"

// External includes
#include "clipper/include/clipper2/clipper.h"

namespace Kratos::BrepClipperUtilities
{

/**
 * @brief Converts all trimming loops of a BRep surface into Clipper paths.
 *
 * The outer trimming loop and all inner loops are tessellated in parameter
 * space and converted to integer coordinates using the specified scaling
 * factor. The outer loop is stored first, followed by the inner loops.
 *
 * @param rBrepSurface BRep surface containing the trimming loops.
 * @param Factor Scaling factor used to convert parameter-space coordinates to integers.
 * @param TessellationTolerance Tolerance used to tessellate the trimming curves.
 * @return Paths representing the outer and inner trimming loops.
 */
Clipper2Lib::Paths64 KRATOS_API(MAPPING_APPLICATION) CreateAllLoops(
    const BrepSurface<PointerVector<Node>, false, PointerVector<Point>>& rBrepSurface,
    const double Factor,
    const double TessellationTolerance = 1.0e-3);

/**
 * @brief Separates the outer trimming path from the inner trimming paths.
 *
 * The first non-empty path is treated as the outer boundary, while the
 * remaining non-empty paths are treated as inner holes.
 *
 * @param rAllLoops Trimming paths with the outer loop stored first.
 * @param rOuterPaths Output container for the outer trimming path.
 * @param rInnerPaths Output container for the inner trimming paths.
 */
void KRATOS_API(MAPPING_APPLICATION) SplitOuterAndInnerPaths(
    const Clipper2Lib::Paths64& rAllLoops,
    Clipper2Lib::Paths64& rOuterPaths,
    Clipper2Lib::Paths64& rInnerPaths);

/**
 * @brief Clips subject paths against a trimmed parameter-space domain.
 *
 * The subject paths are first intersected with the outer boundary. Any regions
 * covered by the inner paths are then removed from the intersection result.
 *
 * @param rSubjectPaths Paths to clip.
 * @param rOuterPaths Paths defining the retained outer domain.
 * @param rInnerPaths Paths defining holes to remove from the retained domain.
 * @return Portions of the subject paths contained in the trimmed domain.
 */
Clipper2Lib::Paths64 KRATOS_API(MAPPING_APPLICATION) ClipPathsWithTrimmedDomain(
    const Clipper2Lib::Paths64& rSubjectPaths,
    const Clipper2Lib::Paths64& rOuterPaths,
    const Clipper2Lib::Paths64& rInnerPaths);

} // namespace Kratos::BrepClipperUtilities
