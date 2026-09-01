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

Clipper2Lib::Paths64 KRATOS_API(MAPPING_APPLICATION) CreateAllLoops(
    const BrepSurface<PointerVector<Node>, false, PointerVector<Point>>& rBrepSurface,
    const double Factor,
    const double TessellationTolerance = 1.0e-3);

void KRATOS_API(MAPPING_APPLICATION) SplitOuterAndInnerPaths(
    const Clipper2Lib::Paths64& rAllLoops,
    Clipper2Lib::Paths64& rOuterPaths,
    Clipper2Lib::Paths64& rInnerPaths);

Clipper2Lib::Paths64 KRATOS_API(MAPPING_APPLICATION) ClipPathsWithTrimmedDomain(
    const Clipper2Lib::Paths64& rSubjectPaths,
    const Clipper2Lib::Paths64& rOuterPaths,
    const Clipper2Lib::Paths64& rInnerPaths);

} // namespace Kratos::BrepClipperUtilities
