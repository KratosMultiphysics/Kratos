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

// System includes
#include <array>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "mapping_application.h"

namespace Kratos
{
namespace MappingTriangulationUtilities
{

using GeometryType = Geometry<Node>;
using GeometryPointerType = GeometryType::Pointer;
using CoordinatesArrayType = GeometryType::CoordinatesArrayType;
using TriangleType = std::array<CoordinatesArrayType, 3>;

/**
 * @brief Intersects a triangle with every overlapping knot-span rectangle.
 * @details Each final intersection polygon is triangulated exactly once. This
 * avoids the artificial diagonals introduced by recursively splitting already
 * triangulated intermediate polygons.
 */
void KRATOS_API(MAPPING_APPLICATION) Triangulation(
    const TriangleType& rOriginalTriangleCoordinates,
    GeometryPointerType pMasterGeometry,
    std::vector<TriangleType>& rNewTriangles);

} // namespace MappingTriangulationUtilities
} // namespace Kratos
