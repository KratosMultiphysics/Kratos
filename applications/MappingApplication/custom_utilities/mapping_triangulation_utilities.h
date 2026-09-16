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
#include "geometries/brep_surface.h"
#include "geometries/local_refined_brep_surface.h"
#include "geometries/thb_surface_geometry.h"
#include "mapping_application.h"

namespace Kratos
{
namespace MappingTriangulationUtilities
{

using GeometryType = Geometry<Node>;
using GeometryPointerType = GeometryType::Pointer;
using CoordinatesArrayType = GeometryType::CoordinatesArrayType;
using TriangleType = std::array<CoordinatesArrayType, 3>;
using PolygonType = std::vector<CoordinatesArrayType>;
using BrepSurfaceType =
    BrepSurface<PointerVector<Node>, false, PointerVector<Point>>;
using BrepLoopArrayType = BrepSurfaceType::BrepCurveOnSurfaceLoopArrayType;
using LocalRefinedBrepSurfaceThbType =
    LocalRefinedBrepSurface<PointerVector<Node>, THBSurfaceGeometry<3, PointerVector<Node>>, false, PointerVector<Point>>;

/**
 * @brief Sorts a parametric polygon counter-clockwise and triangulates it as a fan.
 * @details The first sorted vertex is connected to every non-adjacent pair of
 * vertices. The input polygon is reordered in place and is assumed to be convex.
 * Polygons with fewer than three vertices or negligible area are ignored.
 * @param rPolygon Parametric polygon to sort and triangulate.
 * @param rTriangles Generated three-vertex polygons, appended to the container.
 */
void KRATOS_API(MAPPING_APPLICATION) TriangulatePolygonFan(
    PolygonType& rPolygon,
    std::vector<PolygonType>& rTriangles);

/**
 * @brief Triangulates a clipped region, including regions containing holes.
 * @details The Clipper region is decomposed into vertical strips so that outer
 * boundaries and holes can be triangulated without filling excluded regions.
 * @param rTrimmedRegion Integer-coordinate paths describing the retained region.
 * @param Factor Conversion factor between parametric and Clipper coordinates.
 * @param rTriangles Generated triangles in parametric coordinates.
 */
void KRATOS_API(MAPPING_APPLICATION) TriangulateTrimmedRegion(
    const Clipper2Lib::Paths64& rTrimmedRegion,
    const double Factor,
    std::vector<Matrix>& rTriangles);

/**
 * @brief Clips parametric triangles against outer trimming loops and inner holes.
 * @details Trimming curves are tessellated into Clipper paths. Each candidate is
 * intersected with the outer region, inner regions are subtracted, and partially
 * retained polygons are triangulated.
 * @param rCandidateTriangles Parametric triangles to clip.
 * @param rBrepSurface BRep surface providing its outer and inner trimming loops.
 * @param Factor Conversion factor between parametric and Clipper coordinates.
 * @return Triangles covering the retained part of the trimmed surface.
 */
std::vector<PolygonType> KRATOS_API(MAPPING_APPLICATION)
ClipTrianglesWithTrimmingLoops(
    const std::vector<PolygonType>& rCandidateTriangles,
    const BrepSurfaceType& rBrepSurface,
    const double Factor);

/// THB overload of @ref ClipTrianglesWithTrimmingLoops. Same contract, applied to
/// the trimming loops carried by a LocalRefinedBrepSurface (BrepCurveOnLocalRefinedSurface).
std::vector<PolygonType> KRATOS_API(MAPPING_APPLICATION)
ClipTrianglesWithTrimmingLoops(
    const std::vector<PolygonType>& rCandidateTriangles,
    const LocalRefinedBrepSurfaceThbType& rLocalRefinedBrepSurface,
    const double Factor);

/**
 * @brief Intersects a triangle with every overlapping knot-span rectangle.
 * @details Each final intersection polygon is triangulated exactly once. This
 * avoids the artificial diagonals introduced by recursively splitting already
 * triangulated intermediate polygons.
 * @param rOriginalTriangleCoordinates Coordinates of the triangle to subdivide.
 * @param pMasterGeometry Geometry providing the knot spans in local space.
 * @param rNewTriangles Resulting triangles, replacing any existing contents.
 */
void KRATOS_API(MAPPING_APPLICATION) Triangulation(
    const TriangleType& rOriginalTriangleCoordinates,
    GeometryPointerType pMasterGeometry,
    std::vector<TriangleType>& rNewTriangles);

/**
 * @brief Triangulate a parametric triangle against an explicit set of active
 * cell rectangles.
 * @details Same clipping and area-consistency logic as @ref Triangulation, but
 * the cell rectangles are provided directly (typically the return value of
 * LocalRefinedBrepSurface::GetActiveCells() for THB patches).
 * @param rOriginalTriangleCoordinates Coordinates of the triangle to subdivide.
 * @param rActiveCells Active cell rectangles as @c (u0, u1, v0, v1) tuples.
 * @param rNewTriangles Resulting triangles, replacing any existing contents.
 */
void KRATOS_API(MAPPING_APPLICATION) TriangulationAgainstActiveCells(
    const TriangleType& rOriginalTriangleCoordinates,
    const std::vector<std::array<double, 4>>& rActiveCells,
    std::vector<TriangleType>& rNewTriangles);

} // namespace MappingTriangulationUtilities
} // namespace Kratos
