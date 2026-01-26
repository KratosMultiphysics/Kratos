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
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

#include "geometries/coupling_geometry.h"
#include "geometries/brep_surface.h"
#include "spatial_containers/bins_dynamic.h"

#include "integration/integration_point_utilities.h"
#include "utilities/quadrature_points_utility.h"

namespace Kratos
{
namespace IgaMappingIntersectionUtilities
{
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using GeometryPointerType = GeometryType::Pointer;

    using GeometriesArrayType = typename GeometryType::GeometriesArrayType;
    using CoordinatesArrayType = typename GeometryType::CoordinatesArrayType;
    using IntegrationPointsArrayType = typename GeometryType::IntegrationPointsArrayType;

    using PointType = Node;
    using PointTypePointer = Node::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;

    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;

    using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    using PointerType = DynamicBins::PointerType;
    struct PatchSearchCache
    {
        std::vector<PointTypePointer> points;
        std::unique_ptr<DynamicBins> p_bins;

        // Regular sampling info (needed to reconstruct xi/eta from Point Id)
        IndexType number_pts = 0;   // = n_div+1
        double u_0 = 0.0, v_0 = 0.0;
        double delta_u = 0.0, delta_v = 0.0;
    };

    using PatchCacheMap = std::unordered_map<IndexType, PatchSearchCache>;
    using cInt = int64_t;

        /**
     * @brief Builds a spatial cache (sampling + search bins) for a set of IGA BRep patches.
     *
     * This function samples each IGA patch surface in physical space using a structured grid in its
     * parametric space (u,v). The sampled points are stored together with a spatial search structure
     * (DynamicBins), enabling fast radius queries.
     *
     * The cache is typically used to:
     *  - quickly determine which patches are near a given FEM point/element/condition
     *  - provide good initial guesses for Newton-Raphson surface projections
     *
     * @param patches_id      List of patch geometry IDs (must exist in @p rModelPartIga as geometries)
     * @param rModelPartIga   ModelPart containing the IGA patch geometries (BrepSurface)
     * @param n_div           Number of divisions per parametric direction (sampling resolution = (n_div+1)^2)
     * @return PatchCacheMap  Map patch_id -> PatchSearchCache (sampled points + bins + metadata)
     */
    PatchCacheMap KRATOS_API(MAPPING_APPLICATION) BuildPatchCaches(
        const std::vector<IndexType>& patches_id,
        const ModelPart& rModelPartIga,
        const IndexType n_div);


    /**
     * @brief Creates line coupling geometries between an IGA trimming curve and a FEM curve interface.
     *
     * This utility creates a set of @ref CouplingGeometry objects for line-to-line coupling.
     * Each coupling geometry is composed of:
     *  - an IGA curve on surface (BrepCurveOnSurface) representing the IGA interface
     *  - a FEM NurbsCurveGeometry representing the FEM interface
     *
     * The output coupling geometries are stored inside @p rModelPartResult and can later be used
     * to generate quadrature points (integration entities) along the shared interface.
     *
     * Typical use-case: line coupling between an IGA boundary curve and a FEM curve for mortar-style mapping.
     *
     * @param rModelPartDomainA   ModelPart containing the origin interface
     * @param rModelPartDomainB   ModelPart containing the destination interface
     * @param rIsOriginIga        If true: DomainA is treated as IGA and DomainB as FEM, otherwise swapped
     * @param rModelPartResult    ModelPart where the generated coupling geometries will be stored
     * @param Tolerance           Geometric tolerance used for intersection/projection decisions
     */
    void KRATOS_API(MAPPING_APPLICATION) CreateIgaFEMCouplingGeometriesOnCurve(
        ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        const bool& rIsOriginIga,
        ModelPart& rModelPartResult,
        double Tolerance = 1.0e-6);


    /**
     * @brief Creates quadrature-point coupling geometries along the IGA-FEM coupling interface curve.
     *
     * This utility takes a ModelPart containing coupling geometries (line coupling) and generates
     * quadrature point entities along the intersection/interface curve. These quadrature points are
     * intended for integration-based coupling (e.g. mortar coupling, consistent force transfer).
     *
     * The procedure typically:
     *  1. Computes integration points along the coupling curve(s)
     *  2. Creates quadrature point geometries on both sides (IGA + FEM)
     *  3. Stores the new quadrature coupling entities (usually as conditions) into the same ModelPart
     *
     * @param rModelPartCoupling   ModelPart containing the previously created line coupling geometries
     * @param Tolerance            Geometric tolerance used during evaluation/projection
     */
    void KRATOS_API(MAPPING_APPLICATION) CreateIgaFEMQuadraturePointsOnCurve(
        ModelPart& rModelPartCoupling,
        double Tolerance);

    
        /**
     * @brief Creates surface-to-surface coupling geometries between an IGA interface and a FEM interface.
     *
     * This utility builds a set of @ref CouplingGeometry objects, where each coupling geometry contains:
     *  - The IGA patch surface (represented by a @ref BrepSurface / NURBS surface)
     *  - A FEM interface condition geometry
     *
     * For efficiency, each FEM condition is only paired with a subset of candidate IGA patches, determined
     * through a spatial pre-filtering (radius-based search on a cached patch sampling).
     *
     * @param rModelPartDomainA       ModelPart containing the origin interface (IGA side if @p is_origin_iga = true)
     * @param rModelPartDomainB       ModelPart containing the destination interface (FEM side if @p is_origin_iga = true)
     * @param rModelPartResult        ModelPart where the generated coupling geometries will be stored
     * @param is_origin_iga           If true, DomainA is treated as IGA and DomainB as FEM. Otherwise the opposite.
     * @param search_radius           Radius used for spatial pre-filtering of candidate IGA patches
     * @param rPatchCache             Patch cache (sampling + bins). It is rebuilt internally for the involved patches.
     */
    void KRATOS_API(MAPPING_APPLICATION) CreateIgaFEMCouplingGeometriesOnSurface(
        ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        ModelPart& rModelPartResult,
        bool is_origin_iga,
        const double search_radius,
        PatchCacheMap& rPatchCache);


    /**
     * @brief Creates quadrature-point coupling geometries from previously generated surface coupling geometries.
     *
     * This method takes a ModelPart containing coupling geometries (IGA surface + FEM geometry) and builds
     * quadrature point geometries for integration-based coupling (mortar-like approach).
     *
     * For each coupling geometry:
     *  1. FEM nodes are projected onto the IGA patch in parametric space (u,v)
     *  2. The projected FEM triangle (or clipped triangles, in case of trimming) is triangulated in parameter space
     *  3. Integration points are generated for each parametric triangle
     *  4. The integration points are mapped to the FEM side by projecting the master quadrature-point coordinates
     *     into the FEM local space
     *  5. Quadrature-point geometries are created on both sides and stored as new coupling conditions
     *
     * @param rModelPartCoupling      ModelPart containing the surface coupling geometries (IGA patch + FEM condition)
     * @param origin_is_iga           If true, master = IGA patch and slave = FEM geometry. Otherwise swapped.
     * @param rPatchCache             Patch cache used for fast initial guess generation for surface projections
     * @param search_radius           Radius used in the cached spatial bins search for initial guess estimation
     */
    void KRATOS_API(MAPPING_APPLICATION) CreateIgaFEMQuadraturePointsOnSurface(
        ModelPart& rModelPartCoupling,
        bool origin_is_iga,
        const PatchCacheMap& rPatchCache,
        const double search_radius);


    /**
     * @brief Filters the set of candidate IGA patches for a given FEM entity location.
     *
     * Given a list of patch IDs and a query point (typically the FEM element/condition center), this function
     * performs a radius-based search on each patch sampling cache and returns only the patches for which at least
     * one sampled point lies within @p search_radius of the query point.
     *
     * If no patch is found within the radius, the function falls back to returning the full input list.
     * (This avoids accidentally discarding valid patches due to too small radius.)
     *
     * @param patches_id              List of patch IDs to check
     * @param rCache                  Cache containing sampled points and spatial bins per patch
     * @param element_center          Query point in physical space (x,y,z)
     * @param search_radius           Search radius in physical space
     * @return std::vector<IndexType> List of patch IDs likely containing a valid projection
     */
    std::vector<IndexType> KRATOS_API(MAPPING_APPLICATION) GetPatchesWithProbableProjection(
        const std::vector<IndexType>& patches_id,
        const PatchCacheMap& rCache,
        const CoordinatesArrayType& element_center,
        const double search_radius);


    /**
     * @brief Provides an initial guess for a Newton-Raphson projection of a FEM point onto an IGA patch.
     *
     * This method finds the closest cached sample point (within a radius) on the patch bins structure and converts
     * it back to the corresponding parametric coordinates (u,v). The guess is returned in @p rInitialGuess and can
     * be passed to ProjectionPointGlobalToLocalSpace / Newton projection.
     *
     * @param rSlaveElementNode       Point in physical space to be projected (x,y,z)
     * @param pMasterGeometry         IGA patch geometry (must have a corresponding entry in @p rPatchCache)
     * @param rPatchCache             Patch cache containing sampled points and mapping to parametric indices
     * @param rInitialGuess           Output initial guess in parametric coordinates (u,v,0)
     * @param SearchRadius            Radius used in the spatial bins query
     * @return true                   If a valid guess was found
     * @return false                  If no cached point was found within the search radius
     */
    bool KRATOS_API(MAPPING_APPLICATION) FindInitialGuessNewtonRaphsonProjection(
        const CoordinatesArrayType& rSlaveElementNode,
        GeometryPointerType pMasterGeometry,
        const PatchCacheMap& rPatchCache,
        CoordinatesArrayType& rInitialGuess,
        const double SearchRadius);


    /**
     * @brief Checks whether the projected points lie on the boundary of the NURBS patch parameter domain.
     *
     * This is typically used as a safeguard to avoid degenerate/unstable triangulations when the projected points
     * lie directly on patch boundaries (u = u_min/u_max or v = v_min/v_max).
     *
     * If one point is provided: returns true if it lies on any boundary line.
     * If two or more points are provided: returns true if the first two points lie on the same boundary line.
     *
     * @param points_to_triangulate   Projected points in parameter space (u,v,0)
     * @param nurbs_surface           The NURBS surface defining the parametric domain limits
     * @return true                   If boundary conditions are detected
     * @return false                  Otherwise
     */
    bool KRATOS_API(MAPPING_APPLICATION) AreProjectionsOnParameterSpaceBoundary(
        const std::vector<CoordinatesArrayType>& points_to_triangulate,
        const NurbsSurfaceGeometry<3, PointerVector<Node>>& nurbs_surface);


    /**
     * @brief Computes the intersection between a line segment and the IGA patch using bisection.
     *
     * Given two physical-space points defining a segment:
     *  - r_point_inside  : point known to be projectable onto the patch
     *  - r_point_outside : point known to fail projection onto the patch
     *
     * The method repeatedly bisects the segment and tests the midpoint by attempting a projection onto the patch.
     * The boundary intersection point is approximated once the segment becomes sufficiently small.
     *
     * @param p_geom_master           Patch geometry onto which the projection is attempted
     * @param r_point_inside          Segment endpoint known to be inside (projectable)
     * @param r_point_outside         Segment endpoint known to be outside (not projectable)
     * @param r_initial_guess         Initial guess for local coordinates (u,v), typically from the inside node
     * @param r_intersection_point    Output intersection point in local patch coordinates (u,v,0)
     * @return true                   Always true (current implementation). Can be extended to signal failure.
     */
    bool KRATOS_API(MAPPING_APPLICATION) FindTriangleSegmentSurfaceIntersectionWithBisection(
        GeometryPointerType p_geom_master,
        const CoordinatesArrayType& r_point_inside,
        const CoordinatesArrayType& r_point_outside,
        const CoordinatesArrayType& r_initial_guess,
        CoordinatesArrayType& r_intersection_point);


    /**
     * @brief Sorts polygon vertices counter-clockwise in (u,v) parametric space.
     *
     * This utility ensures a consistent vertex ordering for downstream operations such as triangulation
     * and signed-area computations.
     *
     * The vertices are sorted by their polar angle around the centroid of the polygon.
     *
     * @param r_triangle_vertices     Vertices in parameter space (u,v,0), modified in-place
     */
    void KRATOS_API(MAPPING_APPLICATION) SortVerticesCounterClockwise(
        std::vector<CoordinatesArrayType>& r_triangle_vertices);


}  // namespace IgaMappingIntersectionUtilities.

}  // namespace Kratos.
