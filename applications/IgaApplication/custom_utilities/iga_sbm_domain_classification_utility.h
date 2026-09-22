//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "containers/pointer_vector.h"
#include "geometries/geometry.h"
#include "geometries/brep_surface.h"
#include "geometries/point.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Classifies the knot spans of a trimmed NURBS/Brep surface for domain-level
/// SBM (excluding cut/outside elements from assembly), and extracts the
/// resulting surrogate boundary Gamma_tilde_h.
/**
*   General classification rule (arbitrary trim curve, not specialized to
*   straight/axis-aligned trims): for each knot span S (a rectangle in
*   parametric space), clip the trim boundary polygon against S using the
*   SAME Clipper2 rectangle-clip primitive Kratos's own trimmed-integration
*   pipeline uses internally (kratos/utilities/geometry_utilities/
*   brep_trimming_utilities.cpp), then classify by area ratio:
*     clip_area/span_area < 1e-6       -> OUTSIDE (no material at all)
*     clip_area/span_area > 1-1e-6     -> ACTIVE  (untouched by the trim)
*     otherwise                         -> CUT     (trim boundary passes through S)
*/
class KRATOS_API(IGA_APPLICATION) IgaSbmDomainClassificationUtility
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;
    using SizeType = std::size_t;

    using BrepSurfaceType = BrepSurface<PointerVector<Node>, false, PointerVector<Point>>;

    enum class KnotSpanClassification
    {
        Active = 0,
        Cut = 1,
        Outside = 2
    };

    struct SurrogateBoundarySegment
    {
        array_1d<double, 3> Start;
        array_1d<double, 3> End;
        IndexType ActiveSpanIndexU;
        IndexType ActiveSpanIndexV;
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Classifies every knot span of rBrepSurface as Active/Cut/Outside.
    * @param rBrepSurface The (possibly trimmed) patch surface
    * @param rSpansU Output: knot span breakpoints in U (size = n_spans_u + 1), from rBrepSurface.KnotsU()
    * @param rSpansV Output: knot span breakpoints in V (size = n_spans_v + 1), from rBrepSurface.KnotsV()
    * @param rClassification Output: (n_spans_u x n_spans_v) matrix of KnotSpanClassification (stored as int)
    */
    static void ClassifyKnotSpans(
        const BrepSurfaceType& rBrepSurface,
        std::vector<double>& rSpansU,
        std::vector<double>& rSpansV,
        DenseMatrix<int>& rClassification);

    /**
    * @brief Terminates the simulation with a clear error if no ACTIVE span
    *        exists 
    * @param rClassification As returned by ClassifyKnotSpans
    * @param rPatchName Human-readable name of the patch, for the error message
    */
    static void ValidateNonEmptyActiveDomain(
        const DenseMatrix<int>& rClassification,
        const std::string& rPatchName);

    /**
    * @brief Extracts surrogate boundary: interior knot-line edges bordering exactly
    *        one ACTIVE span 
    * @param rSpansU Knot span breakpoints in U, as returned by ClassifyKnotSpans
    * @param rSpansV Knot span breakpoints in V, as returned by ClassifyKnotSpans
    * @param rClassification As returned by ClassifyKnotSpans
    * @return The surrogate boundary segments, in the patch's own PARAMETRIC space
    */
    static std::vector<SurrogateBoundarySegment> ComputeSurrogateBoundary(
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        const DenseMatrix<int>& rClassification);

    struct ClassificationResult
    {
        std::vector<double> SpansU;
        std::vector<double> SpansV;
        Matrix Classification;
    };

    /**
    * @brief Python-facing entry point: classifies the knot spans of rGeometry.
    * @param rGeometry The candidate geometry 
    * @return SpansU/SpansV/Classification
    */
    static ClassificationResult ClassifyKnotSpansFromGeometry(const Geometry<Node>& rGeometry);

    /**
    * @brief Python-facing entry point: extracts surrogate boundary from a
    *        ClassificationResult 
    * @param rSpansU As in ClassificationResult::SpansU
    * @param rSpansV As in ClassificationResult::SpansV
    * @param rClassification As in ClassificationResult::Classification 
    * @return The surrogate boundary segments
    */
    static std::vector<SurrogateBoundarySegment> ComputeSurrogateBoundaryFromClassification(
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        const Matrix& rClassification);


    static void ValidateNonEmptyActiveDomainFromClassification(
        const Matrix& rClassification,
        const std::string& rPatchName);

    /**
    * @brief The element's own single integration point
    * @param rElement Element with exactly one integration point 
    */
    static array_1d<double, 3> GetElementParametricPosition(const Element& rElement);

    /// One surrogate boundary quadrature point
    struct SurrogateBoundaryQuadraturePointInfo
    {
        array_1d<double, 3> ParametricPosition;
        array_1d<double, 3> PhysicalPosition;
        std::vector<IndexType> NodeIds;
        Vector ShapeFunctionValues;
    };

    /**
    * @brief Builds real, quadrature-bearing geometries along the given
    *        surrogate boundary segments 
    * @param rSegments surrogate boundary segments, as returned by ComputeSurrogateBoundary
    * @param rGeometry The patch's own BrepSurface
    * @param ShapeFunctionDerivativesOrder 
    * @return One SurrogateBoundaryQuadraturePointInfo per quadrature point, across all segments
    */
    static std::vector<SurrogateBoundaryQuadraturePointInfo> CreateSurrogateBoundaryQuadraturePoints(
        const std::vector<SurrogateBoundarySegment>& rSegments,
        const Geometry<Node>& rGeometry,
        const SizeType ShapeFunctionDerivativesOrder);

    /**
    * @brief Creates and adds real Condition instances
    * @param rSegments surrogate boundary segments, 
    * @param rSpansU Knot span breakpoints in U
    * @param rSpansV Knot span breakpoints in V
    * @param rGeometry The patch's own BrepSurface
    * @param ShapeFunctionDerivativesOrder
    * @param rConditionName Registered condition prototype name
    * @param rTargetModelPart Where the new conditions are added
    * @param pProperties Properties assigned to each new condition
    * @param StartId First id to use 
    * @return Number of conditions created
    */
    static SizeType CreateSurrogateBoundaryConditions(
        const std::vector<SurrogateBoundarySegment>& rSegments,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        const Geometry<Node>& rGeometry,
        const SizeType ShapeFunctionDerivativesOrder,
        const std::string& rConditionName,
        ModelPart& rTargetModelPart,
        Properties::Pointer pProperties,
        const IndexType StartId);

    ///@}

}; // Class IgaSbmDomainClassificationUtility

///@}

}  // namespace Kratos.
