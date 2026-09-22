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

// System includes
#include <cmath>
#include <limits>

// External includes
#include "clipper/include/clipper2/clipper.h"

// Project includes
#include "custom_utilities/iga_sbm_domain_classification_utility.h"
#include "utilities/tessellation_utilities/curve_tessellation.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/brep_curve_on_surface.h"
#include "iga_application_variables.h"

namespace Kratos
{

namespace {

using cInt = signed long long;
constexpr double ClipperScaleFactor = 1e-10;

Clipper2Lib::Point64 ToIntPoint(double X, double Y)
{
    Clipper2Lib::Point64 p;
    p.x = static_cast<cInt>(std::round(X / ClipperScaleFactor));
    p.y = static_cast<cInt>(std::round(Y / ClipperScaleFactor));
    return p;
}

// Tessellates one boundary loop (a sequence of BrepCurveOnSurface segments)
// into a single Clipper2 polygon
template<class TLoopType>
Clipper2Lib::Path64 TessellateLoop(const TLoopType& rLoop)
{
    Clipper2Lib::Path64 path;
    Clipper2Lib::Point64 last_point;
    last_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
    last_point.y = static_cast<cInt>(std::numeric_limits<int>::min());

    for (std::size_t j = 0; j < rLoop.size(); ++j) {
        CurveTessellation<PointerVector<Node>> curve_tessellation;
        const auto& r_curve = *(rLoop[j].get());
        curve_tessellation.Tessellate(r_curve, 0.001, 1, true);
        const auto& tessellation = curve_tessellation.GetTessellation();
        for (std::size_t u = 0; u < tessellation.size(); ++u) {
            const double x = std::get<1>(tessellation[u])[0];
            const double y = std::get<1>(tessellation[u])[1];
            const auto new_point = ToIntPoint(x, y);
            if (!(last_point.x == new_point.x && last_point.y == new_point.y)) {
                path.push_back(new_point);
                last_point = new_point;
            }
        }
    }
    return path;
}

double UnionArea(const Clipper2Lib::Paths64& rPaths)
{
    if (rPaths.empty()) {
        return 0.0;
    }
    double area = std::abs(Clipper2Lib::Area(rPaths[0]));
    for (std::size_t k = 1; k < rPaths.size(); ++k) {
        area -= std::abs(Clipper2Lib::Area(rPaths[k]));
    }
    return area;
}

} 

void IgaSbmDomainClassificationUtility::ClassifyKnotSpans(
    const BrepSurfaceType& rBrepSurface,
    std::vector<double>& rSpansU,
    std::vector<double>& rSpansV,
    DenseMatrix<int>& rClassification)
{
    KRATOS_TRY

    rSpansU = rBrepSurface.KnotsU();
    rSpansV = rBrepSurface.KnotsV();

    KRATOS_ERROR_IF(rSpansU.size() < 2 || rSpansV.size() < 2)
        << "IgaSbmDomainClassificationUtility::ClassifyKnotSpans: degenerate knot vector." << std::endl;

    const SizeType n_spans_u = rSpansU.size() - 1;
    const SizeType n_spans_v = rSpansV.size() - 1;

    rClassification.resize(n_spans_u, n_spans_v, false);

    if (!rBrepSurface.IsTrimmed()) {
        for (IndexType i = 0; i < n_spans_u; ++i) {
            for (IndexType j = 0; j < n_spans_v; ++j) {
                rClassification(i, j) = static_cast<int>(KnotSpanClassification::Active);
            }
        }
        return;
    }

    const auto& r_outer_loops = rBrepSurface.GetOuterLoops();
    const auto& r_inner_loops = rBrepSurface.GetInnerLoops();

    Clipper2Lib::Paths64 outer_paths(r_outer_loops.size());
    for (IndexType i = 0; i < r_outer_loops.size(); ++i) {
        outer_paths[i] = TessellateLoop(r_outer_loops[i]);
    }
    Clipper2Lib::Paths64 inner_paths(r_inner_loops.size());
    for (IndexType i = 0; i < r_inner_loops.size(); ++i) {
        inner_paths[i] = TessellateLoop(r_inner_loops[i]);
    }

    for (IndexType i = 0; i < n_spans_u; ++i) {
        for (IndexType j = 0; j < n_spans_v; ++j) {
            Clipper2Lib::Rect64 rectangle(
                static_cast<cInt>(std::round(rSpansU[i] / ClipperScaleFactor)),
                static_cast<cInt>(std::round(rSpansV[j] / ClipperScaleFactor)),
                static_cast<cInt>(std::round(rSpansU[i + 1] / ClipperScaleFactor)),
                static_cast<cInt>(std::round(rSpansV[j + 1] / ClipperScaleFactor)));

            const double span_area = std::abs(Clipper2Lib::Area(rectangle.AsPath()));

            double clip_area = 0.0;
            if (!outer_paths.empty()) {
                const Clipper2Lib::Paths64 solution_outer = Clipper2Lib::RectClip(rectangle, outer_paths);
                clip_area = UnionArea(solution_outer);

                if (!inner_paths.empty() && clip_area > 0.0) {
                    const Clipper2Lib::Paths64 solution_inner = Clipper2Lib::RectClip(rectangle, inner_paths);
                    clip_area -= UnionArea(solution_inner);
                }
            }

            KnotSpanClassification cls;
            if (span_area <= 0.0 || clip_area / span_area < 1e-6) {
                cls = KnotSpanClassification::Outside;
            } else if (std::abs(1.0 - clip_area / span_area) < 1e-6) {
                cls = KnotSpanClassification::Active;
            } else {
                cls = KnotSpanClassification::Cut;
            }
            rClassification(i, j) = static_cast<int>(cls);
        }
    }

    KRATOS_CATCH("")
}

void IgaSbmDomainClassificationUtility::ValidateNonEmptyActiveDomain(
    const DenseMatrix<int>& rClassification,
    const std::string& rPatchName)
{
    KRATOS_TRY

    for (IndexType i = 0; i < rClassification.size1(); ++i) {
        for (IndexType j = 0; j < rClassification.size2(); ++j) {
            if (rClassification(i, j) == static_cast<int>(KnotSpanClassification::Active)) {
                return;
            }
        }
    }

    KRATOS_ERROR
        << "IgaSbmDomainClassificationUtility: patch \"" << rPatchName << "\" has ZERO fully-active knot "
        << "spans (Omega_tilde_h would be empty) -- domain-level SBM is structurally impossible on this "
        << "mesh. Refine (insert_nb_per_span_u/v) until this patch has at least one knot span that does "
        << "not touch the trim boundary." << std::endl;

    KRATOS_CATCH("")
}

std::vector<IgaSbmDomainClassificationUtility::SurrogateBoundarySegment>
IgaSbmDomainClassificationUtility::ComputeSurrogateBoundary(
    const std::vector<double>& rSpansU,
    const std::vector<double>& rSpansV,
    const DenseMatrix<int>& rClassification)
{
    KRATOS_TRY

    std::vector<SurrogateBoundarySegment> result;

    const SizeType n_spans_u = rClassification.size1();
    const SizeType n_spans_v = rClassification.size2();

    const auto is_active = [&](IndexType i, IndexType j) {
        return rClassification(i, j) == static_cast<int>(KnotSpanClassification::Active);
    };

    // U-direction interior edges: between span (i,j) and (i+1,j), at u=rSpansU[i+1].
    for (IndexType i = 0; i + 1 < n_spans_u; ++i) {
        for (IndexType j = 0; j < n_spans_v; ++j) {
            const bool active_left = is_active(i, j);
            const bool active_right = is_active(i + 1, j);
            if (active_left != active_right) {
                SurrogateBoundarySegment segment;
                segment.Start = ZeroVector(3);
                segment.End = ZeroVector(3);
                segment.Start[0] = rSpansU[i + 1];
                segment.Start[1] = rSpansV[j];
                segment.End[0] = rSpansU[i + 1];
                segment.End[1] = rSpansV[j + 1];
                segment.ActiveSpanIndexU = active_left ? i : (i + 1);
                segment.ActiveSpanIndexV = j;
                result.push_back(segment);
            }
        }
    }

    // V-direction interior edges: between span (i,j) and (i,j+1), at v=rSpansV[j+1].
    for (IndexType i = 0; i < n_spans_u; ++i) {
        for (IndexType j = 0; j + 1 < n_spans_v; ++j) {
            const bool active_bottom = is_active(i, j);
            const bool active_top = is_active(i, j + 1);
            if (active_bottom != active_top) {
                SurrogateBoundarySegment segment;
                segment.Start = ZeroVector(3);
                segment.End = ZeroVector(3);
                segment.Start[0] = rSpansU[i];
                segment.Start[1] = rSpansV[j + 1];
                segment.End[0] = rSpansU[i + 1];
                segment.End[1] = rSpansV[j + 1];
                segment.ActiveSpanIndexU = i;
                segment.ActiveSpanIndexV = active_bottom ? j : (j + 1);
                result.push_back(segment);
            }
        }
    }

    return result;

    KRATOS_CATCH("")
}

IgaSbmDomainClassificationUtility::ClassificationResult
IgaSbmDomainClassificationUtility::ClassifyKnotSpansFromGeometry(const Geometry<Node>& rGeometry)
{
    KRATOS_TRY

    const auto* p_brep_surface = dynamic_cast<const BrepSurfaceType*>(&rGeometry);
    KRATOS_ERROR_IF(p_brep_surface == nullptr)
        << "IgaSbmDomainClassificationUtility::ClassifyKnotSpansFromGeometry: the given geometry (id "
        << rGeometry.Id() << ") is not a BrepSurface -- pass one of ModelPart.Geometries's own top-level "
        << "entries (e.g. root_model_part.Geometries[brep_id]), not a domain element's own GetGeometry() "
        << "(that is a lightweight single-quadrature-point wrapper, not the Brep surface itself)." << std::endl;

    ClassificationResult result;
    DenseMatrix<int> classification_int;
    ClassifyKnotSpans(*p_brep_surface, result.SpansU, result.SpansV, classification_int);

    result.Classification.resize(classification_int.size1(), classification_int.size2());
    for (IndexType i = 0; i < classification_int.size1(); ++i) {
        for (IndexType j = 0; j < classification_int.size2(); ++j) {
            result.Classification(i, j) = static_cast<double>(classification_int(i, j));
        }
    }
    return result;

    KRATOS_CATCH("")
}

namespace {

DenseMatrix<int> ToIntClassification(const Matrix& rClassification)
{
    DenseMatrix<int> classification_int(rClassification.size1(), rClassification.size2());
    for (std::size_t i = 0; i < rClassification.size1(); ++i) {
        for (std::size_t j = 0; j < rClassification.size2(); ++j) {
            classification_int(i, j) = static_cast<int>(std::round(rClassification(i, j)));
        }
    }
    return classification_int;
}

} 

std::vector<IgaSbmDomainClassificationUtility::SurrogateBoundarySegment>
IgaSbmDomainClassificationUtility::ComputeSurrogateBoundaryFromClassification(
    const std::vector<double>& rSpansU,
    const std::vector<double>& rSpansV,
    const Matrix& rClassification)
{
    KRATOS_TRY
    return ComputeSurrogateBoundary(rSpansU, rSpansV, ToIntClassification(rClassification));
    KRATOS_CATCH("")
}

void IgaSbmDomainClassificationUtility::ValidateNonEmptyActiveDomainFromClassification(
    const Matrix& rClassification,
    const std::string& rPatchName)
{
    KRATOS_TRY
    ValidateNonEmptyActiveDomain(ToIntClassification(rClassification), rPatchName);
    KRATOS_CATCH("")
}

array_1d<double, 3> IgaSbmDomainClassificationUtility::GetElementParametricPosition(const Element& rElement)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();

    KRATOS_ERROR_IF(r_integration_points.size() != 1)
        << "IgaSbmDomainClassificationUtility::GetElementParametricPosition: element " << rElement.Id()
        << " has " << r_integration_points.size() << " integration points, expected exactly 1." << std::endl;

    array_1d<double, 3> position = ZeroVector(3);
    position[0] = r_integration_points[0][0];
    position[1] = r_integration_points[0][1];
    return position;

    KRATOS_CATCH("")
}

namespace {

struct BuiltSurrogateBoundaryGeometries
{
    std::vector<std::size_t> SegmentIndex;
    Geometry<Node>::GeometriesArrayType Geometries;
};

BuiltSurrogateBoundaryGeometries BuildSurrogateBoundaryGeometries(
    const std::vector<IgaSbmDomainClassificationUtility::SurrogateBoundarySegment>& rSegments,
    const Geometry<Node>& rGeometry,
    const std::size_t ShapeFunctionDerivativesOrder)
{
    using IndexType = std::size_t;
    using BrepSurfaceType = IgaSbmDomainClassificationUtility::BrepSurfaceType;

    const auto* p_brep_surface = dynamic_cast<const BrepSurfaceType*>(&rGeometry);
    KRATOS_ERROR_IF(p_brep_surface == nullptr)
        << "IgaSbmDomainClassificationUtility: the given geometry (id " << rGeometry.Id()
        << ") is not a BrepSurface." << std::endl;

    // NurbsCurveOnSurfaceGeometry::CreateQuadraturePointGeometriesSBM 
    KRATOS_ERROR_IF(ShapeFunctionDerivativesOrder == 0)
        << "IgaSbmDomainClassificationUtility: ShapeFunctionDerivativesOrder must be >= 1 "
        << "(1 = values only; the underlying SBM quadrature-point construction indexes it "
        << "as ShapeFunctionDerivativesOrder-1)." << std::endl;

    using NurbsSurfaceGeometryType = NurbsSurfaceGeometry<3, PointerVector<Node>>;
    using BrepCurveOnSurfaceType = BrepCurveOnSurface<PointerVector<Node>, true, PointerVector<Point>>;

    // BrepSurface's own OuterLoops/InnerLoops store TRIM curves
    auto* p_non_const_brep = const_cast<BrepSurfaceType*>(p_brep_surface);
    auto p_background_geometry = p_non_const_brep->pGetGeometryPart(Geometry<Node>::BACKGROUND_GEOMETRY_INDEX);
    auto p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceGeometryType>(p_background_geometry);
    KRATOS_ERROR_IF(p_nurbs_surface == nullptr)
        << "IgaSbmDomainClassificationUtility: could not obtain the underlying NurbsSurfaceGeometry "
        << "from the given BrepSurface (id " << rGeometry.Id() << ")." << std::endl;

    BuiltSurrogateBoundaryGeometries built;

    for (IndexType seg_idx = 0; seg_idx < rSegments.size(); ++seg_idx) {
        const auto& r_segment = rSegments[seg_idx];

        PointerVector<Point> control_points;
        control_points.push_back(Kratos::make_shared<Point>(r_segment.Start[0], r_segment.Start[1], 0.0));
        control_points.push_back(Kratos::make_shared<Point>(r_segment.End[0], r_segment.End[1], 0.0));
        Vector line_knot_vector(4);
        line_knot_vector[0] = 0.0; line_knot_vector[1] = 0.0;
        line_knot_vector[2] = 1.0; line_knot_vector[3] = 1.0;
        auto p_line_curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(
            control_points, 1, line_knot_vector);

        auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(p_nurbs_surface, p_line_curve);

        Geometry<Node>::GeometriesArrayType segment_quadrature_points;
        IntegrationInfo integration_info = p_brep_curve_on_surface->GetDefaultIntegrationInfo();
        Geometry<Node>::IntegrationPointsArrayType integration_points;
        p_brep_curve_on_surface->CreateIntegrationPoints(integration_points, integration_info);
        p_brep_curve_on_surface->CreateQuadraturePointGeometries(
            segment_quadrature_points, ShapeFunctionDerivativesOrder, integration_points, integration_info);

        for (IndexType k = 0; k < segment_quadrature_points.size(); ++k) {
            built.Geometries.push_back(segment_quadrature_points(k));
            built.SegmentIndex.push_back(seg_idx);
        }
    }

    return built;
}

// The 2D PARAMETRIC (du, dv) direction from a segment toward its OWN active
// span's center 
array_1d<double, 2> ComputeReferenceDirection(
    const IgaSbmDomainClassificationUtility::SurrogateBoundarySegment& rSegment,
    const std::vector<double>& rSpansU,
    const std::vector<double>& rSpansV)
{
    array_1d<double, 2> direction = ZeroVector(2);
    const bool is_u_edge = std::abs(rSegment.Start[0] - rSegment.End[0]) < 1e-12;
    if (is_u_edge) {
        const double u_edge = rSegment.Start[0];
        const double u_active_center = 0.5 * (rSpansU[rSegment.ActiveSpanIndexU] + rSpansU[rSegment.ActiveSpanIndexU + 1]);
        direction[0] = (u_active_center > u_edge) ? 1.0 : -1.0;
    } else {
        const double v_edge = rSegment.Start[1];
        const double v_active_center = 0.5 * (rSpansV[rSegment.ActiveSpanIndexV] + rSpansV[rSegment.ActiveSpanIndexV + 1]);
        direction[1] = (v_active_center > v_edge) ? 1.0 : -1.0;
    }
    return direction;
}

} 

std::vector<IgaSbmDomainClassificationUtility::SurrogateBoundaryQuadraturePointInfo>
IgaSbmDomainClassificationUtility::CreateSurrogateBoundaryQuadraturePoints(
    const std::vector<SurrogateBoundarySegment>& rSegments,
    const Geometry<Node>& rGeometry,
    const SizeType ShapeFunctionDerivativesOrder)
{
    KRATOS_TRY

    const auto built = BuildSurrogateBoundaryGeometries(rSegments, rGeometry, ShapeFunctionDerivativesOrder);

    std::vector<SurrogateBoundaryQuadraturePointInfo> result;

    for (IndexType k = 0; k < built.Geometries.size(); ++k) {
        const auto& r_qp_geometry = *built.Geometries(k);
        const SizeType n_gp = r_qp_geometry.IntegrationPointsNumber();
        const Matrix& r_N_all = r_qp_geometry.ShapeFunctionsValues();
        const SizeType n_cp = r_qp_geometry.PointsNumber();

        for (IndexType gp = 0; gp < n_gp; ++gp) {
            SurrogateBoundaryQuadraturePointInfo info;
            info.ParametricPosition = ZeroVector(3);
            info.ParametricPosition[0] = r_qp_geometry.IntegrationPoints()[gp][0];
            info.ParametricPosition[1] = r_qp_geometry.IntegrationPoints()[gp][1];

            array_1d<double, 3> physical_position = ZeroVector(3);
            for (IndexType i = 0; i < n_cp; ++i) {
                physical_position += r_N_all(gp, i) * r_qp_geometry[i].GetInitialPosition();
            }
            info.PhysicalPosition = physical_position;

            info.NodeIds.resize(n_cp);
            for (IndexType i = 0; i < n_cp; ++i) {
                info.NodeIds[i] = r_qp_geometry[i].Id();
            }
            info.ShapeFunctionValues = row(r_N_all, gp);

            result.push_back(info);
        }
    }

    return result;

    KRATOS_CATCH("")
}

IgaSbmDomainClassificationUtility::SizeType IgaSbmDomainClassificationUtility::CreateSurrogateBoundaryConditions(
    const std::vector<SurrogateBoundarySegment>& rSegments,
    const std::vector<double>& rSpansU,
    const std::vector<double>& rSpansV,
    const Geometry<Node>& rGeometry,
    const SizeType ShapeFunctionDerivativesOrder,
    const std::string& rConditionName,
    ModelPart& rTargetModelPart,
    Properties::Pointer pProperties,
    const IndexType StartId)
{
    KRATOS_TRY

    const auto built = BuildSurrogateBoundaryGeometries(rSegments, rGeometry, ShapeFunctionDerivativesOrder);

    IndexType id = StartId;
    for (IndexType k = 0; k < built.Geometries.size(); ++k) {
        auto p_condition = rTargetModelPart.CreateNewCondition(
            rConditionName, id, built.Geometries(k), pProperties);

        const auto& r_segment = rSegments[built.SegmentIndex[k]];
        const array_1d<double, 2> reference_direction = ComputeReferenceDirection(r_segment, rSpansU, rSpansV);
        Vector reference_direction_vec(2);
        reference_direction_vec[0] = reference_direction[0];
        reference_direction_vec[1] = reference_direction[1];
        p_condition->SetValue(SBM_SURROGATE_REFERENCE_DIRECTION, reference_direction_vec);

        ++id;
    }

    return built.Geometries.size();

    KRATOS_CATCH("")
}

}  // namespace Kratos.
