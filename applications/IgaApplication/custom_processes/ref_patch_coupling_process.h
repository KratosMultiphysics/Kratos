//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli

#pragma once

#include "processes/process.h"
#include "includes/model_part.h"
#include "geometries/geometry.h"
#include "integration/integration_info.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/brep_surface.h"
#include "geometries/nurbs_surface_geometry.h"
#include "includes/properties.h"

namespace Kratos {

// Specialized coupling builder for a base patch and a refinement patch.
// Creates coupling quadrature-point conditions along the refinement patch boundary
// against the base patch surface, segmenting by the union of knot-span breaks
// on both patches to ensure matching integration.
class RefPatchCouplingProcess : public Process {
public:
    KRATOS_CLASS_POINTER_DEFINITION(RefPatchCouplingProcess);

    using IndexType = std::size_t;
    using SizeType  = std::size_t;
    using NodeType  = Node;

    using GeometryType         = Geometry<NodeType>;
    using GeometryPointerType  = typename GeometryType::Pointer;
    using PropertiesPointerType= typename Properties::Pointer;

    using ParamCurveType         = NurbsCurveGeometry<2, PointerVector<Node>>;
    using BrepCurveOnSurfaceSbmType = BrepCurveOnSurface<PointerVector<Node>, true, PointerVector<Node>>;
    using NurbsSurfaceGeometryType        = NurbsSurfaceGeometry<3, PointerVector<Node>>;
    using NurbsSurfaceGeometryPointerType = typename NurbsSurfaceGeometryType::Pointer;
    using IPArray = GeometryType::IntegrationPointsArrayType;

    RefPatchCouplingProcess(
        ModelPart& rModelPart,
        int echo_level,
        double tolerance,
        std::string patch_prefix,
        int base_patch_index,
        int ref_patch_index,
        std::string coupling_condition_name);

    void Execute() override;

private:
    ModelPart&   mrModelPart;
    int          mEchoLevel{0};
    double       mTol{1e-12};
    std::string  mPatchPrefix{"Patch"};
    int          mBaseIndex{1};
    int          mRefIndex{2};
    std::string  mCouplingConditionName{"LaplacianCouplingCondition"};

    struct PatchBox {
        ModelPart* pPatch{nullptr};
        std::string name;
        double u0{}, u1{}, v0{}, v1{};
        int    nu{1}, nv{1}; // spans
        double du{0.0}, dv{0.0}; // span sizes
    };

    static bool IsClose(double a, double b, double tol) {
        return std::abs(a - b) <= tol * (1.0 + std::max(std::abs(a), std::abs(b)));
    }

    static std::vector<double> MakeBreaks(double a, double b, double s, int n_hint);
    static std::vector<double> UnionBreaks(const std::vector<double>& A, const std::vector<double>& B);

    static IndexType NextGeometryId(ModelPart& rModelPart);

    static void CreateAndAddBrepCurve(
        const NurbsSurfaceGeometryPointerType pSurfaceGeometry,
        const Point& rCoordsA,
        const Point& rCoordsB,
        IndexType& rLastGeometryId,
        ModelPart& rModelPart,
        bool MustBeFlipped = false);

    static void CreateConditionsFromBrepCurvesWithMirroredNeighbours(
        const std::vector<GeometryPointerType>& rBrepCurvesPrimary,
        const std::vector<GeometryPointerType>& rBrepCurvesMirror,
        ModelPart& rTargetSubModelPart,
        const std::string& rConditionName,
        SizeType ip_per_span,
        SizeType deriv_order,
        const Vector& rEffectiveKnotSpanSizes,
        int echo_level);

    static Vector ComputeEffectiveKnotSpans(ModelPart& rCouplingMP, ModelPart* pA, ModelPart* pB);

    // Build coupling segments along one edge and create mirrored coupling conditions
    void BuildCouplingSegmentsOnEdge(
        bool vertical_edge,
        double fixed,
        double tmin,
        double tmax,
        bool flipA,
        bool flipB,
        const PatchBox& A,
        const PatchBox& B,
        ModelPart& rCouplingSubModelPart,
        const NurbsSurfaceGeometryPointerType p_surf_A,
        const NurbsSurfaceGeometryPointerType p_surf_B,
        SizeType pA_u,
        SizeType pA_v,
        SizeType pB_u,
        SizeType pB_v);
};

} // namespace Kratos
