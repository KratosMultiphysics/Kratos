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

#include "iga_application_variables.h"
#include "custom_processes/ref_patch_coupling_process.h"

#include "includes/kratos_components.h"
#include "includes/variables.h"
#include "integration/integration_point_utilities.h"

namespace Kratos {

RefPatchCouplingProcess::RefPatchCouplingProcess(
    ModelPart& rModelPart,
    int echo_level,
    double tolerance,
    std::string patch_prefix,
    int base_patch_index,
    int ref_patch_index,
    std::string coupling_condition_name)
    : mrModelPart(rModelPart)
    , mEchoLevel(echo_level)
    , mTol(tolerance)
    , mPatchPrefix(std::move(patch_prefix))
    , mBaseIndex(base_patch_index)
    , mRefIndex(ref_patch_index)
    , mCouplingConditionName(std::move(coupling_condition_name))
{}

// --- helpers ---
std::vector<double> RefPatchCouplingProcess::MakeBreaks(double a, double b, double s, int n_hint)
{
    std::vector<double> r;
    if (s <= 0.0) { r.push_back(a); r.push_back(b); return r; }
    const int n = std::max(1, n_hint);
    r.reserve(n + 1);
    r.push_back(a);
    for (int i = 1; i < n; ++i) {
        const double x = a + i * s;
        if (x < b) r.push_back(x);
    }
    r.push_back(b);
    // unique with tolerance
    std::vector<double> u; u.reserve(r.size());
    for (double t : r) {
        if (u.empty() || !IsClose(u.back(), t, 1e-10)) u.push_back(t);
    }
    return u;
}

std::vector<double> RefPatchCouplingProcess::UnionBreaks(const std::vector<double>& A, const std::vector<double>& B)
{
    std::vector<double> out; out.reserve(A.size() + B.size());
    std::size_t i = 0, j = 0;
    while (i < A.size() || j < B.size()) {
        const double a = (i < A.size() ? A[i] : std::numeric_limits<double>::infinity());
        const double b = (j < B.size() ? B[j] : std::numeric_limits<double>::infinity());
        double pick;
        if (a < b && (j == B.size() || !IsClose(a, b, 1e-10))) { pick = a; ++i; }
        else if (j < B.size() && (i == A.size() || !IsClose(a, b, 1e-10))) { pick = b; ++j; }
        else { pick = a; ++i; ++j; } // equal within tol
        if (out.empty() || !IsClose(out.back(), pick, 1e-10)) out.push_back(pick);
    }
    return out;
}

RefPatchCouplingProcess::IndexType RefPatchCouplingProcess::NextGeometryId(ModelPart& rModelPart)
{
    IndexType max_id = 0;
    for (const auto& g : rModelPart.GetRootModelPart().Geometries()) {
        if (g.Id() > max_id) max_id = g.Id();
    }
    return max_id + 1;
}

void RefPatchCouplingProcess::CreateAndAddBrepCurve(
    const NurbsSurfaceGeometryPointerType pSurfaceGeometry,
    const Point& rCoordsA,
    const Point& rCoordsB,
    IndexType& rLastGeometryId,
    ModelPart& rModelPart,
    bool MustBeFlipped)
{
    Point first_point = rCoordsA;
    Point second_point = rCoordsB;

    Vector active_range_knot_vector = ZeroVector(2);
    if (first_point[0] == second_point[0]) {
        active_range_knot_vector[0] = first_point[1];
        active_range_knot_vector[1] = second_point[1];
    } else {
        active_range_knot_vector[0] = first_point[0];
        active_range_knot_vector[1] = second_point[0];
    }
    std::sort(active_range_knot_vector.begin(), active_range_knot_vector.end());

    if (MustBeFlipped) {
        std::swap(first_point, second_point);
    }

    NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

    // Build degree-1 param-space NURBS curve with two control points
    Vector knot_vector(2);
    knot_vector[0] = active_range_knot_vector[0];
    knot_vector[1] = active_range_knot_vector[1];
    const int p = 1;

    IndexType next_node_id = rModelPart.Nodes().empty() ? 1 : rModelPart.NumberOfNodes() + 1;
    auto pN1 = rModelPart.CreateNewNode(next_node_id++, first_point[0], first_point[1], 0.0);
    auto pN2 = rModelPart.CreateNewNode(next_node_id++, second_point[0], second_point[1], 0.0);

    PointerVector<Node> ctrl_pts;
    ctrl_pts.push_back(pN1);
    ctrl_pts.push_back(pN2);

    using NurbsCurveType = NurbsCurveGeometry<2, PointerVector<Node>>;
    auto p_curve = Kratos::make_shared<NurbsCurveType>(ctrl_pts, p, knot_vector);

    auto p_brep = Kratos::make_shared<BrepCurveOnSurfaceSbmType>(pSurfaceGeometry, p_curve, brep_active_range, true);
    p_brep->SetId(++rLastGeometryId);
    rModelPart.AddGeometry(p_brep);
}

Vector RefPatchCouplingProcess::ComputeEffectiveKnotSpans(ModelPart& rCouplingMP, ModelPart* pA, ModelPart* pB)
{
    Vector fallback;
    if (rCouplingMP.Has(KNOT_SPAN_SIZES)) {
        fallback = rCouplingMP.GetValue(KNOT_SPAN_SIZES);
    } else if (rCouplingMP.GetParentModelPart().Has(KNOT_SPAN_SIZES)) {
        fallback = rCouplingMP.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
    }
    if (fallback.size() == 0) fallback = ZeroVector(2);

    auto fetch = [](ModelPart* p)->Vector {
        if (p && p->Has(KNOT_SPAN_SIZES)) return p->GetValue(KNOT_SPAN_SIZES);
        return Vector();
    };
    Vector a = fetch(pA), b = fetch(pB);
    const SizeType n = std::min(a.size(), b.size());
    if (n > 0) {
        Vector e = ZeroVector(n);
        for (SizeType i = 0; i < n; ++i) e[i] = std::min(a[i], b[i]);
        return e;
    }
    if (a.size() > 0) return a;
    if (b.size() > 0) return b;
    return fallback;
}

void RefPatchCouplingProcess::CreateConditionsFromBrepCurvesWithMirroredNeighbours(
    const std::vector<GeometryPointerType>& rBrepCurvesPrimary,
    const std::vector<GeometryPointerType>& rBrepCurvesMirror,
    ModelPart& rTargetSubModelPart,
    const std::string& rConditionName,
    SizeType ip_per_span,
    SizeType deriv_order,
    const Vector& rEffectiveKnotSpanSizes,
    int echo_level)
{
    Condition const& rRefCond = KratosComponents<Condition>::Get(rConditionName);
    // Next condition id in root
    SizeType next_id = (rTargetSubModelPart.GetRootModelPart().Conditions().empty())
        ? 1
        : (rTargetSubModelPart.GetRootModelPart().Conditions().back().Id() + 1);
    GeometryType::GeometriesArrayType new_qp_geoms;
    std::vector<Condition::Pointer> new_conditions;
    new_conditions.reserve(rBrepCurvesPrimary.size() * std::max<std::size_t>(1, ip_per_span));

    auto gauss_nodes_weights = [](int n) -> const std::vector<std::array<double, 2>>& {
        const auto& gauss_table = IntegrationPointUtilities::s_gauss_legendre;
        KRATOS_ERROR_IF(n < 1 || n > static_cast<int>(gauss_table.size()))
            << "Gauss order " << n << " not implemented." << std::endl;
        return gauss_table[n - 1];
    };

    auto build_shared_ips = [&](GeometryPointerType pBrep) -> IPArray {
        double t0 = 0.0, t1 = 1.0;
        if (auto p_brep = std::dynamic_pointer_cast<BrepCurveOnSurfaceSbmType>(pBrep)) {
            const auto interval = p_brep->DomainInterval();
            t0 = interval.GetT0();
            t1 = interval.GetT1();
        } else {
            KRATOS_ERROR << "Unsupported curve geometry for shared-IP construction: "
                        << pBrep->Info() << std::endl;
        }
        IPArray ips;
        ips.reserve(ip_per_span);
        const double L = t1 - t0;
        const double aL = std::abs(L);
        KRATOS_ERROR_IF(aL < std::numeric_limits<double>::epsilon())
            << "Degenerate Brep curve interval detected." << std::endl;
        const auto& gw = gauss_nodes_weights(static_cast<int>(ip_per_span));
        for (SizeType i = 0; i < ip_per_span; ++i) {
            const double xi = gw[i][0];
            const double wi = gw[i][1];
            ips.emplace_back(IntegrationPoint<1>(t0 + L * xi, wi * aL));
        }
        return ips;
    };

    for (std::size_t s = 0; s < rBrepCurvesPrimary.size(); ++s) {
        auto pPrim = rBrepCurvesPrimary[s];
        auto pMir  = rBrepCurvesMirror[s];

        IPArray ips_shared = build_shared_ips(pPrim);

        GeometryType::GeometriesArrayType qpPrim, qpMir;
        IntegrationInfo infoPrim = pPrim->GetDefaultIntegrationInfo();
        infoPrim.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
        pPrim->CreateQuadraturePointGeometries(qpPrim, deriv_order, ips_shared, infoPrim);

        IntegrationInfo infoMir  = pMir->GetDefaultIntegrationInfo();
        infoMir.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
        pMir->CreateQuadraturePointGeometries(qpMir , deriv_order, ips_shared, infoMir);

        KRATOS_ERROR_IF(qpPrim.size() != qpMir.size())
            << "RefPatchCouplingProcess: mismatched GP counts on segment " << s
            << " (" << qpPrim.size() << " vs " << qpMir.size() << ")." << std::endl;

        for (std::size_t q = 0; q < qpPrim.size(); ++q) {
            auto pGeomPrim = qpPrim(q);
            auto pGeomMir  = qpMir(qpPrim.size() - q - 1); // opposite orientation

            // Compute normals like in PatchIntersectionProcess and enforce opposite directions
            const auto integration_method = pGeomPrim->GetDefaultIntegrationMethod();

            
            const array_1d<double, 3> normal_primary = pGeomPrim->Normal(0, integration_method);
            const array_1d<double, 3> normal_mirror  = pGeomMir ->Normal(0, integration_method);
            const double dot_nm = normal_primary[0]*normal_mirror[0]
                                + normal_primary[1]*normal_mirror[1]
                                + normal_primary[2]*normal_mirror[2];
            if (dot_nm > 0.0) {
                KRATOS_ERROR << "normals are not perpendicular: n*n = " << dot_nm << std::endl;
            }

            Condition::Pointer pCond = rRefCond.Create(next_id++, pGeomPrim, PropertiesPointerType());
            pCond->SetValue(KNOT_SPAN_SIZES, rEffectiveKnotSpanSizes);

            std::vector<Geometry<Node>::Pointer> neigh(1);
            neigh[0] = pGeomMir;
            pCond->SetValue(NEIGHBOUR_GEOMETRIES, neigh);

            // Ensure nodes exist in the target submodelpart
            for (SizeType i = 0; i < pGeomPrim->size(); ++i) {
                rTargetSubModelPart.Nodes().push_back(pGeomPrim->pGetPoint(i));
            }
            new_conditions.push_back(pCond);
        }
    }

    if (!new_conditions.empty()) {
        rTargetSubModelPart.AddConditions(new_conditions.begin(), new_conditions.end());
        KRATOS_INFO_IF("RefPatchCouplingProcess", echo_level > 0)
            << "Added " << new_conditions.size() << " coupling conditions along ref-patch boundary." << std::endl;
    }
}

void RefPatchCouplingProcess::Execute()
{
    const std::string base_name = mPatchPrefix + std::to_string(mBaseIndex);
    const std::string ref_name  = mPatchPrefix + std::to_string(mRefIndex);

    KRATOS_ERROR_IF_NOT(mrModelPart.HasSubModelPart(base_name))
        << "RefPatchCouplingProcess: base patch '" << base_name << "' not found." << std::endl;
    KRATOS_ERROR_IF_NOT(mrModelPart.HasSubModelPart(ref_name))
        << "RefPatchCouplingProcess: ref patch '" << ref_name << "' not found." << std::endl;

    ModelPart& r_base = mrModelPart.GetSubModelPart(base_name);
    ModelPart& r_ref  = mrModelPart.GetSubModelPart(ref_name);

    const Matrix& base_c = r_base.GetValue(PATCH_PARAMETER_SPACE_CORNERS);
    const Matrix& ref_c  = r_ref .GetValue(PATCH_PARAMETER_SPACE_CORNERS);

    PatchBox A, B;
    A.pPatch = &r_base; 
    A.name = base_name;
    B.pPatch = &r_ref ; 
    B.name = ref_name;
    A.u0 = base_c(0,0); A.u1 = base_c(0,1); A.v0 = base_c(1,0); A.v1 = base_c(1,1);
    B.u0 = ref_c (0,0); B.u1 = ref_c (0,1); B.v0 = ref_c (1,0); B.v1 = ref_c (1,1);

    // default
    A.nu = A.nv = 1; A.du = (A.u1 - A.u0); A.dv = (A.v1 - A.v0);
    B.nu = B.nv = 1; B.du = (B.u1 - B.u0); B.dv = (B.v1 - B.v0);

    if (r_base.Has(KNOT_SPAN_SIZES)) {
        const Vector& sz = r_base.GetValue(KNOT_SPAN_SIZES);
        if (sz.size() >= 1 && sz[0] > 0.0) { A.du = sz[0]; A.nu = std::max(1, (int)std::round((A.u1-A.u0)/A.du)); }
        if (sz.size() >= 2 && sz[1] > 0.0) { A.dv = sz[1]; A.nv = std::max(1, (int)std::round((A.v1-A.v0)/A.dv)); }
    }
    if (r_ref.Has(KNOT_SPAN_SIZES)) {
        const Vector& sz = r_ref.GetValue(KNOT_SPAN_SIZES);
        if (sz.size() >= 1 && sz[0] > 0.0) { B.du = sz[0]; B.nu = std::max(1, (int)std::round((B.u1-B.u0)/B.du)); }
        if (sz.size() >= 2 && sz[1] > 0.0) { B.dv = sz[1]; B.nv = std::max(1, (int)std::round((B.v1-B.v0)/B.dv)); }
    }

    ModelPart& coupling = mrModelPart.HasSubModelPart("MultipatchCouplingConditions")
        ? mrModelPart.GetSubModelPart("MultipatchCouplingConditions")
        : mrModelPart.CreateSubModelPart("MultipatchCouplingConditions");

    // Surfaces
    auto p_brep_surf_A = A.pPatch->pGetGeometry(1);
    GeometryType::Pointer p_base_A = p_brep_surf_A->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
    auto p_surf_A = std::dynamic_pointer_cast<NurbsSurfaceGeometryType>(p_base_A);
    KRATOS_ERROR_IF(!p_surf_A) << A.name << ": BACKGROUND_GEOMETRY not a NurbsSurfaceGeometry." << std::endl;

    auto p_brep_surf_B = B.pPatch->pGetGeometry(1);
    GeometryType::Pointer p_base_B = p_brep_surf_B->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
    auto p_surf_B = std::dynamic_pointer_cast<NurbsSurfaceGeometryType>(p_base_B);
    KRATOS_ERROR_IF(!p_surf_B) << B.name << ": BACKGROUND_GEOMETRY not a NurbsSurfaceGeometry." << std::endl;

    // ip per span: use max degree in the swept direction per edge
    const SizeType pA_u = p_surf_A->PolynomialDegree(0);
    const SizeType pA_v = p_surf_A->PolynomialDegree(1);
    const SizeType pB_u = p_surf_B->PolynomialDegree(0);
    const SizeType pB_v = p_surf_B->PolynomialDegree(1);

    // Intersect ref boundary with base box (clip to base domain)
    const double u0 = std::max(B.u0, A.u0);
    const double u1 = std::min(B.u1, A.u1);
    const double v0 = std::max(B.v0, A.v0);
    const double v1 = std::min(B.v1, A.v1);
    if (!(u1 > u0 + mTol && v1 > v0 + mTol)) {
        KRATOS_INFO("RefPatchCouplingProcess") << "Ref patch does not overlap base domain. No coupling created." << std::endl;
        return;
    }

    // Left edge (u = u0)
    BuildCouplingSegmentsOnEdge(true /*vertical*/, B.u0, v0, v1, /*flipA*/false, /*flipB*/true,
        A, B, coupling, p_surf_A, p_surf_B, pA_u, pA_v, pB_u, pB_v);
    // Right edge (u = u1)
    BuildCouplingSegmentsOnEdge(true /*vertical*/, B.u1, v0, v1, /*flipA*/true,  /*flipB*/false,
        A, B, coupling, p_surf_A, p_surf_B, pA_u, pA_v, pB_u, pB_v);
    // Bottom edge (v = v0)
    BuildCouplingSegmentsOnEdge(false /*horizontal*/, B.v0, u0, u1, /*flipA*/true, /*flipB*/false,
        A, B, coupling, p_surf_A, p_surf_B, pA_u, pA_v, pB_u, pB_v);
    // Top edge (v = v1)
    BuildCouplingSegmentsOnEdge(false /*horizontal*/, B.v1, u0, u1, /*flipA*/false,  /*flipB*/true,
        A, B, coupling, p_surf_A, p_surf_B, pA_u, pA_v, pB_u, pB_v);
    // Left edge (u = u0)

}

void RefPatchCouplingProcess::BuildCouplingSegmentsOnEdge(
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
    SizeType pB_v)
{
    // Build breaks for the sweep direction
    std::vector<double> Av = vertical_edge ? MakeBreaks(A.v0, A.v1, A.dv, A.nv)
                                           : MakeBreaks(A.u0, A.u1, A.du, A.nu);
    std::vector<double> Bv = vertical_edge ? MakeBreaks(B.v0, B.v1, B.dv, B.nv)
                                           : MakeBreaks(B.u0, B.u1, B.du, B.nu);
    // Clip to intersection with [tmin, tmax]
    auto clip = [&](std::vector<double>& X){
        std::vector<double> Y; Y.reserve(X.size());
        for (double t : X) {
            if (t >= tmin - 1e-10 && t <= tmax + 1e-10)
                Y.push_back(std::clamp(t, tmin, tmax));
        }
        X.swap(Y);
    };
    clip(Av); clip(Bv);
    const auto breaks = UnionBreaks(Av, Bv);

    // Degree along sweep
    const SizeType p_sweep = vertical_edge ? std::max(pA_v, pB_v) : std::max(pA_u, pB_u);
    const SizeType ip_per_span = 2 * p_sweep + 1;
    const SizeType deriv_order = 2;
    IndexType next_id = NextGeometryId(mrModelPart);

    std::vector<GeometryPointerType> breps_A, breps_B;
    for (std::size_t k = 0; k + 1 < breaks.size(); ++k) {
        const double a = breaks[k];
        const double b = breaks[k+1];
        if (b <= a + mTol) continue;
        Point A0(vertical_edge ? fixed : a, vertical_edge ? a : fixed, 0.0);
        Point A1(vertical_edge ? fixed : b, vertical_edge ? b : fixed, 0.0);
        Point B0 = A0; Point B1 = A1;

        CreateAndAddBrepCurve(p_surf_A, A0, A1, next_id, mrModelPart, flipA);
        breps_A.push_back(mrModelPart.pGetGeometry(next_id));
        CreateAndAddBrepCurve(p_surf_B, B0, B1, next_id, mrModelPart, flipB);
        breps_B.push_back(mrModelPart.pGetGeometry(next_id));
    }

    const Vector eff_spans = ComputeEffectiveKnotSpans(rCouplingSubModelPart, A.pPatch, B.pPatch);
    CreateConditionsFromBrepCurvesWithMirroredNeighbours(
        breps_A, breps_B, rCouplingSubModelPart, mCouplingConditionName,
        ip_per_span, deriv_order, eff_spans, mEchoLevel);
}

} // namespace Kratos
