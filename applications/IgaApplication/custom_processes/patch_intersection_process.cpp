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
#include "custom_processes/patch_intersection_process.h"
#include "includes/variables.h"
#include "integration/integration_point_utilities.h"

namespace Kratos {

// Kratos-style private static helpers
bool PatchIntersectionProcess::IsClose(double a, double b, double tol)
{
    return std::abs(a - b) <= tol * (1.0 + std::max(std::abs(a), std::abs(b)));
}

std::vector<double> PatchIntersectionProcess::MakeBreaks(double a, double b, double s, int nHint)
{
    std::vector<double> r;
    if (s <= 0.0) { r.push_back(a); r.push_back(b); return r; }
    const int n = std::max(1, nHint);
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

std::vector<double> PatchIntersectionProcess::UnionBreaks(const std::vector<double>& rA, const std::vector<double>& rB)
{
    std::vector<double> out; out.reserve(rA.size() + rB.size());
    std::size_t i = 0, j = 0;
    while (i < rA.size() || j < rB.size()) {
        const double a = (i < rA.size() ? rA[i] : std::numeric_limits<double>::infinity());
        const double b = (j < rB.size() ? rB[j] : std::numeric_limits<double>::infinity());
        double pick;
        if (a < b && (j == rB.size() || !IsClose(a, b, 1e-10))) { pick = a; ++i; }
        else if (j < rB.size() && (i == rA.size() || !IsClose(a, b, 1e-10))) { pick = b; ++j; }
        else { pick = a; ++i; ++j; } // equal within tol
        if (out.empty() || !IsClose(out.back(), pick, 1e-10)) out.push_back(pick);
    }
    return out;
}


PatchIntersectionProcess::PatchIntersectionProcess(
    ModelPart& rModelPart,
    int EchoLevel,
    double Tolerance,
    std::string PatchPrefix,
    std::string CouplingConditionName)
    : mrModelPart(rModelPart)
    , mEchoLevel(EchoLevel)
    , mTol(Tolerance)
    , mPatchPrefix(std::move(PatchPrefix))
    , mCouplingConditionName(std::move(CouplingConditionName))
{
    KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2)
        << "[PatchIntersectionProcess]\n"
        << "  patch_prefix     : " << mPatchPrefix << "\n"
        << "  tol              : " << mTol << "\n"
        << "  echo_level       : " << mEchoLevel << "\n"
        << "  coupling_cond    : " << mCouplingConditionName << std::endl;
}

void PatchIntersectionProcess::Execute() {
    KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2)
        << "[PatchIntersectionProcess] Execute()" << std::endl;
    ComputeIntersections();
}

// --------- core logic ---------

void PatchIntersectionProcess::ComputeIntersections()
{
    // 1) Collect all patches + their param boxes and span info.
    std::vector<PatchBox> patches;
    patches.reserve(mrModelPart.NumberOfSubModelParts());

    ModelPart& coupling_sub_model_part = mrModelPart.HasSubModelPart("MultipatchCouplingConditions")
        ? mrModelPart.GetSubModelPart("MultipatchCouplingConditions")
        : mrModelPart.CreateSubModelPart("MultipatchCouplingConditions");

    auto compute_effective_knot_spans = [&](ModelPart* pPatchA, ModelPart* pPatchB) {
        Vector fallback;
        if (coupling_sub_model_part.Has(KNOT_SPAN_SIZES)) {
            fallback = coupling_sub_model_part.GetValue(KNOT_SPAN_SIZES);
        } else if (coupling_sub_model_part.GetParentModelPart().Has(KNOT_SPAN_SIZES)) {
            fallback = coupling_sub_model_part.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
        }
        if (fallback.size() == 0) {
            fallback = ZeroVector(2);
        }

        auto fetch_spans = [&](ModelPart* pPatch) -> Vector {
            if (pPatch && pPatch->Has(KNOT_SPAN_SIZES)) {
                return pPatch->GetValue(KNOT_SPAN_SIZES);
            }
            return Vector();
        };

        Vector span_a = fetch_spans(pPatchA);
        Vector span_b = fetch_spans(pPatchB);

        const SizeType common_size = std::min(span_a.size(), span_b.size());
        if (common_size > 0) {
            Vector effective = ZeroVector(common_size);
            for (SizeType idx = 0; idx < common_size; ++idx) {
                effective[idx] = std::min(span_a[idx], span_b[idx]);
            }
            return effective;
        }

        if (span_a.size() > 0) {
            return span_a;
        }
        if (span_b.size() > 0) {
            return span_b;
        }

        return fallback;
    };

    for (auto& patch : mrModelPart.SubModelParts()) {
        const std::string& name = patch.Name();
        KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2) << "Processing patch: " << name << std::endl;
        if (name.find(mPatchPrefix) != 0) continue;

        const Matrix& parameters_space_corners = patch.GetValue(PATCH_PARAMETER_SPACE_CORNERS);
        KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2) << "  Parameters space corners: " << parameters_space_corners << std::endl;

        // Create patch box
        PatchBox pb;
        pb.pPatch = &patch;
        pb.name = name;
        pb.u0 = parameters_space_corners(0,0); pb.u1 = parameters_space_corners(0,1);
        pb.v0 = parameters_space_corners(1,0); pb.v1 = parameters_space_corners(1,1);

        pb.du = pb.u1 - pb.u0;
        pb.dv = pb.v1 - pb.v0;
        pb.nu = 1; pb.nv = 1;

        if (patch.Has(KNOT_SPAN_SIZES)) {
            const Vector& sz = patch.GetValue(KNOT_SPAN_SIZES);
            if (sz.size() >= 1 && sz[0] > 0.0) { pb.du = sz[0]; pb.nu = std::max(1, static_cast<int>(std::round((pb.u1 - pb.u0)/pb.du))); }
            if (sz.size() >= 2 && sz[1] > 0.0) { pb.dv = sz[1]; pb.nv = std::max(1, static_cast<int>(std::round((pb.v1 - pb.v0)/pb.dv))); }
        }

        patches.push_back(pb);
    }
    KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2) << "Number of patches: " << patches.size() << std::endl;



    // 2) For every adjacent pair, create coupling conditions
    const std::string coupling_condition_name = mCouplingConditionName;

    for (std::size_t i = 0; i < patches.size(); ++i) {
        const auto& A = patches[i];
        KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2) << "patch A: " << A.name << std::endl;
        for (std::size_t j = i + 1; j < patches.size(); ++j) {
            const auto& B = patches[j];
            KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2) << "patch B: " << B.name << std::endl;


            // Vertical adjacency: A.u1 == B.u0 or A.u0 == B.u1, with V-overlap
            const bool a_right_b_left = std::abs(A.u1 - B.u0) <= mTol;
            const bool a_left_b_right = std::abs(A.u0 - B.u1) <= mTol;
            const bool v_overlap      = (std::max(A.v0, B.v0) < std::min(A.v1, B.v1) - mTol);

            if ((a_right_b_left || a_left_b_right) && v_overlap) {
                KRATOS_INFO_IF("\n PatchIntersectionProcess", mEchoLevel > 2)
                    << "  Vertical overlap " << v_overlap << std::endl;

                // Intersected V-range
                const double vmin = std::max(A.v0, B.v0);
                const double vmax = std::min(A.v1, B.v1);

                // Build breaks along v on both sides and take the union
                auto Av = MakeBreaks(A.v0, A.v1, A.dv, A.nv);
                auto Bv = MakeBreaks(B.v0, B.v1, B.dv, B.nv);

                // Clip to [vmin, vmax]
                auto clip_v = [&](std::vector<double>& X){
                    std::vector<double> Y; Y.reserve(X.size());
                    for (double t : X) {
                        if (t >= vmin - 1e-10 && t <= vmax + 1e-10)
                            Y.push_back(std::clamp(t, vmin, vmax));
                    }
                    X.swap(Y);
                };
                clip_v(Av); clip_v(Bv);

                const auto union_breaks_v = UnionBreaks(Av, Bv);
                KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2)
                    << "    Union of breaks (V): " << union_breaks_v << std::endl;

                // Fixed u on both patches along the common vertical edge
                const double u_fixed_a = a_right_b_left ? A.u1 : A.u0;
                const double u_fixed_b = a_right_b_left ? B.u0 : B.u1;

                // -- surfaces on each patch (geom id 1 is the surface)
                auto p_brep_surf_A = A.pPatch->pGetGeometry(1);
                GeometryType::Pointer p_base_A =
                    p_brep_surf_A->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
                auto p_surf_A = std::dynamic_pointer_cast<NurbsSurfaceGeometryType>(p_base_A);
                KRATOS_ERROR_IF(!p_surf_A)
                    << A.name << ": BACKGROUND_GEOMETRY is not a NurbsSurfaceGeometry." << std::endl;

                auto p_brep_surf_B = B.pPatch->pGetGeometry(1);
                GeometryType::Pointer p_base_B =
                    p_brep_surf_B->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
                auto p_surf_B = std::dynamic_pointer_cast<NurbsSurfaceGeometryType>(p_base_B);
                KRATOS_ERROR_IF(!p_surf_B)
                    << B.name << ": BACKGROUND_GEOMETRY is not a NurbsSurfaceGeometry." << std::endl;

                // Choose IPs per span based on the higher polynomial order in v
                const SizeType pA_v = p_surf_A->PolynomialDegree(1);
                const SizeType pB_v = p_surf_B->PolynomialDegree(1);
                const SizeType ip_per_span =
                    2 * static_cast<SizeType>(std::max(pA_v, pB_v)) + 1;  // odd, like elsewhere
                const SizeType deriv_order = 2;

                // fresh geometry id from the root model part
                auto next_geom_id = [](ModelPart& mp)->IndexType {
                    IndexType max_id = 0;
                    for (const auto& g : mp.GetRootModelPart().Geometries()) {
                        if (g.Id() > max_id) max_id = g.Id();
                    }
                    return max_id + 1;
                };
                IndexType next_id = next_geom_id(mrModelPart);

                // Build paired Breps (A-side, B-side) for each V-segment
                std::vector<GeometryPointerType> breps_A;
                std::vector<GeometryPointerType> breps_B;
                breps_A.reserve(union_breaks_v.size());
                breps_B.reserve(union_breaks_v.size());

                for (std::size_t k = 0; k + 1 < union_breaks_v.size(); ++k) {
                    const double v0 = union_breaks_v[k];
                    const double v1 = union_breaks_v[k + 1];
                    if (v1 <= v0 + mTol) continue; // skip degenerate

                    // Endpoints on A and B (same v-interval, different fixed u)
                    Point A0(u_fixed_a, v0, 0.0);
                    Point A1(u_fixed_a, v1, 0.0);
                    Point B0(u_fixed_b, v0, 0.0);
                    Point B1(u_fixed_b, v1, 0.0);

                    const bool flip_b = a_right_b_left;
                    const bool flip_a = a_left_b_right;

                    // Create Brep on A
                    CreateAndAddBrepCurve(p_surf_A, A0, A1, next_id, mrModelPart, flip_a);
                    breps_A.push_back(mrModelPart.pGetGeometry(next_id));

                    // Create Brep on B
                    CreateAndAddBrepCurve(p_surf_B, B0, B1, next_id, mrModelPart, flip_b);
                    breps_B.push_back(mrModelPart.pGetGeometry(next_id));
                }

                const Vector effective_span_sizes = compute_effective_knot_spans(A.pPatch, B.pPatch);

                // Create interface conditions on A and attach mirrored QP geometry from B
                CreateConditionsFromBrepCurvesWithMirroredNeighbours(
                    breps_A, breps_B,
                    coupling_sub_model_part,      // target sub-model part
                    coupling_condition_name,      // condition name
                    ip_per_span,                  // identical IP layout on both sides
                    deriv_order,                  // shape func derivative order
                    effective_span_sizes
                );
            }


            // Horizontal adjacency: A.v1 == B.v0 or A.v0 == B.v1, with U-overlap
            const bool a_top_b_bot = std::abs(A.v1 - B.v0) <= mTol;
            const bool a_bot_b_top = std::abs(A.v0 - B.v1) <= mTol;
            const bool u_overlap   = (std::max(A.u0, B.u0) < std::min(A.u1, B.u1) - mTol);

            if ((a_top_b_bot || a_bot_b_top) && u_overlap) {
                KRATOS_INFO_IF("\n PatchIntersectionProcess", mEchoLevel > 2)
                    << "  Horizontal overlap " << u_overlap << std::endl;

                // Intersected U-range
                const double umin = std::max(A.u0, B.u0);
                const double umax = std::min(A.u1, B.u1);

                // Build breaks along u on both sides and take the union
                auto Au = MakeBreaks(A.u0, A.u1, A.du, A.nu);
                auto Bu = MakeBreaks(B.u0, B.u1, B.du, B.nu);

                // Clip to [umin, umax]
                auto clip_u = [&](std::vector<double>& X){
                    std::vector<double> Y; Y.reserve(X.size());
                    for (double t : X) {
                        if (t >= umin - 1e-10 && t <= umax + 1e-10)
                            Y.push_back(std::clamp(t, umin, umax));
                    }
                    X.swap(Y);
                };
                clip_u(Au); clip_u(Bu);

                const auto union_breaks_u = UnionBreaks(Au, Bu);
                KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2)
                    << "    Union of breaks (U): " << union_breaks_u << std::endl;

                // Fixed v on both patches along the common horizontal edge
                const double v_fixed_a = a_top_b_bot ? A.v1 : A.v0;
                const double v_fixed_b = a_top_b_bot ? B.v0 : B.v1; // opposite side of B

                // -- surfaces on each patch (geom id 1 is the surface)
                auto p_brep_surf_A = A.pPatch->pGetGeometry(1);
                GeometryType::Pointer p_base_A =
                    p_brep_surf_A->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
                auto p_surf_A = std::dynamic_pointer_cast<NurbsSurfaceGeometryType>(p_base_A);
                KRATOS_ERROR_IF(!p_surf_A)
                    << A.name << ": BACKGROUND_GEOMETRY is not a NurbsSurfaceGeometry." << std::endl;

                auto p_brep_surf_B = B.pPatch->pGetGeometry(1);
                GeometryType::Pointer p_base_B =
                    p_brep_surf_B->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
                auto p_surf_B = std::dynamic_pointer_cast<NurbsSurfaceGeometryType>(p_base_B);
                KRATOS_ERROR_IF(!p_surf_B)
                    << B.name << ": BACKGROUND_GEOMETRY is not a NurbsSurfaceGeometry." << std::endl;

                // Choose IPs per span based on the higher polynomial order
                const SizeType pA = p_surf_A->PolynomialDegree(0);
                const SizeType pB = p_surf_B->PolynomialDegree(0);
                const SizeType ip_per_span =
                    2 * static_cast<SizeType>(std::max(pA, pB)) + 1;
                const SizeType deriv_order = 2;

                // fresh geometry id from the root model part
                auto next_geom_id = [](ModelPart& mp)->IndexType {
                    IndexType max_id = 0;
                    for (const auto& g : mp.GetRootModelPart().Geometries()) {
                        if (g.Id() > max_id) max_id = g.Id();
                    }
                    return max_id + 1;
                };
                IndexType next_id = next_geom_id(mrModelPart);

                // Build paired Breps (A-side, B-side) for each U-segment
                std::vector<GeometryPointerType> breps_A;
                std::vector<GeometryPointerType> breps_B;
                breps_A.reserve(union_breaks_u.size());
                breps_B.reserve(union_breaks_u.size());

                for (std::size_t k = 0; k + 1 < union_breaks_u.size(); ++k) {
                    const double u0 = union_breaks_u[k];
                    const double u1 = union_breaks_u[k + 1];
                    if (u1 <= u0 + mTol) continue; // skip degenerate

                    // Endpoints on A and B (same u-interval, different fixed v)
                    Point A0(u0, v_fixed_a, 0.0);
                    Point A1(u1, v_fixed_a, 0.0);
                    Point B0(u0, v_fixed_b, 0.0);
                    Point B1(u1, v_fixed_b, 0.0);

                    const bool flip_b = a_bot_b_top;
                    const bool flip_a = a_top_b_bot;

                    // Create Brep on A
                    CreateAndAddBrepCurve(p_surf_A, A0, A1, next_id, mrModelPart, flip_a);
                    breps_A.push_back(mrModelPart.pGetGeometry(next_id));

                    // Create Brep on B
                    CreateAndAddBrepCurve(p_surf_B, B0, B1, next_id, mrModelPart, flip_b);
                    breps_B.push_back(mrModelPart.pGetGeometry(next_id));
                }

                // Create interface conditions on A and attach mirrored QP geometry from B
                const Vector effective_span_sizes = compute_effective_knot_spans(A.pPatch, B.pPatch);

                CreateConditionsFromBrepCurvesWithMirroredNeighbours(
                    breps_A, breps_B,
                    coupling_sub_model_part,      // target sub-model part
                    coupling_condition_name,      // condition name
                    ip_per_span,                  // identical IP layout on both sides
                    deriv_order,                  // shape func derivative order
                    effective_span_sizes
                );
            } // end horizontal adjacency

        }
    }
}


void PatchIntersectionProcess::CreateAndAddBrepCurve(
    const NurbsSurfaceGeometryPointerType pSurfaceGeometry, // NurbsSurfaceGeometry<3, PointerVector<Node>>::Pointer
    const Point& rCoordsA,                                   // (uA, vA, 0)
    const Point& rCoordsB,                                   // (uB, vB, 0)
    IndexType& rLastGeometryId,                              // consumed & incremented
    ModelPart& rModelPart,
    const bool MustBeFlipped)
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

    // --- build a degree-1 (linear) NURBS curve in UV with Node control points
    Vector knot_vector(2);
    knot_vector[0] = active_range_knot_vector[0];
    knot_vector[1] = active_range_knot_vector[1];
    const int p = 1;

    // Allocate node IDs safely from the root model part
    IndexType next_node_id = rModelPart.Nodes().empty()
        ? 1
        : rModelPart.NumberOfNodes() + 1;

    // Create two parametric nodes (u, v, 0)
    auto pN1 = rModelPart.CreateNewNode(next_node_id++, first_point[0], first_point[1], 0.0);
    auto pN2 = rModelPart.CreateNewNode(next_node_id++, second_point[0], second_point[1], 0.0);

    PointerVector<Node> ctrl_pts;
    ctrl_pts.push_back(pN1);
    ctrl_pts.push_back(pN2);

    using NurbsCurveType = NurbsCurveGeometry<2, PointerVector<Node>>;
    auto p_curve = Kratos::make_shared<NurbsCurveType>(ctrl_pts, p, knot_vector);

    // --- wrap as BrepCurveOnSurface
    auto p_brep = Kratos::make_shared<BrepCurveOnSurfaceType>(pSurfaceGeometry, p_curve, brep_active_range, true);
    
    // assign an id and register it in the model part
    KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2) << "    Assigning ID " << rLastGeometryId << " to BrepCurveOnSurface" << std::endl;
    p_brep->SetId(++rLastGeometryId);
    rModelPart.AddGeometry(p_brep);
}

void PatchIntersectionProcess::CreateConditionsFromBrepCurvesWithMirroredNeighbours(
    const std::vector<GeometryPointerType>& rBrepCurvesPrimary,
    const std::vector<GeometryPointerType>& rBrepCurvesMirror,
    ModelPart& rTargetSubModelPart,
    const std::string& rConditionName,
    SizeType IpPerSpan,     // <- interpreted as points per SEGMENT here
    SizeType DerivOrder,
    const Vector& rEffectiveKnotSpanSizes)
{
    KRATOS_ERROR_IF(rBrepCurvesPrimary.size() != rBrepCurvesMirror.size())
        << "CreateConditionsFromBrepCurvesWithMirroredNeighbours: curve lists have different sizes ("
        << rBrepCurvesPrimary.size() << " vs " << rBrepCurvesMirror.size() << ")." << std::endl;

    KRATOS_ERROR_IF(!KratosComponents<Condition>::Has(rConditionName))
        << rConditionName << " not registered." << std::endl;

    const Condition& rRefCond = KratosComponents<Condition>::Get(rConditionName);

    // Next condition id in the root modelpart
    SizeType next_id = (rTargetSubModelPart.GetRootModelPart().Conditions().empty())
        ? 1
        : (rTargetSubModelPart.GetRootModelPart().Conditions().back().Id() + 1);

    ModelPart::ConditionsContainerType new_conditions;

    // --- small helper: Gauss nodes/weights on [0,1]
    auto gauss_nodes_weights = [](int n) -> const std::vector<std::array<double, 2>>&
    {
        const auto& gauss_table = IntegrationPointUtilities::s_gauss_legendre;
        KRATOS_ERROR_IF(n < 1 || n > static_cast<int>(gauss_table.size()))
            << "Gauss order " << n << " not implemented." << std::endl;
        return gauss_table[n - 1];
    };

    auto build_shared_ips = [&](GeometryPointerType pBrep) -> IPArray
    {
        double t0 = 0.0, t1 = 1.0;

        if (auto p_brep = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(pBrep)) {
            const auto interval = p_brep->DomainInterval();
            t0 = interval.GetT0();
            t1 = interval.GetT1();
        } else {
            KRATOS_ERROR << "Unsupported curve geometry for shared-IP construction: "
                        << pBrep->Info() << " (type has no DomainInterval())." << std::endl;
        }

        // Build Gauss nodes on [t0,t1]
        IPArray ips;
        ips.reserve(IpPerSpan);
        const double length = t1 - t0;
        const double abs_length = std::abs(length);
        KRATOS_ERROR_IF(abs_length < std::numeric_limits<double>::epsilon())
            << "Degenerate Brep curve interval detected (length ~ 0)." << std::endl;
        const auto& gauss_entries = gauss_nodes_weights(static_cast<int>(IpPerSpan));

        for (SizeType i = 0; i < IpPerSpan; ++i) {
            const double xi = gauss_entries[i][0];
            const double wi = gauss_entries[i][1];
            const double t = t0 + length * xi;
            ips.emplace_back(IntegrationPoint<1>(t, wi * abs_length));
        }
        return ips;
    };


    for (std::size_t s = 0; s < rBrepCurvesPrimary.size(); ++s) {
        auto pPrim = rBrepCurvesPrimary[s];
        auto pMir  = rBrepCurvesMirror[s];

        // Build ONE shared Gauss grid (in curve parameter) for THIS segment…
        IPArray ips_shared = build_shared_ips(pPrim);

        // …and reuse it on BOTH sides to create the QP geometries
        GeometryType::GeometriesArrayType qpPrim, qpMir;

        IntegrationInfo infoPrim = pPrim->GetDefaultIntegrationInfo();
        infoPrim.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
        pPrim->CreateQuadraturePointGeometries(qpPrim, DerivOrder, ips_shared, infoPrim);

        IntegrationInfo infoMir  = pMir->GetDefaultIntegrationInfo();
        infoMir.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
        pMir->CreateQuadraturePointGeometries(qpMir , DerivOrder, ips_shared, infoMir);

        // Now counts MUST match (we generated from the same ips_shared)
        KRATOS_ERROR_IF(qpPrim.size() != qpMir.size())
            << "CreateConditionsFromBrepCurvesWithMirroredNeighbours: mismatched GP counts on segment "
            << s << " even with shared ips (" << qpPrim.size() << " vs " << qpMir.size() << ")." << std::endl;

        // Create primary-side conditions and attach the mirrored QP-geometry as neighbour
        for (std::size_t q = 0; q < qpPrim.size(); ++q) {
            auto pGeomPrim = qpPrim(q);
            // the brep are always created in opposite directions
            auto pGeomMir  = qpMir(qpPrim.size()-q-1); 

            const auto integration_method = pGeomPrim->GetDefaultIntegrationMethod();
            const array_1d<double, 3> normal_primary = pGeomPrim->Normal(0, integration_method);
            array_1d<double, 3> normal_mirror = pGeomMir->Normal(0, integration_method);
            
            // Be sure they have opposite orientation (flip mirrored if needed)
            KRATOS_ERROR_IF(MathUtils<double>::Dot(normal_primary, normal_mirror) > 0.0)
                << "normals are not perpendicular: n*n = "
                << MathUtils<double>::Dot(normal_primary, normal_mirror) << std::endl;

            Condition::Pointer pCond = rRefCond.Create(next_id++, pGeomPrim, PropertiesPointerType());
            pCond->SetValue(KNOT_SPAN_SIZES, rEffectiveKnotSpanSizes);

            std::vector<Geometry<Node>::Pointer> neigh(1);
            neigh[0] = pGeomMir; // cross-link with exact matching location
            pCond->SetValue(NEIGHBOUR_GEOMETRIES, neigh);

            // Ensure nodes of the QP-geometry are present
            for (SizeType i = 0; i < pGeomPrim->size(); ++i) {
                rTargetSubModelPart.Nodes().push_back(pGeomPrim->pGetPoint(i));
            }

            new_conditions.push_back(pCond);
        }
    }

    if (!new_conditions.empty()) {
        rTargetSubModelPart.AddConditions(new_conditions.begin(), new_conditions.end());
        KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 1)
            << "[MirroredNeighbour] Added " << new_conditions.size()
            << " interface conditions with mirrored neighbour geometries.\n";
    }
}


} // namespace Kratos
