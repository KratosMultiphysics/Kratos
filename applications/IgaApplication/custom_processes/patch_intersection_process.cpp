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

#include "iga_application_variables.h" // KNOT_SPAN_SIZES, PATCH_PARAMETER_SPACE_CORNERS
#include "custom_processes/patch_intersection_process.h"

#include "includes/kratos_components.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <iostream>

namespace Kratos {

// --------- small local helpers ---------
namespace {

struct PatchBox {
    ModelPart* pPatch{nullptr};
    std::string name;
    double u0{}, u1{}, v0{}, v1{};
    int    nu{1}, nv{1};   // estimated number of spans
    double du{0.0}, dv{0.0}; // span sizes
};

inline bool IsClose(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) <= tol * (1.0 + std::max(std::abs(a), std::abs(b)));
}

// Build [a,b] with step s (inclusive ends). Uses n_hint to avoid drift.
inline std::vector<double> MakeBreaks(double a, double b, double s, int n_hint) {
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

inline std::vector<double> UnionBreaks(const std::vector<double>& A, const std::vector<double>& B) {
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

} // namespace


PatchIntersectionProcess::PatchIntersectionProcess(
    ModelPart& rModelPart,
    int echo_level,
    double tolerance,
    std::string patch_prefix,
    std::string internal_sub,
    std::string external_sub)
    : mrModelPart(rModelPart)
    , mEchoLevel(echo_level)
    , mTol(tolerance)
    , mPatchPrefix(std::move(patch_prefix))
    , mInternalSubName(std::move(internal_sub))
    , mExternalSubName(std::move(external_sub))
{
    if (mEchoLevel > 2) {
        std::cout << "[PatchIntersectionProcess] ctor\n"
                  << "  patch_prefix     : " << mPatchPrefix     << "\n"
                  << "  internal_sub     : " << mInternalSubName << "\n"
                  << "  external_sub     : " << mExternalSubName << "\n"
                  << "  tol              : " << mTol             << "\n"
                  << "  echo_level       : " << mEchoLevel       << std::endl;
    }
}

void PatchIntersectionProcess::Execute() {
    mEchoLevel = 1;
    if (mEchoLevel > 2) std::cout << "[PatchIntersectionProcess] Execute()" << std::endl;
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



    // 2) For every adjacent pair, create body-fitted conditions
    const std::string cond_name_bodyfitted = "CutSbmLaplacianInterfaceCondition"; // <- change if needed // TODO:

    for (std::size_t i = 0; i < patches.size(); ++i) {
        const auto& A = patches[i];
        KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2) << "patch A: " << A.name << std::endl;
        for (std::size_t j = i + 1; j < patches.size(); ++j) {
            const auto& B = patches[j];
            KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2) << "patch B: " << B.name << std::endl;


            // Vertical adjacency: A.u1 == B.u0 or A.u0 == B.u1, with V-overlap
            const bool A_right_B_left = std::abs(A.u1 - B.u0) <= mTol;
            const bool A_left_B_right = std::abs(A.u0 - B.u1) <= mTol;
            const bool v_overlap      = (std::max(A.v0, B.v0) < std::min(A.v1, B.v1) - mTol);

            if ((A_right_B_left || A_left_B_right) && v_overlap) {
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
                const double u_fixed_A = A_right_B_left ? A.u1 : A.u0;
                const double u_fixed_B = A_right_B_left ? B.u0 : B.u1;

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
                    const Point A0(u_fixed_A, v0, 0.0);
                    const Point A1(u_fixed_A, v1, 0.0);
                    const Point B0(u_fixed_B, v0, 0.0);
                    const Point B1(u_fixed_B, v1, 0.0);

                    // Create Brep on A
                    CreateAndAddBrepCurve(p_surf_A, A0, A1, next_id, mrModelPart);
                    breps_A.push_back(mrModelPart.pGetGeometry(next_id));

                    // Create Brep on B
                    CreateAndAddBrepCurve(p_surf_B, B0, B1, next_id, mrModelPart);
                    breps_B.push_back(mrModelPart.pGetGeometry(next_id));
                }

                // Create interface conditions on A and attach mirrored QP geometry from B
                CreateConditionsFromBrepCurvesWithMirroredNeighbours(
                    breps_A, breps_B,
                    coupling_sub_model_part,      // target sub-model part
                    cond_name_bodyfitted,         // condition name
                    ip_per_span,                  // identical IP layout on both sides
                    deriv_order                   // shape func derivative order
                );
            }





            // Horizontal adjacency: A.v1 == B.v0 or A.v0 == B.v1, with U-overlap
            const bool A_top_B_bot = std::abs(A.v1 - B.v0) <= mTol;
            const bool A_bot_B_top = std::abs(A.v0 - B.v1) <= mTol;
            const bool u_overlap   = (std::max(A.u0, B.u0) < std::min(A.u1, B.u1) - mTol);

            if ((A_top_B_bot || A_bot_B_top) && u_overlap) {
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
                const double v_fixed_A = A_top_B_bot ? A.v1 : A.v0;
                const double v_fixed_B = A_top_B_bot ? B.v0 : B.v1; // opposite side of B

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
                    2 * static_cast<SizeType>(std::max(pA, pB)) + 1;  // odd, like elsewhere
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
                    const Point A0(u0, v_fixed_A, 0.0);
                    const Point A1(u1, v_fixed_A, 0.0);
                    const Point B0(u0, v_fixed_B, 0.0);
                    const Point B1(u1, v_fixed_B, 0.0);

                    // Create Brep on A
                    CreateAndAddBrepCurve(p_surf_A, A0, A1, next_id, mrModelPart);
                    breps_A.push_back(mrModelPart.pGetGeometry(next_id));

                    // Create Brep on B
                    CreateAndAddBrepCurve(p_surf_B, B0, B1, next_id, mrModelPart);
                    breps_B.push_back(mrModelPart.pGetGeometry(next_id));
                }

                // Create interface conditions on A and attach mirrored QP geometry from B
                CreateConditionsFromBrepCurvesWithMirroredNeighbours(
                    breps_A, breps_B,
                    coupling_sub_model_part,      // target sub-model part
                    cond_name_bodyfitted,         // condition name
                    ip_per_span,                  // identical IP layout on both sides
                    deriv_order                   // shape func derivative order
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
    ModelPart& rModelPart)
{
    // define active range: the first point and second point are in order
    Vector active_range_knot_vector = ZeroVector(2);
    // Compute the knot vector needed
    if (rCoordsA[0] == rCoordsB[0]) {
        // the brep curve is vertical
        active_range_knot_vector[0] = rCoordsA[1];
        active_range_knot_vector[1] = rCoordsB[1];
    } else {
        // the brep curve is horizontal
        active_range_knot_vector[0] = rCoordsA[0];
        active_range_knot_vector[1] = rCoordsB[0];
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
    auto pN1 = rModelPart.CreateNewNode(next_node_id++, rCoordsA[0], rCoordsA[1], 0.0);
    auto pN2 = rModelPart.CreateNewNode(next_node_id++, rCoordsB[0], rCoordsB[1], 0.0);

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


void PatchIntersectionProcess::CreateConditionsFromBrepCurves(
    const std::vector<GeometryPointerType>& rBrepCurves,
    ModelPart& rTargetSubModelPart,
    const std::string& rConditionName,
    SizeType ip_per_span)
{
    KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 3) << "CreateConditionsFromBrepCurves: number of Brep curves = " << rBrepCurves.size() << std::endl;
    if (rBrepCurves.empty()) return;

    KRATOS_ERROR_IF(!KratosComponents<Condition>::Has(rConditionName))
        << rConditionName << " not registered." << std::endl;

    const Condition& rRefCond = KratosComponents<Condition>::Get(rConditionName);

    // Pick KNOT_SPAN_SIZES from the parent (same as IgaModelerSbm)
    const Vector& knot_span_sizes = rTargetSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
    
    KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 3) << "  Knot span sizes: " << knot_span_sizes << std::endl;
    
    // Next condition id in the root model part
    SizeType next_id = (rTargetSubModelPart.GetRootModelPart().Conditions().empty())
        ? 1
        : (rTargetSubModelPart.GetRootModelPart().Conditions().back().Id() + 1);

    ModelPart::ConditionsContainerType new_conditions;

    for (const auto& pCurve : rBrepCurves) {
        
        KRATOS_ERROR_IF(pCurve.get() == nullptr) << "Null Brep curve geometry." << std::endl;

        // Generate quadrature-point geometries on the BREP curve
        GeometryType::GeometriesArrayType qp_geoms;

        IntegrationInfo info = pCurve->GetDefaultIntegrationInfo();
        info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);

        // Derivative order = 2 (same default as IgaModelerSbm)
        pCurve->CreateQuadraturePointGeometries(qp_geoms, 2, info);

        // Create one condition per QP geometry
        for (auto it = qp_geoms.ptr_begin(); it != qp_geoms.ptr_end(); ++it) {
            Condition::Pointer pNew = rRefCond.Create(next_id, (*it), PropertiesPointerType());

            // Propagate mesh sizes
            pNew->SetValue(KNOT_SPAN_SIZES, knot_span_sizes);

            // Ensure the nodes of the QP geometry are present in the target submodel part
            for (SizeType i = 0; i < (*it)->size(); ++i) {
                rTargetSubModelPart.Nodes().push_back((*it)->pGetPoint(i));
            }

            new_conditions.push_back(pNew);
            ++next_id;
        }
    }

    if (!new_conditions.empty()) {
        rTargetSubModelPart.AddConditions(new_conditions.begin(), new_conditions.end());
        KRATOS_INFO_IF("PatchIntersectionProcess", mEchoLevel > 2) << "Added " << new_conditions.size() << " conditions. \n" << std::endl;
    }
}

void PatchIntersectionProcess::CreateConditionsFromBrepCurvesWithMirroredNeighbours(
    const std::vector<GeometryPointerType>& rBrepCurvesPrimary,
    const std::vector<GeometryPointerType>& rBrepCurvesMirror,
    ModelPart& rTargetSubModelPart,
    const std::string& rConditionName,
    SizeType ip_per_span,     // <- interpreted as points per SEGMENT here
    SizeType deriv_order)
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

    const Vector& knot_span_sizes = rTargetSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    ModelPart::ConditionsContainerType new_conditions;

    // --- small helper: Gauss nodes/weights on [-1,1] for n=1..10
    auto gauss_nodes_weights = [](int n, std::vector<double>& xi, std::vector<double>& w)
    {
        KRATOS_ERROR_IF(n < 1 || n > 10) << "Gauss order " << n << " not implemented (1…10)." << std::endl;
        static const double X[10][10] = {
            { 0 },
            { -0.5773502691896257,  0.5773502691896257 },
            { 0.0, -0.7745966692414834,  0.7745966692414834 },
            { -0.3399810435848563,  0.3399810435848563, -0.8611363115940526,  0.8611363115940526 },
            { 0.0, -0.5384693101056831,  0.5384693101056831, -0.9061798459386640,  0.9061798459386640 },
            { -0.2386191860831969,  0.2386191860831969, -0.6612093864662645,  0.6612093864662645, -0.9324695142031521,  0.9324695142031521 },
            { 0.0, -0.4058451513773972,  0.4058451513773972, -0.7415311855993945,  0.7415311855993945, -0.9491079123427585,  0.9491079123427585 },
            { -0.1834346424956498,  0.1834346424956498, -0.5255324099163290,  0.5255324099163290, -0.7966664774136267,  0.7966664774136267, -0.9602898564975363,  0.9602898564975363 },
            { 0.0, -0.3242534234038089,  0.3242534234038089, -0.6133714327005904,  0.6133714327005904, -0.8360311073266358,  0.8360311073266358, -0.9681602395076261,  0.9681602395076261 },
            { -0.1488743389816312,  0.1488743389816312, -0.4333953941292472,  0.4333953941292472, -0.6794095682990244,  0.6794095682990244, -0.8650633666889845,  0.8650633666889845, -0.9739065285171717,  0.9739065285171717 }
        };
        static const double W[10][10] = {
            { 2.0 },
            { 1.0, 1.0 },
            { 0.8888888888888888, 0.5555555555555556, 0.5555555555555556 },
            { 0.6521451548625461, 0.6521451548625461, 0.3478548451374539, 0.3478548451374539 },
            { 0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891 },
            { 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.3607615730481386, 0.1713244923791704, 0.1713244923791704 },
            { 0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766, 0.1294849661688697, 0.1294849661688697 },
            { 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.3137066458778873, 0.2223810344533745, 0.2223810344533745, 0.1012285362903763, 0.1012285362903763 },
            { 0.3302393550012598, 0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744 },
            { 0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2692667193099963, 0.2190863625159820, 0.2190863625159820, 0.1494513491505806, 0.1494513491505806, 0.0666713443086881, 0.0666713443086881 }
        };
        xi.assign(X[n-1], X[n-1] + n);
        w .assign(W[n-1], W[n-1] + n);
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
        std::vector<double> xi, w;
        gauss_nodes_weights(static_cast<int>(ip_per_span), xi, w);

        IPArray ips;
        ips.reserve(ip_per_span);
        const double half = 0.5 * (t1 - t0);
        const double mid  = 0.5 * (t1 + t0);
        for (SizeType i = 0; i < ip_per_span; ++i) {
            const double t = mid + half * xi[i];
            ips.emplace_back(IntegrationPoint<1>(t, w[i]));
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
        pPrim->CreateQuadraturePointGeometries(qpPrim, deriv_order, ips_shared, infoPrim);

        IntegrationInfo infoMir  = pMir->GetDefaultIntegrationInfo();
        infoMir.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
        pMir->CreateQuadraturePointGeometries(qpMir , deriv_order, ips_shared, infoMir);

        // Now counts MUST match (we generated from the same ips_shared)
        KRATOS_ERROR_IF(qpPrim.size() != qpMir.size())
            << "CreateConditionsFromBrepCurvesWithMirroredNeighbours: mismatched GP counts on segment "
            << s << " even with shared ips (" << qpPrim.size() << " vs " << qpMir.size() << ")." << std::endl;

        // Create primary-side conditions and attach the mirrored QP-geometry as neighbour
        for (std::size_t q = 0; q < qpPrim.size(); ++q) {
            const auto& pGeomPrim = qpPrim(q);
            const auto& pGeomMir  = qpMir(q);

            Condition::Pointer pCond = rRefCond.Create(next_id++, pGeomPrim, PropertiesPointerType());
            pCond->SetValue(KNOT_SPAN_SIZES, knot_span_sizes);

            std::vector<Geometry<Node>::Pointer> neigh(1);
            neigh[0] = pGeomMir; // cross-link with exact matching location
            pCond->SetValue(NEIGHBOUR_GEOMETRIES, neigh);

            // Ensure nodes of the QP-geometry are present
            for (SizeType i = 0; i < pGeomPrim->size(); ++i) {
                rTargetSubModelPart.Nodes().push_back(pGeomPrim->pGetPoint(i));
                // TODO: maybe also add the nodes of the mirror geometry?
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
