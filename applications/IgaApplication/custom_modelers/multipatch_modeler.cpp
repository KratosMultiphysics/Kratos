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

// System includes

// External includes

// Project includes
#include "custom_modelers/multipatch_modeler.h"
#include "iga_application_variables.h"
#include "custom_processes/patch_intersection_process.h"
#include "custom_processes/ref_patch_coupling_process.h"

namespace Kratos
{

namespace
{
using RectangleType = MultipatchModeler::RectangleType;

constexpr std::size_t U_MIN = 0;
constexpr std::size_t U_MAX = 1;
constexpr std::size_t V_MIN = 2;
constexpr std::size_t V_MAX = 3;

RectangleType ClipRectangle(const RectangleType& rRectangle, const RectangleType& rBounds)
{
    RectangleType result;
    result[U_MIN] = std::max(rRectangle[U_MIN], rBounds[U_MIN]);
    result[U_MAX] = std::min(rRectangle[U_MAX], rBounds[U_MAX]);
    result[V_MIN] = std::max(rRectangle[V_MIN], rBounds[V_MIN]);
    result[V_MAX] = std::min(rRectangle[V_MAX], rBounds[V_MAX]);
    return result;
}

bool HasPositiveArea(const RectangleType& rRectangle)
{
    return (rRectangle[U_MAX] - rRectangle[U_MIN] > 0.0) &&
           (rRectangle[V_MAX] - rRectangle[V_MIN] > 0.0);
}

std::string RectangleToString(const RectangleType& rRectangle)
{
    std::ostringstream buffer;
    buffer << "[" << rRectangle[U_MIN] << ", " << rRectangle[U_MAX]
           << "] x [" << rRectangle[V_MIN] << ", " << rRectangle[V_MAX] << "]";
    return buffer.str();
}

} // namespace

ModelPart& MultipatchModeler::CreateOrResetModelPart_(const std::string& rName) const
{
    ModelPart& r_mp = mpModel->HasModelPart(rName)
        ? mpModel->GetModelPart(rName)
        : mpModel->CreateModelPart(rName);

    r_mp.Nodes().clear();
    r_mp.Elements().clear();
    r_mp.Conditions().clear();
    r_mp.Geometries().clear();
    if (!r_mp.HasProperties(0)) {
        r_mp.CreateNewProperties(0);
    }
    return r_mp;
}

void MultipatchModeler::AppendRectangleSkinLoop_(ModelPart& rModelPart, const RectType& rect) const
{
    auto& r_root = rModelPart.GetRootModelPart();
    const int next_node_id = static_cast<int>(r_root.NumberOfNodes()) + 1;
    const int next_cond_id = static_cast<int>(r_root.NumberOfConditions()) + 1;

    const double u0 = rect[U_MIN];
    const double u1 = rect[U_MAX];
    const double v0 = rect[V_MIN];
    const double v1 = rect[V_MAX];

    auto p1 = rModelPart.CreateNewNode(next_node_id + 0, u0, v0, 0.0);
    auto p2 = rModelPart.CreateNewNode(next_node_id + 1, u1, v0, 0.0);
    auto p3 = rModelPart.CreateNewNode(next_node_id + 2, u1, v1, 0.0);
    auto p4 = rModelPart.CreateNewNode(next_node_id + 3, u0, v1, 0.0);

    std::vector<ModelPart::IndexType> conn(2);
    auto p_props = rModelPart.pGetProperties(0);

    conn[0] = p1->Id(); conn[1] = p2->Id();
    rModelPart.CreateNewCondition("LineCondition2D2N", next_cond_id + 0, conn, p_props);

    conn[0] = p2->Id(); conn[1] = p3->Id();
    rModelPart.CreateNewCondition("LineCondition2D2N", next_cond_id + 1, conn, p_props);

    conn[0] = p3->Id(); conn[1] = p4->Id();
    rModelPart.CreateNewCondition("LineCondition2D2N", next_cond_id + 2, conn, p_props);

    conn[0] = p4->Id(); conn[1] = p1->Id();
    rModelPart.CreateNewCondition("LineCondition2D2N", next_cond_id + 3, conn, p_props);
}

ModelPart& MultipatchModeler::CreateSkinCouplingModelPartForRefinements_(const std::string& rSkinModelPartName) const
{
    ModelPart& r_skin = CreateOrResetModelPart_(rSkinModelPartName);
    // Create one rectangular loop per refinement region
    for (const auto& r_reg : mRefinementRegions) {
        AppendRectangleSkinLoop_(r_skin, r_reg.Rectangle);
    }
    return r_skin;
}

MultipatchModeler::MultipatchModeler(Model& rModel, const Parameters ModelParameters)
    : Modeler(rModel, ModelParameters)
    , mpModel(&rModel)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

Modeler::Pointer MultipatchModeler::Create(Model& rModel, const Parameters ModelParameters) const
{
    return Kratos::make_shared<MultipatchModeler>(rModel, ModelParameters);
}

const Parameters MultipatchModeler::GetDefaultParameters() const
{
    return Parameters(R"({
        "model_part_name" : "",
        "base_domain" : {
            "lower_point_uvw" : [0.0, 0.0, 0.0],
            "upper_point_uvw" : [1.0, 1.0, 0.0]
        },
        "refinement_regions" : [],
        "child_patch_prefix" : "Patch",
        "geometry_parameters" : {},
        "analysis_parameters" : {},
        "coupling_conditions_name": "",
        "echo_level" : 0
    })");
}

const std::vector<MultipatchModeler::RectangleType>& MultipatchModeler::GetSubdomains() const
{
    return mSubdomains;
}


void MultipatchModeler::SetupModelPart()
{
    auto& r_parameters = mParameters;

    Parameters geometry_base;
    Parameters analysis_base;
    std::vector<std::string> analysis_target_submodel_parts;
    std::unordered_map<std::string, std::vector<std::string>> condition_name_to_model_parts;

    // Centralized validation, logging, subdivision, and global corners write.
    ModelPart& r_model_part = InitializeSetup_(
        r_parameters,
        geometry_base,
        analysis_base,
        analysis_target_submodel_parts,
        condition_name_to_model_parts);

    // If no geometry parameters were provided, the original code returns early.
    const bool has_geometry_params =
        r_parameters.Has("geometry_parameters") &&
        r_parameters["geometry_parameters"].Has("model_part_name");
    if (!has_geometry_params) {
        return;
    }

    const std::string prefix = r_parameters["child_patch_prefix"].GetString();
    ModelPart& r_parent_model_part = r_model_part;

    // Process only the base domain patch with SBM and the just-created inner skin

    // Create skin_coupling_model_part initial loops:
    // - inner: union of refinement rectangles
    // - outer: base rectangle
    const std::string inner_initial_name = "skin_coupling_model_part_initial";
    const std::string skin_name          = "skin_coupling_model_part";

    // Build inner (from refinement regions)
    ModelPart& r_inner_initial = CreateSkinCouplingModelPartForRefinements_(inner_initial_name);
    KRATOS_INFO_IF("MultipatchModeler", mEchoLevel > 1)
        << "Built inner initial skin '" << inner_initial_name << "' with "
        << r_inner_initial.NumberOfConditions() << " conditions." << std::endl;

    // Note: do NOT set coupling skin names on geometry_base.
    // Base patch will use these names directly in its local patch_geometry.


    /// Process base domain Patch

    const std::string patch_suffix = prefix + std::to_string(1); // Patch naming 
    ModelPart& r_patch_model_part = r_parent_model_part.HasSubModelPart(patch_suffix)
        ? r_parent_model_part.GetSubModelPart(patch_suffix)
        : r_parent_model_part.CreateSubModelPart(patch_suffix);
    const std::string patch_full_name = r_patch_model_part.FullName();

    if (mEchoLevel > 2) {
        KRATOS_INFO("MultipatchModeler")
            << "-- Created " << patch_full_name << std::endl;
    }

    const bool is_base_patch = true;

    // Clone base geometry parameters, scope to this patch
    Parameters patch_geometry = geometry_base.Clone();
    patch_geometry["model_part_name"].SetString(patch_full_name);

    // For base patch, enforce use of the coupling skin we just created
    if (patch_geometry.Has("skin_model_part_inner_initial_name")) patch_geometry.RemoveValue("skin_model_part_inner_initial_name");
    if (patch_geometry.Has("skin_model_part_name")) patch_geometry.RemoveValue("skin_model_part_name");
    patch_geometry.AddEmptyValue("skin_model_part_inner_initial_name").SetString(inner_initial_name);
    patch_geometry.AddEmptyValue("skin_model_part_name").SetString(skin_name);

    // KRATOS_WATCH(patch_geometry)

    // Skin model parts are prepared in SetupModelPart; read skin name if present
    std::string sbm_skin_name_for_base;
    if (is_base_patch && patch_geometry.Has("skin_model_part_name")) {
        sbm_skin_name_for_base = patch_geometry["skin_model_part_name"].GetString();
        KRATOS_INFO_IF("MultipatchModeler", mEchoLevel > 1)
            << "SBM skin for base patch: '" << sbm_skin_name_for_base << "'" << std::endl;
    }

    // Decide SBM/body-fitted policy:
    // - Base (external) patch: use SBM if skin is provided and fully inside
    const bool has_skin_data = patch_geometry.Has("skin_model_part_name") &&
                               (patch_geometry.Has("skin_model_part_outer_initial_name") ||
                                patch_geometry.Has("skin_model_part_inner_initial_name"));

    bool skin_fully_inside = true;
    // For the moment only for inner boundaries
    if (patch_geometry.Has("skin_model_part_inner_initial_name")) {
        const std::string inner_skin_name = patch_geometry["skin_model_part_inner_initial_name"].GetString();
        // For base domain containment check use the base rectangle
        skin_fully_inside = IsSkinFullyInsidePatch_(inner_skin_name, mBaseRect, geometry_base);
        KRATOS_INFO_IF("MultipatchModeler", mEchoLevel > 2)
            << "  Skin containment (inner='" << inner_skin_name << "'): "
            << (skin_fully_inside ? "inside" : "outside") << std::endl;
    }


    // enable SBM on base patch only
    bool use_sbm_for_this_patch = is_base_patch && has_skin_data && skin_fully_inside;
    if (mEchoLevel > 2) {
        KRATOS_INFO("MultipatchModeler")
            << "Creating/Updating patch '" << patch_full_name << "' from subdomain "
            << RectangleToString(mBaseRect) << std::endl;
        KRATOS_INFO("MultipatchModeler")
            << "  Patch classification ....: " << (use_sbm_for_this_patch ? "SBM" : "body-fitted") << std::endl;
        if (!has_skin_data && is_base_patch) {
            KRATOS_INFO("MultipatchModeler")
                << "  Note: missing skin_model_part definitions -> surrogate workflow disabled for base patch." << std::endl;
        }
    }

    // Snapshot geometry IDs before generation (to detect new ones)
    std::vector<ModelPart::IndexType> patch_geometry_ids_before;
    patch_geometry_ids_before.reserve(r_patch_model_part.NumberOfGeometries());
    for (const auto& r_geom : r_patch_model_part.Geometries()) {
        patch_geometry_ids_before.push_back(r_geom.Id());
    }

    // patch_geometry is already prepared above

    // Assign UVW bounds for the patch
    Vector patch_lower_uvw = patch_geometry["lower_point_uvw"].GetVector();
    Vector patch_upper_uvw = patch_geometry["upper_point_uvw"].GetVector();
    patch_lower_uvw[0] = mBaseRect[U_MIN]; patch_upper_uvw[0] = mBaseRect[U_MAX];
    patch_lower_uvw[1] = mBaseRect[V_MIN]; patch_upper_uvw[1] = mBaseRect[V_MAX];
    patch_geometry["lower_point_uvw"].SetVector(patch_lower_uvw);
    patch_geometry["upper_point_uvw"].SetVector(patch_upper_uvw);

    // Map UVW fractions into XYZ using base box and base UVW
    const auto& base_lower_xyz = geometry_base["lower_point_xyz"].GetVector();
    const auto& base_upper_xyz = geometry_base["upper_point_xyz"].GetVector();
    const auto& base_lower_uvw = geometry_base["lower_point_uvw"].GetVector();
    const auto& base_upper_uvw = geometry_base["upper_point_uvw"].GetVector();
    const auto& base_spans     = geometry_base["number_of_knot_spans"].GetVector();

    const double base_u_length = base_upper_uvw[0] - base_lower_uvw[0];
    const double base_v_length = base_upper_uvw[1] - base_lower_uvw[1];

    Vector patch_lower_xyz = patch_geometry["lower_point_xyz"].GetVector();
    Vector patch_upper_xyz = patch_geometry["upper_point_xyz"].GetVector();

    const double u_fraction_min = (mBaseRect[U_MIN] - base_lower_uvw[0]) / base_u_length;
    const double u_fraction_max = (mBaseRect[U_MAX] - base_lower_uvw[0]) / base_u_length;
    const double v_fraction_min = (mBaseRect[V_MIN] - base_lower_uvw[1]) / base_v_length;
    const double v_fraction_max = (mBaseRect[V_MAX] - base_lower_uvw[1]) / base_v_length;

    // Keep Z as-is, interpolate X and Y
    patch_lower_xyz[0] = base_lower_xyz[0] + u_fraction_min * (base_upper_xyz[0] - base_lower_xyz[0]);
    patch_upper_xyz[0] = base_lower_xyz[0] + u_fraction_max * (base_upper_xyz[0] - base_lower_xyz[0]);
    patch_lower_xyz[1] = base_lower_xyz[1] + v_fraction_min * (base_upper_xyz[1] - base_lower_xyz[1]);
    patch_upper_xyz[1] = base_lower_xyz[1] + v_fraction_max * (base_upper_xyz[1] - base_lower_xyz[1]);
    patch_lower_xyz[2] = base_lower_xyz[2];
    patch_upper_xyz[2] = base_upper_xyz[2];

    patch_geometry["lower_point_xyz"].SetVector(patch_lower_xyz);
    patch_geometry["upper_point_xyz"].SetVector(patch_upper_xyz);

    // If SBM is not applicable to this patch (either no region matched or containment failed),
    // remove SBM-related keys so NurbsGeometryModelerSbm goes body-fitted for this patch.
    if (!use_sbm_for_this_patch) {
        if (patch_geometry.Has("skin_model_part_inner_initial_name")) patch_geometry.RemoveValue("skin_model_part_inner_initial_name");
        if (patch_geometry.Has("skin_model_part_outer_initial_name")) patch_geometry.RemoveValue("skin_model_part_outer_initial_name");
        if (patch_geometry.Has("skin_model_part_name")) patch_geometry.RemoveValue("skin_model_part_name");
        if (mEchoLevel > 1) {
            KRATOS_INFO("MultipatchModeler")
                << "  Disabling SBM for patch '" << patch_full_name << "' ("
                << (is_base_patch ? "base patch (no/invalid skin)" : "refinement patch") << ")." << std::endl;
        }
    }

    // Compute per-patch knot spans proportionally, then apply refinement overrides if present
    Vector patch_spans = patch_geometry["number_of_knot_spans"].GetVector();
    const int base_span_u = static_cast<int>(base_spans[0]);
    const int base_span_v = static_cast<int>(base_spans[1]);
    const double u_ratio  = (mBaseRect[U_MAX] - mBaseRect[U_MIN]) / base_u_length;
    const double v_ratio  = (mBaseRect[V_MAX] - mBaseRect[V_MIN]) / base_v_length;
    patch_spans[0] = std::max(1, static_cast<int>(std::round(u_ratio * base_span_u)));
    patch_spans[1] = std::max(1, static_cast<int>(std::round(v_ratio * base_span_v)));

    // Store span sizes as KNOT_SPAN_SIZES on the patch model part
    if (patch_spans.size() >= 2) {
        Vector patch_span_sizes(patch_spans.size());
        for (IndexType span_index = 0; span_index < patch_spans.size(); ++span_index) {
            const double span_count = patch_spans[span_index];
            double span_length = 0.0;
            // First two dimensions are U and V
            span_length = (span_index == 0)
                ? (patch_upper_uvw[0] - patch_lower_uvw[0])
                : (span_index == 1)
                    ? (patch_upper_uvw[1] - patch_lower_uvw[1])
                    : (patch_upper_uvw[span_index] - patch_lower_uvw[span_index]);
            patch_span_sizes[span_index] = span_count > 0.0 ? span_length / span_count : 0.0;
        }
        r_patch_model_part.SetValue(KNOT_SPAN_SIZES, patch_span_sizes);
    }

    // Write PARAMETER_SPACE_CORNERS for this patch (vector and matrix forms)
    std::vector<Vector> patch_parameter_corners(2);
    patch_parameter_corners[0].resize(2);
    patch_parameter_corners[1].resize(2);
    patch_parameter_corners[0][0] = patch_lower_uvw[0];
    patch_parameter_corners[0][1] = patch_upper_uvw[0];
    patch_parameter_corners[1][0] = patch_lower_uvw[1];
    patch_parameter_corners[1][1] = patch_upper_uvw[1];
    r_patch_model_part.SetValue(PARAMETER_SPACE_CORNERS, patch_parameter_corners);

    Matrix patch_parameter_corners_matrix(2, 2, 0.0);
    patch_parameter_corners_matrix(0,0) = patch_lower_uvw[0];
    patch_parameter_corners_matrix(0,1) = patch_upper_uvw[0];
    patch_parameter_corners_matrix(1,0) = patch_lower_uvw[1];
    patch_parameter_corners_matrix(1,1) = patch_upper_uvw[1];
    r_patch_model_part.SetValue(PATCH_PARAMETER_SPACE_CORNERS, patch_parameter_corners_matrix);

    if (mEchoLevel > 0) {
        KRATOS_INFO("MultipatchModeler")
            << "  Divisions (u, v) .......: ["
            << static_cast<int>(patch_spans[0]) << ", "
            << static_cast<int>(patch_spans[1]) << "]"
            << (use_sbm_for_this_patch ? " [refined]" : "") << std::endl;
    }

    // Generate geometry for this patch
    {
        NurbsGeometryModelerSbm geometry_modeler(*mpModel, patch_geometry);
        geometry_modeler.SetupGeometryModel();
        geometry_modeler.PrepareGeometryModel();
        geometry_modeler.SetupModelPart();
    }

    // Detect newly created geometries for BREP classification
    std::vector<ModelPart::IndexType> new_patch_geometry_ids;
    {
        std::vector<ModelPart::IndexType> patch_geometry_ids_after;
        patch_geometry_ids_after.reserve(r_patch_model_part.NumberOfGeometries());
        for (const auto& r_geom : r_patch_model_part.Geometries()) {
            patch_geometry_ids_after.push_back(r_geom.Id());
        }
        for (const auto id : patch_geometry_ids_after) {
            if (std::find(patch_geometry_ids_before.begin(), patch_geometry_ids_before.end(), id) ==
                patch_geometry_ids_before.end()) {
                new_patch_geometry_ids.push_back(id);
            }
        }
        std::sort(patch_geometry_ids_after.begin(), patch_geometry_ids_after.end());
    }

    KRATOS_ERROR_IF(!analysis_base.Has("analysis_model_part_name"))
        << "MultipatchModeler: missing analysis parameters in input." << std::endl;

    Parameters patch_analysis = analysis_base.Clone();
    patch_analysis["analysis_model_part_name"].SetString(patch_full_name);

    // Rebuild element_condition_list for this patch:
    // - normalize iga_model_part names to patch scope
    // - for body-fitted patches: set brep_ids internal/external split

    if (patch_analysis.Has("element_condition_list") && patch_analysis["element_condition_list"].IsArray()) {
        Parameters condition_list = patch_analysis["element_condition_list"];
        Parameters reduced("[]");

        // --- keep elements as-is (ensure geometry_type)
        for (IndexType i = 0; i < condition_list.size(); ++i) {
            const Parameters e = condition_list[i];
            if (e.Has("type") && e["type"].IsString() && e["type"].GetString() == "element") {
                Parameters elem = e.Clone();
                if (!elem.Has("geometry_type")) elem.AddEmptyValue("geometry_type").SetString("GeometrySurface");
                reduced.Append(elem);
            }
        }

        // --- pick ONE external and ONE internal template from the original list
        auto is_edge_cond = [](const Parameters& c)->bool {
            return c.Has("type") && c["type"].IsString() && c["type"].GetString() == "condition";
        };
        auto has_mp = [](const Parameters& c)->bool {
            return c.Has("iga_model_part") && c["iga_model_part"].IsString();
        };
        auto to_lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(),
                                    [](unsigned char ch){ return static_cast<char>(std::tolower(ch)); }); return s; };

        Parameters ext_tmpl, int_tmpl;
        bool has_ext = false, has_int = false;

        // heuristic by name contains "outer/external" vs "inner/internal", else by first/second occurrence
        for (IndexType i = 0; i < condition_list.size(); ++i) {
            const Parameters c = condition_list[i];
            if (!is_edge_cond(c) || !has_mp(c)) continue;
            const std::string mp = c["iga_model_part"].GetString();
            const std::string l  = to_lower(mp);
            if (!has_ext && (l.find("outer") != std::string::npos || l.find("external") != std::string::npos)) {
                ext_tmpl = c.Clone(); has_ext = true; continue;
            }
            if (!has_int && (l.find("inner") != std::string::npos || l.find("internal") != std::string::npos)) {
                int_tmpl = c.Clone(); has_int = true; continue;
            }
        }
        // positional fallback if keywords not found
        for (IndexType i = 0; (i < condition_list.size()) && (!has_ext || !has_int); ++i) {
            const Parameters c = condition_list[i];
            if (!is_edge_cond(c) || !has_mp(c)) continue;
            if (!has_ext) { ext_tmpl = c.Clone(); has_ext = true; }
            else if (!has_int) { int_tmpl = c.Clone(); has_int = true; }
        }
        if (!has_ext && !has_int) { patch_analysis["element_condition_list"] = reduced; /* nothing to add */ return; }
        if (!has_ext && has_int)  { ext_tmpl = int_tmpl.Clone(); has_ext = true; }
        if (!has_int && has_ext)  { int_tmpl = ext_tmpl.Clone(); has_int = true; }

        // --- classify edges once per patch (do NOT sort IDs, keep creation order)
        const auto is_close = [](double a, double b) {
            constexpr double tol = 1e-10;
            return std::abs(a - b) <= tol * (1.0 + std::max(std::abs(a), std::abs(b)));
        };
        const bool is_bottom_external = is_close(patch_lower_uvw[1], base_lower_uvw[1]);
        const bool is_top_external    = is_close(patch_upper_uvw[1], base_upper_uvw[1]);
        const bool is_right_external  = is_close(patch_upper_uvw[0], base_upper_uvw[0]);
        const bool is_left_external   = is_close(patch_lower_uvw[0], base_lower_uvw[0]);

        ModelPart::IndexType surface_id = 0;
        std::vector<ModelPart::IndexType> brep_ids;
        brep_ids.reserve(r_patch_model_part.NumberOfGeometries());
        for (const auto& r_geom : r_patch_model_part.Geometries()) {
            const int ws_dim = r_geom.WorkingSpaceDimension(); // 2 for surface, 1 for edge
            if (ws_dim >= 2 && surface_id == 0) { surface_id = r_geom.Id(); continue; }
            brep_ids.push_back(r_geom.Id());
        }

        std::vector<ModelPart::IndexType> external_ids, internal_ids;
        for (std::size_t edge_idx = 0; edge_idx < brep_ids.size(); ++edge_idx) {
            const auto id = brep_ids[edge_idx];
            bool ext = false;
            switch (edge_idx) {
                case 0: ext = is_bottom_external; break;
                case 1: ext = is_right_external;  break;
                case 2: ext = is_top_external;    break;
                case 3: ext = is_left_external;   break;
                default: ext = false;             break;
            }
            (ext ? external_ids : internal_ids).push_back(id);
        }

        const auto set_brep_ids = [](Parameters& p, const std::vector<ModelPart::IndexType>& ids){
            if (p.Has("brep_ids")) p.RemoveValue("brep_ids");
            if (!ids.empty()) {
                Vector v(static_cast<int>(ids.size()));
                for (std::size_t i = 0; i < ids.size(); ++i) v[i] = static_cast<double>(ids[i]);
                p.AddEmptyValue("brep_ids").SetVector(v);
            }
        };

        if (!external_ids.empty()) {
            Parameters ext = ext_tmpl.Clone();
            if (!ext.Has("geometry_type")) ext.AddEmptyValue("geometry_type").SetString("SurfaceEdge");
            // keep whatever condition "name" and order the external template specified
            set_brep_ids(ext, external_ids);
            // optionally normalize model part name:
            // ext["iga_model_part"].SetString("external_boundaries");
            reduced.Append(ext);
        }

        if (!internal_ids.empty()) {
            Parameters in = int_tmpl.Clone();
            if (!in.Has("geometry_type")) in.AddEmptyValue("geometry_type").SetString("SurfaceEdge");
            // keep internal template "name" and order (can differ from external)
            set_brep_ids(in, internal_ids);
        }

        patch_analysis["element_condition_list"] = reduced;
        
    }

    // Run IGA modeler for this patch
    {
        // For SBM-base-patch, also inform IgaModelerSbm which skin to use
        if (use_sbm_for_this_patch) {
            if (!patch_analysis.Has("skin_model_part_name"))
                patch_analysis.AddEmptyValue("skin_model_part_name");
            if (!sbm_skin_name_for_base.empty())
                patch_analysis["skin_model_part_name"].SetString(sbm_skin_name_for_base);
        }

        KRATOS_INFO_IF("MultipatchModeler::  \n", mEchoLevel > 3)
            << patch_analysis << std::endl;
        IgaModelerSbm iga_modeler(*mpModel, patch_analysis);
        iga_modeler.SetupModelPart();
        KRATOS_INFO_IF("MultipatchModeler:: called IgaModelerSbm", mEchoLevel > 1) << std::endl;
    }

    // Create one simple refinement patch (body-fitted) if a refinement region exists
    if (mSubdomains.size() > 1) {
        const auto& ref_rect = mSubdomains[1];
        ProcessRefPatch_(ref_rect, geometry_base, analysis_base, r_parent_model_part, prefix);
    }

    // Call the intersection process with matching echo level and requested coupling condition
    std::string coupling_condition_name = "";
    if (mParameters.Has("coupling_conditions_name") && mParameters["coupling_conditions_name"].IsString()) {
        coupling_condition_name = mParameters["coupling_conditions_name"].GetString();
    } else {
        KRATOS_ERROR << "MultipatchModeler: no coupling_conditions_name defined" << std::endl;
    }

    // TODO: extend to more than 1 ref region
    // For the multipatch layout (base + one refinement), build coupling only around the refinement patch
    if (mSubdomains.size() > 1) {
        RefPatchCouplingProcess ref_coupling_process(
            r_parent_model_part,
            static_cast<int>(mEchoLevel),
            1e-12,
            prefix,
            /*base_patch_index*/1,
            /*ref_patch_index*/2,
            coupling_condition_name);
        ref_coupling_process.Execute();
    }

    // Collect analysis targets from parameters and build global submodel parts
    std::vector<std::string> targets;
    if (mParameters.Has("analysis_parameters")) {
        const auto& analysis = mParameters["analysis_parameters"];
        if (analysis.Has("element_condition_list") && analysis["element_condition_list"].IsArray()) {
            const auto& ecl = analysis["element_condition_list"];
            for (IndexType i = 0; i < ecl.size(); ++i) {
                if (ecl[i].Has("iga_model_part") && ecl[i]["iga_model_part"].IsString()) {
                    const std::string mp = ecl[i]["iga_model_part"].GetString();
                    if (!mp.empty()) targets.push_back(mp);
                }
            }
        }
    }
    std::sort(targets.begin(), targets.end());
    targets.erase(std::unique(targets.begin(), targets.end()), targets.end());


}

void MultipatchModeler::ProcessRefPatch_(
    const RectType& rect,
    const Parameters& geometry_base,
    const Parameters& analysis_base,
    ModelPart& r_parent_model_part,
    const std::string& prefix)
{
    const std::string patch_suffix = prefix + std::to_string(2);
    ModelPart& r_patch_model_part = r_parent_model_part.HasSubModelPart(patch_suffix)
        ? r_parent_model_part.GetSubModelPart(patch_suffix)
        : r_parent_model_part.CreateSubModelPart(patch_suffix);
    const std::string patch_full_name = r_patch_model_part.FullName();

    // Geometry parameters cloned and set to rect (body-fitted)
    Parameters patch_geometry = geometry_base.Clone();
    patch_geometry["model_part_name"].SetString(patch_full_name);

    Vector patch_lower_uvw = patch_geometry["lower_point_uvw"].GetVector();
    Vector patch_upper_uvw = patch_geometry["upper_point_uvw"].GetVector();
    patch_lower_uvw[0] = rect[U_MIN]; patch_upper_uvw[0] = rect[U_MAX];
    patch_lower_uvw[1] = rect[V_MIN]; patch_upper_uvw[1] = rect[V_MAX];
    patch_geometry["lower_point_uvw"].SetVector(patch_lower_uvw);
    patch_geometry["upper_point_uvw"].SetVector(patch_upper_uvw);

    const auto& base_lower_xyz = geometry_base["lower_point_xyz"].GetVector();
    const auto& base_upper_xyz = geometry_base["upper_point_xyz"].GetVector();
    const auto& base_lower_uvw = geometry_base["lower_point_uvw"].GetVector();
    const auto& base_upper_uvw = geometry_base["upper_point_uvw"].GetVector();
    const auto& base_spans    = geometry_base["number_of_knot_spans"].GetVector();
    const double base_u_length = base_upper_uvw[0] - base_lower_uvw[0];
    const double base_v_length = base_upper_uvw[1] - base_lower_uvw[1];

    Vector patch_lower_xyz = patch_geometry["lower_point_xyz"].GetVector();
    Vector patch_upper_xyz = patch_geometry["upper_point_xyz"].GetVector();
    const double u_fraction_min = (rect[U_MIN] - base_lower_uvw[0]) / base_u_length;
    const double u_fraction_max = (rect[U_MAX] - base_lower_uvw[0]) / base_u_length;
    const double v_fraction_min = (rect[V_MIN] - base_lower_uvw[1]) / base_v_length;
    const double v_fraction_max = (rect[V_MAX] - base_lower_uvw[1]) / base_v_length;
    patch_lower_xyz[0] = base_lower_xyz[0] + u_fraction_min * (base_upper_xyz[0] - base_lower_xyz[0]);
    patch_upper_xyz[0] = base_lower_xyz[0] + u_fraction_max * (base_upper_xyz[0] - base_lower_xyz[0]);
    patch_lower_xyz[1] = base_lower_xyz[1] + v_fraction_min * (base_upper_xyz[1] - base_lower_xyz[1]);
    patch_upper_xyz[1] = base_lower_xyz[1] + v_fraction_max * (base_upper_xyz[1] - base_lower_xyz[1]);
    patch_lower_xyz[2] = base_lower_xyz[2];
    patch_upper_xyz[2] = base_upper_xyz[2];
    patch_geometry["lower_point_xyz"].SetVector(patch_lower_xyz);
    patch_geometry["upper_point_xyz"].SetVector(patch_upper_xyz);

    // Determine if this ref patch should use SBM (inner skin fully inside rect)
    const bool has_skin_data = geometry_base.Has("skin_model_part_name") &&
                               (geometry_base.Has("skin_model_part_outer_initial_name") ||
                                geometry_base.Has("skin_model_part_inner_initial_name"));
    bool skin_fully_inside_ref = true;
    if (geometry_base.Has("skin_model_part_inner_initial_name")) {
        const std::string inner_skin_name = geometry_base["skin_model_part_inner_initial_name"].GetString();
        skin_fully_inside_ref = IsSkinFullyInsidePatch_(inner_skin_name, rect, geometry_base);
        KRATOS_INFO_IF("MultipatchModeler", mEchoLevel > 2)
            << "[RefPatch] Skin containment (inner='" << inner_skin_name << "'): "
            << (skin_fully_inside_ref ? "inside" : "outside") << std::endl;
    }
    const bool use_sbm_for_ref_patch = has_skin_data && skin_fully_inside_ref;

    // If not using SBM, strip skin keys so body-fitted path is used
    if (!use_sbm_for_ref_patch) {
        if (patch_geometry.Has("skin_model_part_inner_initial_name")) patch_geometry.RemoveValue("skin_model_part_inner_initial_name");
        if (patch_geometry.Has("skin_model_part_outer_initial_name")) patch_geometry.RemoveValue("skin_model_part_outer_initial_name");
        if (patch_geometry.Has("skin_model_part_name")) patch_geometry.RemoveValue("skin_model_part_name");
    }

    // Proportional spans
    Vector patch_spans = patch_geometry["number_of_knot_spans"].GetVector();
    const int base_span_u = static_cast<int>(base_spans[0]);
    const int base_span_v = static_cast<int>(base_spans[1]);
    const double u_ratio  = (rect[U_MAX] - rect[U_MIN]) / base_u_length;
    const double v_ratio  = (rect[V_MAX] - rect[V_MIN]) / base_v_length;
    patch_spans[0] = std::max(1, static_cast<int>(std::round(u_ratio * base_span_u)));
    patch_spans[1] = std::max(1, static_cast<int>(std::round(v_ratio * base_span_v)));
    patch_geometry["number_of_knot_spans"].SetVector(patch_spans);

    // Store span sizes
    if (patch_spans.size() >= 2) {
        Vector patch_span_sizes(patch_spans.size());
        double span_length_u = (patch_upper_uvw[0] - patch_lower_uvw[0]);
        double span_length_v = (patch_upper_uvw[1] - patch_lower_uvw[1]);
        patch_span_sizes[0] = patch_spans[0] > 0.0 ? span_length_u / patch_spans[0] : 0.0;
        patch_span_sizes[1] = patch_spans[1] > 0.0 ? span_length_v / patch_spans[1] : 0.0;
        r_patch_model_part.SetValue(KNOT_SPAN_SIZES, patch_span_sizes);
    }

    // Write PARAMETER_SPACE_CORNERS
    {
        std::vector<Vector> patch_parameter_corners(2);
        patch_parameter_corners[0].resize(2);
        patch_parameter_corners[1].resize(2);
        patch_parameter_corners[0][0] = patch_lower_uvw[0];
        patch_parameter_corners[0][1] = patch_upper_uvw[0];
        patch_parameter_corners[1][0] = patch_lower_uvw[1];
        patch_parameter_corners[1][1] = patch_upper_uvw[1];
        r_patch_model_part.SetValue(PARAMETER_SPACE_CORNERS, patch_parameter_corners);

        Matrix patch_parameter_corners_matrix(2, 2, 0.0);
        patch_parameter_corners_matrix(0,0) = patch_lower_uvw[0];
        patch_parameter_corners_matrix(0,1) = patch_upper_uvw[0];
        patch_parameter_corners_matrix(1,0) = patch_lower_uvw[1];
        patch_parameter_corners_matrix(1,1) = patch_upper_uvw[1];
        r_patch_model_part.SetValue(PATCH_PARAMETER_SPACE_CORNERS, patch_parameter_corners_matrix);
    }

    // Create geometry
    {
        NurbsGeometryModelerSbm geometry_modeler(*mpModel, patch_geometry);
        geometry_modeler.SetupGeometryModel();
        geometry_modeler.PrepareGeometryModel();
        geometry_modeler.SetupModelPart();
    }

    // Analysis: elements only by default; if SBM, append sbm-conditions for inner hole
    if (analysis_base.Has("analysis_model_part_name")) {
        Parameters patch_analysis = analysis_base.Clone();
        patch_analysis["analysis_model_part_name"].SetString(patch_full_name);

        if (patch_analysis.Has("element_condition_list") && patch_analysis["element_condition_list"].IsArray()) {
            Parameters condition_list = patch_analysis["element_condition_list"];
            Parameters reduced_elements("[]");
            for (IndexType i = 0; i < condition_list.size(); ++i) {
                const Parameters e = condition_list[i];
                if (e.Has("type") && e["type"].IsString() && e["type"].GetString() == "element") {
                    Parameters elem = e.Clone();
                    if (!elem.Has("geometry_type")) elem.AddEmptyValue("geometry_type").SetString("GeometrySurface");
                    reduced_elements.Append(elem);
                }
            }
            if (use_sbm_for_ref_patch) {
                for (IndexType i = 0; i < condition_list.size(); ++i) {
                    Parameters c = condition_list[i].Clone();
                    if (!(c.Has("type") && c["type"].IsString() && c["type"].GetString() == "condition")) continue;
                    if (!c.Has("sbm_parameters")) continue;
                    if (!c.Has("geometry_type")) c.AddEmptyValue("geometry_type").SetString("SurfaceEdge");
                    reduced_elements.Append(c);
                }
            }
            patch_analysis["element_condition_list"] = reduced_elements;
        }

        if (use_sbm_for_ref_patch && geometry_base.Has("skin_model_part_name")) {
            if (!patch_analysis.Has("skin_model_part_name"))
                patch_analysis.AddEmptyValue("skin_model_part_name");
            patch_analysis["skin_model_part_name"].SetString(geometry_base["skin_model_part_name"].GetString());
        }

        KRATOS_WATCH(patch_analysis)

        IgaModelerSbm iga_modeler(*mpModel, patch_analysis);
        iga_modeler.SetupModelPart();
    } else {
        KRATOS_ERROR << "MultipatchModeler: missing analysis parameters in input." << std::endl;
    }
}

// Validates inputs, prepares defaults, gathers analysis targets,
// logs, generates subdivision, and writes PARAMETER_SPACE_CORNERS on the main ModelPart.
ModelPart& MultipatchModeler::InitializeSetup_(
    Parameters& r_parameters,
    Parameters& rGeometryBaseOut,
    Parameters& rAnalysisBaseOut,
    std::vector<std::string>& rAnalysisTargetsOut,
    std::unordered_map<std::string, std::vector<std::string>>& rCondNameToPartsOut)
{
    // Required: model part name
    const std::string& model_part_name = r_parameters["model_part_name"].GetString();
    KRATOS_ERROR_IF(model_part_name.empty())
        << "MultipatchModeler: \"model_part_name\" must be specified." << std::endl;

    ModelPart& r_model_part = mpModel->GetModelPart(model_part_name);

    // Ensure DOMAIN_SIZE is set on the root process info (some solvers set it later).
    // For surface IGA in this modeler, default to 2 if unset/zero.
    {
        auto& r_process_info = r_model_part.GetProcessInfo();
        int current_domain_size = 0;
        if (r_process_info.Has(DOMAIN_SIZE)) {
            current_domain_size = static_cast<int>(r_process_info[DOMAIN_SIZE]);
        }
        if (current_domain_size == 0) {
            r_process_info.SetValue(DOMAIN_SIZE, 2);
        }
    }

    // Echo level
    const int echo_level_input = r_parameters["echo_level"].GetInt();
    mEchoLevel = static_cast<SizeType>(echo_level_input);

    // Presence flags for geometry/analysis blocks
    const bool has_geometry_params =
        r_parameters.Has("geometry_parameters") &&
        r_parameters["geometry_parameters"].Has("model_part_name");

    const bool has_analysis_params =
        r_parameters.Has("analysis_parameters") &&
        r_parameters["analysis_parameters"].Has("analysis_model_part_name");

    // Copy base parameter blocks as working templates
    if (has_geometry_params) rGeometryBaseOut = r_parameters["geometry_parameters"];
    if (has_analysis_params) rAnalysisBaseOut = r_parameters["analysis_parameters"];

    // Small utility: set default vector if missing or empty
    const auto assign_default_vector = [](Parameters& p, const std::string& key, const Vector& def) {
        if (!p.Has(key)) { p.AddEmptyValue(key).SetVector(def); return; }
        const auto cur = p[key].GetVector();
        if (cur.size() == 0) p[key].SetVector(def);
    };

    // Ensure geometry base has consistent UVW/XYZ bounds defaults
    if (has_geometry_params) {
        const Vector base_lower_uvw = r_parameters["base_domain"]["lower_point_uvw"].GetVector();
        const Vector base_upper_uvw = r_parameters["base_domain"]["upper_point_uvw"].GetVector();

        if (!rGeometryBaseOut.Has("model_part_name") ||
            rGeometryBaseOut["model_part_name"].GetString().empty()) {
            rGeometryBaseOut.AddEmptyValue("model_part_name").SetString(model_part_name);
        }

        // Keep behavior: use UVW defaults also for XYZ if missing
        assign_default_vector(rGeometryBaseOut, "lower_point_uvw", base_lower_uvw);
        assign_default_vector(rGeometryBaseOut, "upper_point_uvw", base_upper_uvw);
        assign_default_vector(rGeometryBaseOut, "lower_point_xyz", base_lower_uvw);
        assign_default_vector(rGeometryBaseOut, "upper_point_xyz", base_upper_uvw);
    }

    // Collect unique analysis target submodel parts and map condition-name -> parts
    if (has_analysis_params) {
        const auto& elem_cond_list = rAnalysisBaseOut["element_condition_list"];
        std::vector<std::string> targets;
        for (IndexType i = 0; i < elem_cond_list.size(); ++i) {
            if (elem_cond_list[i].Has("iga_model_part")) {
                targets.push_back(elem_cond_list[i]["iga_model_part"].GetString());
            }
            if (elem_cond_list[i].Has("name")) {
                const std::string cond_name = elem_cond_list[i]["name"].GetString();
                if (!cond_name.empty()) {
                    rCondNameToPartsOut[cond_name].push_back(
                        elem_cond_list[i].Has("iga_model_part")
                            ? elem_cond_list[i]["iga_model_part"].GetString()
                            : std::string());
                }
            }
        }
        std::sort(targets.begin(), targets.end());
        targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
        rAnalysisTargetsOut.swap(targets);
    }

    // Verbose logging (unchanged semantics)
    if (mEchoLevel > 2) {
        KRATOS_INFO("MultipatchModeler")
            << "SetupModelPart for model part '" << model_part_name << "'" << std::endl;
        KRATOS_INFO("MultipatchModeler")
            << "  echo_level ............: " << mEchoLevel << std::endl;

        const auto& base_lower = r_parameters["base_domain"]["lower_point_uvw"].GetVector();
        const auto& base_upper = r_parameters["base_domain"]["upper_point_uvw"].GetVector();
        KRATOS_INFO("MultipatchModeler")
            << "  base_domain UV ........: [" << base_lower[0] << ", " << base_upper[0]
            << "] x [" << base_lower[1] << ", " << base_upper[1] << "]" << std::endl;

        const auto& refinement_array = r_parameters["refinement_regions"];
        KRATOS_INFO("MultipatchModeler")
            << "  refinement_regions ....: " << refinement_array.size() << std::endl;
        for (IndexType i = 0; i < refinement_array.size(); ++i) {
            const auto& lower = refinement_array[i]["lower_point_uvw"].GetVector();
            const auto& upper = refinement_array[i]["upper_point_uvw"].GetVector();
            KRATOS_INFO("MultipatchModeler")
                << "    - region " << i + 1 << " -> [" << lower[0] << ", " << upper[0]
                << "] x [" << lower[1] << ", " << upper[1] << "]" << std::endl;
        }

        KRATOS_INFO("MultipatchModeler")
            << "  geometry_parameters ...: " << (has_geometry_params ? "provided" : "not provided") << std::endl;
        if (has_geometry_params)
            KRATOS_INFO("MultipatchModeler") << rGeometryBaseOut.PrettyPrintJsonString() << std::endl;

        KRATOS_INFO("MultipatchModeler")
            << "  analysis_parameters ...: " << (has_analysis_params ? "provided" : "not provided") << std::endl;
        if (has_analysis_params)
            KRATOS_INFO("MultipatchModeler") << rAnalysisBaseOut.PrettyPrintJsonString() << std::endl;
    }

    // Build subdomains based on inputs
    GenerateSubdivision();

    if (mEchoLevel > 0) {
        KRATOS_INFO("MultipatchModeler")
            << "Generated " << mSubdomains.size() << " subdomains" << std::endl;

        if (mEchoLevel > 1) {
            for (std::size_t i = 0; i < mSubdomains.size(); ++i) {
                const auto& rect = mSubdomains[i];
                KRATOS_INFO("MultipatchModeler")
                    << "  Subdomain " << i + 1 << ": " << RectangleToString(rect) << std::endl;
            }
        }
    }

    // Write all rectangles to PARAMETER_SPACE_CORNERS on the main model part
    std::vector<Vector> rectangles_data;
    rectangles_data.reserve(mSubdomains.size());
    for (const auto& rect : mSubdomains) {
        Vector data(4);
        data[0] = rect[U_MIN]; data[1] = rect[U_MAX];
        data[2] = rect[V_MIN]; data[3] = rect[V_MAX];
        rectangles_data.push_back(std::move(data));
    }
    if (mEchoLevel > 2) {
        KRATOS_INFO("MultipatchModeler")
            << "Assigning " << rectangles_data.size()
            << " rectangles to PARAMETER_SPACE_CORNERS" << std::endl;
    }
    r_model_part.SetValue(PARAMETER_SPACE_CORNERS, rectangles_data);

    return r_model_part;
}















void MultipatchModeler::GenerateSubdivision()
{
    // New strategy:
    // - Always create ONE external patch equal to the base domain
    // - For each refinement region, create ONE body-fitted patch equal to the region (clipped to base)
    // No grid subdivision into many small patches.

    mSubdomains.clear();
    mRefinementRegions.clear();

    const auto& r_parameters = mParameters;

    const auto& base_lower = r_parameters["base_domain"]["lower_point_uvw"].GetVector();
    const auto& base_upper = r_parameters["base_domain"]["upper_point_uvw"].GetVector();

    mBaseRect = RectangleType{base_lower[0], base_upper[0], base_lower[1], base_upper[1]};

    if (mEchoLevel > 2) {
        KRATOS_INFO("MultipatchModeler")
            << "GenerateSubdivision (embedded) -> base rectangle: " << RectangleToString(mBaseRect) << std::endl;
    }

    // 1) Add the external/base patch
    if (HasPositiveArea(mBaseRect)) {
        mSubdomains.push_back(mBaseRect);
    } else {
        KRATOS_ERROR << "MultipatchModeler: base_domain has non-positive area." << std::endl;
    }

    // 2) Add one patch per refinement region
    const auto& refinement_array = r_parameters["refinement_regions"];
    for (IndexType i = 0; i < refinement_array.size(); ++i) {
        const auto& ref_lower = refinement_array[i]["lower_point_uvw"].GetVector();
        const auto& ref_upper = refinement_array[i]["upper_point_uvw"].GetVector();

        RectangleType region{ref_lower[0], ref_upper[0], ref_lower[1], ref_upper[1]};
        RectangleType clipped = ClipRectangle(region, mBaseRect);

        if (mEchoLevel > 2) {
            KRATOS_INFO("MultipatchModeler")
                << "  Refinement region " << i + 1 << " requested: " << RectangleToString(region) << std::endl;
            if (!HasPositiveArea(clipped)) {
                KRATOS_INFO("MultipatchModeler")
                    << "    -> skipped (no overlap with base domain)." << std::endl;
            } else if (clipped != region) {
                KRATOS_INFO("MultipatchModeler")
                    << "    -> clipped to: " << RectangleToString(clipped) << std::endl;
            } else {
                KRATOS_INFO("MultipatchModeler")
                    << "    -> used as provided." << std::endl;
            }
        }
        if (!HasPositiveArea(clipped)) {
            continue;
        }

        // Store refinement-region-specific overrides
        RefinementRegionData region_data;
        region_data.Rectangle = clipped;

        if (refinement_array[i].Has("polynomial_order")) {
            Parameters poly_param = refinement_array[i]["polynomial_order"];
            const std::size_t size = poly_param.size();
            KRATOS_ERROR_IF(size == 0)
                << "MultipatchModeler: refinement region polynomial_order array is empty." << std::endl;
            region_data.HasPolynomialOrder = true;
            region_data.PolynomialOrder.resize(size);
            for (IndexType j = 0; j < size; ++j) {
                region_data.PolynomialOrder[j] = poly_param.GetArrayItem(j).GetInt();
            }
        }

        if (refinement_array[i].Has("number_of_knot_spans")) {
            Parameters span_param = refinement_array[i]["number_of_knot_spans"];
            const std::size_t size = span_param.size();
            KRATOS_ERROR_IF(size == 0)
                << "MultipatchModeler: refinement region number_of_knot_spans array is empty." << std::endl;
            region_data.HasNumberOfKnotSpans = true;
            region_data.NumberOfKnotSpans.resize(size);
            for (IndexType j = 0; j < size; ++j) {
                region_data.NumberOfKnotSpans[j] = span_param.GetArrayItem(j).GetInt();
            }
        }

        mRefinementRegions.push_back(region_data);
        mSubdomains.push_back(clipped);
    }

    if (mEchoLevel > 2) {
        KRATOS_INFO("MultipatchModeler")
            << "GenerateSubdivision (embedded) -> total patches: " << mSubdomains.size() << std::endl;
        for (std::size_t i = 0; i < mSubdomains.size(); ++i) {
            KRATOS_INFO("MultipatchModeler")
                << "  Patch[" << i << "] = " << RectangleToString(mSubdomains[i]) << std::endl;
        }
    }
}

bool MultipatchModeler::IsSkinFullyInsidePatch_(
    const std::string& rSkinModelPartName,
    const RectType& rect,
    const Parameters& geometry_base) const
{
    // physical == parameter space: use (u,v) directly
    const double patch_min_u = rect[U_MIN];
    const double patch_max_u = rect[U_MAX];
    const double patch_min_v = rect[V_MIN];
    const double patch_max_v = rect[V_MAX];

    const ModelPart& r_skin_mp = mpModel->GetModelPart(rSkinModelPartName);
    if (r_skin_mp.NumberOfGeometries() == 0 && r_skin_mp.NumberOfConditions() == 0) {
        return false;
    }

    // Sample each curve geometry to build skin bounds in (u,v)
    double skin_min_u = std::numeric_limits<double>::max();
    double skin_max_u = std::numeric_limits<double>::lowest();
    double skin_min_v = std::numeric_limits<double>::max();
    double skin_max_v = std::numeric_limits<double>::lowest();

    if (r_skin_mp.NumberOfGeometries() > 0) {
        const int samples_per_curve = 100; // few points are enough
        for (const auto& rGeom : r_skin_mp.Geometries()) {
            for (int s = 0; s <= samples_per_curve; ++s) {
                const double t = static_cast<double>(s) / samples_per_curve;
                array_1d<double,3> local; local[0] = t; local[1] = 0.0; local[2] = 0.0;
                array_1d<double,3> global;
                rGeom.GlobalCoordinates(global, local); // (x,y) == (u,v)
                const double u = global[0];
                const double v = global[1];
                if (u < skin_min_u) skin_min_u = u;
                if (u > skin_max_u) skin_max_u = u;
                if (v < skin_min_v) skin_min_v = v;
                if (v > skin_max_v) skin_max_v = v;
            }
        }
    }
    if (r_skin_mp.NumberOfConditions() > 0) {
        for (const auto& rCond : r_skin_mp.Conditions()) {
            const auto& g = rCond.GetGeometry();
            for (std::size_t i = 0; i < g.size(); ++i) {
                const auto& P = g[i];
                const double u = P.X();
                const double v = P.Y();
                if (u < skin_min_u) skin_min_u = u;
                if (u > skin_max_u) skin_max_u = u;
                if (v < skin_min_v) skin_min_v = v;
                if (v > skin_max_v) skin_max_v = v;
            }
        }
    }
    // Tolerance
    const double tol = 1e-12;
    const bool inside_u = (skin_min_u >= patch_min_u - tol) && (skin_max_u <= patch_max_u + tol);
    const bool inside_v = (skin_min_v >= patch_min_v - tol) && (skin_max_v <= patch_max_v + tol);

    return inside_u && inside_v;
}

void MultipatchModeler::BuildGlobalSubModelParts(ModelPart& r_parent_model_part, const std::vector<std::string>& rTargetNames) const
{
    ModelPart& r_root = r_parent_model_part; // your main "IgaModelPart"

    for (const std::string& target_name : rTargetNames) {
        // Create/get the global submodelpart for this target
        ModelPart& r_global = r_root.HasSubModelPart(target_name)
            ? r_root.GetSubModelPart(target_name)
            : r_root.CreateSubModelPart(target_name);

        // Clear previous memberships (does not delete entities from the root)
        r_global.Elements().clear();
        r_global.Conditions().clear();
        r_global.Nodes().clear();

        std::size_t added_elems = 0;
        std::size_t added_conds = 0;
        std::size_t added_nodes = 0;

        // Loop all Patch* submodel parts
        for (auto& rPatch : r_root.SubModelParts()) {
            const std::string& patch_name = rPatch.Name();
            if (patch_name.rfind("Patch", 0) != 0) continue; // only Patch*

            // Each patch should have a submodelpart with this target name
            if (!rPatch.HasSubModelPart(target_name)) {
                KRATOS_INFO_IF("MultipatchModeler", mEchoLevel > 2)
                    << "Patch '" << patch_name << "' has no '" << target_name << "' — skipping.\n";
                continue;
            }

            ModelPart& r_src = rPatch.GetSubModelPart(target_name);
            if (r_src.NumberOfGeometries() > 0) {
                r_global.AddGeometries(r_src.GeometriesBegin(), r_src.GeometriesEnd());
            }

            // 1) Nodes (prefer direct copy; otherwise harvest from elems/conds)
            if (r_src.NumberOfNodes() > 0) {
                r_global.AddNodes(r_src.NodesBegin(), r_src.NodesEnd());
                added_nodes += r_src.NumberOfNodes();
            } else {
                // Fallback: gather nodes from geometries (avoid duplicates)
                std::set<ModelPart::IndexType> seen_ids;
                for (const auto& rElem : r_src.Elements()) {
                    const auto& g = rElem.GetGeometry();
                    for (std::size_t i = 0; i < g.size(); ++i) {
                        const auto pN = g.pGetPoint(i);
                        if (pN && seen_ids.insert(pN->Id()).second) {
                            r_global.Nodes().push_back(pN);
                        }
                    }
                }
                for (const auto& rCond : r_src.Conditions()) {
                    const auto& g = rCond.GetGeometry();
                    for (std::size_t i = 0; i < g.size(); ++i) {
                        const auto pN = g.pGetPoint(i);
                        if (pN && seen_ids.insert(pN->Id()).second) {
                            r_global.Nodes().push_back(pN);
                        }
                    }
                }
                added_nodes += seen_ids.size();
            }
            // 2) Elements
            if (r_src.NumberOfElements() > 0) {
                r_global.AddElements(r_src.ElementsBegin(), r_src.ElementsEnd());
                added_elems += r_src.NumberOfElements();
            }
            // 3) Conditions
            if (r_src.NumberOfConditions() > 0) {
                r_global.AddConditions(r_src.ConditionsBegin(), r_src.ConditionsEnd());
                added_conds += r_src.NumberOfConditions();
            }


            KRATOS_INFO_IF("MultipatchModeler", mEchoLevel > 2)
                << "Collected from '" << rPatch.FullName() << "' into global '"
                << r_global.FullName() << "': elems=" << r_src.NumberOfElements()
                << ", conds=" << r_src.NumberOfConditions()
                << ", nodes~=" << r_src.NumberOfNodes() << "\n";
        }

        KRATOS_INFO_IF("MultipatchModeler", mEchoLevel > 0)
            << "[Global '" << target_name << "'] SubModelPart '" << r_global.FullName()
            << "' now has " << r_global.NumberOfElements() << " elements, "
            << r_global.NumberOfConditions() << " conditions and "
            << r_global.NumberOfNodes() << " nodes (inserted ~"
            << added_elems << " elems, ~" << added_conds << " conds, ~"
            << added_nodes << " nodes).\n";

    }

}

} // namespace Kratos
