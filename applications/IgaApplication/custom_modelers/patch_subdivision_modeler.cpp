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
#include <limits>

// External includes

// Project includes
#include "patch_subdivision_modeler.h"
#include "iga_application_variables.h"
#include "custom_processes/patch_intersection_process.h"

namespace Kratos
{

namespace
{
using RectangleType = PatchSubdivisionModeler::RectangleType;

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

PatchSubdivisionModeler::PatchSubdivisionModeler(Model& rModel, const Parameters ModelParameters)
    : Modeler(rModel, ModelParameters)
    , mpModel(&rModel)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

Modeler::Pointer PatchSubdivisionModeler::Create(Model& rModel, const Parameters ModelParameters) const
{
    return Kratos::make_shared<PatchSubdivisionModeler>(rModel, ModelParameters);
}

const Parameters PatchSubdivisionModeler::GetDefaultParameters() const
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
        "echo_level" : 0
    })");
}

const std::vector<PatchSubdivisionModeler::RectangleType>& PatchSubdivisionModeler::GetSubdomains() const
{
    return mSubdomains;
}

// PatchSubdivisionModeler.cpp

void PatchSubdivisionModeler::SetupModelPart()
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

    // One call per subdomain
    for (std::size_t i_patch = 0; i_patch < mSubdomains.size(); ++i_patch) {
        const auto& rect = mSubdomains[i_patch];
        ProcessPatch_(
            i_patch,
            rect,
            geometry_base,
            analysis_base,
            analysis_target_submodel_parts,
            condition_name_to_model_parts,
            r_parent_model_part,
            prefix);
    }

    // Call the intersection process
    PatchIntersectionProcess intersection_process(r_parent_model_part);
    intersection_process.Execute(); // HEREEEEE


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

    // Create unique sub_model_parts for each analysis target
    // this->BuildGlobalSubModelParts(r_model_part, targets);

}

// Validates inputs, prepares defaults, gathers analysis targets,
// logs, generates subdivision, and writes PARAMETER_SPACE_CORNERS on the main ModelPart.
ModelPart& PatchSubdivisionModeler::InitializeSetup_(
    Parameters& r_parameters,
    Parameters& rGeometryBaseOut,
    Parameters& rAnalysisBaseOut,
    std::vector<std::string>& rAnalysisTargetsOut,
    std::unordered_map<std::string, std::vector<std::string>>& rCondNameToPartsOut)
{
    // Required: model part name
    const std::string& model_part_name = r_parameters["model_part_name"].GetString();
    KRATOS_ERROR_IF(model_part_name.empty())
        << "PatchSubdivisionModeler: \"model_part_name\" must be specified." << std::endl;

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
        KRATOS_INFO("PatchSubdivisionModeler")
            << "SetupModelPart for model part '" << model_part_name << "'" << std::endl;
        KRATOS_INFO("PatchSubdivisionModeler")
            << "  echo_level ............: " << mEchoLevel << std::endl;

        const auto& base_lower = r_parameters["base_domain"]["lower_point_uvw"].GetVector();
        const auto& base_upper = r_parameters["base_domain"]["upper_point_uvw"].GetVector();
        KRATOS_INFO("PatchSubdivisionModeler")
            << "  base_domain UV ........: [" << base_lower[0] << ", " << base_upper[0]
            << "] x [" << base_lower[1] << ", " << base_upper[1] << "]" << std::endl;

        const auto& refinement_array = r_parameters["refinement_regions"];
        KRATOS_INFO("PatchSubdivisionModeler")
            << "  refinement_regions ....: " << refinement_array.size() << std::endl;
        for (IndexType i = 0; i < refinement_array.size(); ++i) {
            const auto& lower = refinement_array[i]["lower_point_uvw"].GetVector();
            const auto& upper = refinement_array[i]["upper_point_uvw"].GetVector();
            KRATOS_INFO("PatchSubdivisionModeler")
                << "    - region " << i + 1 << " -> [" << lower[0] << ", " << upper[0]
                << "] x [" << lower[1] << ", " << upper[1] << "]" << std::endl;
        }

        KRATOS_INFO("PatchSubdivisionModeler")
            << "  geometry_parameters ...: " << (has_geometry_params ? "provided" : "not provided") << std::endl;
        if (has_geometry_params)
            KRATOS_INFO("PatchSubdivisionModeler") << rGeometryBaseOut.PrettyPrintJsonString() << std::endl;

        KRATOS_INFO("PatchSubdivisionModeler")
            << "  analysis_parameters ...: " << (has_analysis_params ? "provided" : "not provided") << std::endl;
        if (has_analysis_params)
            KRATOS_INFO("PatchSubdivisionModeler") << rAnalysisBaseOut.PrettyPrintJsonString() << std::endl;
    }

    // Build subdomains based on inputs
    GenerateSubdivision();

    if (mEchoLevel > 0) {
        KRATOS_INFO("PatchSubdivisionModeler")
            << "Generated " << mSubdomains.size() << " subdomains" << std::endl;

        if (mEchoLevel > 1) {
            for (std::size_t i = 0; i < mSubdomains.size(); ++i) {
                const auto& rect = mSubdomains[i];
                KRATOS_INFO("PatchSubdivisionModeler")
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
        KRATOS_INFO("PatchSubdivisionModeler")
            << "Assigning " << rectangles_data.size()
            << " rectangles to PARAMETER_SPACE_CORNERS" << std::endl;
    }
    r_model_part.SetValue(PARAMETER_SPACE_CORNERS, rectangles_data);

    return r_model_part;
}

// Handles one subdomain: builds geometry and analysis for a single patch,
// applies overrides (knot spans, polynomial order), classifies BREP edges,
void PatchSubdivisionModeler::ProcessPatch_(
    std::size_t i_patch,
    const RectType& rect,
    const Parameters& geometry_base,
    const Parameters& analysis_base,
    const std::vector<std::string>& analysis_target_submodel_parts,
    const std::unordered_map<std::string, std::vector<std::string>>& condition_name_to_model_parts,
    ModelPart& r_parent_model_part,
    const std::string& prefix)
{
    // Patch naming
    const std::string patch_suffix = prefix + std::to_string(i_patch + 1);
    ModelPart& r_patch_model_part = r_parent_model_part.HasSubModelPart(patch_suffix)
        ? r_parent_model_part.GetSubModelPart(patch_suffix)
        : r_parent_model_part.CreateSubModelPart(patch_suffix);
    const std::string patch_full_name = r_patch_model_part.FullName(); // HEREEEEE

    if (mEchoLevel > 2) {
        KRATOS_INFO("PatchSubdivisionModeler")
            << "-- Debug: entering patch loop for " << patch_full_name << std::endl;
    }

    // Identify matching refinement region (if any)
    const RefinementRegionData* p_refinement_region = nullptr;
    {
        const double tol = 1e-12;
        for (const auto& r_region : mRefinementRegions) {
            if ((rect[U_MIN] >= r_region.Rectangle[U_MIN] - tol) &&
                (rect[U_MAX] <= r_region.Rectangle[U_MAX] + tol) &&
                (rect[V_MIN] >= r_region.Rectangle[V_MIN] - tol) &&
                (rect[V_MAX] <= r_region.Rectangle[V_MAX] + tol)) {
                p_refinement_region = &r_region;
                break;
            }
        }
    }

    // SBM is applicable only if a refinement region matched AND skin data is configured
    const bool has_skin_data = geometry_base.Has("skin_model_part_name") &&
                               (geometry_base.Has("skin_model_part_outer_initial_name") ||
                                geometry_base.Has("skin_model_part_inner_initial_name"));
    
    // Check if the inner skin (if provided) is fully inside this patch
    bool skin_fully_inside = true;
    if (geometry_base.Has("skin_model_part_inner_initial_name")) {
        const std::string inner_skin_name = geometry_base["skin_model_part_inner_initial_name"].GetString();
        skin_fully_inside = IsSkinFullyInsidePatch_(inner_skin_name, rect, geometry_base);
        if (mEchoLevel > 2) {
            KRATOS_INFO("PatchSubdivisionModeler")
                << "  Skin containment (inner='" << inner_skin_name << "'): "
                << (skin_fully_inside ? "inside" : "outside") << std::endl;
        }
    }

    bool use_sbm_for_this_patch = (p_refinement_region != nullptr) && has_skin_data && skin_fully_inside;

    if (mEchoLevel > 2) {
        KRATOS_INFO("PatchSubdivisionModeler")
            << "Creating/Updating patch '" << patch_full_name << "' from subdomain "
            << RectangleToString(rect) << std::endl;
        KRATOS_INFO("PatchSubdivisionModeler")
            << "  Patch classification ....: " << (use_sbm_for_this_patch ? "refined" : "body-fitted") << std::endl;
        if (!has_skin_data && p_refinement_region != nullptr) {
            KRATOS_INFO("PatchSubdivisionModeler")
                << "  Note: missing skin_model_part definitions -> surrogate workflow disabled." << std::endl;
        }
    }

    // Snapshot geometry IDs before generation (to detect new ones)
    std::vector<ModelPart::IndexType> patch_geometry_ids_before;
    patch_geometry_ids_before.reserve(r_patch_model_part.NumberOfGeometries());
    for (const auto& r_geom : r_patch_model_part.Geometries()) {
        patch_geometry_ids_before.push_back(r_geom.Id());
    }

    // Clone base geometry parameters, scope to this patch
    Parameters patch_geometry = geometry_base.Clone();
    patch_geometry["model_part_name"].SetString(patch_full_name);

    // Assign UVW bounds for the patch
    Vector patch_lower_uvw = patch_geometry["lower_point_uvw"].GetVector();
    Vector patch_upper_uvw = patch_geometry["upper_point_uvw"].GetVector();
    patch_lower_uvw[0] = rect[U_MIN]; patch_upper_uvw[0] = rect[U_MAX];
    patch_lower_uvw[1] = rect[V_MIN]; patch_upper_uvw[1] = rect[V_MAX];
    patch_geometry["lower_point_uvw"].SetVector(patch_lower_uvw);
    patch_geometry["upper_point_uvw"].SetVector(patch_upper_uvw);

    // Map UVW fractions into XYZ using base box and base UVW
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
            KRATOS_INFO("PatchSubdivisionModeler")
                << "  Disabling SBM for patch '" << patch_full_name << "' ("
                << (p_refinement_region ? "skin outside" : "not a refinement patch") << ")." << std::endl;
        }
    }

    // Compute per-patch knot spans proportionally, then apply refinement overrides if present
    Vector patch_spans = patch_geometry["number_of_knot_spans"].GetVector();
    const int base_span_u = static_cast<int>(base_spans[0]);
    const int base_span_v = static_cast<int>(base_spans[1]);
    const double u_ratio  = (rect[U_MAX] - rect[U_MIN]) / base_u_length;
    const double v_ratio  = (rect[V_MAX] - rect[V_MIN]) / base_v_length;
    patch_spans[0] = std::max(1, static_cast<int>(std::round(u_ratio * base_span_u)));
    patch_spans[1] = std::max(1, static_cast<int>(std::round(v_ratio * base_span_v)));

    if (p_refinement_region && p_refinement_region->HasNumberOfKnotSpans) {
        // Strict size match with patch dimension
        KRATOS_ERROR_IF(p_refinement_region->NumberOfKnotSpans.size() != patch_spans.size())
            << "PatchSubdivisionModeler: refinement region number_of_knot_spans size mismatch." << std::endl;
        for (IndexType i_dir = 0; i_dir < patch_spans.size(); ++i_dir) {
            patch_spans[i_dir] = static_cast<double>(p_refinement_region->NumberOfKnotSpans[i_dir]);
        }
    }
    patch_geometry["number_of_knot_spans"].SetVector(patch_spans);

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

    // Polynomial order override from refinement region, if any
    if (p_refinement_region && p_refinement_region->HasPolynomialOrder) {
        Vector patch_order = patch_geometry["polynomial_order"].GetVector();
        KRATOS_ERROR_IF(p_refinement_region->PolynomialOrder.size() != patch_order.size())
            << "PatchSubdivisionModeler: refinement region polynomial_order size mismatch." << std::endl;
        for (IndexType i_dir = 0; i_dir < patch_order.size(); ++i_dir) {
            patch_order[i_dir] = static_cast<double>(p_refinement_region->PolynomialOrder[i_dir]);
        }
        patch_geometry["polynomial_order"].SetVector(patch_order);
    }

    if (mEchoLevel > 0) {
        KRATOS_INFO("PatchSubdivisionModeler")
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

    // If analysis parameters exist, build per-patch analysis parameters and run IGA modeler
    if (analysis_base.Has("analysis_model_part_name")) {
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

            
            if (use_sbm_for_this_patch) {
                // For SBM patches: start from the computed body-fitted list (reduced)
                // and append only the original conditions that have sbm_parameters.
                Parameters reduced_sbm_first_stage("[]");
                for (IndexType i = 0; i < reduced.size(); ++i) {
                    reduced_sbm_first_stage.Append(reduced[i]);
                }

                auto is_cond_1 = [](const Parameters& c)->bool {
                    return c.Has("type") && c["type"].IsString() && c["type"].GetString() == "condition";
                };

                // condition_list contains the original per-analysis list from inputs
                for (IndexType i = 0; i < condition_list.size(); ++i) {
                    Parameters e = condition_list[i].Clone();
                    if (!is_cond_1(e)) continue;
                    if (!e.Has("sbm_parameters")) continue;
                    if (!e.Has("geometry_type")) e.AddEmptyValue("geometry_type").SetString("SurfaceEdge");
                    reduced_sbm_first_stage.Append(e);
                }
                patch_analysis["element_condition_list"] = reduced_sbm_first_stage;
                
                Parameters full_list = reduced_sbm_first_stage;

                // SBM patch: keep elements and only conditions with sbm_parameters or non-empty brep_ids
                Parameters reduced_sbm("[]");

                auto is_cond = [](const Parameters& c)->bool {
                    return c.Has("type") && c["type"].IsString() && c["type"].GetString() == "condition";
                };

                for (IndexType i = 0; i < full_list.size(); ++i) {
                    Parameters e = full_list[i].Clone();
                    if (e.Has("type") && e["type"].IsString() && e["type"].GetString() == "element") {
                        if (!e.Has("geometry_type")) e.AddEmptyValue("geometry_type").SetString("GeometrySurface");
                        reduced_sbm.Append(e);
                        continue;
                    }
                    if (!is_cond(e)) continue;
                    const bool has_sbm = e.Has("sbm_parameters");
                    bool has_brep = false;
                    if (e.Has("brep_ids") && e["brep_ids"].IsArray()) {
                        has_brep = (e["brep_ids"].size() > 0);
                    }
                    if (has_sbm || has_brep) {
                        if (!e.Has("geometry_type")) e.AddEmptyValue("geometry_type").SetString("SurfaceEdge");
                        reduced_sbm.Append(e);
                    }
                }
                patch_analysis["element_condition_list"] = reduced_sbm;
            } else {
                // body-fitted patch
                patch_analysis["element_condition_list"] = reduced;
            }
            
        }

        // Run IGA modeler for this patch
        {
            KRATOS_INFO_IF("PatchSubdivisionModeler::  \n", mEchoLevel > 1)
                << patch_analysis << std::endl;
            IgaModelerSbm iga_modeler(*mpModel, patch_analysis);
            iga_modeler.SetupModelPart();
        }
        
    } else {
        // Original behavior: analysis parameters missing is an error at this point
        KRATOS_ERROR << "PatchSubdivisionModeler: missing analysis parameters in input." << std::endl;
    }

}



void PatchSubdivisionModeler::GenerateSubdivision()
{
    mSubdomains.clear();
    mRefinementRegions.clear();

    const auto& r_parameters = mParameters;

    const auto& base_lower = r_parameters["base_domain"]["lower_point_uvw"].GetVector();
    const auto& base_upper = r_parameters["base_domain"]["upper_point_uvw"].GetVector();

    RectangleType base_rect{base_lower[0], base_upper[0], base_lower[1], base_upper[1]};

    if (mEchoLevel > 2) {
        KRATOS_INFO("PatchSubdivisionModeler")
            << "GenerateSubdivision -> base rectangle: " << RectangleToString(base_rect) << std::endl;
    }

    std::set<double> u_breaks{base_rect[U_MIN], base_rect[U_MAX]};
    std::set<double> v_breaks{base_rect[V_MIN], base_rect[V_MAX]};

    const auto& refinement_array = r_parameters["refinement_regions"];
    for (IndexType i = 0; i < refinement_array.size(); ++i) {
        const auto& ref_lower = refinement_array[i]["lower_point_uvw"].GetVector();
        const auto& ref_upper = refinement_array[i]["upper_point_uvw"].GetVector();

        RectangleType region{ref_lower[0], ref_upper[0], ref_lower[1], ref_upper[1]};
        RectangleType clipped = ClipRectangle(region, base_rect);
        if (mEchoLevel > 2) {
            KRATOS_INFO("PatchSubdivisionModeler")
                << "  Refinement region " << i + 1 << " requested: " << RectangleToString(region) << std::endl;
            if (!HasPositiveArea(clipped)) {
                KRATOS_INFO("PatchSubdivisionModeler")
                    << "    -> skipped (no overlap with base domain)." << std::endl;
            } else if (clipped != region) {
                KRATOS_INFO("PatchSubdivisionModeler")
                    << "    -> clipped to: " << RectangleToString(clipped) << std::endl;
            } else {
                KRATOS_INFO("PatchSubdivisionModeler")
                    << "    -> used as provided." << std::endl;
            }
        }
        if (!HasPositiveArea(clipped)) {
            continue;
        }

        u_breaks.insert(clipped[U_MIN]);
        u_breaks.insert(clipped[U_MAX]);
        v_breaks.insert(clipped[V_MIN]);
        v_breaks.insert(clipped[V_MAX]);

        RefinementRegionData region_data;
        region_data.Rectangle = clipped;

        if (refinement_array[i].Has("polynomial_order")) {
            Parameters poly_param = refinement_array[i]["polynomial_order"];
            const std::size_t size = poly_param.size();
            KRATOS_ERROR_IF(size == 0)
                << "PatchSubdivisionModeler: refinement region polynomial_order array is empty." << std::endl;
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
                << "PatchSubdivisionModeler: refinement region number_of_knot_spans array is empty." << std::endl;
            region_data.HasNumberOfKnotSpans = true;
            region_data.NumberOfKnotSpans.resize(size);
            for (IndexType j = 0; j < size; ++j) {
                region_data.NumberOfKnotSpans[j] = span_param.GetArrayItem(j).GetInt();
            }
        }

        mRefinementRegions.push_back(std::move(region_data));
    }

    std::vector<double> u_values(u_breaks.begin(), u_breaks.end());
    std::vector<double> v_values(v_breaks.begin(), v_breaks.end());

    if (mEchoLevel > 2) {
        std::ostringstream u_stream;
        u_stream << "  u-breaks: ";
        for (const double value : u_values) {
            u_stream << value << " ";
        }
        KRATOS_INFO("PatchSubdivisionModeler") << u_stream.str() << std::endl;

        std::ostringstream v_stream;
        v_stream << "  v-breaks: ";
        for (const double value : v_values) {
            v_stream << value << " ";
        }
        KRATOS_INFO("PatchSubdivisionModeler") << v_stream.str() << std::endl;
    }

    for (std::size_t i = 0; i + 1 < u_values.size(); ++i) {
        for (std::size_t j = 0; j + 1 < v_values.size(); ++j) {
            RectangleType rect{
                u_values[i], u_values[i + 1],
                v_values[j], v_values[j + 1]
            };

            if (mEchoLevel > 2) {
                KRATOS_INFO("PatchSubdivisionModeler")
                    << "  Candidate cell: " << RectangleToString(rect) << std::endl;
            }
            if (!HasPositiveArea(rect)) {
                if (mEchoLevel > 2) {
                    KRATOS_INFO("PatchSubdivisionModeler")
                        << "    -> discarded (zero area)." << std::endl;
                }
                continue;
            }

            rect = ClipRectangle(rect, base_rect);
            if (!HasPositiveArea(rect)) {
                if (mEchoLevel > 2) {
                    KRATOS_INFO("PatchSubdivisionModeler")
                        << "    -> discarded after clipping (outside base domain)." << std::endl;
                }
                continue;
            }

            if (mEchoLevel > 2) {
                KRATOS_INFO("PatchSubdivisionModeler")
                    << "    -> accepted: " << RectangleToString(rect) << std::endl;
            }
            mSubdomains.push_back(rect);
        }
    }

    if (mEchoLevel > 2) {
        KRATOS_INFO("PatchSubdivisionModeler")
            << "GenerateSubdivision -> total accepted subdomains: " << mSubdomains.size() << std::endl;
    }
}

bool PatchSubdivisionModeler::IsSkinFullyInsidePatch_(
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
    if (r_skin_mp.NumberOfGeometries() == 0) {
        return false;
    }

    // Sample each curve geometry to build skin bounds in (u,v)
    double skin_min_u = std::numeric_limits<double>::max();
    double skin_max_u = std::numeric_limits<double>::lowest();
    double skin_min_v = std::numeric_limits<double>::max();
    double skin_max_v = std::numeric_limits<double>::lowest();

    const int samples_per_curve = 100; // few points are enough
    for (const auto& rGeom : r_skin_mp.Geometries()) {
        // Assume curve-like geometry with local dimension 1, parameter t in [0,1]
        for (int s = 0; s <= samples_per_curve; ++s) {
            const double t = static_cast<double>(s) / samples_per_curve;

            array_1d<double,3> local;
            local[0] = t; local[1] = 0.0; local[2] = 0.0;
            array_1d<double,3> global;
            rGeom.GlobalCoordinates(global, local); // returns (x,y,*) == (u,v,*)

            const double u = global[0];
            const double v = global[1];

            if (u < skin_min_u) skin_min_u = u;
            if (u > skin_max_u) skin_max_u = u;
            if (v < skin_min_v) skin_min_v = v;
            if (v > skin_max_v) skin_max_v = v;
        }
    }
    // Tolerance
    const double tol = 1e-12;
    const bool inside_u = (skin_min_u >= patch_min_u - tol) && (skin_max_u <= patch_max_u + tol);
    const bool inside_v = (skin_min_v >= patch_min_v - tol) && (skin_max_v <= patch_max_v + tol);

    return inside_u && inside_v;
}

void PatchSubdivisionModeler::BuildGlobalSubModelParts(ModelPart& r_parent_model_part, const std::vector<std::string>& rTargetNames) const
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
                KRATOS_INFO_IF("PatchSubdivisionModeler", mEchoLevel > 2)
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


            KRATOS_INFO_IF("PatchSubdivisionModeler", mEchoLevel > 2)
                << "Collected from '" << rPatch.FullName() << "' into global '"
                << r_global.FullName() << "': elems=" << r_src.NumberOfElements()
                << ", conds=" << r_src.NumberOfConditions()
                << ", nodes~=" << r_src.NumberOfNodes() << "\n";
        }

        KRATOS_INFO_IF("PatchSubdivisionModeler", mEchoLevel > 0)
            << "[Global '" << target_name << "'] SubModelPart '" << r_global.FullName()
            << "' now has " << r_global.NumberOfElements() << " elements, "
            << r_global.NumberOfConditions() << " conditions and "
            << r_global.NumberOfNodes() << " nodes (inserted ~"
            << added_elems << " elems, ~" << added_conds << " conds, ~"
            << added_nodes << " nodes).\n";

        // // Now detach entities from each Patch.* target submodelpart (membership only).
        // for (auto& rPatch : r_root.SubModelParts()) {
        //     const std::string& patch_name = rPatch.Name();
        //     if (patch_name.rfind("Patch", 0) != 0) continue;
        //     if (!rPatch.HasSubModelPart(target_name)) continue;

        //     ModelPart& r_src = rPatch.GetSubModelPart(target_name);

        //     // Drop memberships (does NOT delete from the root model part).
        //     r_src.Elements().clear();
        //     r_src.Conditions().clear();
        //     r_src.Nodes().clear();

        //     // Also clear any geometry memberships (if any)
        //     r_src.Geometries().clear();
        //     // Remove the now-empty child submodel part from this Patch
        //     rPatch.RemoveSubModelPart(target_name);
        // }
    }



    // // ---- Remove all Patch* submodel parts from the root ----
    // std::vector<std::string> patches_to_delete;
    // // Phase 1: clear contents (and children) of each Patch*
    // for (auto& rPatch : r_root.SubModelParts()) {
    //     const std::string& patch_name = rPatch.Name();
    //     if (patch_name.rfind("Patch", 0) != 0) continue; // only Patch*
    //     // Clear child submodel parts first
    //     std::vector<std::string> child_names;
    //     child_names.reserve(rPatch.NumberOfSubModelParts());
    //     for (auto& rChild : rPatch.SubModelParts()) {
    //         rChild.Elements().clear();
    //         rChild.Conditions().clear();
    //         rChild.Nodes().clear();
    //         rChild.Geometries().clear();
    //         child_names.push_back(rChild.Name());
    //     }
    //     for (const auto& child_name : child_names) {
    //         rPatch.RemoveSubModelPart(child_name);
    //     }
    //     // Clear the Patch itself
    //     rPatch.Elements().clear();
    //     rPatch.Conditions().clear();
    //     rPatch.Nodes().clear();
    //     rPatch.Geometries().clear();
    //     patches_to_delete.push_back(patch_name);
    // }



    // // Print info for all nodes in the root model part
    // int count = 0;
    // int count_dof = 0;
    // for (auto& r_node : r_root.Nodes()) {
    //     KRATOS_WATCH(r_node.Id())
    //     KRATOS_WATCH(r_node.Coordinates())
    //     KRATOS_WATCH(r_node)
    //     count++;
    //     // Print DoFs info
    //     for (auto it_dof = r_node.GetDofs().begin(); it_dof != r_node.GetDofs().end(); ++it_dof) {
    //         KRATOS_WATCH((*it_dof)->GetVariable().Name())
    //         KRATOS_WATCH((*it_dof)->GetSolutionStepValue())
    //         count_dof++;
    //     }
    // }
    // KRATOS_INFO_IF("PatchSubdivisionModeler", mEchoLevel > 0)
    //     << "Root model part '" << r_root.FullName() << "' has "
    //     << r_root.NumberOfNodes() << " nodes (printed " << count << "), "
    //     << count_dof << " DoFs." << std::endl;
    // exit(0);





    // // Phase 2: actually remove Patch* submodel parts from the root
    // for (const auto& name : patches_to_delete) {
    //     if (r_root.HasSubModelPart(name)) {
    //         r_root.RemoveSubModelPart(name);
    //     }
    // }

    // KRATOS_WATCH(r_parent_model_part)
    // exit(0);

}

} // namespace Kratos
