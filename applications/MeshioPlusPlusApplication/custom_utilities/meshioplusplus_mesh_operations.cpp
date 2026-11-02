//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                        Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <algorithm>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <set>

// External includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/skin.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/boundary_layer.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"
#include "meshioplusplus/operations/error_size.hpp"
#include "meshioplusplus/operations/extrude.hpp"
#include "meshioplusplus/operations/find_intersections.hpp"
#include "meshioplusplus/operations/mirror.hpp"
#include "meshioplusplus/operations/mmg.hpp"
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/refine_patches.hpp"
#include "meshioplusplus/operations/select.hpp"
#include "meshioplusplus/operations/tet_poly.hpp"
#include "meshioplusplus/operations/tetrahedralize.hpp"
#include "meshioplusplus/operations/thickness.hpp"
#include "meshioplusplus/operations/blend.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/operations/feature_edges.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/periodic.hpp"
#include "meshioplusplus/operations/quality_gate.hpp"
#include "meshioplusplus/operations/region_ops.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/conservative_interpolate.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/crop.hpp"
#include "meshioplusplus/operations/curvature.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_calc.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/data_condition.hpp"
#include "meshioplusplus/operations/data_info.hpp"
#include "meshioplusplus/operations/data_integrate.hpp"
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/decimate_volume.hpp"
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/operations/error.hpp"
#include "meshioplusplus/operations/gradient.hpp"
#include "meshioplusplus/operations/hessian.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/neighbors.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/remesh.hpp"
#include "meshioplusplus/operations/remesh_volume.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/repair.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/operations/shrinkwrap.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/sobolev_deform.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/subdivide.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/tensor_invariants.hpp"
#include "meshioplusplus/operations/tessellate.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

// Project includes
#include "custom_utilities/meshioplusplus_mesh_operations.h"
#include "custom_utilities/meshioplusplus_conversion_utilities.h"
#include "includes/global_variables.h"
#include "containers/model.h"

namespace mio = meshioplusplus;

namespace Kratos
{
namespace
{

const std::vector<std::string>& OperationNames()
{
    static const std::vector<std::string> names = [] {
        std::vector<std::string> sorted_names = {
            "agglomerate", "attach_quality", "cell_data_to_point_data", "clean", "compute_normals",
            "compute_sdf", "convert_cells", "crop_bbox", "crop_halfspace", "crop_predicate",
            "curvature", "data_calc", "data_condition", "data_info", "data_integrate", "data_manage",
            "decimate", "decimate_volume", "edit_regions", "estimate_error", "extract_skin",
            "extract_surface", "feature_edges", "fit_cad", "fit_primitives", "gradient", "hessian", "isosurface", "optimize_volume", "partition",
            "point_data_to_cell_data", "quality", "check_quality", "refine", "remesh", "remesh_volume", "reorder",
            "repair", "segment", "slice", "smooth", "sobolev_deform", "split", "stats", "subdivide",
            "tensor_invariants", "tessellate", "transform", "voxelize",
            "adapt", "boundary_layer", "check_conformity", "discretize_levelset", "error_size",
            "extrude", "find_intersections", "mesh_domain", "metric_from_hessian", "mirror", "move",
            "octree", "octree_mesh", "refine_patches", "revolve", "select", "tet_to_poly",
            "tetrahedralize", "thickness"
        };
        std::sort(sorted_names.begin(), sorted_names.end());
        return sorted_names;
    }();
    return names;
}

/***********************************************************************************/
/***********************************************************************************/

/**
 * @brief The meshio++ primitive-type bitmask of a "fit_types" list.
 * @details meshio++ has no public name-to-type helper for fit_primitives, so the five names its
 * own command line and pipeline accept are mapped here; the bit of a type is 1 << its enum value.
 */
std::uint32_t ReadFitTypes(Parameters Settings, const std::string& rKey)
{
    KRATOS_ERROR_IF_NOT(Settings.Has(rKey)) << "Missing \"" << rKey << "\" setting" << std::endl;
    KRATOS_ERROR_IF_NOT(Settings[rKey].IsArray() && Settings[rKey].size() > 0)
        << "\"" << rKey << "\" must be a non-empty list of plane, sphere, cylinder, cone or torus" << std::endl;
    std::uint32_t mask = 0;
    for (std::size_t i = 0; i < Settings[rKey].size(); ++i) {
        const std::string name = Settings[rKey][i].GetString();
        if (name == "plane") {
            mask |= 1u << static_cast<std::uint32_t>(mio::FitPrimitiveType::Plane);
        } else if (name == "sphere") {
            mask |= 1u << static_cast<std::uint32_t>(mio::FitPrimitiveType::Sphere);
        } else if (name == "cylinder") {
            mask |= 1u << static_cast<std::uint32_t>(mio::FitPrimitiveType::Cylinder);
        } else if (name == "cone") {
            mask |= 1u << static_cast<std::uint32_t>(mio::FitPrimitiveType::Cone);
        } else if (name == "torus") {
            mask |= 1u << static_cast<std::uint32_t>(mio::FitPrimitiveType::Torus);
        } else {
            KRATOS_ERROR << "Unknown primitive type \"" << name << "\" in \"" << rKey
                         << "\"; expected plane, sphere, cylinder, cone or torus" << std::endl;
        }
    }
    return mask;
}

/***********************************************************************************/
/***********************************************************************************/

/**
 * @brief The enum value a setting names, among the (name, value) pairs it accepts.
 * @details For the meshio++ enums that have no public `*_from_name` parser of their own; the
 * names are the ones meshio++'s command line and pipeline accept.
 */
template <class TEnum>
TEnum ReadEnum(
    const Parameters& rSettings,
    const std::string& rKey,
    std::initializer_list<std::pair<const char*, TEnum>> Values
    )
{
    const std::string name = rSettings[rKey].GetString();
    std::string accepted;
    for (const auto& r_value : Values) {
        if (name == r_value.first) {
            return r_value.second;
        }
        accepted += (accepted.empty() ? "" : ", ") + std::string(r_value.first);
    }
    KRATOS_ERROR << "Unknown \"" << rKey << "\" value \"" << name << "\"; expected one of "
                 << accepted << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

/// A meshio++ point as a Kratos Vector, for a report.
Vector ToVector(const std::array<double, 3>& rPoint)
{
    Vector result(3);
    for (std::size_t i = 0; i < 3; ++i) {
        result[i] = rPoint[i];
    }
    return result;
}

/***********************************************************************************/
/***********************************************************************************/

/// Adds a double to a report only when it is finite: meshio++ reports "not measured" as NaN,
/// which a JSON document cannot hold.
void AddFiniteDouble(Parameters& rReport, const std::string& rKey, const double Value)
{
    if (std::isfinite(Value)) {
        rReport.AddDouble(rKey, Value);
    }
}

/***********************************************************************************/
/***********************************************************************************/

/// The plane of a 2-D skin for "voxelize" and "octree" ("plane": auto picks the axis of zero
/// extent, z first); ignored by their 3-D variants.
mio::detail::PlanarPlane ReadPlanarPlane(const Parameters& rSettings)
{
    return ReadEnum<mio::detail::PlanarPlane>(rSettings, "plane",
        {{"auto", mio::detail::PlanarPlane::Auto}, {"xy", mio::detail::PlanarPlane::XY},
         {"xz", mio::detail::PlanarPlane::XZ}, {"yz", mio::detail::PlanarPlane::YZ}});
}

/***********************************************************************************/
/***********************************************************************************/

/**
 * @brief The segment options of a "segment" or "fit_primitives" operation.
 * @details "feature_angle" and "region" are the settings feature_edges and refine already use.
 */
mio::SegmentOptions ReadSegmentOptions(Parameters Settings)
{
    mio::SegmentOptions options;
    options.mFeatureAngleDeg = Settings["feature_angle"].GetDouble();
    options.mNormalToleranceDeg = Settings["normal_tolerance"].GetDouble();
    options.mCurvatureTolerance = Settings["curvature_tolerance"].GetDouble();
    options.mMinPatchCells = Settings["min_patch_cells"].GetInt();
    options.mRegion = Settings["region"].GetString();
    options.mPrefix = Settings["patch_prefix"].GetString();
    options.mMethod = ReadEnum<mio::SegmentMethod>(Settings, "segment_method",
        {{"angle", mio::SegmentMethod::Angle}, {"primitive", mio::SegmentMethod::Primitive}});
    options.mTolerance = Settings["segment_tolerance"].GetDouble();
    return options;
}

/***********************************************************************************/
/***********************************************************************************/

/**
 * @brief The options of "refine_patches", also embedded in "fit_cad"'s
 * "cad_refine_segmentation".
 * @details "feature_angle", "fit_types", "fit_relative" and "patch_prefix" are the settings
 * segment and fit_primitives already use; a negative "refine_patches_seed" keeps meshio++'s own.
 */
mio::RefinePatchesOptions ReadRefinePatchesOptions(Parameters Settings)
{
    mio::RefinePatchesOptions options;
    options.mTolerance = Settings["refine_patches_tolerance"].GetDouble();
    options.mNormalWeight = Settings["refine_patches_normal_weight"].GetDouble();
    options.mEpsilon = Settings["refine_patches_epsilon"].GetDouble();
    options.mFeatureAngleDeg = Settings["feature_angle"].GetDouble();
    options.mMaxIterations = Settings["refine_patches_max_iterations"].GetInt();
    options.mMinPatchCells = Settings["refine_patches_min_patch_cells"].GetInt();
    if (Settings["refine_patches_seed"].GetInt() >= 0) {
        options.mRansacSeed = static_cast<std::uint64_t>(Settings["refine_patches_seed"].GetInt());
    }
    options.mTypes = ReadFitTypes(Settings, "fit_types");
    options.mRelative = Settings["fit_relative"].GetBool();
    options.mAllowSplit = Settings["refine_patches_allow_split"].GetBool();
    options.mAllowMerge = Settings["refine_patches_allow_merge"].GetBool();
    options.mPrefix = Settings["patch_prefix"].GetString();
    return options;
}

/***********************************************************************************/
/***********************************************************************************/

std::array<double, 3> ReadVector3(Parameters Settings, const std::string& rKey)
{
    std::array<double, 3> value{{0.0, 0.0, 0.0}};
    KRATOS_ERROR_IF_NOT(Settings.Has(rKey)) << "Missing \"" << rKey << "\" setting" << std::endl;
    const Vector read = Settings[rKey].GetVector();
    for (std::size_t i = 0; i < std::min<std::size_t>(3, read.size()); ++i) {
        value[i] = read[i];
    }
    return value;
}

/***********************************************************************************/
/***********************************************************************************/

/// Builds the shared @ref Internals::FieldDataSelection from the settings @ref Execute /
/// @ref Interpolate / @ref MeshioPlusPlusMeshOperations::Merge have in common with
/// @ref MeshioPlusPlusIO - empty/false by default, so an operation run without them stages
/// nothing (byte-identical to the pre-field-data behavior).
Internals::FieldDataSelection BuildFieldDataSelection(const Parameters& rSettings)
{
    Internals::FieldDataSelection selection;
    selection.NodalSolutionStepVariables = rSettings["nodal_solution_step_data_variables"].GetStringArray();
    selection.NodalDataValueVariables = rSettings["nodal_data_value_variables"].GetStringArray();
    selection.NodalFlags = rSettings["nodal_flags"].GetStringArray();
    selection.ElementDataValueVariables = rSettings["element_data_value_variables"].GetStringArray();
    selection.ElementFlags = rSettings["element_flags"].GetStringArray();
    selection.ConditionDataValueVariables = rSettings["condition_data_value_variables"].GetStringArray();
    selection.ConditionFlags = rSettings["condition_flags"].GetStringArray();
    selection.GaussPointVariables = rSettings["gauss_point_variables_in_elements"].GetStringArray();
    selection.WriteIds = rSettings["write_ids"].GetBool();
    selection.WriteMdpaIds = rSettings["write_mdpa_ids"].GetBool();
    return selection;
}

/***********************************************************************************/
/***********************************************************************************/

/// Reads the surface-distance settings meshio++ embeds by value in both VoxelOptions and
/// SdfOptions, so "voxelize" and "compute_sdf" cannot drift on the shared half.
mio::SurfaceDistanceOptions BuildSurfaceDistanceOptions(const Parameters& rSettings)
{
    mio::SurfaceDistanceOptions options;
    options.mSign = mio::sdf_sign_from_name(rSettings["sdf_sign"].GetString());
    options.mWeight = mio::sdf_weight_from_name(rSettings["sdf_weight"].GetString());
    options.mLocation = mio::sdf_location_from_name(rSettings["sdf_location"].GetString());
    options.mBand = rSettings["band"].GetDouble();
    options.mRecordClosestCell = rSettings["record_closest_cell"].GetBool();
    options.mRecordInside = rSettings["record_inside"].GetBool();
    options.mWatertightCheck = mio::sdf_watertight_check_from_name(rSettings["watertight_check"].GetString());
    options.mGridCellSize = rSettings["grid_cell_size"].GetDouble();
    options.mMaxWindingWork = rSettings["max_winding_work"].GetDouble();
    return options;
}

/***********************************************************************************/
/***********************************************************************************/

/// The "bounds" setting (xlo, ylo, zlo, xhi, yhi, zhi): an empty array leaves it unset, the
/// input's bounding box then being used.
std::optional<std::array<double, 6>> ReadBounds(const Parameters& rSettings)
{
    const Vector bounds = rSettings["bounds"].GetVector();
    if (bounds.size() == 0) {
        return std::nullopt;
    }
    KRATOS_ERROR_IF(bounds.size() != 6)
        << "\"bounds\" needs six entries (xlo, ylo, zlo, xhi, yhi, zhi), got "
        << bounds.size() << std::endl;
    std::array<double, 6> value{};
    for (std::size_t i = 0; i < 6; ++i) {
        value[i] = bounds[i];
    }
    return value;
}

/***********************************************************************************/
/***********************************************************************************/

/// Reads the lattice settings "voxelize" and "compute_sdf" share. Both spell "unset" as an
/// empty array / a non-positive scalar, since meshio++ takes std::optional there and setting
/// neither (or both) of resolution/cell size is its own named error.
///
/// One caveat the shared key buys: upstream's own default for "padding_relative" is not the same
/// for every operation taking it (voxelize 0.0, compute_sdf and remesh_volume 0.1). A single
/// Kratos setting can only carry one default, and it carries voxelize's 0.0 - so compute_sdf and
/// remesh_volume are padded less here than a direct meshio++ call would pad them unless the
/// setting is given explicitly. Kept at 0.0 rather than "corrected" because changing it would
/// silently change what every existing voxelize call produces.
template <class TOptions>
void ReadLatticeSettings(const Parameters& rSettings, TOptions& rOptions)
{
    const Vector resolution = rSettings["resolution"].GetVector();
    if (resolution.size() > 0) {
        KRATOS_ERROR_IF(resolution.size() != 3)
            << "\"resolution\" needs three entries, got " << resolution.size() << std::endl;
        rOptions.mResolution = std::array<std::int64_t, 3>{{static_cast<std::int64_t>(resolution[0]),
                                                            static_cast<std::int64_t>(resolution[1]),
                                                            static_cast<std::int64_t>(resolution[2])}};
    }
    const double cell_size = rSettings["cell_size"].GetDouble();
    if (cell_size > 0.0) {
        rOptions.mCellSize = cell_size;
    }
    rOptions.mBounds = ReadBounds(rSettings);
    rOptions.mPadding = rSettings["padding"].GetDouble();
    rOptions.mPaddingRelative = rSettings["padding_relative"].GetDouble();
    rOptions.mMaxCells = rSettings["max_cells"].GetInt();
    rOptions.mDistance = BuildSurfaceDistanceOptions(rSettings);
}

/***********************************************************************************/
/***********************************************************************************/

/// Renames a result array produced under one of meshio++'s own colon-namespaced names to the
/// caller's "output". Without this the sdf/voxelize results are unreachable from Kratos: a
/// Kratos entity can only hold data under a registered Variable, and "sdf:distance" never is.
void RenameResultArray(
    mio::Mesh& rMesh,
    const mio::DataLocation Location,
    const std::string& rFrom,
    const std::string& rTo
    )
{
    if (rTo.empty() || rTo == rFrom || !mio::data_has(rMesh, Location, rFrom)) {
        return;
    }
    rMesh = mio::data_rename(rMesh, Location, rFrom, rTo);
}

/***********************************************************************************/
/***********************************************************************************/

/// Appends a lattice description (dims / origin / spacing) to an operation's report.
void AppendLatticeReport(
    Parameters& rReport,
    const std::array<std::int64_t, 3>& rDims,
    const std::array<double, 3>& rOrigin,
    const std::array<double, 3>& rSpacing
    )
{
    Parameters dims(R"([])");
    Parameters origin(R"([])");
    Parameters spacing(R"([])");
    for (std::size_t i = 0; i < 3; ++i) {
        dims.Append(Parameters(std::to_string(rDims[i])));
        origin.Append(Parameters(std::to_string(rOrigin[i])));
        spacing.Append(Parameters(std::to_string(rSpacing[i])));
    }
    rReport.AddValue("dims", dims);
    rReport.AddValue("origin", origin);
    rReport.AddValue("spacing", spacing);
}

/***********************************************************************************/
/***********************************************************************************/

/// Flattens a meshio++ surface verdict into a report.
Parameters SurfaceQualityReport(const mio::SurfaceQuality& rQuality)
{
    Parameters report(R"({})");
    report.AddBool("watertight", rQuality.mWatertight);
    report.AddInt("boundary_edges", static_cast<int>(rQuality.mBoundaryEdges));
    report.AddInt("non_manifold_edges", static_cast<int>(rQuality.mNonManifoldEdges));
    report.AddInt("inconsistent_pairs", static_cast<int>(rQuality.mInconsistentPairs));
    report.AddInt("degenerate_triangles", static_cast<int>(rQuality.mDegenerateTriangles));
    return report;
}

/***********************************************************************************/
/***********************************************************************************/

/// Reads a "location" setting ("point"/"cell") into the meshio++ enum.
mio::DataLocation ReadLocation(const Parameters& rSettings, const std::string& rKey = "location")
{
    return mio::data_location_from_name(rSettings[rKey].GetString());
}

/// Reads a JSON array of {"location", "name"} objects into meshio++ DataKeys (data_manage's
/// "drop"/"keep" settings).
std::vector<mio::DataKey> ReadDataKeys(Parameters Settings)
{
    std::vector<mio::DataKey> keys;
    keys.reserve(Settings.size());
    for (std::size_t i = 0; i < Settings.size(); ++i) {
        mio::DataKey key;
        key.location = ReadLocation(Settings[i]);
        key.name = Settings[i]["name"].GetString();
        keys.push_back(key);
    }
    return keys;
}

/// Reads a JSON array of {"location", "from", "to"} objects into meshio++ DataRenames
/// (data_manage's "rename" setting).
std::vector<mio::DataRename> ReadDataRenames(Parameters Settings)
{
    std::vector<mio::DataRename> renames;
    renames.reserve(Settings.size());
    for (std::size_t i = 0; i < Settings.size(); ++i) {
        mio::DataRename rename;
        rename.location = ReadLocation(Settings[i]);
        rename.from = Settings[i]["from"].GetString();
        rename.to = Settings[i]["to"].GetString();
        renames.push_back(rename);
    }
    return renames;
}

/// Reads a JSON array of integers into meshio++'s block-major global cell indices (refine's
/// "cells" selector). Parameters has no int64 accessor, so this widens the same way GetInt()
/// itself narrows the underlying JSON number - as everywhere else in this file (mSeed,
/// mGhostLayers, ... are all read the same way).
std::vector<std::int64_t> ReadCellIndices(Parameters Settings)
{
    std::vector<std::int64_t> cells;
    cells.reserve(Settings.size());
    for (std::size_t i = 0; i < Settings.size(); ++i) {
        cells.push_back(static_cast<std::int64_t>(Settings[i].GetInt()));
    }
    return cells;
}

/***********************************************************************************/
/***********************************************************************************/

/// Reads a pin mask ("frozen") into meshio++'s per-point flag vector. Empty means "no extra
/// pins", which is what every operation taking one spells as an empty vector - so an unset
/// setting reaches meshio++ unchanged rather than as a zero-filled mask of the wrong length.
std::vector<std::uint8_t> ReadFrozenFlags(Parameters Settings)
{
    std::vector<std::uint8_t> frozen;
    frozen.reserve(Settings.size());
    for (std::size_t i = 0; i < Settings.size(); ++i) {
        frozen.push_back(Settings[i].GetBool() ? std::uint8_t(1) : std::uint8_t(0));
    }
    return frozen;
}

/***********************************************************************************/
/***********************************************************************************/

/// The criteria of "select"; each "select_*" switch turns one on, and the ones on intersect.
mio::SelectOptions ReadSelectOptions(Parameters Settings)
{
    mio::SelectOptions options;
    options.mName = Settings["select_name"].GetString();
    KRATOS_ERROR_IF(options.mName.empty())
        << "The \"select\" operation needs a \"select_name\" (the sub model part it writes)" << std::endl;
    options.mKind = mio::region_kind_from_name(Settings["select_kind"].GetString());
    options.mTag = Settings["select_tag"].GetInt();
    options.mMode = mio::select_mode_from_name(Settings["select_mode"].GetString());
    options.mBoundaryOnly = Settings["select_boundary_only"].GetBool();
    options.mBox = Settings["select_box"].GetBool();
    options.mBoxLo = ReadVector3(Settings, "box_min");
    options.mBoxHi = ReadVector3(Settings, "box_max");
    options.mPlane = Settings["select_plane"].GetBool();
    options.mPlanePoint = ReadVector3(Settings, "origin");
    options.mPlaneNormal = ReadVector3(Settings, "normal");
    options.mSphere = Settings["select_sphere"].GetBool();
    options.mSphereCenter = ReadVector3(Settings, "select_sphere_center");
    options.mSphereRadius = Settings["select_sphere_radius"].GetDouble();
    options.mWhereArray = Settings["select_where_array"].GetString();
    options.mWhereOp = mio::refine_compare_from_name(Settings["select_where_op"].GetString());
    options.mWhereValue = Settings["select_where_value"].GetDouble();
    options.mQualityMetric = Settings["select_quality_metric"].GetString();
    options.mQualityOp = mio::refine_compare_from_name(Settings["select_quality_op"].GetString());
    options.mQualityValue = Settings["select_quality_value"].GetDouble();
    options.mComponent = Settings["select_component"].GetInt();
    options.mSeed = Settings["select_seed"].GetBool();
    options.mSeedPoint = ReadVector3(Settings, "select_seed_point");
    options.mFeatureAngleDeg = Settings["feature_angle"].GetDouble();
    options.mNormalToleranceDeg = Settings["normal_tolerance"].GetDouble();
    return options;
}

/***********************************************************************************/
/***********************************************************************************/

/// The options of "find_intersections"; "region" restricts the first surface (a sub model part).
mio::FindIntersectionsOptions ReadFindIntersectionsOptions(Parameters Settings)
{
    mio::FindIntersectionsOptions options;
    options.mRegion = Settings["region"].GetString();
    options.mRegionB = Settings["intersections_region_b"].GetString();
    options.mAttachFlags = false;
    options.mSegments = true;
    options.mMaxPairs = Settings["intersections_max_pairs"].GetInt();
    return options;
}

/***********************************************************************************/
/***********************************************************************************/

/// The report of "find_intersections": the counts, and the intersecting cell pairs as
/// container positions of the cells written to meshio++ (elements, then conditions).
Parameters FindIntersectionsReport(const mio::FindIntersectionsResult& rResult)
{
    Parameters report(R"({})");
    report.AddBool("clean", rResult.mClean);
    report.AddInt("number_of_pairs", static_cast<int>(rResult.mNumPairs));
    report.AddInt("number_of_triangle_pairs", static_cast<int>(rResult.mNumTrianglePairs));
    report.AddInt("number_of_coplanar", static_cast<int>(rResult.mNumCoplanar));
    report.AddInt("number_of_segments", static_cast<int>(rResult.mNumSegments));
    report.AddInt("number_of_cells_first", static_cast<int>(rResult.mNumCellsA));
    report.AddInt("number_of_cells_second", static_cast<int>(rResult.mNumCellsB));
    report.AddBool("truncated", rResult.mTruncated);
    Parameters pairs(R"([])");
    for (const auto& r_pair : rResult.mPairs) {
        Parameters entry(R"([])");
        entry.Append(static_cast<int>(r_pair[0]));
        entry.Append(static_cast<int>(r_pair[1]));
        pairs.Append(entry);
    }
    report.AddValue("pairs", pairs);
    report.AddStringArray("warnings", rResult.mWarnings);
    return report;
}

/***********************************************************************************/
/***********************************************************************************/

/// The options of "thickness"; "region" restricts the measured surface (a sub model part).
mio::ThicknessOptions ReadThicknessOptions(Parameters Settings)
{
    mio::ThicknessOptions options;
    options.mDirection = ReadEnum<mio::ThicknessDirection>(Settings, "thickness_direction",
        {{"auto", mio::ThicknessDirection::Auto}, {"inward", mio::ThicknessDirection::Inward},
         {"outward", mio::ThicknessDirection::Outward}});
    options.mMaxDistance = Settings["thickness_max_distance"].GetDouble();
    options.mThreshold = Settings["thickness_threshold"].GetDouble();
    options.mRegion = Settings["region"].GetString();
    options.mRegionTarget = Settings["thickness_region_target"].GetString();
    options.mAttach = true;
    return options;
}

/***********************************************************************************/
/***********************************************************************************/

/// The report of "thickness": the least, largest and mean distance over the rays that hit.
Parameters ThicknessReport(const mio::ThicknessResult& rResult)
{
    Parameters report(R"({})");
    report.AddString("direction", rResult.mDirection == mio::ThicknessDirection::Outward ? "outward" : "inward");
    report.AddInt("number_of_cells", static_cast<int>(rResult.mNumCells));
    report.AddInt("number_of_hits", static_cast<int>(rResult.mNumHit));
    report.AddInt("number_of_misses", static_cast<int>(rResult.mNumMissed));
    report.AddInt("number_below_threshold", static_cast<int>(rResult.mNumBelow));
    report.AddDouble("min", rResult.mMin);
    report.AddDouble("max", rResult.mMax);
    report.AddDouble("mean", rResult.mMean);
    report.AddInt("min_cell", static_cast<int>(rResult.mMinCell));
    report.AddVector("min_point", ToVector(rResult.mMinPoint));
    report.AddDouble("max_distance", rResult.mMaxDistance);
    report.AddStringArray("warnings", rResult.mWarnings);
    return report;
}

/***********************************************************************************/
/***********************************************************************************/

/// The MMG parameters every MMG operation shares ("mmg_*"; -1 keeps MMG's own default).
/// "frozen" lists the points MMG must keep, "mmg_required" the regions, and "mmg_local" holds
/// {"region", "hmin", "hmax", "hausd"} per-region overrides.
mio::MmgParameters ReadMmgParameters(Parameters Settings)
{
    mio::MmgParameters parameters;
    parameters.mHmin = Settings["mmg_hmin"].GetDouble();
    parameters.mHmax = Settings["mmg_hmax"].GetDouble();
    parameters.mHausd = Settings["mmg_hausd"].GetDouble();
    parameters.mHgrad = Settings["mmg_hgrad"].GetDouble();
    parameters.mHgradreq = Settings["mmg_hgradreq"].GetDouble();
    parameters.mAngleDetection = Settings["mmg_angle_detection"].GetDouble();
    parameters.mOptim = Settings["mmg_optim"].GetBool();
    parameters.mNoInsert = Settings["mmg_noinsert"].GetBool();
    parameters.mNoSwap = Settings["mmg_noswap"].GetBool();
    parameters.mNoMove = Settings["mmg_nomove"].GetBool();
    parameters.mNoSurf = Settings["mmg_nosurf"].GetBool();
    parameters.mNreg = Settings["mmg_nreg"].GetBool();
    parameters.mOpnbdy = Settings["mmg_opnbdy"].GetBool();
    parameters.mMem = Settings["mmg_mem"].GetInt();
    parameters.mVerbose = Settings["mmg_verbose"].GetInt();
    parameters.mFrozen = ReadFrozenFlags(Settings["frozen"]);
    parameters.mRequired = Settings["mmg_required"].GetStringArray();
    for (std::size_t i = 0; i < Settings["mmg_local"].size(); ++i) {
        Parameters local = Settings["mmg_local"][i];
        local.AddMissingParameters(Parameters(R"({"region": "", "hmin": -1.0, "hmax": -1.0, "hausd": -1.0})"));
        mio::MmgLocalParameter entry;
        entry.mRegion = local["region"].GetString();
        entry.mHmin = local["hmin"].GetDouble();
        entry.mHmax = local["hmax"].GetDouble();
        entry.mHausd = local["hausd"].GetDouble();
        parameters.mLocal.push_back(entry);
    }
    return parameters;
}

/***********************************************************************************/
/***********************************************************************************/

/// Flattens one "data_integrate" region reduction. The per-component vectors stay arrays even
/// for a scalar array, so a caller reading the report does not need a different shape per
/// component count.
Parameters FieldIntegralRegionReport(const mio::FieldIntegralRegion& rRegion)
{
    Parameters entry(R"({})");
    entry.AddString("name", rRegion.mName);
    entry.AddInt("number_of_cells", static_cast<int>(rRegion.mNumCells));
    entry.AddInt("number_of_skipped", static_cast<int>(rRegion.mNumSkipped));
    const auto to_vector = [](const std::vector<double>& rValues) {
        Vector vector(rValues.size());
        std::copy(rValues.begin(), rValues.end(), vector.begin());
        return vector;
    };
    entry.AddVector("domain_measure", to_vector(rRegion.mDomainMeasurePerComponent));
    entry.AddVector("total", to_vector(rRegion.mTotalPerComponent));
    entry.AddVector("mean", to_vector(rRegion.mMeanPerComponent));
    Parameters number_of_nan(R"([])");
    for (const std::int64_t count : rRegion.mNumNanPerComponent) {
        number_of_nan.Append(Parameters(std::to_string(count)));
    }
    entry.AddValue("number_of_nan", number_of_nan);
    return entry;
}

/***********************************************************************************/
/***********************************************************************************/

/// Stores the result of a *polyhedral* operation ("subdivide", "agglomerate").
///
/// Both emit `polyhedronN` cells by construction - one child per face, one cell per group - and
/// a Kratos ModelPart cannot hold one: an Element's geometry is a registered type with a fixed
/// node count, and there is no polyhedral geometry to map onto. Written straight through, the
/// destination comes back silently *empty*, which is the worst of the available outcomes.
///
/// So "simplexify_result" (default true) runs meshio++'s own remedy - the
/// `convert_cells(Simplexify)` that `refine` and `decimate` already name when they refuse a
/// polyhedron - decomposing the result into tetrahedra before the write-back. The tetrahedra
/// are what Kratos can actually hold, so the operation produces a usable model part in one
/// call rather than requiring a second pass that no ModelPart could carry the input for.
///
/// Setting it false hands back the raw polyhedral mesh, which is only useful when the result is
/// on its way to a meshio++ *file* rather than to a model part - so an empty destination is
/// reported by name instead of silently returned.
void StorePolyhedralResult(
    mio::Mesh& rMesh,
    ModelPart& rDestination,
    const std::string& rOperation,
    const bool Simplexify
    )
{
    if (Simplexify) {
        mio::ConvertCellsOptions options;
        options.mMode = mio::ConvertCellsMode::Simplexify;
        mio::ConvertCellsResult converted = mio::convert_cells(rMesh, options);
        Internals::MeshToModelPart(converted.mMesh, rDestination);
    } else {
        Internals::MeshToModelPart(rMesh, rDestination);
    }

    KRATOS_ERROR_IF(rDestination.NumberOfElements() == 0 && rDestination.NumberOfConditions() == 0)
        << "The \"" << rOperation << "\" operation produced polyhedral cells, which a Kratos "
        << "ModelPart cannot hold - an Element's geometry is a registered type with a fixed node "
        << "count. Set \"simplexify_result\" : true (the default) to decompose the result into "
        << "tetrahedra first." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void StoreResult(mio::Mesh& rMesh, ModelPart& rDestination)
{
    Internals::MeshToModelPart(rMesh, rDestination);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPieces>
Parameters StorePieces(TPieces& rPieces, ModelPart& rDestination, const std::string& rPrefix)
{
    Parameters report(R"({"pieces": []})");
    for (std::size_t i = 0; i < rPieces.size(); ++i) {
        const std::string name = rDestination.Name() + "_" + rPrefix + "_" + std::to_string(i);
        ModelPart& r_piece = rDestination.GetModel().HasModelPart(name)
            ? rDestination.GetModel().GetModelPart(name)
            : rDestination.GetModel().CreateModelPart(name);
        Internals::MeshToModelPart(rPieces[i].mMesh, r_piece);

        Parameters entry(R"({})");
        entry.AddString("name", name);
        entry.AddInt("number_of_nodes", static_cast<int>(r_piece.NumberOfNodes()));
        entry.AddInt("number_of_elements", static_cast<int>(r_piece.NumberOfElements()));
        entry.AddInt("number_of_conditions", static_cast<int>(r_piece.NumberOfConditions()));
        report["pieces"].Append(entry);
    }
    report.AddInt("number_of_pieces", static_cast<int>(rPieces.size()));
    return report;
}

/// Reads a {"name", "kind", "dim", "tag"} region selector; "kind" is "any" (default), "point",
/// "cell" or "side", and an absent "dim"/"tag" matches any.
mio::RegionSelector ReadRegionSelector(Parameters Selector, const std::string& rWhere)
{
    KRATOS_ERROR_IF_NOT(Selector.Has("name") && Selector["name"].IsString())
        << rWhere << ": a region selector needs a \"name\" string" << std::endl;
    mio::RegionSelector selector;
    selector.mName = Selector["name"].GetString();
    const std::string kind = Selector.Has("kind") ? Selector["kind"].GetString() : "any";
    if (kind == "any") {
        selector.mKind = -1;
    } else if (kind == "point") {
        selector.mKind = 0;
    } else if (kind == "cell") {
        selector.mKind = 1;
    } else if (kind == "side") {
        selector.mKind = 2;
    } else {
        KRATOS_ERROR << rWhere << ": unknown region \"kind\" \"" << kind
                     << "\" (use \"any\", \"point\", \"cell\" or \"side\")" << std::endl;
    }
    if (Selector.Has("dim")) {
        selector.mDim = Selector["dim"].GetInt();
    }
    if (Selector.Has("tag")) {
        selector.mTag = Selector["tag"].GetInt();
    }
    return selector;
}

/// Builds the affine transform of the "translation", "rotation_axis"/"rotation_angle" (radians)
/// or "matrix" (16 values, row-major) settings; the matrix wins when present.
mio::AffineTransform ReadAffineTransform(Parameters Settings)
{
    if (Settings["matrix"].size() > 0) {
        const Vector matrix = Settings["matrix"].GetVector();
        KRATOS_ERROR_IF_NOT(matrix.size() == 16)
            << "\"matrix\" needs 16 row-major values, got " << matrix.size() << std::endl;
        return mio::transform_from_matrix(&matrix[0]);
    }
    const auto axis = ReadVector3(Settings, "rotation_axis");
    const auto translation = ReadVector3(Settings, "translation");
    mio::AffineTransform xform;
    const double angle = Settings["rotation_angle"].GetDouble();
    if (std::abs(angle) > 0.0) {
        xform = mio::transform_rotation(axis[0], axis[1], axis[2], angle);
    }
    return mio::transform_compose(
        mio::transform_translation(translation[0], translation[1], translation[2]), xform);
}

} // namespace

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> MeshioPlusPlusMeshOperations::GetSupportedOperations()
{
    return OperationNames();
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::GetCapabilities()
{
    Parameters capabilities(R"({})");
    capabilities.AddBool("occt", mio::fit_cad_has_occt());
    capabilities.AddBool("netgen", mio::tetrahedralize_has_netgen());
    capabilities.AddBool("mmg", mio::mmg_available());
    capabilities.AddBool("mmg_elas", mio::mmg_has_elas());
    return capabilities;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::GetDefaultParameters()
{
    Parameters params(R"({
        "operation"                                   : "clean",
        "record_provenance"                           : true,
        "entity_type"                                 : "automatic",
        "use_deformed_configuration"                  : false,
        "weld"                                         : false,
        "tolerance"                                    : 1e-8,
        "remove_orphans"                               : true,
        "drop_degenerate"                              : true,
        "drop_duplicate_cells"                         : true,
        "translation"                                  : [0.0, 0.0, 0.0],
        "scale"                                        : [1.0, 1.0, 1.0],
        "rotation_axis"                                : [0.0, 0.0, 1.0],
        "rotation_angle"                               : 0.0,
        "rotate_vector_data"                           : false,
        "mode"                                         : "linearize",
        "levels"                                       : 1,
        "target_ratio"                                 : -1.0,
        "target_faces"                                 : -1,
        "max_error"                                    : -1.0,
        "preserve_boundary"                            : true,
        "preserve_features"                            : true,
        "feature_angle"                                : 30.0,
        "method"                                       : "taubin",
        "iterations"                                   : 10,
        "lambda"                                       : -1.0,
        "mu"                                           : -0.34,
        "fix_boundary"                                 : true,
        "guard_inversion"                              : true,

        "fix_orientation"                              : true,
        "orient_outward"                               : true,
        "fill_holes"                                   : true,
        "split_non_manifold"                           : true,
        "repair_record_provenance"                     : false,
        "max_hole_edges"                               : 10,
        "weld_tolerance"                               : 0.0,

        "length_scale"                                 : 0.0,
        "fixed_points_array"                           : "",
        "sobolev_fix_boundary"                         : false,
        "record_filtered"                              : false,
        "sobolev_max_iterations"                       : 128,
        "sobolev_tolerance"                            : 1e-10,

        "number_of_parts"                              : 2,
        "imbalance"                                    : 0.03,
        "seed"                                         : 0,
        "ghost_layers"                                 : 0,
        "weights_variable"                             : "",
        "split_by"                                     : "type",
        "tag_name"                                     : "",
        "box_min"                                      : [0.0, 0.0, 0.0],
        "box_max"                                      : [1.0, 1.0, 1.0],
        "origin"                                       : [0.0, 0.0, 0.0],
        "normal"                                       : [0.0, 0.0, 1.0],
        "keep_partial_cells"                           : false,
        "array_name"                                   : "",
        "isovalues"                                    : [0.0],
        "linearize"                                    : false,
        "record_parent_ids"                            : false,

        "nodal_solution_step_data_variables"           : [],
        "nodal_data_value_variables"                   : [],
        "nodal_flags"                                  : [],
        "element_data_value_variables"                 : [],
        "element_flags"                                : [],
        "condition_data_value_variables"               : [],
        "condition_flags"                              : [],
        "gauss_point_variables_in_elements"            : [],
        "write_ids"                                    : false,

        "expression"                                   : "",
        "location"                                     : "point",
        "output"                                        : "",
        "output_overwrite"                             : false,

        "names"                                        : [],
        "scope"                                        : "component",
        "lo"                                           : 0.0,
        "hi"                                           : 1.0,
        "nan_policy"                                   : "ignore",
        "nan_replacement"                              : 0.0,
        "output_suffix"                                : "",
        "preserve_dtype"                               : true,

        "drop"                                         : [],
        "keep"                                         : [],
        "rename"                                       : [],
        "ignore_missing"                               : false,

        "weight"                                       : "uniform",
        "prefix"                                       : "",
        "overwrite"                                    : true,

        "extrapolate"                                  : false,
        "default_value"                                : 0.0,
        "on_conflict"                                  : "error",

        "cells"                                        : [],
        "region"                                       : "",
        "predicate_array"                              : "",
        "predicate_op"                                 : "<",
        "predicate_value"                              : 0.0,
        "closure"                                      : "redgreen",
        "record_levels"                                : false,
        "record_hierarchy"                             : false,

        "gradient_operator"                            : "gradient",
        "gradient_method"                              : "green_gauss",
        "component"                                    : -1,

        "error_method"                                 : "zz",
        "marking"                                      : "none",
        "marking_value"                                : 0.0,
        "marked_name"                                  : "",

        "curvature_mean"                               : true,
        "curvature_gaussian"                           : true,
        "dual_area"                                    : "mixed-voronoi",
        "include_boundary"                             : false,
        "record_area"                                  : false,
        "record_principal"                             : false,

        "point_normals"                                : true,
        "cell_normals"                                 : false,
        "normal_weight"                                : "angle",
        "split"                                        : false,
        "split_angle"                                  : 30.0,

        "invariants"                                   : ["all"],

        "resolution"                                   : [],
        "cell_size"                                    : -1.0,
        "bounds"                                       : [],
        "padding"                                      : 0.0,
        "padding_relative"                             : 0.0,
        "max_cells"                                    : 20000000,
        "fill"                                         : "all",
        "attach_occupancy"                             : false,

        "structure"                                    : "voxel",
        "root_resolution"                              : 8,
        "max_depth"                                    : 4,
        "band_cells"                                   : 1.0,

        "sdf_sign"                                     : "pseudonormal",
        "sdf_weight"                                   : "angle",
        "sdf_location"                                 : "corner",
        "band"                                         : 0.0,
        "record_closest_cell"                          : false,
        "record_inside"                                : false,
        "watertight_check"                             : "warn",
        "grid_cell_size"                               : 0.0,
        "max_winding_work"                             : 2.0e9,

        "offset"                                       : 0.0,
        "max_distance"                                 : 0.0,
        "weights"                                      : "",
        "target_region"                                : "",
        "record_distance"                              : false,

        "dims"                                         : [1, 1, 1],
        "spacing"                                      : [1.0, 1.0, 1.0],

        "target_group_size"                            : 8,
        "simplexify_result"                            : true,
        "target_cells"                                 : -1,
        "placement"                                    : "optimal",
        "frozen"                                       : [],

        "num_clusters"                                 : 0,
        "subdivide_level"                              : -1,
        "subsample_ratio"                              : 10.0,
        "max_subdivide"                                : 4,
        "max_iterations"                               : 100,
        "max_repair_passes"                            : 10,
        "metric"                                       : "isotropic",
        "gradation"                                    : 0.0,
        "max_anisotropy"                               : 4.0,

        "max_tets"                                     : 20000000,
        "warp_fraction"                                : 0.35,

        "optimize_iterations"                          : 10,
        "relocate"                                     : true,
        "flip"                                         : true,
        "min_improvement"                              : 1.0e-6,

        "merge_coplanar_faces"                         : false,
        "coplanar_angle"                               : 1.0,
        "min_sphericity"                               : 0.0,

        "include_feature_edges"                        : true,
        "include_boundary_edges"                       : true,
        "include_non_manifold_edges"                   : true,
        "include_inconsistent_edges"                   : true,

        "normal_tolerance"                             : 0.0,
        "curvature_tolerance"                          : 0.0,
        "min_patch_cells"                              : 1,
        "patch_prefix"                                 : "patch",
        "fit_tolerance"                                : 1.0e-3,
        "fit_relative"                                 : true,
        "min_inlier_fraction"                          : 0.9,
        "fit_types"                                    : ["plane", "sphere", "cylinder", "cone", "torus"],
        "fit_max_iterations"                           : 50,

        "cad_step_path"                                : "",
        "cad_json_path"                                : "",
        "cad_tolerance"                                : 0.0,
        "cad_mesh_deflection"                          : 0.0,
        "cad_exact_edges"                              : true,
        "cad_freeform"                                 : true,
        "cad_freeform_degree"                          : 3,
        "cad_freeform_max_control"                     : 64,
        "cad_freeform_split"                           : true,
        "cad_fillets"                                  : true,
        "cad_threads"                                  : true,
        "cad_max_demotion_rounds"                      : 8,

        "tessellate_levels"                            : 2,
        "tessellate_curved"                            : true,
        "tessellate_fields"                            : true,
        "tessellate_record_stencil"                    : false,

        "edits"                                        : [],

        "require"                                      : "",
        "max_inverted"                                 : 0,
        "max_degenerate"                               : 0,

        "write_mdpa_ids"                               : false,

        "segment_method"                               : "angle",
        "segment_tolerance"                            : 1.0e-3,

        "refine_patches_tolerance"                     : 1.0e-3,
        "refine_patches_normal_weight"                 : 1.0,
        "refine_patches_epsilon"                       : 1.0e-3,
        "refine_patches_max_iterations"                : 10,
        "refine_patches_min_patch_cells"               : 3,
        "refine_patches_seed"                          : -1,
        "refine_patches_allow_split"                   : true,
        "refine_patches_allow_merge"                   : true,

        "cad_backend"                                  : "auto",
        "cad_merge_faces"                              : true,
        "cad_refine_segmentation"                      : false,
        "cad_freeform_general"                         : false,
        "cad_freeform_tolerance"                       : 0.0,
        "cad_freeform_smoothing"                       : 1.0e-12,
        "cad_max_freeform_poles"                       : 40,
        "cad_max_freeform_iterations"                  : 30,
        "cad_max_freeform_splits"                      : 2,
        "cad_sew_tolerance"                            : 0.0,
        "cad_extend_surfaces"                          : false,
        "cad_snap_junctions"                           : false,
        "cad_face_deviation"                           : false,
        "cad_freeform_mode"                            : "legacy",
        "cad_freeform_parameterization"                : "auto",
        "cad_extraordinary_vertex"                     : "g1x",
        "cad_local_refine"                             : "knots",
        "cad_freeform_g1_tolerance"                    : 1.0,
        "cad_freeform_max_control_points"              : 400,
        "cad_layout_target_patch_size"                 : 0.1,
        "cad_layout_max_valence_deviation"             : 2,
        "cad_resegment_rounds"                         : 2,
        "cad_layout_snap_to_features"                  : true,
        "cad_validate_by_sampling"                     : true,
        "cad_constrain_primitives"                     : false,

        "optimize_objective"                           : "worst",
        "max_ring_size"                                : 7,

        "plane"                                        : "auto",

        "tet_backend"                                  : "auto",
        "tet_max_size"                                 : 0.0,
        "tet_max_volume"                               : 0.0,
        "tet_grading"                                  : 0.3,
        "tet_optimize_steps"                           : 3,
        "tet_max_cells"                                : 50000000,
        "tet_face_array"                               : "cad:face",

        "bl_first_height"                              : 0.0,
        "bl_growth"                                    : 1.2,
        "bl_number_of_layers"                          : 3,
        "bl_reverse"                                   : false,
        "bl_collision_fraction"                        : 0.45,
        "bl_smooth_iterations"                         : 5,
        "bl_max_shrink"                                : 6,
        "bl_collapse"                                  : "points",
        "bl_normal_smooth_iterations"                  : 0,
        "bl_front_check"                               : true,

        "extrude_direction"                            : [0.0, 0.0, 1.0],
        "extrude_layers"                               : 1,
        "extrude_ratio"                                : 1.0,
        "extrude_offsets"                              : [],
        "revolve_axis_point"                           : [0.0, 0.0, 0.0],
        "revolve_axis"                                 : [0.0, 1.0, 0.0],
        "revolve_angle"                                : 360.0,
        "revolve_layers"                               : 16,
        "revolve_tolerance"                            : 1.0e-9,
        "record_layers"                                : false,

        "mirror_normal"                                : [1.0, 0.0, 0.0],
        "mirror_keep_original"                         : false,
        "mirror_weld"                                  : true,
        "mirror_tolerance"                             : 1.0e-9,
        "reflect_vector_data"                          : false,

        "octree_root_cells"                            : [1, 1, 1],
        "octree_min_level"                             : 0,
        "octree_max_level"                             : 4,
        "octree_dimension"                             : 0,
        "octree_sample_distance"                       : 0.0,
        "octree_balance"                               : "face",
        "octree_attach_distance"                       : false,
        "octree_even_odd"                              : false,
        "octree_transition"                            : "hanging",
        "octree_snap"                                  : false,
        "octree_surface_level"                         : -1,
        "octree_curvature_threshold"                   : 0.0,
        "octree_feature_angle"                         : 0.0,
        "octree_size_array"                            : "",

        "octree_mesh_output"                           : "tetra",
        "min_size"                                     : 0.0,
        "max_size"                                     : 0.0,
        "octree_mesh_max_level"                        : -1,
        "octree_mesh_padding_relative"                 : 0.1,
        "octree_mesh_curvature_angle"                  : 15.0,
        "octree_mesh_corner_angle"                     : 30.0,
        "octree_mesh_feature_size"                     : 0.0,
        "octree_mesh_gap_cells"                        : 2.0,
        "octree_mesh_gap_outside"                      : false,
        "octree_mesh_size_array"                       : "",
        "octree_mesh_surface_region"                   : "",
        "octree_mesh_grading_layers"                   : 1,
        "octree_mesh_skin_grading"                     : "graded",
        "octree_mesh_surface_size"                     : 0.0,
        "octree_mesh_winding_threshold"                : 0.5,
        "octree_mesh_winding_beta"                     : 2.0,
        "octree_mesh_conformity"                       : "stuffing",
        "octree_mesh_alpha_long"                       : 0.24799,
        "octree_mesh_alpha_short"                      : 0.40173,
        "octree_mesh_fit_iterations"                   : 5,
        "octree_mesh_dihedral_floor_guard"             : true,
        "octree_mesh_fill"                             : "inside",
        "octree_mesh_keep_outside"                     : false,
        "octree_mesh_attach_distance"                  : true,
        "octree_mesh_attach_quality"                   : false,
        "octree_mesh_pillow"                           : true,
        "octree_mesh_hex_grading"                      : "uniform",
        "octree_mesh_untangle_iterations"              : 50,
        "octree_mesh_max_tets"                         : 100000000,

        "select_name"                                  : "",
        "select_kind"                                  : "cell",
        "select_tag"                                   : -1,
        "select_mode"                                  : "all",
        "select_boundary_only"                         : true,
        "select_box"                                   : false,
        "select_plane"                                 : false,
        "select_sphere"                                : false,
        "select_sphere_center"                         : [0.0, 0.0, 0.0],
        "select_sphere_radius"                         : 0.0,
        "select_where_array"                           : "",
        "select_where_op"                              : "<",
        "select_where_value"                           : 0.0,
        "select_quality_metric"                        : "",
        "select_quality_op"                            : "<",
        "select_quality_value"                         : 0.0,
        "select_component"                             : -1,
        "select_seed"                                  : false,
        "select_seed_point"                            : [0.0, 0.0, 0.0],

        "conformity_tolerance"                         : 1.0e-9,

        "intersections_region_b"                       : "",
        "intersections_max_pairs"                      : 0,

        "thickness_direction"                          : "auto",
        "thickness_max_distance"                       : 0.0,
        "thickness_threshold"                          : 0.0,
        "thickness_region_target"                      : "",

        "error_array"                                  : "",
        "target_error"                                 : 0.0,
        "error_order"                                  : 1.0,

        "output_file"                                  : "",
        "output_format"                                : "",
        "tet_poly_feature_angle"                       : 1.0,
        "point_data_to_cell_data"                      : false,
        "record_vertex_ids"                            : false,

        "mmg_hmin"                                     : -1.0,
        "mmg_hmax"                                     : -1.0,
        "mmg_hausd"                                    : -1.0,
        "mmg_hgrad"                                    : -1.0,
        "mmg_hgradreq"                                 : -1.0,
        "mmg_angle_detection"                          : 45.0,
        "mmg_optim"                                    : false,
        "mmg_noinsert"                                 : false,
        "mmg_noswap"                                   : false,
        "mmg_nomove"                                   : false,
        "mmg_nosurf"                                   : false,
        "mmg_nreg"                                     : false,
        "mmg_opnbdy"                                   : false,
        "mmg_mem"                                      : -1,
        "mmg_verbose"                                  : -1,
        "mmg_required"                                 : [],
        "mmg_local"                                    : [],
        "mmg_hsiz"                                     : -1.0,
        "size_array"                                   : "",
        "metric_array"                                 : "",
        "level_set_array"                              : "",
        "isovalue"                                     : 0.0,
        "remove_small_components"                      : -1.0,
        "mmg_materials"                                : [],
        "mmg_holes"                                    : [],
        "displacement_array"                           : "",
        "moving_regions"                               : [],
        "lagrangian"                                   : 1,
        "hessian_array"                                : "",
        "metric_error"                                 : 1.0e-2,
        "metric_hmin"                                  : 1.0e-3,
        "metric_hmax"                                  : 1.0
    })");

    return params;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::Execute(
    const ModelPart& rSource,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSource.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    Settings.AddMissingParameters(GetDefaultParameters());
    const std::string operation = Settings["operation"].GetString();

    // "record_provenance" rather than the IO's "provenance": that one selects a *mode*
    // (off/best_effort/required) for a write, this one only says whether this operation joins
    // the chain a later write records. Different questions, so deliberately different names.
    //
    // Append to meshio++'s thread-local provenance chain rather than opening a scope: the
    // operation's result usually reaches disk through a *later* MeshioPlusPlusIO write, and a
    // scope closed here would take the record with it. A caller wanting the chain bounded opens
    // its own scope around the operation-then-write sequence.
    if (Settings["record_provenance"].GetBool()) {
        mio::detail::provenance_add_operation(operation);
    }

    const std::string entity_type = Settings["entity_type"].GetString();
    const bool write_elements = entity_type != "conditions";
    const bool write_conditions = entity_type != "elements";

    mio::Mesh mesh = Internals::ModelPartToMeshWithData(
        rSource, write_elements, write_conditions, Settings["use_deformed_configuration"].GetBool(),
        BuildFieldDataSelection(Settings));

    const auto crop_mode = Settings["keep_partial_cells"].GetBool() ? mio::CropMode::Any
                                                                   : mio::CropMode::All;
    Parameters report(R"({})");

    if (operation == "clean") {
        mio::CleanOptions options;
        options.weld = Settings["weld"].GetBool();
        options.atol = Settings["tolerance"].GetDouble();
        options.remove_orphans = Settings["remove_orphans"].GetBool();
        options.drop_degenerate = Settings["drop_degenerate"].GetBool();
        options.drop_duplicate_cells = Settings["drop_duplicate_cells"].GetBool();
        mio::CleanResult result = mio::clean(mesh, options);
        report.AddInt("points_welded", static_cast<int>(result.mPointsWelded));
        report.AddInt("points_removed_orphan", static_cast<int>(result.mPointsRemovedOrphan));
        report.AddInt("cells_dropped_degenerate", static_cast<int>(result.mCellsDroppedDegenerate));
        report.AddInt("cells_dropped_duplicate", static_cast<int>(result.mCellsDroppedDuplicate));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "transform") {
        // Composed in the usual order: scale, then rotate, then translate.
        const auto scale = ReadVector3(Settings, "scale");
        const auto axis = ReadVector3(Settings, "rotation_axis");
        const auto translation = ReadVector3(Settings, "translation");
        mio::AffineTransform xform = mio::transform_scale(scale[0], scale[1], scale[2]);
        const double angle = Settings["rotation_angle"].GetDouble();
        if (std::abs(angle) > 0.0) {
            xform = mio::transform_compose(mio::transform_rotation(axis[0], axis[1], axis[2], angle), xform);
        }
        xform = mio::transform_compose(
            mio::transform_translation(translation[0], translation[1], translation[2]), xform);
        mio::Mesh result = mio::transform(mesh, xform, Settings["rotate_vector_data"].GetBool());
        StoreResult(result, rDestination);

    } else if (operation == "convert_cells") {
        mio::ConvertCellsOptions options;
        options.mMode = mio::convert_cells_mode_from_name(Settings["mode"].GetString());
        options.mRecordParentIds = Settings["record_parent_ids"].GetBool();
        mio::ConvertCellsResult result = mio::convert_cells(mesh, options);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "tessellate") {
        // The isoparametric subdivision of triangle6, quad8, quad9, tetra10 and hexahedron27
        // cells into linear triangles and tetrahedra: the visualization/export route for Kratos
        // quadratic geometries. Every other cell passes through unchanged.
        mio::TessellateOptions options;
        options.mLevels = Settings["tessellate_levels"].GetInt();
        options.mCurved = Settings["tessellate_curved"].GetBool();
        options.mFields = Settings["tessellate_fields"].GetBool();
        options.mRecordStencil = Settings["tessellate_record_stencil"].GetBool();
        mio::TessellateResult result = mio::tessellate(mesh, options);
        report.AddInt("source_points", static_cast<int>(result.mNumSourcePoints));
        report.AddInt("source_cells", static_cast<int>(result.mNumSourceCells));
        report.AddInt("curved_source_cells", static_cast<int>(result.mNumCurvedSourceCells));
        report.AddInt("pass_through_source_cells", static_cast<int>(result.mNumPassThroughSourceCells));
        report.AddInt("points", static_cast<int>(result.mNumPoints));
        report.AddInt("cells", static_cast<int>(result.mNumCells));
        report.AddInt("boundary_facets", static_cast<int>(result.mNumBoundaryFacets));
        report.AddInt("interior_facets", static_cast<int>(result.mNumInteriorFacets));
        report.AddInt("facets_with_bad_count", static_cast<int>(result.mNumFacetsWithBadCount));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "refine") {
        mio::RefineOptions options;
        options.mLevels = Settings["levels"].GetInt();
        options.mRecordParentIds = Settings["record_parent_ids"].GetBool();
        // At most one of "cells"/"region"/"predicate_array" may be set; mio::refine itself
        // rejects two. None set is uniform refinement.
        options.mCells = ReadCellIndices(Settings["cells"]);
        options.mRegion = Settings["region"].GetString();
        options.mPredicateArray = Settings["predicate_array"].GetString();
        options.mPredicateOp = mio::refine_compare_from_name(Settings["predicate_op"].GetString());
        options.mPredicateValue = Settings["predicate_value"].GetDouble();
        options.mClosure = mio::refine_closure_from_name(Settings["closure"].GetString());
        options.mRecordLevels = Settings["record_levels"].GetBool();
        // Attaches "refine:cell_id"/"refine:parent_id" (and forces "refine:entity"), the
        // persistent parent/child hierarchy a multigrid caller resolves across the sequence of
        // meshes it keeps - and what "undo_green" needs to identify a green group later.
        options.mRecordHierarchy = Settings["record_hierarchy"].GetBool();
        mio::RefineResult result = mio::refine(mesh, options);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "subdivide") {
        // Polyhedral refinement: one polyhedral child per face, so unlike "refine" it handles
        // an arbitrary polyhedron and needs no per-type template - and needs no closure either,
        // since a shared face is never touched and no hanging node can appear. The result is
        // simplexified before the write-back, see StorePolyhedralResult.
        mio::SubdivideOptions options;
        options.mRecordParentIds = Settings["record_parent_ids"].GetBool();
        mio::SubdivideResult result = mio::subdivide(mesh, options);
        StorePolyhedralResult(result.mMesh, rDestination, operation,
                              Settings["simplexify_result"].GetBool());

    } else if (operation == "agglomerate") {
        mio::AgglomerateOptions options;
        options.mTargetGroupSize = static_cast<std::size_t>(Settings["target_group_size"].GetInt());
        options.mMergeCoplanarFaces = Settings["merge_coplanar_faces"].GetBool();
        options.mCoplanarAngleDeg = Settings["coplanar_angle"].GetDouble();
        options.mMinSphericity = Settings["min_sphericity"].GetDouble();
        mio::AgglomerateResult result = mio::agglomerate(mesh, options);
        report.AddInt("faces_merged", static_cast<int>(result.mNumFacesMerged));
        report.AddInt("candidates_rejected", static_cast<int>(result.mNumRejected));
        StorePolyhedralResult(result.mMesh, rDestination, operation,
                              Settings["simplexify_result"].GetBool());

    } else if (operation == "decimate") {
        mio::DecimateOptions options;
        options.mTargetRatio = Settings["target_ratio"].GetDouble();
        options.mTargetFaces = Settings["target_faces"].GetInt();
        options.mMaxError = Settings["max_error"].GetDouble();
        options.mPreserveBoundary = Settings["preserve_boundary"].GetBool();
        options.mPreserveFeatures = Settings["preserve_features"].GetBool();
        options.mFeatureAngleDeg = Settings["feature_angle"].GetDouble();
        mio::DecimateResult result = mio::decimate(mesh, options);
        report.AddInt("faces_removed", static_cast<int>(result.mFacesRemoved));
        report.AddInt("points_removed", static_cast<int>(result.mPointsRemoved));
        report.AddInt("collapses_rejected", static_cast<int>(result.mCollapsesRejected));
        report.AddDouble("max_error_applied", result.mMaxErrorApplied);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "decimate_volume") {
        // The tetrahedral counterpart of "decimate": the same quadric-collapse idea, but the
        // budget is a cell count rather than a face count, hence its own "target_cells".
        mio::DecimateVolumeOptions options;
        options.mTargetRatio = Settings["target_ratio"].GetDouble();
        options.mTargetCells = static_cast<std::int64_t>(Settings["target_cells"].GetInt());
        options.mMaxError = Settings["max_error"].GetDouble();
        options.mPlacement = mio::decimate_placement_from_name(Settings["placement"].GetString());
        options.mPreserveBoundary = Settings["preserve_boundary"].GetBool();
        options.mPreserveFeatures = Settings["preserve_features"].GetBool();
        options.mFeatureAngleDeg = Settings["feature_angle"].GetDouble();
        options.mFrozen = ReadFrozenFlags(Settings["frozen"]);
        mio::DecimateVolumeResult result = mio::decimate_volume(mesh, options);
        report.AddInt("tets_removed", static_cast<int>(result.mTetsRemoved));
        report.AddInt("points_removed", static_cast<int>(result.mPointsRemoved));
        report.AddInt("collapses_rejected", static_cast<int>(result.mCollapsesRejected));
        report.AddDouble("max_error_applied", result.mMaxErrorApplied);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "remesh") {
        // Surface remeshing (ACVD clustering): a triangle mesh in, a better-shaped triangle
        // mesh out. The volume counterpart is "remesh_volume", a separate operation rather
        // than a mode of this one.
        mio::RemeshOptions options;
        options.mNumClusters = static_cast<std::int64_t>(Settings["num_clusters"].GetInt());
        options.mSubdivide = Settings["subdivide_level"].GetInt();
        options.mSubsampleRatio = Settings["subsample_ratio"].GetDouble();
        options.mMaxSubdivide = Settings["max_subdivide"].GetInt();
        options.mMaxIterations = Settings["max_iterations"].GetInt();
        options.mMaxRepairPasses = Settings["max_repair_passes"].GetInt();
        options.mMetric = mio::remesh_metric_from_name(Settings["metric"].GetString());
        options.mGradation = Settings["gradation"].GetDouble();
        options.mPreserveBoundary = Settings["preserve_boundary"].GetBool();
        options.mMaxAnisotropy = Settings["max_anisotropy"].GetDouble();
        mio::RemeshResult result = mio::remesh(mesh, options);
        report.AddInt("number_of_clusters", static_cast<int>(result.mNumClusters));
        report.AddInt("number_of_iterations", static_cast<int>(result.mNumIterations));
        report.AddInt("subdivide_applied", result.mSubdivideApplied);
        report.AddInt("number_of_isolated_clusters", static_cast<int>(result.mNumIsolatedClusters));
        report.AddInt("number_of_non_manifold_vertices", static_cast<int>(result.mNumNonManifoldVertices));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "remesh_volume") {
        // A closed surface in, a tetrahedral volume out: a lattice is cut against the surface's
        // signed distance, so it shares the whole lattice/SDF settings block with "voxelize"
        // and "compute_sdf" rather than declaring its own.
        mio::RemeshVolumeOptions options;
        ReadLatticeSettings(Settings, options);
        options.mMaxTets = static_cast<std::int64_t>(Settings["max_tets"].GetInt());
        options.mWarpFraction = Settings["warp_fraction"].GetDouble();
        mio::RemeshVolumeResult result = mio::remesh_volume(mesh, options);
        report.AddInt("number_of_tets", static_cast<int>(result.mNumTets));
        report.AddInt("number_of_vertices_warped", static_cast<int>(result.mNumVerticesWarped));
        report.AddInt("number_of_tets_rejected", static_cast<int>(result.mNumTetsRejected));
        report.AddInt("number_of_non_manifold_edges", static_cast<int>(result.mNumNonManifoldEdges));
        report.AddValue("surface_quality", SurfaceQualityReport(result.mQuality));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "optimize_volume") {
        // Improves a tetrahedral mesh without changing its cell budget (vertex relocation plus
        // 2-3/3-2 flips) - the natural follow-up to "remesh_volume". "optimize_iterations"
        // rather than the shared "max_iterations": upstream's default is 10 here against
        // "remesh"'s 100, and one key cannot carry two defaults.
        mio::OptimizeVolumeOptions options;
        options.mMaxIterations = Settings["optimize_iterations"].GetInt();
        options.mRelocate = Settings["relocate"].GetBool();
        options.mFlip = Settings["flip"].GetBool();
        options.mPreserveBoundary = Settings["preserve_boundary"].GetBool();
        options.mMinImprovement = Settings["min_improvement"].GetDouble();
        options.mFrozen = ReadFrozenFlags(Settings["frozen"]);
        // "worst" (meshio++'s default since the v18.0.0 tag) raises the least dihedral angle and
        // never lowers it; "energy" is the former behaviour, byte for byte.
        options.mObjective = mio::optimize_volume_objective_from_name(Settings["optimize_objective"].GetString());
        options.mMaxRingSize = Settings["max_ring_size"].GetInt();
        mio::OptimizeVolumeResult result = mio::optimize_volume(mesh, options);
        report.AddInt("number_of_flips", static_cast<int>(result.mNumFlips));
        report.AddInt("number_of_2_3_flips", static_cast<int>(result.mNum23Flips));
        report.AddInt("number_of_3_2_flips", static_cast<int>(result.mNum32Flips));
        report.AddInt("number_of_edge_removals", static_cast<int>(result.mNumEdgeRemovals));
        report.AddInt("number_of_multi_face_removals", static_cast<int>(result.mNumMultiFaceRemovals));
        report.AddInt("number_of_worst_smoothed", static_cast<int>(result.mNumWorstSmoothed));
        report.AddInt("number_of_vertices_moved", static_cast<int>(result.mNumVerticesMoved));
        report.AddInt("number_of_tets", static_cast<int>(result.mNumTets));
        report.AddDouble("min_quality_before", result.mMinQualityBefore);
        report.AddDouble("min_quality_after", result.mMinQualityAfter);
        report.AddDouble("min_dihedral_before", result.mMinDihedralBefore);
        report.AddDouble("min_dihedral_after", result.mMinDihedralAfter);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "smooth") {
        mio::SmoothOptions options;
        options.mMethod = mio::smooth_method_from_name(Settings["method"].GetString());
        options.mIterations = Settings["iterations"].GetInt();
        options.mLambda = Settings["lambda"].GetDouble();
        options.mMu = Settings["mu"].GetDouble();
        options.mFixBoundary = Settings["fix_boundary"].GetBool();
        options.mPreserveFeatures = Settings["preserve_features"].GetBool();
        options.mFeatureAngleDeg = Settings["feature_angle"].GetDouble();
        options.mGuardInversion = Settings["guard_inversion"].GetBool();
        mio::SmoothResult result = mio::smooth(mesh, options);
        report.AddInt("nodes_moved", static_cast<int>(result.mNumNodesMoved));
        report.AddDouble("max_displacement", result.mMaxDisplacement);
        report.AddInt("skipped_inversion", static_cast<int>(result.mNumSkippedInversion));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "repair") {
        // A surface topology fix (orientation, holes, pinched vertices) - "repair_record_provenance"
        // rather than the shared "record_provenance": that one appends to meshio++'s provenance
        // chain for a later write, this one is RepairOptions::mRecordProvenance, attaching
        // "repair:parent_point"/"repair:hole" onto the *mesh itself* - a different question with
        // a different default, so it cannot share the key.
        mio::RepairOptions options;
        options.mFixOrientation = Settings["fix_orientation"].GetBool();
        options.mOrientOutward = Settings["orient_outward"].GetBool();
        options.mFillHoles = Settings["fill_holes"].GetBool();
        options.mSplitNonManifold = Settings["split_non_manifold"].GetBool();
        options.mRecordProvenance = Settings["repair_record_provenance"].GetBool();
        options.mMaxHoleEdges = static_cast<std::int64_t>(Settings["max_hole_edges"].GetInt());
        options.mWeldTolerance = Settings["weld_tolerance"].GetDouble();
        mio::RepairResult result = mio::repair(mesh, options);
        report.AddValue("surface_quality_before", SurfaceQualityReport(result.mQualityBefore));
        report.AddValue("surface_quality_after", SurfaceQualityReport(result.mQualityAfter));
        report.AddInt("number_of_flipped", static_cast<int>(result.mNumFlipped));
        report.AddInt("number_of_components", static_cast<int>(result.mNumComponents));
        report.AddInt("largest_component", static_cast<int>(result.mLargestComponent));
        report.AddInt("number_of_oriented_outward", static_cast<int>(result.mNumOrientedOutward));
        report.AddInt("number_of_unorientable", static_cast<int>(result.mNumUnorientable));
        report.AddInt("number_of_vertices_split", static_cast<int>(result.mNumVerticesSplit));
        report.AddInt("number_of_holes_detected", static_cast<int>(result.mNumHolesDetected));
        report.AddInt("number_of_holes_filled", static_cast<int>(result.mNumHolesFilled));
        report.AddInt("number_of_holes_skipped", static_cast<int>(result.mNumHolesSkipped));
        report.AddInt("number_of_faces_added", static_cast<int>(result.mNumFacesAdded));
        report.AddInt("number_of_points_added", static_cast<int>(result.mNumPointsAdded));
        report.AddInt("points_welded", static_cast<int>(result.mPointsWelded));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "sobolev_deform") {
        // Sobolev-filters a raw displacement field, then moves the mesh's points by the
        // filtered result - "sobolev_fix_boundary"/"sobolev_max_iterations"/"sobolev_tolerance"
        // rather than the shared "fix_boundary"/"max_iterations"/"tolerance": upstream's own
        // defaults for this operation (false, 128, 1e-10) differ from what those keys already
        // default to for "smooth"/"remesh"/"clean", and a single Kratos key can only carry one
        // default (the "optimize_iterations" precedent, see above).
        mio::SobolevOptions options;
        options.mArrayName = Settings["array_name"].GetString();
        KRATOS_ERROR_IF(options.mArrayName.empty())
            << "The \"sobolev_deform\" operation needs an \"array_name\"" << std::endl;
        options.mLengthScale = Settings["length_scale"].GetDouble();
        options.mFixedPoints = ReadFrozenFlags(Settings["frozen"]);
        options.mFixedPointsArray = Settings["fixed_points_array"].GetString();
        options.mFixBoundary = Settings["sobolev_fix_boundary"].GetBool();
        options.mRecordFiltered = Settings["record_filtered"].GetBool();
        options.mMaxIterations = Settings["sobolev_max_iterations"].GetInt();
        options.mTolerance = Settings["sobolev_tolerance"].GetDouble();
        mio::SobolevResult result = mio::sobolev_deform(mesh, options);
        report.AddInt("number_of_iterations", static_cast<int>(result.mNumIterations));
        report.AddDouble("residual", result.mResidual);
        report.AddBool("converged", result.mConverged);
        report.AddInt("number_of_fixed", static_cast<int>(result.mNumFixed));
        report.AddInt("number_of_isolated", static_cast<int>(result.mNumIsolated));
        report.AddDouble("max_displacement", result.mMaxDisplacement);
        RenameResultArray(result.mMesh, mio::DataLocation::Point, mio::kSobolevDisplacementName,
                          Settings["output"].GetString());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "reorder") {
        const std::int64_t before = mio::compute_bandwidth(mesh);
        mio::ReorderResult result = mio::reorder(mesh, mio::reorder_method_from_name(Settings["method"].GetString()));
        report.AddInt("bandwidth_before", static_cast<int>(before));
        report.AddInt("bandwidth_after", static_cast<int>(mio::compute_bandwidth(result.mMesh)));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "extract_surface") {
        mio::Mesh result = mio::extract_surface(mesh, Settings["record_parent_ids"].GetBool());
        StoreResult(result, rDestination);

    } else if (operation == "extract_skin") {
        mio::Mesh result = mio::extract_skin(mesh, Settings["linearize"].GetBool());
        StoreResult(result, rDestination);

    } else if (operation == "feature_edges") {
        mio::FeatureEdgeOptions options;
        options.mFeatureAngleDeg = Settings["feature_angle"].GetDouble();
        options.mFeature = Settings["include_feature_edges"].GetBool();
        options.mBoundary = Settings["include_boundary_edges"].GetBool();
        options.mNonManifold = Settings["include_non_manifold_edges"].GetBool();
        options.mInconsistent = Settings["include_inconsistent_edges"].GetBool();
        options.mRegion = Settings["region"].GetString();
        mio::FeatureEdgeResult result = mio::feature_edges(mesh, options);
        report.AddInt("feature_edges", static_cast<int>(result.mNumFeature));
        report.AddInt("boundary_edges", static_cast<int>(result.mNumBoundary));
        report.AddInt("non_manifold_edges", static_cast<int>(result.mNumNonManifold));
        report.AddInt("inconsistent_edges", static_cast<int>(result.mNumInconsistent));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "segment") {
        const mio::SegmentOptions options = ReadSegmentOptions(Settings);
        mio::SegmentResult result = mio::segment(mesh, options);
        report.AddInt("number_of_patches", static_cast<int>(result.mNumPatches));
        report.AddInt("merged_patches", static_cast<int>(result.mNumMerged));
        report.AddInt("feature_edges", static_cast<int>(result.mNumFeatureEdges));
        report.AddInt("boundary_edges", static_cast<int>(result.mNumBoundaryEdges));
        report.AddInt("non_manifold_edges", static_cast<int>(result.mNumNonManifoldEdges));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "fit_cad") {
        // fit_cad needs the closed skin that segment and fit_primitives annotated; their
        // colon-namespaced arrays cannot ride on a model part, so all three run back to back
        // here, with the segment and fit_* settings.
        //
        // "cad_backend" : "auto" uses OpenCASCADE when meshio++ has it and its native backend
        // otherwise (since v17.9.0); only an explicit "occt" needs the kernel.
        const auto backend = ReadEnum<mio::FitCadBackend>(Settings, "cad_backend",
            {{"auto", mio::FitCadBackend::Auto}, {"occt", mio::FitCadBackend::Occt},
             {"native", mio::FitCadBackend::Native}});
        KRATOS_ERROR_IF(backend == mio::FitCadBackend::Occt && !mio::fit_cad_has_occt())
            << "The \"fit_cad\" operation with \"cad_backend\" : \"occt\" needs a meshio++ built "
            << "with OpenCASCADE (-DMESHIOPLUSPLUS_WITH_OCCT=ON); \"auto\" or \"native\" use the "
            << "native backend instead" << std::endl;
        const mio::SegmentResult segmented = mio::segment(mesh, ReadSegmentOptions(Settings));
        mio::FitPrimitivesOptions fit_options;
        fit_options.mPatchData = mio::kSegmentPatchName;
        fit_options.mTolerance = Settings["fit_tolerance"].GetDouble();
        fit_options.mRelative = Settings["fit_relative"].GetBool();
        fit_options.mMinInlierFraction = Settings["min_inlier_fraction"].GetDouble();
        fit_options.mTypes = ReadFitTypes(Settings, "fit_types");
        fit_options.mMaxIterations = Settings["fit_max_iterations"].GetInt();
        const mio::FitPrimitivesResult fitted = mio::fit_primitives(segmented.mMesh, fit_options);
        mio::FitCadOptions options;
        options.mTolerance = Settings["cad_tolerance"].GetDouble();
        options.mExactEdges = Settings["cad_exact_edges"].GetBool();
        options.mTypes = ReadFitTypes(Settings, "fit_types");
        options.mMeshDeflection = Settings["cad_mesh_deflection"].GetDouble();
        options.mMaxDemotionRounds = Settings["cad_max_demotion_rounds"].GetInt();
        options.mFreeform = Settings["cad_freeform"].GetBool();
        options.mFreeformDegree = Settings["cad_freeform_degree"].GetInt();
        options.mFreeformMaxControl = Settings["cad_freeform_max_control"].GetInt();
        options.mFreeformSplit = Settings["cad_freeform_split"].GetBool();
        options.mFillets = Settings["cad_fillets"].GetBool();
        options.mThreads = Settings["cad_threads"].GetBool();
        options.mMergeFaces = Settings["cad_merge_faces"].GetBool();
        options.mBackend = backend;
        options.mRefineSegmentation = Settings["cad_refine_segmentation"].GetBool();
        options.mRefine = ReadRefinePatchesOptions(Settings);
        options.mFreeformGeneral = Settings["cad_freeform_general"].GetBool();
        options.mFreeformTolerance = Settings["cad_freeform_tolerance"].GetDouble();
        options.mFreeformSmoothing = Settings["cad_freeform_smoothing"].GetDouble();
        options.mMaxFreeformPoles = Settings["cad_max_freeform_poles"].GetInt();
        options.mMaxFreeformIterations = Settings["cad_max_freeform_iterations"].GetInt();
        options.mMaxFreeformSplits = Settings["cad_max_freeform_splits"].GetInt();
        options.mSewTolerance = Settings["cad_sew_tolerance"].GetDouble();
        options.mExtendSurfaces = Settings["cad_extend_surfaces"].GetBool();
        options.mSnapJunctions = Settings["cad_snap_junctions"].GetBool();
        options.mFaceDeviation = Settings["cad_face_deviation"].GetBool();
        options.mFreeformMode = ReadEnum<mio::FitCadFreeformMode>(Settings, "cad_freeform_mode",
            {{"legacy", mio::FitCadFreeformMode::Legacy}, {"adaptive", mio::FitCadFreeformMode::Adaptive},
             {"network", mio::FitCadFreeformMode::Network}, {"auto", mio::FitCadFreeformMode::Auto}});
        options.mFreeformParameterization = ReadEnum<mio::FitCadParameterization>(Settings, "cad_freeform_parameterization",
            {{"auto", mio::FitCadParameterization::Auto}, {"projection", mio::FitCadParameterization::Projection},
             {"lscm", mio::FitCadParameterization::Lscm}, {"harmonic", mio::FitCadParameterization::Harmonic}});
        options.mExtraordinaryVertex = ReadEnum<mio::FitCadExtraordinaryVertex>(Settings, "cad_extraordinary_vertex",
            {{"g0", mio::FitCadExtraordinaryVertex::G0}, {"g1", mio::FitCadExtraordinaryVertex::G1},
             {"g1x", mio::FitCadExtraordinaryVertex::G1x}, {"g1xx", mio::FitCadExtraordinaryVertex::G1xx}});
        options.mLocalRefine = ReadEnum<mio::FitCadLocalRefine>(Settings, "cad_local_refine",
            {{"knots", mio::FitCadLocalRefine::Knots}, {"thb", mio::FitCadLocalRefine::Thb}});
        options.mFreeformG1ToleranceDeg = Settings["cad_freeform_g1_tolerance"].GetDouble();
        options.mFreeformMaxControlPoints = Settings["cad_freeform_max_control_points"].GetInt();
        options.mLayoutTargetPatchSize = Settings["cad_layout_target_patch_size"].GetDouble();
        options.mLayoutMaxValenceDeviation = Settings["cad_layout_max_valence_deviation"].GetInt();
        options.mResegmentRounds = Settings["cad_resegment_rounds"].GetInt();
        options.mLayoutSnapToFeatures = Settings["cad_layout_snap_to_features"].GetBool();
        options.mValidateBySurfaceSampling = Settings["cad_validate_by_sampling"].GetBool();
        options.mConstrainPrimitives = Settings["cad_constrain_primitives"].GetBool();
        mio::FitCadResult result = mio::fit_cad(fitted.mMesh, Settings["cad_step_path"].GetString(), options);
        const std::string cad_json_path = Settings["cad_json_path"].GetString();
        if (!cad_json_path.empty()) {
            mio::write_cad_json(cad_json_path, result.mMesh, result.mBrep);
        }
        report.AddBool("valid", result.mValid);
        report.AddDouble("tolerance", result.mTolerance);
        report.AddDouble("volume", result.mVolume);
        report.AddDouble("skin_volume", result.mSkinVolume);
        report.AddDouble("volume_ratio", result.mVolumeRatio);
        report.AddDouble("max_deviation", result.mMaxDeviation);
        report.AddInt("number_of_faces", static_cast<int>(result.mNumFaces));
        report.AddInt("planes", static_cast<int>(result.mNumPlane));
        report.AddInt("spheres", static_cast<int>(result.mNumSphere));
        report.AddInt("cylinders", static_cast<int>(result.mNumCylinder));
        report.AddInt("cones", static_cast<int>(result.mNumCone));
        report.AddInt("tori", static_cast<int>(result.mNumTorus));
        report.AddInt("faceted_faces", static_cast<int>(result.mNumFacetedFaces));
        report.AddInt("number_of_patches", static_cast<int>(result.mNumPatches));
        report.AddInt("demoted_patches", static_cast<int>(result.mNumDemoted));
        report.AddInt("exact_edges", static_cast<int>(result.mNumExactEdges));
        report.AddInt("approximated_edges", static_cast<int>(result.mNumApproxEdges));
        report.AddInt("solids", static_cast<int>(result.mNumSolids));
        report.AddInt("bspline_faces", static_cast<int>(result.mNumBspline));
        report.AddInt("split_patches", static_cast<int>(result.mNumSplitPatches));
        report.AddInt("fillets", static_cast<int>(result.mNumFillet));
        report.AddInt("variable_fillets", static_cast<int>(result.mNumVariableFillet));
        report.AddInt("blends", static_cast<int>(result.mNumBlend));
        report.AddInt("conic_blends", static_cast<int>(result.mNumConicBlend));
        report.AddInt("vertex_blends", static_cast<int>(result.mNumVertexBlend));
        report.AddInt("chamfers", static_cast<int>(result.mNumChamfer));
        report.AddInt("threads", static_cast<int>(result.mNumHelical));
        report.AddInt("rounded_threads", static_cast<int>(result.mNumRoundedThread));
        report.AddInt("tapered_threads", static_cast<int>(result.mNumTaperedThread));
        report.AddInt("variable_pitch_threads", static_cast<int>(result.mNumVariablePitchThread));
        report.AddString("backend", result.mBackend == mio::FitCadBackend::Native ? "native" : "occt");
        report.AddBool("closed", result.mClosed);
        report.AddBool("brep_checked", result.mBrepChecked);
        AddFiniteDouble(report, "hausdorff_max", result.mHausdorffMax);
        AddFiniteDouble(report, "hausdorff_rms", result.mHausdorffRms);
        report.AddBool("refined", result.mRefined);
        report.AddBool("refinement_kept", result.mRefinementKept);
        report.AddInt("refine_iterations", static_cast<int>(result.mRefineIterations));
        report.AddInt("refine_splits", static_cast<int>(result.mNumRefineSplits));
        report.AddInt("refine_merges", static_cast<int>(result.mNumRefineMerges));
        report.AddDouble("initial_explained_area_fraction", result.mInitialExplainedAreaFraction);
        report.AddDouble("explained_area_fraction", result.mExplainedAreaFraction);
        report.AddInt("freeform_general_faces", static_cast<int>(result.mNumFreeformGeneral));
        report.AddInt("freeform_refinements", static_cast<int>(result.mNumFreeformRefinements));
        report.AddDouble("average_freeform_poles", result.mAvgFreeformPoles);
        report.AddInt("extended_edges", static_cast<int>(result.mNumExtendedEdges));
        report.AddInt("snapped_junctions", static_cast<int>(result.mNumSnappedJunctions));
        report.AddInt("freeform_adaptive_faces", static_cast<int>(result.mNumFreeformAdaptive));
        report.AddInt("network_patches", static_cast<int>(result.mNumNetworkPatches));
        report.AddInt("network_edges", static_cast<int>(result.mNumNetworkEdges));
        report.AddInt("extraordinary_vertices", static_cast<int>(result.mNumExtraordinaryVertices));
        report.AddInt("resegmented_patches", static_cast<int>(result.mNumResegmented));
        report.AddInt("escalations", static_cast<int>(result.mNumEscalations));
        report.AddDouble("facet_area_fraction", result.mFacetAreaFraction);
        report.AddDouble("max_g1_angle", result.mMaxG1AngleDeg);
        report.AddDouble("max_freeform_distance", result.mMaxFreeformDistance);
        report.AddDouble("rms_freeform_distance", result.mRmsFreeformDistance);
        report.AddInt("constrained_primitives", static_cast<int>(result.mNumConstrainedPrimitives));
        report.AddInt("constraint_groups", static_cast<int>(result.mNumConstraintGroups));
        // One entry per face of the solid: the patch it came from, its kind and how well it fits.
        Parameters faces(R"([])");
        for (const auto& r_face : result.mFaces) {
            Parameters face(R"({})");
            face.AddInt("patch", static_cast<int>(r_face.mPatch));
            face.AddInt("face", static_cast<int>(r_face.mFace));
            face.AddInt("kind", static_cast<int>(r_face.mKind));
            face.AddInt("demotion_reason", static_cast<int>(r_face.mDemotionReason));
            face.AddInt("number_of_control_points", static_cast<int>(r_face.mNumControl));
            AddFiniteDouble(face, "max_distance", r_face.mMaxDist);
            AddFiniteDouble(face, "rms_distance", r_face.mRmsDist);
            AddFiniteDouble(face, "max_g1_angle", r_face.mMaxG1AngleDeg);
            faces.Append(face);
        }
        report.AddValue("faces", faces);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "fit_primitives") {
        // fit_primitives reads the patches from the "segment:patch" cell data a segment run
        // attaches. A colon-namespaced name never survives the trip onto a model part, so the
        // two run back to back here, with the segment settings, instead of through the model part.
        const mio::SegmentResult segmented = mio::segment(mesh, ReadSegmentOptions(Settings));
        mio::FitPrimitivesOptions options;
        options.mPatchData = mio::kSegmentPatchName;
        options.mTolerance = Settings["fit_tolerance"].GetDouble();
        options.mRelative = Settings["fit_relative"].GetBool();
        options.mMinInlierFraction = Settings["min_inlier_fraction"].GetDouble();
        options.mTypes = ReadFitTypes(Settings, "fit_types");
        options.mMaxIterations = Settings["fit_max_iterations"].GetInt();
        mio::FitPrimitivesResult result = mio::fit_primitives(segmented.mMesh, options);
        report.AddInt("number_of_patches", static_cast<int>(result.mNumPatches));
        report.AddInt("planes", static_cast<int>(result.mNumPlane));
        report.AddInt("spheres", static_cast<int>(result.mNumSphere));
        report.AddInt("cylinders", static_cast<int>(result.mNumCylinder));
        report.AddInt("cones", static_cast<int>(result.mNumCone));
        report.AddInt("tori", static_cast<int>(result.mNumTorus));
        report.AddInt("unfitted_patches", static_cast<int>(result.mNumUnfitted));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "edit_regions") {
        // Each edit: {"operation": union|intersection|difference|rename|retag|delete,
        // "inputs": [selectors], "output_name", "output_dim", "output_tag", "keep_inputs"}.
        std::vector<mio::RegionEdit> edits;
        for (std::size_t i_edit = 0; i_edit < Settings["edits"].size(); ++i_edit) {
            Parameters r_edit = Settings["edits"][i_edit];
            mio::RegionEdit edit;
            KRATOS_ERROR_IF_NOT(r_edit.Has("operation"))
                << "edit_regions: every edit needs an \"operation\"" << std::endl;
            edit.mOp = mio::region_op_from_name(r_edit["operation"].GetString());
            if (r_edit.Has("inputs")) {
                for (std::size_t i_input = 0; i_input < r_edit["inputs"].size(); ++i_input) {
                    edit.mInputs.push_back(ReadRegionSelector(r_edit["inputs"][i_input], "edit_regions"));
                }
            }
            if (r_edit.Has("output_name")) {
                edit.mOutputName = r_edit["output_name"].GetString();
            }
            if (r_edit.Has("output_dim")) {
                edit.mOutputDim = r_edit["output_dim"].GetInt();
            }
            if (r_edit.Has("output_tag")) {
                edit.mOutputTag = r_edit["output_tag"].GetInt();
            }
            if (r_edit.Has("keep_inputs")) {
                edit.mKeepInputs = r_edit["keep_inputs"].GetBool();
            }
            edits.push_back(std::move(edit));
        }
        mio::Mesh result = mio::edit_regions(mesh, edits);
        report.AddInt("number_of_edits", static_cast<int>(edits.size()));
        StoreResult(result, rDestination);

    } else if (operation == "crop_bbox") {
        const auto lo = ReadVector3(Settings, "box_min");
        const auto hi = ReadVector3(Settings, "box_max");
        mio::CropResult result = mio::crop_bbox(mesh, lo.data(), hi.data(), crop_mode,
                                                Settings["record_parent_ids"].GetBool());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "crop_halfspace") {
        const auto origin = ReadVector3(Settings, "origin");
        const auto normal = ReadVector3(Settings, "normal");
        mio::CropResult result = mio::crop_halfspace(mesh, origin.data(), normal.data(), crop_mode,
                                                     Settings["record_parent_ids"].GetBool());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "crop_predicate") {
        // Deliberately no "keep_partial_cells": bbox and half-space test *points* and then need
        // a rule reducing a cell's several nodes to one verdict, whereas a cell_data predicate
        // is already one value per cell. meshio++ rejects the mode outright for the same reason.
        const std::string array_name = Settings["predicate_array"].GetString();
        KRATOS_ERROR_IF(array_name.empty())
            << "The \"crop_predicate\" operation needs a \"predicate_array\"" << std::endl;
        mio::CropResult result = mio::crop_predicate(
            mesh, array_name, mio::refine_compare_from_name(Settings["predicate_op"].GetString()),
            Settings["predicate_value"].GetDouble(), Settings["record_parent_ids"].GetBool());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "gradient") {
        mio::GradientOptions options;
        options.mArrayName = Settings["array_name"].GetString();
        KRATOS_ERROR_IF(options.mArrayName.empty())
            << "The \"gradient\" operation needs an \"array_name\"" << std::endl;
        options.mOperator = mio::gradient_operator_from_name(Settings["gradient_operator"].GetString());
        options.mMethod = mio::gradient_method_from_name(Settings["gradient_method"].GetString());
        options.mLocation = ReadLocation(Settings);
        options.mOutputName = Settings["output"].GetString();
        options.mOverwrite = Settings["output_overwrite"].GetBool();
        // Negative is the "every component" sentinel every meshio++ binding boundary uses.
        const int component = Settings["component"].GetInt();
        if (component >= 0) {
            options.mComponent = component;
        }
        mio::GradientResult result = mio::gradient(mesh, options);
        report.AddInt("number_of_skipped", static_cast<int>(result.mNumSkipped));
        report.AddInt("number_of_fallback", static_cast<int>(result.mNumFallback));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "hessian") {
        mio::HessianOptions options;
        options.mArrayName = Settings["array_name"].GetString();
        KRATOS_ERROR_IF(options.mArrayName.empty())
            << "The \"hessian\" operation needs an \"array_name\"" << std::endl;
        options.mMethod = mio::gradient_method_from_name(Settings["gradient_method"].GetString());
        options.mLocation = ReadLocation(Settings);
        options.mOutputName = Settings["output"].GetString();
        options.mOverwrite = Settings["output_overwrite"].GetBool();
        mio::HessianResult result = mio::hessian(mesh, options);
        report.AddInt("number_of_skipped", static_cast<int>(result.mNumSkipped));
        report.AddInt("number_of_fallback", static_cast<int>(result.mNumFallback));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "estimate_error") {
        // The error indicator "refine"'s predicate selector was built for: run this, then
        // "refine" with "predicate_array" pointed at the marked array. Both of its outputs are
        // colon-namespaced ("error:zz", "error:marked"), so neither is reachable from Kratos
        // without "output"/"marked_name" naming a registered Variable - see RenameResultArray.
        mio::ErrorOptions options;
        options.mArrayName = Settings["array_name"].GetString();
        KRATOS_ERROR_IF(options.mArrayName.empty())
            << "The \"estimate_error\" operation needs an \"array_name\"" << std::endl;
        options.mMethod = mio::error_method_from_name(Settings["error_method"].GetString());
        options.mMarking = mio::error_marking_from_name(Settings["marking"].GetString());
        options.mMarkingValue = Settings["marking_value"].GetDouble();
        options.mOutputName = Settings["output"].GetString();
        options.mMarkedName = Settings["marked_name"].GetString();
        options.mOverwrite = Settings["output_overwrite"].GetBool();
        mio::ErrorResult result = mio::estimate_error(mesh, options);
        report.AddDouble("global_error", result.mGlobalError);
        report.AddInt("number_of_skipped", static_cast<int>(result.mNumSkipped));
        report.AddInt("number_of_marked", static_cast<int>(result.mNumMarked));
        // mOutputName/mMarkedName are honored by meshio++ itself, so these only fire for the
        // default-named results.
        RenameResultArray(result.mMesh, mio::DataLocation::Cell, mio::kErrorZzName,
                          Settings["output"].GetString());
        RenameResultArray(result.mMesh, mio::DataLocation::Cell, mio::kErrorMarkedName,
                          Settings["marked_name"].GetString());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "curvature") {
        // Per-vertex mean/Gaussian curvature of a surface, over the standard discrete
        // estimators (angle defect for K, cotangent Laplace-Beltrami for H). Reuses "region"
        // (already the "restrict to this named Cell region" key for "refine") and "output" (the
        // primary result's rename target, the same single-array precedent "voxelize"/
        // "compute_sdf" use) - "curvature:gaussian"/"curvature:area"/"curvature:principal" stay
        // colon-namespaced and unreachable from Kratos, the same documented gap
        // "attach_quality"'s own invented names have.
        mio::CurvatureOptions options;
        options.mMean = Settings["curvature_mean"].GetBool();
        options.mGaussian = Settings["curvature_gaussian"].GetBool();
        options.mDualArea = mio::curvature_dual_area_from_name(Settings["dual_area"].GetString());
        options.mIncludeBoundary = Settings["include_boundary"].GetBool();
        options.mRecordArea = Settings["record_area"].GetBool();
        options.mRecordPrincipal = Settings["record_principal"].GetBool();
        options.mRegion = Settings["region"].GetString();
        mio::CurvatureResult result = mio::compute_curvature(mesh, options);
        report.AddValue("surface_quality", SurfaceQualityReport(result.mQuality));
        report.AddInt("number_of_boundary", static_cast<int>(result.mNumBoundary));
        report.AddInt("number_of_isolated", static_cast<int>(result.mNumIsolated));
        report.AddInt("number_of_degenerate", static_cast<int>(result.mNumDegenerate));
        report.AddDouble("total_angle_defect", result.mTotalAngleDefect);
        RenameResultArray(result.mMesh, mio::DataLocation::Point, mio::kCurvatureMeanName,
                          Settings["output"].GetString());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "compute_normals") {
        // Unit normals of a surface: per point (the angle/area/uniform-weighted average of the
        // incident faces) and optionally per cell. With "split" the points where the surface
        // creases by more than "split_angle" degrees are duplicated so that every point carries
        // one normal - the added nodes follow the originals and "record_parent_ids" names the
        // node each copy came from. A volume mesh is refused by name upstream (run
        // "extract_surface" first). "output" renames the "normals" array(s), e.g. to NORMAL.
        mio::NormalsOptions options;
        options.mPointNormals = Settings["point_normals"].GetBool();
        options.mCellNormals = Settings["cell_normals"].GetBool();
        options.mWeight = mio::sdf_weight_from_name(Settings["normal_weight"].GetString());
        options.mSplit = Settings["split"].GetBool();
        options.mSplitAngle = Settings["split_angle"].GetDouble();
        options.mRecordParentIds = Settings["record_parent_ids"].GetBool();
        options.mRegion = Settings["region"].GetString();
        mio::NormalsResult result = mio::compute_normals(mesh, options);
        report.AddValue("surface_quality", SurfaceQualityReport(result.mQuality));
        report.AddInt("number_of_isolated", static_cast<int>(result.mNumIsolated));
        report.AddInt("number_of_undefined", static_cast<int>(result.mNumUndefined));
        report.AddInt("number_of_degenerate", static_cast<int>(result.mNumDegenerate));
        report.AddInt("number_of_split_points", static_cast<int>(result.mNumSplitPoints));
        report.AddInt("number_of_added_points", static_cast<int>(result.mNumAddedPoints));
        RenameResultArray(result.mMesh, mio::DataLocation::Point, mio::kNormalsName,
                          Settings["output"].GetString());
        RenameResultArray(result.mMesh, mio::DataLocation::Cell, mio::kNormalsName,
                          Settings["output"].GetString());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "tensor_invariants") {
        // Von Mises, principal values (ascending), hydrostatic and deviatoric parts of each
        // six- (xx yy zz xy yz zx) or nine-component (row-major 3x3) array in "names", at the
        // "location". "invariants" picks the outputs ("mises", "principal", "hydrostatic",
        // "deviatoric" or "all"); the result arrays are named <prefix><name>_<invariant>
        // <output_suffix>, with the same "prefix"/"output_suffix"/"overwrite" keys the
        // data-averaging operations use. None of those names is a Kratos Variable, so
        // "rename" (data_manage's {"location", "from", "to"} list) is applied to the result -
        // the way to land e.g. "CAUCHY_STRESS_VECTOR_mises" on VON_MISES_STRESS.
        const auto invariant_names = Settings["invariants"].GetStringArray();
        KRATOS_ERROR_IF(invariant_names.empty())
            << "\"tensor_invariants\" needs at least one entry in \"invariants\"" << std::endl;
        unsigned outputs = 0;
        for (const auto& r_name : invariant_names) {
            outputs |= static_cast<unsigned>(mio::tensor_invariant_from_name(r_name));
        }
        mio::TensorInvariantsOptions options;
        options.location = mio::data_location_from_name(Settings["location"].GetString());
        options.names = Settings["names"].GetStringArray();
        options.outputs = static_cast<mio::TensorInvariant>(outputs);
        options.prefix = Settings["prefix"].GetString();
        options.suffix = Settings["output_suffix"].GetString();
        options.overwrite = Settings["overwrite"].GetBool();
        mio::Mesh result = mio::tensor_invariants(mesh, options);
        if (Settings["rename"].size() > 0) {
            mio::DataManageOptions rename_options;
            rename_options.rename = ReadDataRenames(Settings["rename"]);
            result = mio::data_manage(result, rename_options).mMesh;
        }
        StoreResult(result, rDestination);

    } else if (operation == "voxelize") {
        mio::VoxelOptions options;
        ReadLatticeSettings(Settings, options);
        options.mFill = mio::voxel_fill_from_name(Settings["fill"].GetString());
        options.mAttachOccupancy = Settings["attach_occupancy"].GetBool();
        options.mPlane = ReadPlanarPlane(Settings);
        mio::VoxelResult result = mio::voxelize(mesh, options);
        report.AddInt("number_of_occupied", static_cast<int>(result.mNumOccupied));
        AppendLatticeReport(report, result.mDims, result.mOrigin, result.mSpacing);
        RenameResultArray(result.mMesh, mio::DataLocation::Cell, mio::kVoxelOccupancyName,
                          Settings["output"].GetString());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "compute_sdf") {
        mio::SdfOptions options;
        ReadLatticeSettings(Settings, options);
        options.mStructure = mio::sdf_structure_from_name(Settings["structure"].GetString());
        options.mRootResolution = Settings["root_resolution"].GetInt();
        options.mMaxDepth = Settings["max_depth"].GetInt();
        options.mBandCells = Settings["band_cells"].GetDouble();
        options.mRecordLevels = Settings["record_levels"].GetBool();
        mio::SdfResult result = mio::compute_sdf(mesh, options);
        report.AddInt("max_depth", static_cast<int>(result.mMaxDepth));
        report.AddInt("number_of_banded", static_cast<int>(result.mNumBanded));
        AppendLatticeReport(report, result.mDims, result.mOrigin, result.mSpacing);
        report.AddValue("surface_quality", SurfaceQualityReport(result.mQuality));
        RenameResultArray(result.mMesh, mio::DataLocation::Point, mio::kSdfDistanceName,
                          Settings["output"].GetString());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "slice") {
        mio::SliceOptions options;
        options.mOrigin = ReadVector3(Settings, "origin");
        options.mNormal = ReadVector3(Settings, "normal");
        options.mRecordParentIds = Settings["record_parent_ids"].GetBool();
        mio::Mesh result = mio::slice(mesh, options);
        StoreResult(result, rDestination);

    } else if (operation == "isosurface") {
        mio::IsosurfaceOptions options;
        options.mArrayName = Settings["array_name"].GetString();
        KRATOS_ERROR_IF(options.mArrayName.empty())
            << "The \"isosurface\" operation needs an \"array_name\"" << std::endl;
        const Vector isovalues = Settings["isovalues"].GetVector();
        options.mIsovalues.assign(isovalues.begin(), isovalues.end());
        options.mRecordParentIds = Settings["record_parent_ids"].GetBool();
        mio::Mesh result = mio::isosurface(mesh, options);
        StoreResult(result, rDestination);

    } else if (operation == "attach_quality") {
        mio::Mesh result = mio::attach_quality(mesh);
        StoreResult(result, rDestination);

    } else if (operation == "data_calc") {
        const std::string expression = Settings["expression"].GetString();
        KRATOS_ERROR_IF(expression.empty())
            << "The \"data_calc\" operation needs an \"expression\"" << std::endl;
        mio::DataCalcOptions options;
        options.location = ReadLocation(Settings);
        options.output = Settings["output"].GetString();
        KRATOS_ERROR_IF(options.output.empty())
            << "The \"data_calc\" operation needs an \"output\" array name" << std::endl;
        options.overwrite = Settings["output_overwrite"].GetBool();
        mio::Mesh result = mio::data_calc(mesh, expression, options);
        StoreResult(result, rDestination);

    } else if (operation == "data_condition") {
        mio::DataConditionOptions options;
        options.location = ReadLocation(Settings);
        options.names = Settings["names"].GetStringArray();
        options.mode = mio::condition_mode_from_name(Settings["mode"].GetString());
        options.scope = mio::condition_scope_from_name(Settings["scope"].GetString());
        options.lo = Settings["lo"].GetDouble();
        options.hi = Settings["hi"].GetDouble();
        options.nan_policy = mio::nan_policy_from_name(Settings["nan_policy"].GetString());
        options.nan_replacement = Settings["nan_replacement"].GetDouble();
        options.suffix = Settings["output_suffix"].GetString();
        options.preserve_dtype = Settings["preserve_dtype"].GetBool();
        mio::Mesh result = mio::data_condition(mesh, options);
        StoreResult(result, rDestination);

    } else if (operation == "data_manage") {
        mio::DataManageOptions options;
        options.keep = ReadDataKeys(Settings["keep"]);
        options.drop = ReadDataKeys(Settings["drop"]);
        options.rename = ReadDataRenames(Settings["rename"]);
        options.ignore_missing = Settings["ignore_missing"].GetBool();
        mio::DataManageResult result = mio::data_manage(mesh, options);
        report.AddStringArray("dropped", result.mDropped);
        Parameters renamed(R"([])");
        for (const auto& r_rename : result.mRenamed) {
            Parameters entry(R"({})");
            entry.AddString("from", r_rename.first);
            entry.AddString("to", r_rename.second);
            renamed.Append(entry);
        }
        report.AddValue("renamed", renamed);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "data_info") {
        const mio::DataInfoReport info = mio::data_info(mesh);
        report.AddInt("number_of_point_data", static_cast<int>(info.mNumPointData));
        report.AddInt("number_of_cell_data", static_cast<int>(info.mNumCellData));
        report.AddInt("number_of_field_data", static_cast<int>(info.mNumFieldData));
        Parameters arrays(R"([])");
        for (const auto& r_array : info.mArrays) {
            Parameters entry(R"({})");
            entry.AddString("location", mio::data_location_name(r_array.mLocation));
            entry.AddString("name", r_array.mName);
            entry.AddInt("num_components", static_cast<int>(r_array.mNumComponents));
            entry.AddInt("num_values", static_cast<int>(r_array.mNumValues));
            entry.AddDouble("min", r_array.mMin);
            entry.AddDouble("max", r_array.mMax);
            entry.AddDouble("mean", r_array.mMean);
            entry.AddInt("num_nan", static_cast<int>(r_array.mNumNan));
            entry.AddInt("num_inf", static_cast<int>(r_array.mNumInf));
            arrays.Append(entry);
        }
        report.AddValue("arrays", arrays);
        // Report-only: nothing is written into rDestination.
        return report;

    } else if (operation == "data_integrate") {
        // Cell-measure-weighted reduction: the integral and the measure-weighted mean of each
        // array over the whole mesh and over every named cell region. Report-only, like
        // "data_info" - the mesh is unchanged, the answer is entirely in the report.
        mio::DataIntegrateOptions options;
        options.mArrayNames = Settings["names"].GetStringArray();
        const mio::DataIntegrateReport integral = mio::data_integrate(mesh, options);
        Parameters arrays(R"([])");
        for (const auto& r_array : integral.mArrays) {
            Parameters entry(R"({})");
            entry.AddString("name", r_array.mName);
            entry.AddInt("num_components", static_cast<int>(r_array.mNumComponents));
            entry.AddValue("domain", FieldIntegralRegionReport(r_array.mDomain));
            Parameters regions(R"([])");
            for (const auto& r_region : r_array.mRegions) {
                regions.Append(FieldIntegralRegionReport(r_region));
            }
            entry.AddValue("regions", regions);
            arrays.Append(entry);
        }
        report.AddValue("arrays", arrays);
        // Report-only: nothing is written into rDestination.
        return report;

    } else if (operation == "point_data_to_cell_data" || operation == "cell_data_to_point_data") {
        mio::DataAverageOptions options;
        options.names = Settings["names"].GetStringArray();
        options.weight = mio::cell_point_weight_from_name(Settings["weight"].GetString());
        options.prefix = Settings["prefix"].GetString();
        options.suffix = Settings["output_suffix"].GetString();
        options.overwrite = Settings["overwrite"].GetBool();
        options.nan_policy = mio::nan_policy_from_name(Settings["nan_policy"].GetString());
        options.nan_replacement = Settings["nan_replacement"].GetDouble();
        mio::Mesh result = operation == "point_data_to_cell_data"
            ? mio::point_data_to_cell_data(mesh, options)
            : mio::cell_data_to_point_data(mesh, options);
        StoreResult(result, rDestination);

    } else if (operation == "split") {
        mio::SplitResult result = mio::split(mesh, mio::split_by_from_name(Settings["split_by"].GetString()),
                                             Settings["tag_name"].GetString());
        report = StorePieces(result.mPieces, rDestination, "split");

    } else if (operation == "partition") {
        mio::PartitionOptions options;
        options.mNParts = Settings["number_of_parts"].GetInt();
        options.mImbalance = Settings["imbalance"].GetDouble();
        options.mSeed = Settings["seed"].GetInt();
        options.mGhostLayers = Settings["ghost_layers"].GetInt();
        options.mWeightsKey = Settings["weights_variable"].GetString();
        options.mRecordIds = Settings["record_parent_ids"].GetBool();
        mio::PartitionResult result = mio::partition(mesh, options);
        report = StorePieces(result.mPieces, rDestination, "partition");

    } else if (operation == "stats") {
        return ComputeStatistics(rSource);

    } else if (operation == "quality") {
        return ComputeQuality(rSource);

    } else if (operation == "check_quality") {
        // "require" is meshio++'s one threshold text on every surface, e.g.
        // "scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%".
        mio::QualityGateOptions options;
        options.mThresholds = mio::parse_quality_thresholds(Settings["require"].GetString());
        options.mMaxInverted = Settings["max_inverted"].GetInt();
        options.mMaxDegenerate = Settings["max_degenerate"].GetInt();
        const mio::QualityGateResult result = mio::check_quality(mesh, options);
        report.AddBool("passed", result.mPassed);
        report.AddString("summary", mio::quality_gate_summary(result));
        Parameters checks(R"([])");
        for (const auto& r_check : result.mChecks) {
            Parameters entry(R"({})");
            entry.AddString("name", r_check.mName);
            entry.AddBool("passed", r_check.mPassed);
            entry.AddInt("evaluated", static_cast<int>(r_check.mEvaluated));
            entry.AddInt("violations", static_cast<int>(r_check.mViolations));
            entry.AddDouble("fraction", r_check.mFraction);
            entry.AddDouble("worst", r_check.mWorst);
            entry.AddInt("worst_cell", static_cast<int>(r_check.mWorstCell));
            checks.Append(entry);
        }
        report.AddValue("checks", checks);
        return report;

    } else if (operation == "tetrahedralize") {
        // Fills a closed triangle surface with tetrahedra: Netgen keeps the surface's triangles,
        // the native octree backend fits it ("auto" takes Netgen when meshio++ has it).
        mio::TetrahedralizeOptions options;
        options.mBackend = ReadEnum<mio::TetrahedralizeBackend>(Settings, "tet_backend",
            {{"auto", mio::TetrahedralizeBackend::Auto}, {"netgen", mio::TetrahedralizeBackend::Netgen},
             {"octree", mio::TetrahedralizeBackend::Octree}});
        options.mMaxSize = Settings["tet_max_size"].GetDouble();
        options.mMaxVolume = Settings["tet_max_volume"].GetDouble();
        options.mGrading = Settings["tet_grading"].GetDouble();
        options.mOptimizeSteps = Settings["tet_optimize_steps"].GetInt();
        options.mMaxCells = Settings["tet_max_cells"].GetInt();
        options.mFaceArray = Settings["tet_face_array"].GetString();
        mio::TetrahedralizeResult result = mio::tetrahedralize(mesh, options);
        report.AddString("backend", result.mBackend);
        report.AddInt("number_of_tetrahedra", static_cast<int>(result.mNumTetrahedra));
        report.AddInt("number_of_points", static_cast<int>(result.mNumPoints));
        report.AddInt("number_of_new_points", static_cast<int>(result.mNumNewPoints));
        report.AddInt("number_of_regions", static_cast<int>(result.mNumRegions));
        report.AddInt("number_of_boundary_triangles", static_cast<int>(result.mNumBoundaryTriangles));
        report.AddBool("boundary_preserved", result.mBoundaryPreserved);
        report.AddDouble("hausdorff_to_surface", result.mHausdorffToSurface);
        report.AddDouble("min_dihedral", result.mMinDihedralDeg);
        report.AddDouble("max_dihedral", result.mMaxDihedralDeg);
        report.AddDouble("volume", result.mVolume);
        report.AddStringArray("warnings", result.mWarnings);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "boundary_layer") {
        // Grows prism (wedge) layers from a surface, the near-wall mesh of a CFD model.
        mio::BoundaryLayerOptions options;
        options.mFirstHeight = Settings["bl_first_height"].GetDouble();
        options.mGrowth = Settings["bl_growth"].GetDouble();
        options.mNumLayers = Settings["bl_number_of_layers"].GetInt();
        options.mReverse = Settings["bl_reverse"].GetBool();
        options.mCollisionFraction = Settings["bl_collision_fraction"].GetDouble();
        options.mSmoothIterations = Settings["bl_smooth_iterations"].GetInt();
        options.mMaxShrink = Settings["bl_max_shrink"].GetInt();
        options.mCollapse = ReadEnum<mio::BoundaryLayerCollapse>(Settings, "bl_collapse",
            {{"points", mio::BoundaryLayerCollapse::Points}, {"layers", mio::BoundaryLayerCollapse::Layers}});
        options.mNormalSmoothIterations = Settings["bl_normal_smooth_iterations"].GetInt();
        options.mFrontCheck = Settings["bl_front_check"].GetBool();
        mio::BoundaryLayerResult result = mio::boundary_layer(mesh, options);
        report.AddInt("number_of_cells", static_cast<int>(result.mNumCells));
        report.AddInt("wedges", static_cast<int>(result.mNumWedges));
        report.AddInt("hexahedra", static_cast<int>(result.mNumHexahedra));
        report.AddInt("pyramids", static_cast<int>(result.mNumPyramids));
        report.AddInt("tetrahedra", static_cast<int>(result.mNumTetrahedra));
        report.AddInt("shrunk_points", static_cast<int>(result.mNumShrunk));
        report.AddInt("collapsed_points", static_cast<int>(result.mNumCollapsed));
        report.AddInt("not_visible_points", static_cast<int>(result.mNumNotVisible));
        report.AddInt("fewer_layers_points", static_cast<int>(result.mNumFewerLayers));
        report.AddInt("front_crossings", static_cast<int>(result.mNumFrontCrossings));
        report.AddDouble("min_scaled_jacobian", result.mMinScaledJacobian);
        report.AddDouble("min_thickness", result.mMinThickness);
        report.AddDouble("max_thickness", result.mMaxThickness);
        report.AddBool("cap_closed", result.mCapClosed);
        report.AddStringArray("warnings", result.mWarnings);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "extrude" || operation == "revolve") {
        // Sweeps a surface (or a curve) into a volume (or a surface): along a direction in
        // layers ("extrude_direction" is the whole displacement, split into the layers), or
        // around an axis.
        mio::SweepResult result;
        if (operation == "extrude") {
            mio::ExtrudeOptions options;
            options.mDirection = ReadVector3(Settings, "extrude_direction");
            options.mNumLayers = Settings["extrude_layers"].GetInt();
            options.mRatio = Settings["extrude_ratio"].GetDouble();
            const Vector offsets = Settings["extrude_offsets"].GetVector();
            options.mOffsets.assign(offsets.begin(), offsets.end());
            options.mRecordLayers = Settings["record_layers"].GetBool();
            options.mMaxCells = Settings["max_cells"].GetInt();
            result = mio::extrude(mesh, options);
        } else {
            mio::RevolveOptions options;
            options.mAxisPoint = ReadVector3(Settings, "revolve_axis_point");
            options.mAxis = ReadVector3(Settings, "revolve_axis");
            options.mAngle = Settings["revolve_angle"].GetDouble();
            options.mNumLayers = Settings["revolve_layers"].GetInt();
            options.mTolerance = Settings["revolve_tolerance"].GetDouble();
            options.mRecordLayers = Settings["record_layers"].GetBool();
            options.mMaxCells = Settings["max_cells"].GetInt();
            result = mio::revolve(mesh, options);
        }
        report.AddInt("number_of_cells", static_cast<int>(result.mNumCells));
        report.AddInt("hexahedra", static_cast<int>(result.mNumHexahedra));
        report.AddInt("wedges", static_cast<int>(result.mNumWedges));
        report.AddInt("pyramids", static_cast<int>(result.mNumPyramids));
        report.AddInt("tetrahedra", static_cast<int>(result.mNumTetrahedra));
        report.AddInt("quadrilaterals", static_cast<int>(result.mNumQuads));
        report.AddInt("triangles", static_cast<int>(result.mNumTriangles));
        report.AddInt("flipped", static_cast<int>(result.mNumFlipped));
        report.AddInt("on_axis", static_cast<int>(result.mNumOnAxis));
        report.AddDouble("min_scaled_jacobian", result.mMinScaledJacobian);
        report.AddBool("full_revolution", result.mFull);
        report.AddStringArray("warnings", result.mWarnings);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "mirror") {
        // "mirror_normal" rather than the shared "normal": upstream's default plane is x = 0,
        // where "slice" and "crop_halfspace" default to z, and one key cannot carry both.
        mio::MirrorOptions options;
        options.mOrigin = ReadVector3(Settings, "origin");
        options.mNormal = ReadVector3(Settings, "mirror_normal");
        options.mKeepOriginal = Settings["mirror_keep_original"].GetBool();
        options.mWeld = Settings["mirror_weld"].GetBool();
        options.mTolerance = Settings["mirror_tolerance"].GetDouble();
        options.mReflectVectorData = Settings["reflect_vector_data"].GetBool();
        mio::MirrorResult result = mio::mirror(mesh, options);
        report.AddInt("flipped_cells", static_cast<int>(result.mNumFlippedCells));
        report.AddInt("points_on_plane", static_cast<int>(result.mNumPointsOnPlane));
        report.AddInt("dropped_cells", static_cast<int>(result.mNumDroppedCells));
        report.AddInt("dropped_sides", static_cast<int>(result.mNumDroppedSides));
        report.AddStringArray("warnings", result.mWarnings);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "octree") {
        // An adaptive octree (quadtree for a planar skin) refined toward the input surface. The
        // lattice keys ("bounds", "padding", "padding_relative", "max_cells", "fill", "sdf_*")
        // are voxelize's; its octree:* arrays are colon-namespaced and stay behind.
        mio::OctreeOptions options;
        const Vector root_cells = Settings["octree_root_cells"].GetVector();
        KRATOS_ERROR_IF(root_cells.size() != 3)
            << "\"octree_root_cells\" needs three entries, got " << root_cells.size() << std::endl;
        options.mRootCells = {{static_cast<std::int64_t>(root_cells[0]), static_cast<std::int64_t>(root_cells[1]),
                               static_cast<std::int64_t>(root_cells[2])}};
        options.mMinLevel = Settings["octree_min_level"].GetInt();
        options.mMaxLevel = Settings["octree_max_level"].GetInt();
        options.mBounds = ReadBounds(Settings);
        options.mPadding = Settings["padding"].GetDouble();
        options.mPaddingRelative = Settings["padding_relative"].GetDouble();
        options.mDimension = Settings["octree_dimension"].GetInt();
        options.mSampleDistance = Settings["octree_sample_distance"].GetDouble();
        options.mBalance = mio::octree_balance_from_name(Settings["octree_balance"].GetString());
        options.mFill = mio::voxel_fill_from_name(Settings["fill"].GetString());
        options.mAttachDistance = Settings["octree_attach_distance"].GetBool();
        options.mEvenOdd = Settings["octree_even_odd"].GetBool();
        options.mMaxCells = Settings["max_cells"].GetInt();
        options.mDistance = BuildSurfaceDistanceOptions(Settings);
        options.mTransition = mio::octree_transition_from_name(Settings["octree_transition"].GetString());
        options.mSnap = Settings["octree_snap"].GetBool();
        options.mSurfaceLevel = Settings["octree_surface_level"].GetInt();
        options.mCurvatureThreshold = Settings["octree_curvature_threshold"].GetDouble();
        options.mFeatureAngle = Settings["octree_feature_angle"].GetDouble();
        options.mSizeArray = Settings["octree_size_array"].GetString();
        options.mPlane = ReadPlanarPlane(Settings);
        mio::OctreeResult result = mio::octree(mesh, options);
        report.AddInt("dimension", static_cast<int>(result.mDimension));
        report.AddVector("origin", ToVector(result.mOrigin));
        report.AddDouble("root_size", result.mRootSize);
        report.AddInt("number_of_leaves", static_cast<int>(result.mNumLeaves));
        report.AddInt("hanging", static_cast<int>(result.mNumHanging));
        report.AddInt("finest_level", static_cast<int>(result.mFinestLevel));
        report.AddInt("snapped", static_cast<int>(result.mNumSnapped));
        report.AddInt("refused", static_cast<int>(result.mNumRefused));
        RenameResultArray(result.mMesh, mio::DataLocation::Point, mio::kSdfDistanceName,
                          Settings["output"].GetString());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "octree_mesh") {
        // Meshes the volume a skin encloses from a graded octree: tetrahedra (isosurface stuffing
        // or cleaving), hexahedra, or the castellated cells.
        mio::OctreeMeshOptions options;
        options.mOutput = mio::octree_mesh_output_from_name(Settings["octree_mesh_output"].GetString());
        options.mMinSize = Settings["min_size"].GetDouble();
        options.mMaxSize = Settings["max_size"].GetDouble();
        options.mMaxLevel = Settings["octree_mesh_max_level"].GetInt();
        options.mBounds = ReadBounds(Settings);
        options.mPadding = Settings["padding"].GetDouble();
        options.mPaddingRelative = Settings["octree_mesh_padding_relative"].GetDouble();
        options.mCurvatureAngle = Settings["octree_mesh_curvature_angle"].GetDouble();
        options.mFeatureAngle = Settings["feature_angle"].GetDouble();
        options.mCornerAngle = Settings["octree_mesh_corner_angle"].GetDouble();
        options.mFeatureSize = Settings["octree_mesh_feature_size"].GetDouble();
        options.mGapCells = Settings["octree_mesh_gap_cells"].GetDouble();
        options.mGapOutside = Settings["octree_mesh_gap_outside"].GetBool();
        options.mSizeArray = Settings["octree_mesh_size_array"].GetString();
        options.mSurfaceRegion = Settings["octree_mesh_surface_region"].GetString();
        options.mGradingLayers = Settings["octree_mesh_grading_layers"].GetInt();
        options.mSkinGrading = mio::octree_mesh_skin_grading_from_name(Settings["octree_mesh_skin_grading"].GetString());
        options.mSurfaceSize = Settings["octree_mesh_surface_size"].GetDouble();
        options.mWindingThreshold = Settings["octree_mesh_winding_threshold"].GetDouble();
        options.mWindingBeta = Settings["octree_mesh_winding_beta"].GetDouble();
        options.mConformity = mio::octree_mesh_conformity_from_name(Settings["octree_mesh_conformity"].GetString());
        options.mAlphaLong = Settings["octree_mesh_alpha_long"].GetDouble();
        options.mAlphaShort = Settings["octree_mesh_alpha_short"].GetDouble();
        options.mPreserveFeatures = Settings["preserve_features"].GetBool();
        options.mFitIterations = Settings["octree_mesh_fit_iterations"].GetInt();
        options.mDihedralFloorGuard = Settings["octree_mesh_dihedral_floor_guard"].GetBool();
        options.mFill = mio::octree_mesh_fill_from_name(Settings["octree_mesh_fill"].GetString());
        options.mKeepOutside = Settings["octree_mesh_keep_outside"].GetBool();
        options.mAttachDistance = Settings["octree_mesh_attach_distance"].GetBool();
        options.mAttachQuality = Settings["octree_mesh_attach_quality"].GetBool();
        options.mPillow = Settings["octree_mesh_pillow"].GetBool();
        options.mHexGrading = mio::octree_mesh_hex_grading_from_name(Settings["octree_mesh_hex_grading"].GetString());
        options.mUntangleIterations = Settings["octree_mesh_untangle_iterations"].GetInt();
        options.mMaxCells = Settings["max_cells"].GetInt();
        options.mMaxTets = Settings["octree_mesh_max_tets"].GetInt();
        mio::OctreeMeshResult result = mio::octree_mesh(mesh, options);
        report.AddString("output", mio::octree_mesh_output_name(result.mOutput));
        report.AddVector("origin", ToVector(result.mOrigin));
        report.AddDouble("root_size", result.mRootSize);
        report.AddInt("max_level", static_cast<int>(result.mMaxLevel));
        report.AddInt("finest_level", static_cast<int>(result.mFinestLevel));
        report.AddInt("surface_level", static_cast<int>(result.mSurfaceLevel));
        report.AddDouble("min_size", result.mMinSize);
        report.AddDouble("max_size", result.mMaxSize);
        report.AddInt("number_of_leaves", static_cast<int>(result.mNumLeaves));
        Parameters leaves_per_level(R"([])");
        for (const auto leaves : result.mLeavesPerLevel) {
            leaves_per_level.Append(static_cast<int>(leaves));
        }
        report.AddValue("leaves_per_level", leaves_per_level);
        report.AddInt("number_of_points", static_cast<int>(result.mNumPoints));
        report.AddInt("number_of_cells", static_cast<int>(result.mNumCells));
        report.AddInt("hexahedra", static_cast<int>(result.mNumHexahedra));
        report.AddInt("tetrahedra", static_cast<int>(result.mNumTetrahedra));
        report.AddInt("inside_cells", static_cast<int>(result.mNumInside));
        report.AddInt("boundary_cells", static_cast<int>(result.mNumBoundaryCells));
        report.AddInt("surrogate_faces", static_cast<int>(result.mNumSurrogateFaces));
        report.AddInt("hanging", static_cast<int>(result.mNumHanging));
        report.AddInt("cut_points", static_cast<int>(result.mNumCutPoints));
        report.AddInt("warped", static_cast<int>(result.mNumWarped));
        report.AddInt("warp_refused", static_cast<int>(result.mNumWarpRefused));
        report.AddDouble("min_dihedral", result.mMinDihedral);
        report.AddDouble("max_dihedral", result.mMaxDihedral);
        report.AddDouble("min_scaled_jacobian", result.mMinScaledJacobian);
        report.AddInt("inverted", static_cast<int>(result.mNumInverted));
        report.AddDouble("hausdorff_to_skin", result.mHausdorffToSkin);
        report.AddInt("open_boundary_edges", static_cast<int>(result.mNumBoundaryEdgesOpen));
        report.AddInt("non_manifold_boundary_edges", static_cast<int>(result.mNumNonManifoldBoundaryEdges));
        report.AddInt("corners", static_cast<int>(result.mNumCorners));
        report.AddInt("corners_snapped", static_cast<int>(result.mNumCornersSnapped));
        report.AddInt("feature_chains", static_cast<int>(result.mNumFeatureChains));
        report.AddInt("feature_chains_missed", static_cast<int>(result.mNumFeatureChainsMissed));
        report.AddInt("feature_moves", static_cast<int>(result.mNumFeatureMoves));
        report.AddInt("optimized_moves", static_cast<int>(result.mNumOptimizedMoves));
        report.AddInt("flips", static_cast<int>(result.mNumFlips));
        report.AddBool("guarantee_applies", result.mGuaranteeApplies);
        report.AddInt("winding_rechecks", static_cast<int>(result.mNumWindingRechecks));
        report.AddInt("winding_skipped", static_cast<int>(result.mNumWindingSkipped));
        report.AddStringArray("warnings", result.mWarnings);
        RenameResultArray(result.mMesh, mio::DataLocation::Point, mio::kSdfDistanceName,
                          Settings["output"].GetString());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "refine_patches") {
        // Improves a segmentation by Lloyd iterations over primitive proxies. It reads the
        // "segment:patch" cell data that no model part can carry, so the segmentation runs here
        // first, with the segment settings (the "fit_primitives" shape); the refined patches come
        // back as sub model parts <patch_prefix>:<k>.
        const mio::SegmentResult segmented = mio::segment(mesh, ReadSegmentOptions(Settings));
        mio::RefinePatchesResult result = mio::refine_patches(segmented.mMesh, ReadRefinePatchesOptions(Settings));
        report.AddInt("number_of_patches", static_cast<int>(result.mNumPatches));
        report.AddInt("iterations", static_cast<int>(result.mNumIterations));
        report.AddInt("splits", static_cast<int>(result.mNumSplits));
        report.AddInt("merges", static_cast<int>(result.mNumMerges));
        report.AddInt("absorbed", static_cast<int>(result.mNumAbsorbed));
        report.AddInt("primitive_patches", static_cast<int>(result.mNumPrimitive));
        report.AddInt("freeform_patches", static_cast<int>(result.mNumFreeform));
        report.AddInt("facet_patches", static_cast<int>(result.mNumFacets));
        report.AddDouble("tolerance", result.mTolerance);
        report.AddDouble("explained_area_fraction", result.mExplainedAreaFraction);
        report.AddDouble("initial_explained_area_fraction", result.mInitialExplainedAreaFraction);
        report.AddDouble("hausdorff_max", result.mHausdorffMax);
        report.AddDouble("hausdorff_rms", result.mHausdorffRms);
        report.AddDouble("initial_energy", result.mInitialEnergy);
        report.AddDouble("energy", result.mEnergy);
        report.AddInt("stop", static_cast<int>(result.mStop));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "select") {
        // Writes one named region - a sub model part of the destination - chosen by geometry,
        // data, quality, connectivity or growth from a seed; the criteria intersect.
        mio::SelectResult result = mio::select_region(mesh, ReadSelectOptions(Settings));
        report.AddInt("number_of_candidates", static_cast<int>(result.mNumCandidates));
        report.AddInt("number_of_selected", static_cast<int>(result.mNumSelected));
        report.AddInt("number_grown", static_cast<int>(result.mNumGrown));
        report.AddInt("number_of_components", static_cast<int>(result.mNumComponents));
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "check_conformity") {
        // Report-only, like "check_quality": its per-point flags are colon-namespaced.
        mio::ConformityOptions options;
        options.mTolerance = Settings["conformity_tolerance"].GetDouble();
        options.mAttachFlags = false;
        const mio::ConformityResult result = mio::check_conformity(mesh, options);
        report.AddBool("conforming", result.mConforming);
        report.AddInt("dimension", static_cast<int>(result.mDimension));
        report.AddInt("hanging_nodes", static_cast<int>(result.mNumHangingNodes));
        report.AddInt("t_junctions", static_cast<int>(result.mNumTJunctions));
        report.AddInt("partial_overlaps", static_cast<int>(result.mNumPartialOverlaps));
        report.AddInt("unwelded_facets", static_cast<int>(result.mNumUnweldedFacets));
        report.AddInt("boundary_facets", static_cast<int>(result.mNumBoundaryFacets));
        report.AddInt("non_manifold_facets", static_cast<int>(result.mNumNonManifoldFacets));
        report.AddString("interface_mode", result.mInterfaceMode);
        report.AddStringArray("warnings", result.mWarnings);
        return report;

    } else if (operation == "find_intersections") {
        // Where a surface intersects itself: the cell pairs in the report, the crossing
        // segments (line cells) in the destination.
        mio::FindIntersectionsResult result = mio::find_intersections(mesh, ReadFindIntersectionsOptions(Settings));
        report = FindIntersectionsReport(result);
        StoreResult(result.mSegments, rDestination);

    } else if (operation == "thickness") {
        // Wall thickness per facet, by a ray from each centroid along its normal; "output"
        // renames "thickness:distance" to a registered variable.
        mio::ThicknessResult result = mio::thickness(mesh, ReadThicknessOptions(Settings));
        report = ThicknessReport(result);
        RenameResultArray(result.mMesh, mio::DataLocation::Cell, mio::kThicknessDistanceName,
                          Settings["output"].GetString());
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "error_size") {
        // A target point size from a cell error indicator (the input "adapt" and "remesh" read);
        // "error_array" names a registered variable the field data settings carried in.
        mio::ErrorSizeOptions options;
        options.mErrorArray = Settings["error_array"].GetString();
        KRATOS_ERROR_IF(options.mErrorArray.empty())
            << "The \"error_size\" operation needs an \"error_array\"" << std::endl;
        options.mTargetError = Settings["target_error"].GetDouble();
        options.mOrder = Settings["error_order"].GetDouble();
        options.mMinSize = Settings["min_size"].GetDouble();
        options.mMaxSize = Settings["max_size"].GetDouble();
        if (!Settings["output"].GetString().empty()) {
            options.mOutputName = Settings["output"].GetString();
        }
        options.mOverwrite = Settings["output_overwrite"].GetBool();
        mio::ErrorSizeResult result = mio::error_size(mesh, options);
        report.AddInt("number_of_cells", static_cast<int>(result.mNumCells));
        report.AddInt("skipped", static_cast<int>(result.mNumSkipped));
        report.AddInt("clamped", static_cast<int>(result.mNumClamped));
        report.AddDouble("min_size", result.mMinSize);
        report.AddDouble("max_size", result.mMaxSize);
        StoreResult(result.mMesh, rDestination);

    } else if (operation == "tet_to_poly") {
        // The median dual of a tetrahedral mesh, one polyhedron per vertex. No ModelPart can
        // hold a polyhedron (and simplexifying would undo the conversion), so the result goes to
        // "output_file" and the destination is left untouched.
        const std::string output_file = Settings["output_file"].GetString();
        KRATOS_ERROR_IF(output_file.empty())
            << "The \"tet_to_poly\" operation produces polyhedra, which a Kratos ModelPart cannot "
            << "hold; set \"output_file\" (and optionally \"output_format\") to write them" << std::endl;
        mio::TetToPolyOptions options;
        options.mFeatureAngleDeg = Settings["tet_poly_feature_angle"].GetDouble();
        options.mPointDataToCellData = Settings["point_data_to_cell_data"].GetBool();
        options.mRecordVertexIds = Settings["record_vertex_ids"].GetBool();
        const mio::TetToPolyResult result = mio::tet_to_poly(mesh, options);
        const std::string format = mio::resolve_write_format(output_file, Settings["output_format"].GetString());
        mio::registry_writers().at(format)(output_file, result.mMesh);
        report.AddInt("number_of_new_points", static_cast<int>(result.mNumNewPoints));
        report.AddInt("number_of_boundary_faces", static_cast<int>(result.mNumBoundaryFaces));
        report.AddInt("number_of_split_vertices", static_cast<int>(result.mNumSplitVertices));
        report.AddString("output_file", output_file);
        report.AddString("output_format", format);
        return report;

    } else if (operation == "adapt" || operation == "discretize_levelset" || operation == "mesh_domain"
               || operation == "move") {
        // The MMG remeshers (mmg2d, mmgs, mmg3d chosen by the mesh). The size, metric, level set
        // and displacement arrays are registered variable names the field data settings carried
        // in. Without MMG (or, for "move", without its ELAS add-on) meshio++ refuses by name.
        const mio::MmgParameters parameters = ReadMmgParameters(Settings);
        mio::MmgReport mmg_report;
        mio::Mesh result_mesh;
        if (operation == "adapt") {
            mio::AdaptOptions options;
            options.mParameters = parameters;
            options.mHsiz = Settings["mmg_hsiz"].GetDouble();
            options.mSize = Settings["size_array"].GetString();
            options.mMetric = Settings["metric_array"].GetString();
            mio::AdaptResult result = mio::adapt(mesh, options);
            mmg_report = std::move(result.mReport);
            result_mesh = std::move(result.mMesh);
        } else if (operation == "discretize_levelset") {
            mio::DiscretizeLevelsetOptions options;
            options.mParameters = parameters;
            options.mLevelSet = Settings["level_set_array"].GetString();
            options.mIsovalue = Settings["isovalue"].GetDouble();
            options.mRemoveSmallComponents = Settings["remove_small_components"].GetDouble();
            for (std::size_t i = 0; i < Settings["mmg_materials"].size(); ++i) {
                Parameters material = Settings["mmg_materials"][i];
                material.AddMissingParameters(Parameters(R"({"ref": 0, "split": true, "ref_inside": 3, "ref_outside": 2})"));
                mio::MmgMaterial entry;
                entry.mRef = material["ref"].GetInt();
                entry.mSplit = material["split"].GetBool();
                entry.mRefInside = material["ref_inside"].GetInt();
                entry.mRefOutside = material["ref_outside"].GetInt();
                options.mMaterials.push_back(entry);
            }
            mio::DiscretizeLevelsetResult result = mio::discretize_levelset(mesh, options);
            mmg_report = std::move(result.mReport);
            result_mesh = std::move(result.mMesh);
        } else if (operation == "mesh_domain") {
            mio::MeshDomainOptions options;
            options.mParameters = parameters;
            options.mHsiz = Settings["mmg_hsiz"].GetDouble();
            options.mSize = Settings["size_array"].GetString();
            const Vector holes = Settings["mmg_holes"].GetVector();
            options.mHoles.assign(holes.begin(), holes.end());
            mio::MeshDomainResult result = mio::mesh_domain(mesh, options);
            report.AddInt("number_of_subdomains", static_cast<int>(result.mNumSubdomains));
            report.AddInt("number_of_holes", static_cast<int>(result.mNumHoles));
            mmg_report = std::move(result.mReport);
            result_mesh = std::move(result.mMesh);
        } else {
            mio::MoveOptions options;
            options.mParameters = parameters;
            options.mDisplacement = Settings["displacement_array"].GetString();
            options.mMoving = Settings["moving_regions"].GetStringArray();
            options.mLagrangian = Settings["lagrangian"].GetInt();
            mio::AdaptResult result = mio::move(mesh, options);
            mmg_report = std::move(result.mReport);
            result_mesh = std::move(result.mMesh);
        }
        report.AddBool("success", mmg_report.mStatus == mio::MmgStatus::Success);
        report.AddString("library", mmg_report.mLibrary);
        report.AddInt("number_of_points_in", static_cast<int>(mmg_report.mNumPointsIn));
        report.AddInt("number_of_points_out", static_cast<int>(mmg_report.mNumPointsOut));
        report.AddInt("number_of_elements_in", static_cast<int>(mmg_report.mNumElementsIn));
        report.AddInt("number_of_elements_out", static_cast<int>(mmg_report.mNumElementsOut));
        AddFiniteDouble(report, "min_quality", mmg_report.mQuality.mMin);
        AddFiniteDouble(report, "max_quality", mmg_report.mQuality.mMax);
        AddFiniteDouble(report, "mean_quality", mmg_report.mQuality.mMean);
        report.AddStringArray("warnings", mmg_report.mWarnings);
        StoreResult(result_mesh, rDestination);

    } else if (operation == "metric_from_hessian") {
        // The anisotropic metric "adapt" reads, from a recovered Hessian ("hessian" first, with
        // "output" naming its result); needs no MMG.
        mio::MetricFromHessianOptions options;
        options.mHessian = Settings["hessian_array"].GetString();
        KRATOS_ERROR_IF(options.mHessian.empty())
            << "The \"metric_from_hessian\" operation needs a \"hessian_array\"" << std::endl;
        options.mError = Settings["metric_error"].GetDouble();
        options.mHmin = Settings["metric_hmin"].GetDouble();
        options.mHmax = Settings["metric_hmax"].GetDouble();
        if (!Settings["output"].GetString().empty()) {
            options.mOutput = Settings["output"].GetString();
        }
        mio::Mesh result = mio::metric_from_hessian(mesh, options);
        StoreResult(result, rDestination);

    } else {
        std::ostringstream supported;
        for (const auto& r_name : OperationNames()) {
            supported << " \"" << r_name << "\"";
        }
        KRATOS_ERROR << "Unknown meshio++ operation \"" << operation << "\". Supported operations are:"
                     << supported.str() << std::endl;
    }

    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::Merge(
    const std::vector<const ModelPart*>& rSources,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSources.size() < 2) << "Merging needs at least two model parts, got "
        << rSources.size() << std::endl;

    Settings.AddMissingParameters(GetDefaultParameters());

    // The meshes must outlive the merge call, which takes pointers into them.
    std::vector<mio::Mesh> meshes;
    meshes.reserve(rSources.size());
    for (const ModelPart* p_source : rSources) {
        KRATOS_ERROR_IF(p_source->IsDistributed())
            << "The meshio++ operations do not support distributed model parts" << std::endl;
        meshes.push_back(Internals::ModelPartToMesh(*p_source, true, true,
                                                    Settings["use_deformed_configuration"].GetBool()));
    }

    std::vector<const mio::Mesh*> pointers;
    pointers.reserve(meshes.size());
    for (const mio::Mesh& r_mesh : meshes) {
        pointers.push_back(&r_mesh);
    }

    mio::MergeOptions options;
    options.weld = Settings["weld"].GetBool();
    options.atol = Settings["tolerance"].GetDouble();
    options.drop_duplicate_cells = Settings["drop_duplicate_cells"].GetBool();

    mio::MergeResult result = mio::merge(pointers, options);
    Internals::MeshToModelPart(result.mMesh, rDestination);

    Parameters report(R"({})");
    report.AddInt("number_of_sources", static_cast<int>(rSources.size()));
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::Interpolate(
    const ModelPart& rSource,
    const ModelPart& rTarget,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSource.IsDistributed() || rTarget.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    Settings.AddMissingParameters(GetDefaultParameters());

    const bool deformed = Settings["use_deformed_configuration"].GetBool();
    const Internals::FieldDataSelection selection = BuildFieldDataSelection(Settings);
    mio::Mesh source = Internals::ModelPartToMeshWithData(rSource, true, true, deformed, selection);
    mio::Mesh target = Internals::ModelPartToMeshWithData(rTarget, true, true, deformed, selection);

    mio::InterpolateOptions options;
    options.mMethod = mio::interpolate_method_from_name(Settings["method"].GetString());
    options.mArrays = Settings["names"].GetStringArray();
    options.mExtrapolate = Settings["extrapolate"].GetBool();
    options.mDefaultValue = Settings["default_value"].GetDouble();
    options.mOnConflict = mio::interpolate_conflict_from_name(Settings["on_conflict"].GetString());

    mio::Mesh result = mio::interpolate(source, target, options);
    Internals::MeshToModelPart(result, rDestination);

    Parameters report(R"({})");
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::ConservativeInterpolate(
    const ModelPart& rSource,
    const ModelPart& rTarget,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSource.IsDistributed() || rTarget.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    Settings.AddMissingParameters(GetDefaultParameters());

    const bool deformed = Settings["use_deformed_configuration"].GetBool();
    const Internals::FieldDataSelection selection = BuildFieldDataSelection(Settings);
    mio::Mesh source = Internals::ModelPartToMeshWithData(rSource, true, true, deformed, selection);
    mio::Mesh target = Internals::ModelPartToMeshWithData(rTarget, true, true, deformed, selection);

    mio::ConservativeInterpolateOptions options;
    options.mArrays = Settings["names"].GetStringArray();
    options.mDefaultValue = Settings["default_value"].GetDouble();
    options.mOnConflict = mio::conservative_interpolate_conflict_from_name(
        Settings["on_conflict"].GetString());

    mio::Mesh result = mio::conservative_interpolate(source, target, options);
    Internals::MeshToModelPart(result, rDestination);

    Parameters report(R"({})");
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::Shrinkwrap(
    const ModelPart& rSource,
    const ModelPart& rTarget,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSource.IsDistributed() || rTarget.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    Settings.AddMissingParameters(GetDefaultParameters());

    const bool deformed = Settings["use_deformed_configuration"].GetBool();
    const Internals::FieldDataSelection selection = BuildFieldDataSelection(Settings);
    mio::Mesh source = Internals::ModelPartToMeshWithData(rSource, true, true, deformed, selection);
    mio::Mesh target = Internals::ModelPartToMeshWithData(rTarget, true, true, deformed, selection);

    // "sdf_weight" is the same SdfPseudonormalWeight axis DistanceToSurface/compute_sdf already
    // expose - the hit feature's pseudonormal weighting, not a per-point blend factor (that is
    // "weights", the array name below).
    mio::ShrinkwrapOptions options;
    options.mOffset = Settings["offset"].GetDouble();
    options.mMaxDistance = Settings["max_distance"].GetDouble();
    options.mWeights = Settings["weights"].GetString();
    options.mTargetRegion = Settings["target_region"].GetString();
    options.mNormalWeight = mio::sdf_weight_from_name(Settings["sdf_weight"].GetString());
    options.mRecordDistance = Settings["record_distance"].GetBool();
    options.mRecordClosestCell = Settings["record_closest_cell"].GetBool();
    options.mGridCellSize = Settings["grid_cell_size"].GetDouble();

    mio::ShrinkwrapResult result = mio::shrinkwrap(source, target, options);
    RenameResultArray(result.mMesh, mio::DataLocation::Point, mio::kShrinkwrapDistanceName,
                      Settings["output"].GetString());
    Internals::MeshToModelPart(result.mMesh, rDestination);

    Parameters report(R"({})");
    report.AddValue("target_quality", SurfaceQualityReport(result.mQuality));
    report.AddInt("number_of_projected", static_cast<int>(result.mNumProjected));
    report.AddInt("number_of_missed", static_cast<int>(result.mNumMissed));
    report.AddInt("number_of_skipped", static_cast<int>(result.mNumSkipped));
    report.AddDouble("max_displacement", result.mMaxDisplacement);
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::Grid(
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    Settings.AddMissingParameters(GetDefaultParameters());

    const Vector dims = Settings["dims"].GetVector();
    KRATOS_ERROR_IF(dims.size() != 3)
        << "\"dims\" needs three entries, got " << dims.size() << std::endl;
    const std::array<std::int64_t, 3> counts{{static_cast<std::int64_t>(dims[0]),
                                              static_cast<std::int64_t>(dims[1]),
                                              static_cast<std::int64_t>(dims[2])}};
    const auto origin = ReadVector3(Settings, "origin");
    const auto spacing = ReadVector3(Settings, "spacing");

    mio::Mesh result = mio::grid(counts, origin, spacing,
                                 static_cast<std::int64_t>(Settings["max_cells"].GetInt()));
    Internals::MeshToModelPart(result, rDestination);

    Parameters report(R"({})");
    AppendLatticeReport(report, counts, origin, spacing);
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::FindIntersections(
    const ModelPart& rFirst,
    const ModelPart& rSecond,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rFirst.IsDistributed() || rSecond.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    Settings.ValidateAndAssignDefaults(Parameters(R"({
        "region"                  : "",
        "intersections_region_b"  : "",
        "intersections_max_pairs" : 0
    })"));

    const mio::Mesh first = Internals::ModelPartToMesh(rFirst);
    const mio::Mesh second = Internals::ModelPartToMesh(rSecond);
    mio::FindIntersectionsResult result = mio::find_intersections(first, second, ReadFindIntersectionsOptions(Settings));
    Parameters report = FindIntersectionsReport(result);
    StoreResult(result.mSegments, rDestination);
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::Thickness(
    const ModelPart& rSource,
    const ModelPart& rTarget,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSource.IsDistributed() || rTarget.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    Settings.ValidateAndAssignDefaults(Parameters(R"({
        "thickness_direction"     : "auto",
        "thickness_max_distance"  : 0.0,
        "thickness_threshold"     : 0.0,
        "region"                  : "",
        "thickness_region_target" : "",
        "output"                  : ""
    })"));

    const mio::Mesh source = Internals::ModelPartToMesh(rSource);
    const mio::Mesh target = Internals::ModelPartToMesh(rTarget);
    mio::ThicknessResult result = mio::thickness(source, target, ReadThicknessOptions(Settings));
    Parameters report = ThicknessReport(result);
    RenameResultArray(result.mMesh, mio::DataLocation::Cell, mio::kThicknessDistanceName,
                      Settings["output"].GetString());
    StoreResult(result.mMesh, rDestination);
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::Primitive(
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    // Every generator keeps its own upstream defaults (a disk's "n" is 8, a ball's 4), so a
    // setting left out is not filled from one shared table: it falls back per shape.
    static const std::set<std::string> accepted_keys = {
        "primitive", "dims", "size", "origin", "center", "cell_type", "split", "order", "complete",
        "max_cells", "radius", "inner_radius", "outer_radius", "major_radius", "minor_radius",
        "height", "width", "edge", "n", "nu", "nv", "nr", "nz", "nt", "nxy", "ntheta", "ncap",
        "subdivisions", "sides", "rings", "angle", "pitch", "turns", "num_twists", "ratio", "seed",
        "start", "end", "name", "types", "pairs", "loops", "target_size", "min_angle", "max_points"};
    for (auto it = Settings.begin(); it != Settings.end(); ++it) {
        KRATOS_ERROR_IF(accepted_keys.count(it.name()) == 0)
            << "Unknown \"" << it.name() << "\" setting for MeshioPlusPlusMeshOperations.Primitive" << std::endl;
    }
    KRATOS_ERROR_IF_NOT(Settings.Has("primitive"))
        << "MeshioPlusPlusMeshOperations.Primitive needs a \"primitive\" (box, disk, ball, ...)" << std::endl;

    const auto real = [&Settings](const char* pKey, const double Default) {
        return Settings.Has(pKey) ? Settings[pKey].GetDouble() : Default;
    };
    const auto integer = [&Settings](const char* pKey, const std::int64_t Default) {
        return Settings.Has(pKey) ? static_cast<std::int64_t>(Settings[pKey].GetInt()) : Default;
    };
    const auto text = [&Settings](const char* pKey, const std::string& rDefault) {
        return Settings.Has(pKey) ? Settings[pKey].GetString() : rDefault;
    };
    const auto vector3 = [&Settings](const char* pKey) {
        std::array<double, 3> value{{0.0, 0.0, 0.0}};
        return Settings.Has(pKey) ? ReadVector3(Settings, pKey) : value;
    };
    const auto pair_of = [&Settings](const char* pKey, const std::array<double, 2>& rDefault) {
        if (!Settings.Has(pKey)) {
            return rDefault;
        }
        const Vector value = Settings[pKey].GetVector();
        KRATOS_ERROR_IF(value.size() != 2) << "\"" << pKey << "\" needs two entries" << std::endl;
        return std::array<double, 2>{{value[0], value[1]}};
    };
    const auto triple_of = [&Settings](const char* pKey, const std::array<double, 3>& rDefault) {
        return Settings.Has(pKey) ? ReadVector3(Settings, pKey) : rDefault;
    };
    const auto counts = [&Settings](const char* pKey, const std::size_t Size, const std::int64_t Default) {
        std::array<std::int64_t, 3> value{{Default, Default, Default}};
        if (Settings.Has(pKey)) {
            const Vector read = Settings[pKey].GetVector();
            KRATOS_ERROR_IF(read.size() != Size) << "\"" << pKey << "\" needs " << Size << " entries" << std::endl;
            for (std::size_t i = 0; i < Size; ++i) {
                value[i] = static_cast<std::int64_t>(read[i]);
            }
        }
        return value;
    };

    const std::string name = Settings["primitive"].GetString();
    const int order = static_cast<int>(integer("order", 1));
    const std::int64_t max_cells = integer("max_cells", 20000000);
    const bool complete = Settings.Has("complete") && Settings["complete"].GetBool();
    const auto center = vector3("center");
    const auto origin = vector3("origin");

    mio::Mesh result;
    if (name == "rectangle") {
        const auto dims = counts("dims", 2, 1);
        result = mio::rectangle({{dims[0], dims[1]}}, pair_of("size", {{1.0, 1.0}}), origin,
                                text("cell_type", "quad"), text("split", "zigzag"), order, max_cells, complete);
    } else if (name == "box") {
        result = mio::box(counts("dims", 3, 1), triple_of("size", {{1.0, 1.0, 1.0}}), origin,
                          text("cell_type", "hexahedron"), order, max_cells, complete);
    } else if (name == "box_surface") {
        result = mio::box_surface(counts("dims", 3, 1), triple_of("size", {{1.0, 1.0, 1.0}}), origin,
                                  text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "disk") {
        result = mio::disk(real("radius", 1.0), integer("n", 8), center, text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "uv_sphere") {
        result = mio::uv_sphere(real("radius", 1.0), integer("nu", 16), integer("nv", 8), center,
                                text("cell_type", "triangle"), order, max_cells, complete);
    } else if (name == "icosa_sphere") {
        result = mio::icosa_sphere(real("radius", 1.0), integer("subdivisions", 2), center, order, max_cells, complete);
    } else if (name == "tetra_sphere") {
        result = mio::tetra_sphere(real("radius", 1.0), integer("subdivisions", 3), center, order, max_cells, complete);
    } else if (name == "octa_sphere") {
        result = mio::octa_sphere(real("radius", 1.0), integer("subdivisions", 2), center, order, max_cells, complete);
    } else if (name == "cube_sphere") {
        result = mio::cube_sphere(real("radius", 1.0), integer("n", 4), center, text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "geo_sphere") {
        result = mio::geo_sphere(real("radius", 1.0), integer("n", 4), center, order, max_cells, complete);
    } else if (name == "cylinder") {
        result = mio::cylinder(real("radius", 1.0), real("height", 1.0), integer("nxy", 8), integer("nz", 8), center,
                               text("cell_type", "hexahedron"), order, max_cells, complete);
    } else if (name == "reference_triangle") {
        result = mio::reference_triangle(integer("n", 4), pair_of("size", {{1.0, 1.0}}), origin, order, max_cells, complete);
    } else if (name == "ngon") {
        result = mio::ngon(integer("sides", 6), real("radius", 1.0), integer("rings", 1), center, order, max_cells, complete);
    } else if (name == "annulus") {
        result = mio::annulus(real("inner_radius", 0.5), real("outer_radius", 1.0), integer("ntheta", 16), integer("nr", 4),
                              real("angle", 2.0 * Globals::Pi), center, text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "lshape") {
        result = mio::lshape(integer("n", 4), real("size", 2.0), center, text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "plate_with_hole") {
        result = mio::plate_with_hole(real("radius", 0.25), real("size", 1.0), integer("nt", 8), integer("nr", 4), center,
                                      text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "honeycomb") {
        const auto dims = counts("dims", 2, 4);
        result = mio::honeycomb({{dims[0], dims[1]}}, real("edge", 1.0), origin, max_cells);
    } else if (name == "line") {
        result = mio::line(integer("n", 8), triple_of("start", {{0.0, 0.0, 0.0}}), triple_of("end", {{1.0, 0.0, 0.0}}),
                           order, max_cells, complete);
    } else if (name == "circle") {
        result = mio::circle(real("radius", 1.0), integer("n", 16), center, order, max_cells, complete);
    } else if (name == "helix") {
        result = mio::helix(real("radius", 1.0), real("pitch", 1.0), real("turns", 2.0), integer("n", 32), center,
                            order, max_cells, complete);
    } else if (name == "sphere_points") {
        result = mio::sphere_points(integer("n", 100), real("radius", 1.0), center, max_cells);
    } else if (name == "box_points") {
        result = mio::box_points(integer("n", 100), triple_of("size", {{1.0, 1.0, 1.0}}), origin,
                                 static_cast<std::uint64_t>(integer("seed", 0)), max_cells);
    } else if (name == "torus") {
        result = mio::torus(real("major_radius", 1.0), real("minor_radius", 0.25), integer("nu", 32), integer("nv", 16),
                            center, text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "moebius") {
        result = mio::moebius(real("radius", 1.0), real("width", 0.5), integer("nu", 32), integer("nv", 4),
                              integer("num_twists", 1), center, text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "tube") {
        result = mio::tube(real("radius", 1.0), real("height", 1.0), integer("ntheta", 16), integer("nz", 4), center,
                           text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "cone") {
        result = mio::cone(real("radius", 1.0), real("height", 1.0), integer("ntheta", 16), integer("nz", 4), center,
                           text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "capsule") {
        result = mio::capsule(real("radius", 0.5), real("height", 1.0), integer("ntheta", 16), integer("nz", 4),
                              integer("ncap", 4), center, text("cell_type", "quad"), order, max_cells, complete);
    } else if (name == "platonic") {
        result = mio::platonic(text("name", "icosahedron"), real("radius", 1.0), center, max_cells);
    } else if (name == "ball") {
        result = mio::ball(real("radius", 1.0), integer("n", 4), center, text("cell_type", "hexahedron"), order, max_cells, complete);
    } else if (name == "spherical_shell") {
        result = mio::spherical_shell(real("inner_radius", 0.5), real("outer_radius", 1.0), integer("n", 4), integer("nr", 2),
                                      real("ratio", 1.0), center, order, max_cells, complete);
    } else if (name == "pipe") {
        result = mio::pipe(real("inner_radius", 0.5), real("outer_radius", 1.0), real("height", 1.0), integer("ntheta", 16),
                           integer("nr", 2), integer("nz", 4), real("ratio", 1.0), center, order, max_cells, complete);
    } else if (name == "fichera") {
        result = mio::fichera(integer("n", 2), real("size", 2.0), center, text("cell_type", "hexahedron"), order, max_cells, complete);
    } else if (name == "solid_torus") {
        result = mio::solid_torus(real("major_radius", 1.0), real("minor_radius", 0.25), integer("n", 4), integer("nu", 32),
                                  center, order, max_cells, complete);
    } else if (name == "reference_cells") {
        result = mio::reference_cells(Settings.Has("types") ? Settings["types"].GetStringArray() : std::vector<std::string>{},
                                      Settings.Has("pairs") && Settings["pairs"].GetBool());
    } else if (name == "polygon") {
        // A constrained Delaunay triangulation of planar loops: the first is the outer boundary,
        // the rest holes, each a list of [x, y] points.
        KRATOS_ERROR_IF_NOT(Settings.Has("loops"))
            << "The \"polygon\" primitive needs \"loops\": a list of [[x, y], ...] point lists" << std::endl;
        std::vector<std::vector<std::array<double, 2>>> loops;
        for (std::size_t i_loop = 0; i_loop < Settings["loops"].size(); ++i_loop) {
            std::vector<std::array<double, 2>> loop;
            for (std::size_t i_point = 0; i_point < Settings["loops"][i_loop].size(); ++i_point) {
                const Vector point = Settings["loops"][i_loop][i_point].GetVector();
                KRATOS_ERROR_IF(point.size() != 2) << "Every \"loops\" point needs two coordinates" << std::endl;
                loop.push_back({{point[0], point[1]}});
            }
            loops.push_back(std::move(loop));
        }
        mio::TriangulatePolygon2dOptions options;
        options.mTargetSize = real("target_size", 0.0);
        options.mMinAngleDeg = real("min_angle", 0.0);
        options.mMaxPoints = integer("max_points", 1000000);
        result = mio::triangulate_polygon_2d(loops, options);
    } else {
        KRATOS_ERROR << "Unknown primitive \"" << name << "\". Supported primitives are: annulus, ball, box, "
                     << "box_points, box_surface, capsule, circle, cone, cube_sphere, cylinder, disk, fichera, "
                     << "geo_sphere, helix, honeycomb, icosa_sphere, line, lshape, moebius, ngon, octa_sphere, "
                     << "pipe, plate_with_hole, platonic, polygon, rectangle, reference_cells, "
                     << "reference_triangle, solid_torus, sphere_points, spherical_shell, tetra_sphere, "
                     << "torus, tube, uv_sphere" << std::endl;
    }
    Internals::MeshToModelPart(result, rDestination);

    Parameters report(R"({})");
    report.AddString("primitive", name);
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::DistanceToSurface(
    const ModelPart& rQuery,
    const ModelPart& rSurface,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rQuery.IsDistributed() || rSurface.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    Settings.AddMissingParameters(GetDefaultParameters());

    const bool deformed = Settings["use_deformed_configuration"].GetBool();
    const Internals::FieldDataSelection selection = BuildFieldDataSelection(Settings);
    mio::Mesh query = Internals::ModelPartToMeshWithData(rQuery, true, true, deformed, selection);
    mio::Mesh surface = Internals::ModelPartToMeshWithData(rSurface, true, true, deformed, selection);

    mio::SurfaceDistanceResult result =
        mio::distance_to_surface(query, surface, BuildSurfaceDistanceOptions(Settings));

    RenameResultArray(result.mMesh, mio::DataLocation::Point, mio::kSdfDistanceName,
                      Settings["output"].GetString());
    Internals::MeshToModelPart(result.mMesh, rDestination);

    Parameters report(R"({})");
    report.AddInt("number_of_banded", static_cast<int>(result.mNumBanded));
    report.AddValue("surface_quality", SurfaceQualityReport(result.mQuality));
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::CheckSurfaceWatertight(const ModelPart& rSurface)
{
    KRATOS_TRY

    mio::Mesh surface = Internals::ModelPartToMesh(rSurface);
    return SurfaceQualityReport(mio::surface_watertight_check(surface));

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::Diff(
    const ModelPart& rFirst,
    const ModelPart& rSecond,
    Parameters Settings
    )
{
    KRATOS_TRY

    Settings.AddMissingParameters(Parameters(R"({
        "absolute_tolerance" : 1e-12,
        "relative_tolerance" : 1e-9
    })"));

    mio::Mesh first = Internals::ModelPartToMesh(rFirst);
    mio::Mesh second = Internals::ModelPartToMesh(rSecond);

    const double atol = Settings["absolute_tolerance"].GetDouble();
    const double rtol = Settings["relative_tolerance"].GetDouble();

    Parameters report(R"({})");
    report.AddBool("equal", mio::meshes_equal(first, second, atol, rtol));
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::ComputeStatistics(const ModelPart& rSource)
{
    KRATOS_TRY

    mio::Mesh mesh = Internals::ModelPartToMesh(rSource);
    const mio::StatsReport stats = mio::compute_stats(mesh);

    Parameters report(R"({})");
    report.AddInt("number_of_points", static_cast<int>(stats.mNumPoints));
    report.AddInt("number_of_cells", static_cast<int>(stats.mNumCells));
    report.AddDouble("total_area", stats.mTotalArea);
    report.AddDouble("signed_volume", stats.mSignedVolume);
    report.AddDouble("unsigned_volume", stats.mUnsignedVolume);
    report.AddInt("number_of_inverted", static_cast<int>(stats.mNumInverted));

    Parameters counts(R"({})");
    for (const auto& r_entry : stats.mCellTypeCounts) {
        counts.AddInt(r_entry.first, static_cast<int>(r_entry.second));
    }
    report.AddValue("cell_type_counts", counts);
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::ComputeQuality(const ModelPart& rSource)
{
    KRATOS_TRY

    mio::Mesh mesh = Internals::ModelPartToMesh(rSource);
    const mio::QualityReport quality = mio::compute_quality(mesh);

    Parameters report(R"({})");
    report.AddInt("number_of_cells", static_cast<int>(quality.mNumCells));
    report.AddInt("number_of_inverted", static_cast<int>(quality.mNumInverted));
    report.AddInt("number_of_degenerate", static_cast<int>(quality.mNumDegenerate));

    Parameters metrics(R"({})");
    for (const auto& r_entry : quality.mMetrics) {
        Parameters summary(R"({})");
        summary.AddDouble("min", r_entry.second.mMin);
        summary.AddDouble("max", r_entry.second.mMax);
        summary.AddDouble("mean", r_entry.second.mMean);
        summary.AddInt("count", static_cast<int>(r_entry.second.mCount));
        metrics.AddValue(r_entry.first, summary);
    }
    report.AddValue("metrics", metrics);
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t MeshioPlusPlusMeshOperations::ComputeBandwidth(const ModelPart& rSource)
{
    KRATOS_TRY

    mio::Mesh mesh = Internals::ModelPartToMesh(rSource);
    return static_cast<std::size_t>(mio::compute_bandwidth(mesh));

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::pair<std::size_t, std::size_t>> MeshioPlusPlusMeshOperations::NeighborPairs(
    const ModelPart& rSource,
    Parameters Settings
    )
{
    KRATOS_TRY

    Settings.ValidateAndAssignDefaults(Parameters(R"({
        "method"                     : "radius",
        "radius"                     : 0.0,
        "k"                          : 0,
        "box"                        : [],
        "cell_size"                  : 0.0,
        "use_deformed_configuration" : false
    })"));

    mio::NeighborOptions options;
    const std::string method = Settings["method"].GetString();
    if (method == "radius") {
        options.mMethod = mio::NeighborMethod::Radius;
        options.mRadius = Settings["radius"].GetDouble();
        KRATOS_ERROR_IF_NOT(options.mRadius > 0.0)
            << "NeighborPairs: \"radius\" must be positive for the \"radius\" method, got "
            << options.mRadius << std::endl;
    } else if (method == "k_nearest") {
        options.mMethod = mio::NeighborMethod::KNearest;
        options.mK = Settings["k"].GetInt();
        KRATOS_ERROR_IF_NOT(options.mK > 0)
            << "NeighborPairs: \"k\" must be positive for the \"k_nearest\" method, got "
            << options.mK << std::endl;
    } else {
        KRATOS_ERROR << "NeighborPairs: unknown \"method\" \"" << method
                     << "\" (use \"radius\" or \"k_nearest\")" << std::endl;
    }
    const Vector box = Settings["box"].GetVector();
    options.mBox.assign(box.begin(), box.end());
    options.mCellSize = Settings["cell_size"].GetDouble();

    // The node coordinates as meshio++'s (N, 3) point array, in container order; the returned
    // row indices are mapped back to node Ids through the same order.
    const std::size_t number_of_nodes = rSource.NumberOfNodes();
    if (number_of_nodes == 0) {
        return {};
    }
    const bool deformed = Settings["use_deformed_configuration"].GetBool();
    mio::NDArray points = mio::NDArray::Uninit(mio::DType::Float64, {number_of_nodes, 3});
    double* p_points = points.As<double>();
    std::vector<std::size_t> node_ids(number_of_nodes);
    std::size_t row = 0;
    for (const auto& r_node : rSource.Nodes()) {
        const auto& r_coordinates = deformed ? r_node.Coordinates() : r_node.GetInitialPosition().Coordinates();
        p_points[3 * row] = r_coordinates[0];
        p_points[3 * row + 1] = r_coordinates[1];
        p_points[3 * row + 2] = r_coordinates[2];
        node_ids[row] = r_node.Id();
        ++row;
    }

    const mio::NeighborPairs pairs = mio::neighbor_pairs(points, options);
    const std::size_t number_of_pairs = pairs.mSource.Size();
    const std::int64_t* p_source = pairs.mSource.As<std::int64_t>();
    const std::int64_t* p_target = pairs.mTarget.As<std::int64_t>();
    std::vector<std::pair<std::size_t, std::size_t>> result(number_of_pairs);
    for (std::size_t i = 0; i < number_of_pairs; ++i) {
        result[i] = {node_ids[static_cast<std::size_t>(p_source[i])],
                     node_ids[static_cast<std::size_t>(p_target[i])]};
    }
    return result;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::HausdorffDistance(
    const ModelPart& rFirst,
    const ModelPart& rSecond,
    Parameters Settings
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rFirst.IsDistributed() || rSecond.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    Settings.ValidateAndAssignDefaults(Parameters(R"({
        "face_samples"   : 0,
        "region_first"   : "",
        "region_second"  : "",
        "grid_cell_size" : 0.0
    })"));

    mio::HausdorffOptions options;
    options.mFaceSamples = Settings["face_samples"].GetInt();
    options.mRegionA = Settings["region_first"].GetString();
    options.mRegionB = Settings["region_second"].GetString();
    options.mGridCellSize = Settings["grid_cell_size"].GetDouble();

    const mio::Mesh first = Internals::ModelPartToMesh(rFirst);
    const mio::Mesh second = Internals::ModelPartToMesh(rSecond);
    const mio::HausdorffResult result = mio::hausdorff_distance(first, second, options);

    Parameters report(R"({})");
    report.AddDouble("distance", result.mDistance);
    report.AddDouble("first_to_second", result.mAtoB);
    report.AddDouble("second_to_first", result.mBtoA);
    report.AddDouble("mean_first_to_second", result.mMeanAtoB);
    report.AddDouble("rms_first_to_second", result.mRmsAtoB);
    report.AddDouble("mean_second_to_first", result.mMeanBtoA);
    report.AddDouble("rms_second_to_first", result.mRmsBtoA);
    report.AddInt("samples_first", static_cast<int>(result.mNumSamplesA));
    report.AddInt("samples_second", static_cast<int>(result.mNumSamplesB));
    Vector worst_first(3), worst_second(3);
    for (std::size_t i = 0; i < 3; ++i) {
        worst_first[i] = result.mWorstPointA[i];
        worst_second[i] = result.mWorstPointB[i];
    }
    report.AddVector("worst_point_first", worst_first);
    report.AddVector("worst_point_second", worst_second);
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::pair<std::size_t, std::size_t>> MeshioPlusPlusMeshOperations::MatchPeriodicNodes(
    const ModelPart& rSource,
    Parameters Settings
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSource.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    // The selectors are free-form sub-objects, which a recursive validation against "{}" would
    // reject, so the defaults are merged rather than validated.
    Settings.AddMissingParameters(Parameters(R"({
        "slave"                      : {},
        "master"                     : {},
        "translation"                : [0.0, 0.0, 0.0],
        "rotation_axis"              : [0.0, 0.0, 1.0],
        "rotation_angle"             : 0.0,
        "matrix"                     : [],
        "tolerance"                  : 1e-8,
        "require_complete"           : true,
        "use_deformed_configuration" : false
    })"));

    mio::PeriodicOptions options;
    options.mTransform = ReadAffineTransform(Settings);
    options.mAtol = Settings["tolerance"].GetDouble();
    options.mRequireComplete = Settings["require_complete"].GetBool();
    const mio::RegionSelector slave = ReadRegionSelector(Settings["slave"], "MatchPeriodicNodes");
    const mio::RegionSelector master = ReadRegionSelector(Settings["master"], "MatchPeriodicNodes");

    const mio::Mesh mesh = Internals::ModelPartToMesh(rSource, true, true,
                                                      Settings["use_deformed_configuration"].GetBool());
    const mio::PeriodicPairs pairs = mio::match_periodic_nodes(mesh, slave, master, options);

    // The pair rows are point rows, which follow the node container order of the model part.
    std::vector<std::size_t> node_ids;
    node_ids.reserve(rSource.NumberOfNodes());
    for (const auto& r_node : rSource.Nodes()) {
        node_ids.push_back(r_node.Id());
    }
    const std::int64_t* p_slave = pairs.mSlave.As<std::int64_t>();
    const std::int64_t* p_master = pairs.mMaster.As<std::int64_t>();
    std::vector<std::pair<std::size_t, std::size_t>> result;
    result.reserve(pairs.mSlave.Size());
    for (std::size_t i = 0; i < pairs.mSlave.Size(); ++i) {
        result.emplace_back(node_ids[static_cast<std::size_t>(p_slave[i])],
                            node_ids[static_cast<std::size_t>(p_master[i])]);
    }
    return result;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::BlendSteps(
    const ModelPart& rFirst,
    const ModelPart& rSecond,
    const double Weight,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rFirst.IsDistributed() || rSecond.IsDistributed())
        << "The meshio++ operations do not support distributed model parts" << std::endl;

    Settings.AddMissingParameters(GetDefaultParameters());
    Settings.AddMissingParameters(Parameters(R"({"blend_points" : false})"));

    const bool deformed = Settings["use_deformed_configuration"].GetBool();
    const Internals::FieldDataSelection selection = BuildFieldDataSelection(Settings);
    const mio::Mesh first = Internals::ModelPartToMeshWithData(rFirst, true, true, deformed, selection);
    const mio::Mesh second = Internals::ModelPartToMeshWithData(rSecond, true, true, deformed, selection);

    mio::BlendOptions options;
    options.mBlendPoints = Settings["blend_points"].GetBool();
    mio::Mesh result = mio::blend_steps(first, second, Weight, options);
    Internals::MeshToModelPart(result, rDestination);

    Parameters report(R"({})");
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    return report;

    KRATOS_CATCH("")
}


} // namespace Kratos
