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
#include <array>
#include <cmath>

// External includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/skin.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/crop.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_calc.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/data_condition.hpp"
#include "meshioplusplus/operations/data_info.hpp"
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/transform.hpp"

// Project includes
#include "custom_utilities/meshioplusplus_mesh_operations.h"
#include "custom_utilities/meshioplusplus_conversion_utilities.h"
#include "containers/model.h"

namespace mio = meshioplusplus;

namespace Kratos
{
namespace
{

const std::vector<std::string>& OperationNames()
{
    static const std::vector<std::string> names = {
        "attach_quality", "cell_data_to_point_data", "clean", "convert_cells", "crop_bbox",
        "crop_halfspace", "data_calc", "data_condition", "data_info", "data_manage",
        "decimate", "extract_skin", "extract_surface", "isosurface", "partition",
        "point_data_to_cell_data", "quality", "refine", "reorder", "slice", "smooth",
        "split", "stats", "transform"
    };
    return names;
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
    return selection;
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

} // namespace

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> MeshioPlusPlusMeshOperations::GetSupportedOperations()
{
    return OperationNames();
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusMeshOperations::GetDefaultParameters()
{
    return Parameters(R"({
        "operation"                                   : "clean",
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
        "on_conflict"                                  : "error"
    })");
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

    } else if (operation == "refine") {
        mio::RefineOptions options;
        options.mLevels = Settings["levels"].GetInt();
        options.mRecordParentIds = Settings["record_parent_ids"].GetBool();
        mio::RefineResult result = mio::refine(mesh, options);
        StoreResult(result.mMesh, rDestination);

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

} // namespace Kratos
