//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
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
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/diff.hpp"
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

/// The operations Execute() dispatches on, and whether each produces a mesh.
const std::vector<std::string>& OperationNames()
{
    static const std::vector<std::string> names = {
        "attach_quality", "clean", "convert_cells", "crop_bbox", "crop_halfspace",
        "decimate", "extract_skin", "extract_surface", "isosurface", "partition",
        "quality", "refine", "reorder", "slice", "smooth", "split", "stats", "transform"
    };
    return names;
}

/// Reads a 3-component vector setting, defaulting each missing component to zero.
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

/// Writes a mesh-producing operation's result into the destination model part.
void StoreResult(mio::Mesh& rMesh, ModelPart& rDestination)
{
    Internals::MeshToModelPart(rMesh, rDestination);
}

/// Writes one model part per piece, for the multi-output operations.
/// @note The pieces cannot be sub model parts of a shared destination: a Kratos sub model
/// part is a view into its root's entity containers, and every piece is independently
/// renumbered from 1, so their nodes would collide on insertion. They are therefore created
/// as siblings of the destination, named after it, and the destination itself stays empty.
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
        "operation"                    : "clean",
        "entity_type"                  : "automatic",
        "use_deformed_configuration"   : false,
        "weld"                         : false,
        "tolerance"                    : 1e-8,
        "remove_orphans"               : true,
        "drop_degenerate"              : true,
        "drop_duplicate_cells"         : true,
        "translation"                  : [0.0, 0.0, 0.0],
        "scale"                        : [1.0, 1.0, 1.0],
        "rotation_axis"                : [0.0, 0.0, 1.0],
        "rotation_angle"               : 0.0,
        "rotate_vector_data"           : false,
        "mode"                         : "linearize",
        "levels"                       : 1,
        "target_ratio"                 : -1.0,
        "target_faces"                 : -1,
        "max_error"                    : -1.0,
        "preserve_boundary"            : true,
        "preserve_features"            : true,
        "feature_angle"                : 30.0,
        "method"                       : "taubin",
        "iterations"                   : 10,
        "lambda"                       : -1.0,
        "mu"                           : -0.34,
        "fix_boundary"                 : true,
        "guard_inversion"              : true,
        "number_of_parts"              : 2,
        "imbalance"                    : 0.03,
        "seed"                         : 0,
        "ghost_layers"                 : 0,
        "weights_variable"             : "",
        "split_by"                     : "type",
        "tag_name"                     : "",
        "box_min"                      : [0.0, 0.0, 0.0],
        "box_max"                      : [1.0, 1.0, 1.0],
        "origin"                       : [0.0, 0.0, 0.0],
        "normal"                       : [0.0, 0.0, 1.0],
        "keep_partial_cells"           : false,
        "array_name"                   : "",
        "isovalues"                    : [0.0],
        "linearize"                    : false,
        "record_parent_ids"            : false
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

    mio::Mesh mesh = Internals::ModelPartToMesh(rSource, write_elements, write_conditions,
                                                Settings["use_deformed_configuration"].GetBool());

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
