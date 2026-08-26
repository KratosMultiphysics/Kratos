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

// External includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/skin.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/conservative_interpolate.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/crop.hpp"
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
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/remesh.hpp"
#include "meshioplusplus/operations/remesh_volume.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/subdivide.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

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
        "agglomerate", "attach_quality", "cell_data_to_point_data", "clean", "compute_sdf",
        "convert_cells", "crop_bbox", "crop_halfspace", "crop_predicate", "data_calc",
        "data_condition", "data_info", "data_integrate", "data_manage", "decimate",
        "decimate_volume", "estimate_error", "extract_skin", "extract_surface", "gradient",
        "hessian", "isosurface", "optimize_volume", "partition", "point_data_to_cell_data",
        "quality", "refine", "remesh", "remesh_volume", "reorder", "slice", "smooth", "split",
        "stats", "subdivide", "transform", "voxelize"
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
    const Vector bounds = rSettings["bounds"].GetVector();
    if (bounds.size() > 0) {
        KRATOS_ERROR_IF(bounds.size() != 6)
            << "\"bounds\" needs six entries (xlo, ylo, zlo, xhi, yhi, zhi), got "
            << bounds.size() << std::endl;
        std::array<double, 6> value{};
        for (std::size_t i = 0; i < 6; ++i) {
            value[i] = bounds[i];
        }
        rOptions.mBounds = value;
    }
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

        "write_mdpa_ids"                               : false
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
        mio::AgglomerateResult result = mio::agglomerate(mesh, options);
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
        mio::OptimizeVolumeResult result = mio::optimize_volume(mesh, options);
        report.AddInt("number_of_flips", static_cast<int>(result.mNumFlips));
        report.AddInt("number_of_2_3_flips", static_cast<int>(result.mNum23Flips));
        report.AddInt("number_of_3_2_flips", static_cast<int>(result.mNum32Flips));
        report.AddInt("number_of_vertices_moved", static_cast<int>(result.mNumVerticesMoved));
        report.AddInt("number_of_tets", static_cast<int>(result.mNumTets));
        report.AddDouble("min_quality_before", result.mMinQualityBefore);
        report.AddDouble("min_quality_after", result.mMinQualityAfter);
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

    } else if (operation == "voxelize") {
        mio::VoxelOptions options;
        ReadLatticeSettings(Settings, options);
        options.mFill = mio::voxel_fill_from_name(Settings["fill"].GetString());
        options.mAttachOccupancy = Settings["attach_occupancy"].GetBool();
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

} // namespace Kratos
