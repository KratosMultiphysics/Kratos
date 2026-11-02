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
#include <cstdint>
#include <string>
#include <utility>
#include <variant>
#include <vector>

// External includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/nurbs.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/iges.hpp"
#include "meshioplusplus/formats/mfem.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/operations/iga.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/project_to_cad.hpp"

// Project includes
#include "custom_utilities/meshioplusplus_iga_utilities.h"
#include "custom_utilities/meshioplusplus_conversion_utilities.h"

namespace mio = meshioplusplus;

namespace Kratos
{
namespace
{

/// A model read from a file: a B-rep (cad_json, step, iges) or a NURBS mesh (mfem).
using IgaModel = std::variant<mio::BrepInfo, mio::NurbsInfo>;

/***********************************************************************************/
/***********************************************************************************/

/// The meshio++ format name of a model file ("format" : "auto" resolves it from the name).
std::string ResolveModelFormat(
    const std::string& rFileName,
    const Parameters& rSettings
    )
{
    const std::string format = rSettings["format"].GetString();
    const std::string name = mio::resolve_format(rFileName, format == "auto" ? std::string() : format);
    KRATOS_ERROR_IF(name != "cad_json" && name != "step" && name != "iges" && name != "mfem")
        << "\"" << rFileName << "\" resolves to format \"" << name << "\", which carries no NURBS model "
        << "(expected cad_json, step, iges or mfem)" << std::endl;
    return name;
}

/***********************************************************************************/
/***********************************************************************************/

/// Reads the model and, when pMesh is not null, the mesh it was sampled into by the reader.
IgaModel ReadModel(
    const std::string& rFileName,
    const Parameters& rSettings,
    mio::Mesh* pMesh = nullptr
    )
{
    const std::string format = ResolveModelFormat(rFileName, rSettings);
    mio::Mesh mesh;
    if (format == "mfem") {
        mio::NurbsInfo info;
        mesh = mio::read_mfem(rFileName, {}, mio::ReadOptions{}, info);
        KRATOS_ERROR_IF(info.mPatches.empty())
            << "The MFEM mesh \"" << rFileName << "\" is not a NURBS mesh" << std::endl;
        if (pMesh != nullptr) {
            *pMesh = std::move(mesh);
        }
        return info;
    }
    mio::BrepInfo info;
    if (format == "cad_json") {
        mesh = mio::read_cad_json(rFileName, info);
    } else {
        mio::CadReadOptions options;
        options.mDeflection = rSettings["cad_deflection"].GetDouble();
        options.mAngularDeflection = rSettings["cad_angular_deflection"].GetDouble();
        mesh = format == "step" ? mio::read_step(rFileName, info, options)
                                : mio::read_iges(rFileName, info, options);
    }
    if (pMesh != nullptr) {
        *pMesh = std::move(mesh);
    }
    return info;
}

/***********************************************************************************/
/***********************************************************************************/

/// The patches of either kind of model.
const std::vector<mio::NurbsPatch>& Patches(const IgaModel& rModel)
{
    return std::holds_alternative<mio::BrepInfo>(rModel)
        ? std::get<mio::BrepInfo>(rModel).mPatches : std::get<mio::NurbsInfo>(rModel).mPatches;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters DescribeModel(const IgaModel& rModel)
{
    Parameters result(R"({})");
    const bool is_brep = std::holds_alternative<mio::BrepInfo>(rModel);
    result.AddString("kind", is_brep ? "brep" : "nurbs");
    result.AddInt("number_of_patches", static_cast<int>(Patches(rModel).size()));
    int number_of_trimmed_faces = 0;
    if (is_brep) {
        const auto& r_info = std::get<mio::BrepInfo>(rModel);
        for (const auto& r_face : r_info.mFaces) {
            number_of_trimmed_faces += r_face.mTrimmed ? 1 : 0;
        }
        result.AddInt("number_of_breps", static_cast<int>(r_info.mBreps.size()));
        result.AddInt("number_of_faces", static_cast<int>(r_info.mFaces.size()));
        result.AddInt("number_of_edges", static_cast<int>(r_info.mEdges.size()));
        result.AddInt("number_of_vertices", static_cast<int>(r_info.mVertices.size()));
    } else {
        result.AddInt("number_of_breps", 0);
        result.AddInt("number_of_faces", 0);
        result.AddInt("number_of_edges", 0);
        result.AddInt("number_of_vertices", 0);
    }
    result.AddInt("number_of_trimmed_faces", number_of_trimmed_faces);
    result.AddEmptyArray("patches");
    for (const auto& r_patch : Patches(rModel)) {
        Parameters patch(R"({})");
        std::vector<int> degrees;
        for (const auto degree : r_patch.mDegrees) {
            degrees.push_back(static_cast<int>(degree));
        }
        patch.AddEmptyArray("degrees");
        for (const int degree : degrees) {
            patch["degrees"].Append(degree);
        }
        patch.AddInt("number_of_control_points", static_cast<int>(r_patch.mWeights.Size()));
        bool rational = false;
        const double* p_weights = r_patch.mWeights.As<double>();
        for (std::size_t i = 0; i < r_patch.mWeights.Size(); ++i) {
            rational = rational || p_weights[i] != 1.0;
        }
        patch.AddBool("rational", rational);
        result["patches"].Append(patch);
    }
    return result;
}

/***********************************************************************************/
/***********************************************************************************/

/// Reads a 3-integer list setting.
std::array<std::int32_t, 3> ReadTriple(
    const Parameters& rSettings,
    const std::string& rKey,
    const std::int32_t Minimum
    )
{
    std::array<std::int32_t, 3> value{{Minimum, Minimum, Minimum}};
    KRATOS_ERROR_IF_NOT(rSettings[rKey].IsArray() && rSettings[rKey].size() <= 3)
        << "\"" << rKey << "\" must be a list of up to 3 integers" << std::endl;
    for (std::size_t i = 0; i < rSettings[rKey].size(); ++i) {
        value[i] = rSettings[rKey][i].GetInt();
        KRATOS_ERROR_IF(value[i] < Minimum)
            << "\"" << rKey << "\" entries must be at least " << Minimum << std::endl;
    }
    return value;
}

/***********************************************************************************/
/***********************************************************************************/

/// An N x Columns matrix of the Float64 array.
Matrix ToMatrix(
    const mio::NDArray& rArray,
    const std::size_t Rows,
    const std::size_t Columns
    )
{
    Matrix result(Rows, Columns);
    const double* p_data = rArray.As<double>();
    for (std::size_t i = 0; i < Rows; ++i) {
        for (std::size_t j = 0; j < Columns; ++j) {
            result(i, j) = p_data[i * Columns + j];
        }
    }
    return result;
}

/***********************************************************************************/
/***********************************************************************************/

/// The (face or edge id, target size) pairs of a "face_sizes"/"edge_sizes" list of [id, size].
std::vector<std::pair<std::int64_t, double>> ReadSizes(const Parameters& rSettings, const std::string& rKey)
{
    std::vector<std::pair<std::int64_t, double>> sizes;
    for (std::size_t i = 0; i < rSettings[rKey].size(); ++i) {
        const Vector entry = rSettings[rKey][i].GetVector();
        KRATOS_ERROR_IF(entry.size() != 2)
            << "Every \"" << rKey << "\" entry is [id, size]" << std::endl;
        sizes.emplace_back(static_cast<std::int64_t>(entry[0]), entry[1]);
    }
    return sizes;
}


} // namespace

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIgaUtilities::GetDefaultParameters()
{
    return Parameters(R"({
        "format"                 : "auto",
        "cad_deflection"         : 0.0,
        "cad_angular_deflection" : 0.5,
        "subdivisions"           : [1, 1, 1],
        "elevate"                : [0, 0, 0],
        "sample_subdivisions"    : 0,
        "trim_samples"           : 16,
        "trim"                   : true,
        "edges"                  : true,
        "index"                  : 0,
        "parameters"             : [[0.0, 0.0]],
        "derivatives"            : 0,

        "target_size"            : 0.0,
        "max_chordal_deviation"  : 0.0,
        "min_points_per_edge"    : 0,
        "geometric_tolerance"    : 0.0,
        "strict"                 : false,
        "max_points"             : 50000000,
        "min_angle"              : 0.0,
        "smooth_iterations"      : 0,
        "smooth_method"          : "laplacian",
        "grading"                : 0.0,
        "curvature_deviation"    : 0.0,
        "cells_across_gaps"      : 0.0,
        "metric"                 : "mean",
        "pole_rings"             : 0,
        "reparam_condition_number" : 0.0,
        "order"                  : 1,
        "face_sizes"             : [],
        "edge_sizes"             : [],
        "virtual_edges"          : false,
        "model_tolerances"       : false,
        "virtual_faces"          : false,
        "lop_corners"            : false,
        "untangle"               : false,
        "quads"                  : "none",
        "quad_pairs"             : "matching",
        "face_ids_variable"      : "",
        "edge_ids_variable"      : "",
        "background_size_variable"    : "",
        "background_hessian_variable" : "",
        "metric_error"           : 0.0,

        "max_distance"           : 0.0
    })");
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIgaUtilities::Describe(
    const std::string& rFileName,
    Parameters Settings
    )
{
    KRATOS_TRY

    Settings.AddMissingParameters(GetDefaultParameters());
    return DescribeModel(ReadModel(rFileName, Settings));

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIgaUtilities::Sample(
    const std::string& rFileName,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rDestination.IsDistributed())
        << "MeshioPlusPlusIgaUtilities does not support distributed model parts" << std::endl;
    KRATOS_ERROR_IF(rDestination.NumberOfNodes() > 0)
        << "The destination model part \"" << rDestination.FullName() << "\" must be empty" << std::endl;
    Settings.AddMissingParameters(GetDefaultParameters());

    mio::IgaSampleOptions options;
    options.mSubdivisions = Settings["sample_subdivisions"].GetInt();
    options.mTrimSamples = Settings["trim_samples"].GetInt();
    options.mTrim = Settings["trim"].GetBool();
    options.mEdges = Settings["edges"].GetBool();

    const IgaModel model = ReadModel(rFileName, Settings);
    mio::Mesh mesh = std::holds_alternative<mio::BrepInfo>(model)
        ? mio::iga_sample(std::get<mio::BrepInfo>(model), options)
        : mio::iga_sample(std::get<mio::NurbsInfo>(model), options);
    Internals::MeshToModelPart(mesh, rDestination);

    Parameters report(R"({})");
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    return report;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIgaUtilities::Refine(
    const std::string& rInputFileName,
    const std::string& rOutputFileName,
    Parameters Settings
    )
{
    KRATOS_TRY

    Settings.AddMissingParameters(GetDefaultParameters());

    mio::IgaRefineOptions options;
    options.mSubdivisions = ReadTriple(Settings, "subdivisions", 1);
    options.mElevate = ReadTriple(Settings, "elevate", 0);

    const IgaModel model = ReadModel(rInputFileName, Settings);
    IgaModel refined;
    if (std::holds_alternative<mio::BrepInfo>(model)) {
        mio::BrepInfo info = mio::iga_refine(std::get<mio::BrepInfo>(model), options);
        mio::write_cad_json(rOutputFileName, mio::sample_brep(info), info);
        refined = std::move(info);
    } else {
        mio::NurbsInfo info = mio::iga_refine(std::get<mio::NurbsInfo>(model), options);
        mio::write_mfem(rOutputFileName, mio::sample_mfem_nurbs(info), /*GridFunctions=*/false, info);
        refined = std::move(info);
    }
    return DescribeModel(refined);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIgaUtilities::Evaluate(
    const std::string& rFileName,
    Parameters Settings
    )
{
    KRATOS_TRY

    Settings.AddMissingParameters(GetDefaultParameters());

    const int index = Settings["index"].GetInt();
    const int derivatives = Settings["derivatives"].GetInt();
    KRATOS_ERROR_IF(derivatives < 0 || derivatives > 1) << "\"derivatives\" must be 0 or 1" << std::endl;
    const IgaModel model = ReadModel(rFileName, Settings);
    const bool is_brep = std::holds_alternative<mio::BrepInfo>(model);
    const std::size_t number_of_entities = is_brep ? std::get<mio::BrepInfo>(model).mFaces.size()
                                                   : std::get<mio::NurbsInfo>(model).mPatches.size();
    KRATOS_ERROR_IF(index < 0 || static_cast<std::size_t>(index) >= number_of_entities)
        << "\"index\" " << index << " is out of range: the model has " << number_of_entities
        << (is_brep ? " faces" : " patches") << std::endl;

    // The parameters, as the N x parametric dimension Float64 array meshio++ takes.
    const Matrix parameters = Settings["parameters"].GetMatrix();
    mio::NDArray params(mio::DType::Float64, {parameters.size1(), parameters.size2()});
    double* p_params = params.As<double>();
    for (std::size_t i = 0; i < parameters.size1(); ++i) {
        for (std::size_t j = 0; j < parameters.size2(); ++j) {
            p_params[i * parameters.size2() + j] = parameters(i, j);
        }
    }

    mio::IgaEvaluation evaluation;
    std::vector<std::string> field_names;
    if (is_brep) {
        evaluation = mio::iga_evaluate(std::get<mio::BrepInfo>(model), static_cast<std::size_t>(index), params,
                                       derivatives, Settings["trim_samples"].GetInt());
    } else {
        const auto& r_info = std::get<mio::NurbsInfo>(model);
        evaluation = mio::iga_evaluate(r_info, static_cast<std::size_t>(index), params, derivatives);
        field_names = r_info.mFieldNames;
    }

    const std::size_t n = parameters.size1();
    Parameters result(R"({})");
    const auto& r_point_shape = evaluation.mPoints.Shape();
    result.AddMatrix("points", ToMatrix(evaluation.mPoints, n, r_point_shape.size() > 1 ? r_point_shape[1] : 0));
    if (evaluation.mDerivatives.Size() > 0) {
        result.AddMatrix("derivatives", ToMatrix(evaluation.mDerivatives, n, evaluation.mDerivatives.Size() / std::max<std::size_t>(n, 1)));
    }
    if (evaluation.mNormals.Size() > 0) {
        result.AddMatrix("normals", ToMatrix(evaluation.mNormals, n, 3));
    }
    if (evaluation.mInside.Size() > 0) {
        Vector inside(n);
        const auto* p_inside = evaluation.mInside.As<std::uint8_t>();
        for (std::size_t i = 0; i < n; ++i) {
            inside[i] = static_cast<double>(p_inside[i]);
        }
        result.AddVector("inside", inside);
    }
    for (std::size_t i = 0; i < evaluation.mFields.size() && i < field_names.size(); ++i) {
        const auto& r_field = evaluation.mFields[i];
        result.AddMatrix("field_" + field_names[i], ToMatrix(r_field, n, r_field.Size() / std::max<std::size_t>(n, 1)));
    }
    return result;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

namespace
{

/// MeshCad with an optional background model part: its "background_size_variable" nodal values
/// are the size field (MeshCadSizeSources::pBackground) and its "background_hessian_variable" ones a
/// Hessian whose metric drives the mesh (pMetricBackground, with "metric_error").
Parameters MeshCadImpl(
    const std::string& rFileName,
    Parameters Settings,
    ModelPart& rDestination,
    const ModelPart* pBackground
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rDestination.IsDistributed())
        << "MeshioPlusPlusIgaUtilities does not support distributed model parts" << std::endl;
    KRATOS_ERROR_IF(rDestination.NumberOfNodes() > 0)
        << "The destination model part \"" << rDestination.FullName() << "\" must be empty" << std::endl;
    Settings.AddMissingParameters(MeshioPlusPlusIgaUtilities::GetDefaultParameters());

    const IgaModel model = ReadModel(rFileName, Settings);
    KRATOS_ERROR_IF_NOT(std::holds_alternative<mio::BrepInfo>(model))
        << "MeshCad meshes a B-rep (cad_json, step or iges); \"" << rFileName << "\" is a NURBS mesh" << std::endl;

    mio::MeshCadSurfaceOptions options;
    options.mTargetSize = Settings["target_size"].GetDouble();
    options.mMaxChordalDeviation = Settings["max_chordal_deviation"].GetDouble();
    options.mMinPointsPerEdge = Settings["min_points_per_edge"].GetInt();
    options.mGeometricTolerance = Settings["geometric_tolerance"].GetDouble();
    options.mEdges = Settings["edges"].GetBool();
    options.mStrict = Settings["strict"].GetBool();
    options.mMaxPoints = Settings["max_points"].GetInt();
    options.mMinAngleDeg = Settings["min_angle"].GetDouble();
    options.mSmoothIterations = Settings["smooth_iterations"].GetInt();
    const std::string smooth_method = Settings["smooth_method"].GetString();
    KRATOS_ERROR_IF(smooth_method != "laplacian" && smooth_method != "odt")
        << "\"smooth_method\" is \"laplacian\" or \"odt\", not \"" << smooth_method << "\"" << std::endl;
    options.mSmoothMethod = smooth_method == "odt" ? mio::MeshCadSmoothing::Odt : mio::MeshCadSmoothing::Laplacian;
    options.mGrading = Settings["grading"].GetDouble();
    options.mCurvatureDeviation = Settings["curvature_deviation"].GetDouble();
    options.mCellsAcrossGaps = Settings["cells_across_gaps"].GetDouble();
    const std::string metric = Settings["metric"].GetString();
    KRATOS_ERROR_IF(metric != "mean" && metric != "varying")
        << "\"metric\" is \"mean\" or \"varying\", not \"" << metric << "\"" << std::endl;
    options.mMetric = metric == "varying" ? mio::MeshCadMetric::Varying : mio::MeshCadMetric::Mean;
    options.mPoleRings = Settings["pole_rings"].GetInt();
    options.mReparamConditionNumber = Settings["reparam_condition_number"].GetDouble();
    options.mOrder = Settings["order"].GetInt();

    mio::MeshCadSizeSources sources;
    sources.mFaceSizes = ReadSizes(Settings, "face_sizes");
    sources.mEdgeSizes = ReadSizes(Settings, "edge_sizes");
    sources.mVirtualEdges = Settings["virtual_edges"].GetBool();

    const std::string size_variable = Settings["background_size_variable"].GetString();
    const std::string hessian_variable = Settings["background_hessian_variable"].GetString();
    mio::Mesh background;
    if (pBackground != nullptr) {
        KRATOS_ERROR_IF(size_variable.empty() && hessian_variable.empty())
            << "MeshCad was given a background model part, but neither \"background_size_variable\" "
            << "nor \"background_hessian_variable\" names the nodal data to take from it" << std::endl;
        Internals::FieldDataSelection selection;
        for (const auto& r_name : {size_variable, hessian_variable}) {
            if (!r_name.empty()) {
                selection.NodalDataValueVariables.push_back(r_name);
            }
        }
        background = Internals::ModelPartToMeshWithData(*pBackground, true, true, false, selection);
        if (!size_variable.empty()) {
            sources.pBackground = &background;
            sources.mBackgroundArray = size_variable;
        }
        if (!hessian_variable.empty()) {
            KRATOS_ERROR_IF(Settings["metric_error"].GetDouble() <= 0.0)
                << "\"background_hessian_variable\" needs a positive \"metric_error\"" << std::endl;
            sources.pMetricBackground = &background;
            sources.mMetricArray = hessian_variable;
            sources.mMetricError = Settings["metric_error"].GetDouble();
        }
    } else {
        KRATOS_ERROR_IF(!size_variable.empty() || !hessian_variable.empty())
            << "\"background_size_variable\" and \"background_hessian_variable\" need a background "
            << "model part (the four-argument MeshCad)" << std::endl;
    }

    mio::MeshCadSurfaceExtras extras;
    extras.mModelTolerances = Settings["model_tolerances"].GetBool();
    extras.mVirtualFaces = Settings["virtual_faces"].GetBool();
    extras.mLopCorners = Settings["lop_corners"].GetBool();
    extras.mUntangle = Settings["untangle"].GetBool();
    const std::string quads = Settings["quads"].GetString();
    KRATOS_ERROR_IF(quads != "none" && quads != "dominant" && quads != "all")
        << "\"quads\" is \"none\", \"dominant\" or \"all\", not \"" << quads << "\"" << std::endl;
    extras.mQuads = quads == "all" ? mio::MeshCadQuads::All
                  : quads == "dominant" ? mio::MeshCadQuads::Dominant : mio::MeshCadQuads::None;
    const std::string quad_pairs = Settings["quad_pairs"].GetString();
    KRATOS_ERROR_IF(quad_pairs != "matching" && quad_pairs != "greedy")
        << "\"quad_pairs\" is \"matching\" or \"greedy\", not \"" << quad_pairs << "\"" << std::endl;
    extras.mQuadPairs = quad_pairs == "greedy" ? mio::MeshCadQuadPairs::Greedy : mio::MeshCadQuadPairs::Matching;

    mio::MeshCadSizeReport size_report;
    mio::MeshCadSurfaceExtrasReport extras_report;
    mio::MeshCadSurfaceResult result = mio::mesh_cad_surface(
        std::get<mio::BrepInfo>(model), options, sources, extras, &size_report, &extras_report);

    // The face and edge of every cell are colon-namespaced ("cad:face", "cad:edge"), which no
    // Kratos Variable is; naming a registered int Variable carries them onto the entities, which
    // is what ProjectToCad reads back. The per-face sub model parts "face:<k>" come through anyway.
    const std::string face_ids_variable = Settings["face_ids_variable"].GetString();
    if (!face_ids_variable.empty() && mio::data_has(result.mMesh, mio::DataLocation::Cell, mio::kMeshCadFaceName)) {
        result.mMesh = mio::data_rename(result.mMesh, mio::DataLocation::Cell, mio::kMeshCadFaceName, face_ids_variable);
    }
    const std::string edge_ids_variable = Settings["edge_ids_variable"].GetString();
    if (!edge_ids_variable.empty() && mio::data_has(result.mMesh, mio::DataLocation::Cell, mio::kMeshCadEdgeName)) {
        result.mMesh = mio::data_rename(result.mMesh, mio::DataLocation::Cell, mio::kMeshCadEdgeName, edge_ids_variable);
    }
    Internals::MeshToModelPart(result.mMesh, rDestination);

    Parameters report(R"({})");
    report.AddInt("number_of_nodes", static_cast<int>(rDestination.NumberOfNodes()));
    report.AddInt("number_of_elements", static_cast<int>(rDestination.NumberOfElements()));
    report.AddInt("number_of_conditions", static_cast<int>(rDestination.NumberOfConditions()));
    report.AddInt("number_of_points", static_cast<int>(result.mNumPoints));
    report.AddInt("number_of_triangles", static_cast<int>(result.mNumTriangles));
    report.AddInt("number_of_edge_segments", static_cast<int>(result.mNumEdgeSegments));
    report.AddInt("number_of_faces", static_cast<int>(result.mNumFaces));
    report.AddInt("number_of_failed_faces", static_cast<int>(result.mNumFailed));
    report.AddInt("number_of_edges", static_cast<int>(result.mNumEdges));
    report.AddInt("open_edges", static_cast<int>(result.mNumOpenEdges));
    report.AddInt("non_manifold_edges", static_cast<int>(result.mNumNonManifoldEdges));
    report.AddInt("misoriented_edges", static_cast<int>(result.mNumMisorientedEdges));
    report.AddInt("degenerate_dropped", static_cast<int>(result.mNumDegenerateDropped));
    report.AddInt("refinement_points", static_cast<int>(result.mNumRefinementPoints));
    report.AddInt("faces_below_bound", static_cast<int>(result.mNumFacesBelowBound));
    report.AddInt("smoothed", static_cast<int>(result.mNumSmoothed));
    report.AddInt("curved_straightened", static_cast<int>(result.mNumCurvedStraightened));
    report.AddDouble("tolerance", result.mTolerance);
    report.AddDouble("min_angle", result.mMinAngleDeg);
    report.AddDouble("max_angle", result.mMaxAngleDeg);
    report.AddDouble("min_edge", result.mMinEdge);
    report.AddDouble("max_edge", result.mMaxEdge);
    report.AddDouble("max_edge_gap", result.mMaxEdgeGap);
    Parameters histogram(R"([])");
    for (const auto count : result.mAngleHistogram) {
        histogram.Append(static_cast<int>(count));
    }
    report.AddValue("angle_histogram", histogram);
    Parameters faces(R"([])");
    for (const auto& r_face : result.mFaces) {
        Parameters face(R"({})");
        face.AddInt("face", static_cast<int>(r_face.mFace));
        face.AddInt("id", static_cast<int>(r_face.mId));
        face.AddString("status", mio::mesh_cad_face_status_name(r_face.mStatus));
        face.AddString("reason", r_face.mReason);
        face.AddInt("number_of_triangles", static_cast<int>(r_face.mNumTriangles));
        face.AddDouble("min_angle", r_face.mMinAngleDeg);
        face.AddDouble("max_angle", r_face.mMaxAngleDeg);
        face.AddBool("below_bound", r_face.mBelowBound);
        faces.Append(face);
    }
    report.AddValue("faces", faces);
    report.AddInt("untangled", static_cast<int>(extras_report.mNumUntangled));
    report.AddInt("quadrilaterals", static_cast<int>(extras_report.mNumQuads));
    report.AddInt("leftover_triangles", static_cast<int>(extras_report.mNumLeftoverTriangles));
    report.AddInt("reparametrised_faces", static_cast<int>(extras_report.mReparametrised.size()));
    report.AddInt("merged_faces", static_cast<int>(extras_report.mMergedFaces.size()));
    report.AddInt("lopped_corners", static_cast<int>(extras_report.mLoppedCorners.size()));
    report.AddStringArray("warnings", result.mWarnings);
    return report;

    KRATOS_CATCH("")
}

} // namespace

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIgaUtilities::MeshCad(
    const std::string& rFileName,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    return MeshCadImpl(rFileName, Settings, rDestination, nullptr);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIgaUtilities::MeshCad(
    const std::string& rFileName,
    Parameters Settings,
    ModelPart& rDestination,
    const ModelPart& rBackground
    )
{
    KRATOS_ERROR_IF(rBackground.IsDistributed())
        << "MeshioPlusPlusIgaUtilities does not support distributed model parts" << std::endl;
    return MeshCadImpl(rFileName, Settings, rDestination, &rBackground);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIgaUtilities::ProjectToCad(
    const ModelPart& rSource,
    const std::string& rFileName,
    Parameters Settings,
    ModelPart& rDestination
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSource.IsDistributed() || rDestination.IsDistributed())
        << "MeshioPlusPlusIgaUtilities does not support distributed model parts" << std::endl;
    Settings.AddMissingParameters(GetDefaultParameters());

    const std::string face_ids_variable = Settings["face_ids_variable"].GetString();
    KRATOS_ERROR_IF(face_ids_variable.empty())
        << "ProjectToCad needs \"face_ids_variable\": the registered int Variable holding each "
        << "entity's CAD face, as MeshCad wrote it" << std::endl;
    const std::string edge_ids_variable = Settings["edge_ids_variable"].GetString();

    const IgaModel model = ReadModel(rFileName, Settings);
    KRATOS_ERROR_IF_NOT(std::holds_alternative<mio::BrepInfo>(model))
        << "ProjectToCad projects onto a B-rep (cad_json, step or iges); \"" << rFileName
        << "\" is a NURBS mesh" << std::endl;

    Internals::FieldDataSelection selection;
    selection.ElementDataValueVariables = {face_ids_variable};
    selection.ConditionDataValueVariables = {face_ids_variable};
    if (!edge_ids_variable.empty()) {
        selection.ElementDataValueVariables.push_back(edge_ids_variable);
        selection.ConditionDataValueVariables.push_back(edge_ids_variable);
    }
    const mio::Mesh mesh = Internals::ModelPartToMeshWithData(rSource, true, true, false, selection);

    mio::ProjectToCadOptions options;
    options.mFaceArray = face_ids_variable;
    options.mEdgeArray = edge_ids_variable;
    options.mMaxDistance = Settings["max_distance"].GetDouble();
    options.mTrimSamples = Settings["trim_samples"].GetInt();
    mio::ProjectToCadResult result = mio::project_to_cad(mesh, std::get<mio::BrepInfo>(model), options);
    Internals::MeshToModelPart(result.mMesh, rDestination);

    Parameters report(R"({})");
    report.AddInt("projected", static_cast<int>(result.mNumProjected));
    report.AddInt("on_edges", static_cast<int>(result.mNumOnEdges));
    report.AddInt("on_curves", static_cast<int>(result.mNumOnCurves));
    report.AddInt("corners", static_cast<int>(result.mNumCorners));
    report.AddInt("ambiguous", static_cast<int>(result.mNumAmbiguous));
    report.AddInt("too_far", static_cast<int>(result.mNumTooFar));
    report.AddInt("failed", static_cast<int>(result.mNumFailed));
    report.AddInt("outside_trim", static_cast<int>(result.mNumOutsideTrim));
    report.AddDouble("max_move", result.mMaxMove);
    report.AddStringArray("warnings", result.mWarnings);
    return report;

    KRATOS_CATCH("")
}

} // namespace Kratos
