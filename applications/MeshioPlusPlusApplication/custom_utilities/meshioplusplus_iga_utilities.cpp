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
#include "meshioplusplus/operations/iga.hpp"

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
        "derivatives"            : 0
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

} // namespace Kratos
