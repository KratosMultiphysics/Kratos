//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                        Application
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <cmath>

// External includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"

// Project includes
#include "custom_io/meshioplusplus_io.h"
#include "custom_utilities/meshioplusplus_conversion_utilities.h"
#include "includes/kratos_components.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "processes/integration_values_extrapolation_to_nodes_process.h"
#include "utilities/geometry_utilities.h"

namespace mio = meshioplusplus;

namespace Kratos
{
namespace
{

using DataArray = mio::XdmfTimeSeriesWriter::NamedArray;

/**
 * @brief Collects the listed variables of a container as flat data arrays.
 * @details Supported types (mirroring VtkOutput): double, int, bool (scalar),
 * array_1d<double, 3/4/6/9> and Vector (multi-component, Vector size taken
 * from the first entity). Unknown or unsupported variables are skipped with a
 * warning. rGetValue(entity, variable) provides the value, rValidate(variable)
 * runs an optional per-variable check (e.g. historical availability).
 */
template <class TContainer, class TGetter, class TValidator>
void CollectVariableDataArrays(
    const std::vector<std::string>& rVariableNames,
    const TContainer& rEntities,
    const std::size_t NumberOfEntities,
    const TGetter& rGetValue,
    const TValidator& rValidate,
    std::vector<DataArray>& rOutput
    )
{
    for (const auto& r_variable_name : rVariableNames) {
        DataArray data;
        data.mName = r_variable_name;

        auto collect_scalar = [&](const auto& rVariable) {
            rValidate(rVariable);
            data.mNumComponents = 1;
            data.mValues.reserve(NumberOfEntities);
            for (const auto& r_entity : rEntities) {
                data.mValues.push_back(static_cast<double>(rGetValue(r_entity, rVariable)));
            }
        };
        auto collect_vector = [&](const auto& rVariable, const std::size_t NumberOfComponents) {
            rValidate(rVariable);
            data.mNumComponents = NumberOfComponents;
            data.mValues.reserve(NumberOfEntities * NumberOfComponents);
            for (const auto& r_entity : rEntities) {
                const auto& r_value = rGetValue(r_entity, rVariable);
                for (std::size_t i = 0; i < NumberOfComponents; ++i) {
                    data.mValues.push_back(i < r_value.size() ? static_cast<double>(r_value[i]) : 0.0);
                }
            }
        };

        if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            collect_scalar(KratosComponents<Variable<double>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<int>>::Has(r_variable_name)) {
            collect_scalar(KratosComponents<Variable<int>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<bool>>::Has(r_variable_name)) {
            collect_scalar(KratosComponents<Variable<bool>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name), 3);
        } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 4>>>::Get(r_variable_name), 4);
        } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 6>>>::Get(r_variable_name), 6);
        } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 9>>>::Get(r_variable_name), 9);
        } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
            const auto& r_variable = KratosComponents<Variable<Vector>>::Get(r_variable_name);
            const std::size_t number_of_components =
                NumberOfEntities > 0 ? rGetValue(*rEntities.begin(), r_variable).size() : 0;
            if (number_of_components == 0) {
                KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Vector variable \"" << r_variable_name
                    << "\" has no components on the first entity - skipping it" << std::endl;
                continue;
            }
            collect_vector(r_variable, number_of_components);
        } else {
            KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Variable \"" << r_variable_name
                << "\" is not registered with a type suitable for MeshioPlusPlusIO - skipping it" << std::endl;
            continue;
        }

        rOutput.push_back(std::move(data));
    }
}

/**
 * @brief Collects the listed flags of a container as 1/0/-1 scalar arrays
 * (VtkOutput convention: -1 when the flag is not defined on the entity).
 */
template <class TContainer>
void CollectFlagDataArrays(
    const std::vector<std::string>& rFlagNames,
    const TContainer& rEntities,
    const std::size_t NumberOfEntities,
    std::vector<DataArray>& rOutput
    )
{
    for (const auto& r_flag_name : rFlagNames) {
        if (!KratosComponents<Flags>::Has(r_flag_name)) {
            KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Flag \"" << r_flag_name
                << "\" is not registered - skipping it" << std::endl;
            continue;
        }
        const Flags& r_flag = KratosComponents<Flags>::Get(r_flag_name);

        DataArray data;
        data.mName = r_flag_name;
        data.mNumComponents = 1;
        data.mValues.reserve(NumberOfEntities);
        for (const auto& r_entity : rEntities) {
            data.mValues.push_back(r_entity.IsDefined(r_flag) ? (r_entity.Is(r_flag) ? 1.0 : 0.0) : -1.0);
        }
        rOutput.push_back(std::move(data));
    }
}

/**
 * @brief Collects entity ids as a scalar array named rName.
 */
template <class TContainer>
DataArray CollectIdsArray(
    const TContainer& rEntities,
    const std::size_t NumberOfEntities,
    const std::string& rName
    )
{
    DataArray data;
    data.mName = rName;
    data.mNumComponents = 1;
    data.mValues.reserve(NumberOfEntities);
    for (const auto& r_entity : rEntities) {
        data.mValues.push_back(static_cast<double>(r_entity.Id()));
    }
    return data;
}

/**
 * @brief Collects the entities' properties ids as the PROPERTIES_ID array.
 */
template <class TContainer>
DataArray CollectPropertiesIdsArray(
    const TContainer& rEntities,
    const std::size_t NumberOfEntities
    )
{
    DataArray data;
    data.mName = "PROPERTIES_ID";
    data.mNumComponents = 1;
    data.mValues.reserve(NumberOfEntities);
    for (const auto& r_entity : rEntities) {
        data.mValues.push_back(static_cast<double>(r_entity.GetProperties().Id()));
    }
    return data;
}

/**
 * @brief Collects gauss point results averaged over the integration points
 * (VtkOutput convention). Variables whose CalculateOnIntegrationPoints returns
 * nothing (e.g. the generic core entities) are skipped with a warning.
 * @note CalculateOnIntegrationPoints is non-const, hence the non-const container.
 */
template <class TContainer>
void CollectGaussPointDataArrays(
    const std::vector<std::string>& rVariableNames,
    TContainer& rEntities,
    const std::size_t NumberOfEntities,
    const ProcessInfo& rProcessInfo,
    std::vector<DataArray>& rOutput
    )
{
    for (const auto& r_variable_name : rVariableNames) {
        DataArray data;
        data.mName = r_variable_name;

        // Averages the integration point values of one entity into rAppendTo
        auto collect_scalar = [&]<class TDataType>(const Variable<TDataType>& rVariable) -> bool {
            std::vector<TDataType> gp_values;
            if (NumberOfEntities == 0) {
                return false;
            }
            rEntities.begin()->CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
            if (gp_values.empty()) {
                KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Gauss point variable \"" << r_variable_name
                    << "\" returns no integration point values - skipping it" << std::endl;
                return false;
            }
            data.mNumComponents = 1;
            data.mValues.reserve(NumberOfEntities);
            for (auto& r_entity : rEntities) {
                r_entity.CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
                double average = 0.0;
                for (const auto& r_value : gp_values) {
                    average += static_cast<double>(r_value);
                }
                data.mValues.push_back(gp_values.empty() ? 0.0 : average / gp_values.size());
            }
            return true;
        };
        auto collect_vector = [&]<class TDataType>(const Variable<TDataType>& rVariable) -> bool {
            std::vector<TDataType> gp_values;
            if (NumberOfEntities == 0) {
                return false;
            }
            rEntities.begin()->CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
            if (gp_values.empty() || gp_values[0].size() == 0) {
                KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Gauss point variable \"" << r_variable_name
                    << "\" returns no integration point values - skipping it" << std::endl;
                return false;
            }
            const std::size_t number_of_components = gp_values[0].size();
            data.mNumComponents = number_of_components;
            data.mValues.reserve(NumberOfEntities * number_of_components);
            std::vector<double> average(number_of_components);
            for (auto& r_entity : rEntities) {
                r_entity.CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
                std::fill(average.begin(), average.end(), 0.0);
                for (const auto& r_value : gp_values) {
                    for (std::size_t i = 0; i < number_of_components && i < r_value.size(); ++i) {
                        average[i] += r_value[i];
                    }
                }
                for (std::size_t i = 0; i < number_of_components; ++i) {
                    data.mValues.push_back(gp_values.empty() ? 0.0 : average[i] / gp_values.size());
                }
            }
            return true;
        };

        bool collected = false;
        if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            collected = collect_scalar(KratosComponents<Variable<double>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<int>>::Has(r_variable_name)) {
            collected = collect_scalar(KratosComponents<Variable<int>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<bool>>::Has(r_variable_name)) {
            collected = collect_scalar(KratosComponents<Variable<bool>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
            collected = collect_vector(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(r_variable_name)) {
            collected = collect_vector(KratosComponents<Variable<array_1d<double, 6>>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
            collected = collect_vector(KratosComponents<Variable<Vector>>::Get(r_variable_name));
        } else {
            KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Gauss point variable \"" << r_variable_name
                << "\" is not registered with a type suitable for MeshioPlusPlusIO - skipping it" << std::endl;
        }

        if (collected) {
            rOutput.push_back(std::move(data));
        }
    }
}

/**
 * @brief Merges per-kind cell arrays into combined ones covering the full cell
 * range of the written mesh (element rows first, then condition rows),
 * zero-filling the rows of the entity kind an array does not apply to.
 */
std::vector<DataArray> MergeCellDataParts(
    std::vector<DataArray>&& rElementPart,
    const std::size_t NumberOfElementRows,
    std::vector<DataArray>&& rConditionPart,
    const std::size_t NumberOfConditionRows
    )
{
    std::vector<DataArray> merged;
    merged.reserve(rElementPart.size() + rConditionPart.size());

    for (auto& r_element_data : rElementPart) {
        const auto it_condition = std::find_if(rConditionPart.begin(), rConditionPart.end(),
            [&r_element_data](const DataArray& rData) { return rData.mName == r_element_data.mName; });
        if (it_condition != rConditionPart.end() &&
            it_condition->mNumComponents == r_element_data.mNumComponents) {
            r_element_data.mValues.insert(r_element_data.mValues.end(),
                                         it_condition->mValues.begin(), it_condition->mValues.end());
            rConditionPart.erase(it_condition);
        } else {
            r_element_data.mValues.resize(
                r_element_data.mValues.size() + NumberOfConditionRows * r_element_data.mNumComponents, 0.0);
        }
        merged.push_back(std::move(r_element_data));
    }

    for (auto& r_condition_data : rConditionPart) {
        r_condition_data.mValues.insert(r_condition_data.mValues.begin(),
                                       NumberOfElementRows * r_condition_data.mNumComponents, 0.0);
        merged.push_back(std::move(r_condition_data));
    }

    return merged;
}

/**
 * @brief Builds a meshio++ NDArray ({n} or {n, components} Float64) from a DataArray.
 */
mio::NDArray ToNDArray(const DataArray& rData)
{
    const std::size_t number_of_rows =
        rData.mNumComponents > 0 ? rData.mValues.size() / rData.mNumComponents : 0;
    mio::NDArray array = mio::NDArray::Uninit(
        mio::DType::Float64,
        rData.mNumComponents == 1 ? std::vector<std::size_t>{number_of_rows}
                                      : std::vector<std::size_t>{number_of_rows, rData.mNumComponents});
    std::copy(rData.mValues.begin(), rData.mValues.end(), array.As<double>());
    return array;
}

// Format enum <-> canonical meshio++ format name. Keep in sync with the
// Format enum in meshioplusplus_io.h and meshio++'s registry.cpp.
#define KRATOS_MESHIOPLUSPLUS_FORMATS(X) \
    X(ABAQUS, "abaqus")                  \
    X(ANSYS, "ansys")                    \
    X(ANSYSINP, "ansysinp")              \
    X(AVSUCD, "avsucd")                  \
    X(CGNS, "cgns")                      \
    X(DEX, "dex")                        \
    X(DOLFIN, "dolfin")                  \
    X(ENSIGHT, "ensight")                \
    X(EXODUS, "exodus")                  \
    X(FLAC3D, "flac3d")                  \
    X(FLUX, "flux")                      \
    X(FREEFEM, "freefem")                \
    X(GMSH, "gmsh")                      \
    X(GMSH22, "gmsh22")                  \
    X(H5M, "h5m")                        \
    X(HMF, "hmf")                        \
    X(IP, "ip")                          \
    X(MDPA, "mdpa")                      \
    X(MED, "med")                        \
    X(MEDIT, "medit")                    \
    X(MFF, "mff")                        \
    X(MFM, "mfm")                        \
    X(MPHTXT, "mphtxt")                  \
    X(NASTRAN, "nastran")                \
    X(NETGEN, "netgen")                  \
    X(OBJ, "obj")                        \
    X(OFF, "off")                        \
    X(OPENFOAM, "openfoam")              \
    X(PERMAS, "permas")                  \
    X(PLY, "ply")                        \
    X(STL, "stl")                        \
    X(SU2, "su2")                        \
    X(SVG, "svg")                        \
    X(TECPLOT, "tecplot")                \
    X(TETGEN, "tetgen")                  \
    X(TIKZ, "tikz")                      \
    X(TRIANGLE, "triangle")              \
    X(UGRID, "ugrid")                    \
    X(UNV, "unv")                        \
    X(VTK, "vtk")                        \
    X(VTP, "vtp")                        \
    X(VTU, "vtu")                        \
    X(WKT, "wkt")                        \
    X(XDMF, "xdmf")

const std::unordered_map<std::string, MeshioPlusPlusIO::Format>& GetFormatNameMap()
{
    static const std::unordered_map<std::string, MeshioPlusPlusIO::Format> format_map = {
#define KRATOS_MESHIOPLUSPLUS_FORMAT_LOOKUP(EnumName, Name) {Name, MeshioPlusPlusIO::Format::EnumName},
        KRATOS_MESHIOPLUSPLUS_FORMATS(KRATOS_MESHIOPLUSPLUS_FORMAT_LOOKUP)
#undef KRATOS_MESHIOPLUSPLUS_FORMAT_LOOKUP
    };
    return format_map;
}

/**
 * @brief Writes rMesh honoring an ascii/binary "file_format" override for the
 * formats exposing such a flag, bypassing the meshio++ registry defaults.
 * @param Skin Whether stl/ply extract and write the boundary skin of volume
 *        meshes (the "skin" setting); ignored by the other formats.
 * @return false when the format has no ascii/binary variant (caller falls
 *         back to the registry writer).
 */
bool WriteWithFileFormatOverride(
    const std::string& rFormatName,
    const bool Binary,
    const bool Skin,
    const std::filesystem::path& rPath,
    const mio::Mesh& rMesh
    )
{
    const std::string path = rPath.string();
    if (rFormatName == "vtu") {
        mio::write_vtu(path, rMesh, Binary, /*zlib=*/Binary);
        return true;
    }
    if (rFormatName == "vtp") {
        mio::write_vtp(path, rMesh, Binary, /*zlib=*/Binary);
        return true;
    }
    if (rFormatName == "vtk") {
        mio::write_vtk(path, rMesh, Binary, /*v51=*/true);
        return true;
    }
    if (rFormatName == "gmsh") {
        mio::write_gmsh41(path, rMesh, Binary);
        return true;
    }
    if (rFormatName == "stl") {
        mio::write_stl(path, rMesh, Binary, Skin);
        return true;
    }
    if (rFormatName == "ply") {
        mio::write_ply(path, rMesh, Binary, Skin);
        return true;
    }
    if (rFormatName == "ansys") {
        mio::write_ansys(path, rMesh, Binary);
        return true;
    }
    if (rFormatName == "flac3d") {
        mio::write_flac3d(path, rMesh, ".16e", Binary);
        return true;
    }
    if (rFormatName == "ensight") {
        mio::write_ensight(path, rMesh, Binary);
        return true;
    }
    return false;
}

} // namespace


/***********************************************************************************/
/***********************************************************************************/

MeshioPlusPlusIO::MeshioPlusPlusIO(
    const std::filesystem::path& rFileName,
    Parameters ThisParameters
    ) : mFileName(rFileName),
        mParameters(ThisParameters)
{
    KRATOS_TRY

    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Eagerly validate the enumerated settings so misconfiguration fails at
    // construction and not in the middle of a simulation.
    const std::string time_series = mParameters["time_series"].GetString();
    KRATOS_ERROR_IF(time_series != "automatic" && time_series != "file_series" && time_series != "single_file")
        << "Unknown \"time_series\" setting \"" << time_series
        << "\" (use \"automatic\", \"file_series\" or \"single_file\")" << std::endl;

    const std::string output_control = mParameters["output_control_type"].GetString();
    KRATOS_ERROR_IF(output_control != "step" && output_control != "time")
        << "Unknown \"output_control_type\" setting \"" << output_control
        << "\" (use \"step\" or \"time\")" << std::endl;

    const std::string xdmf_data_format = mParameters["xdmf_data_format"].GetString();
    KRATOS_ERROR_IF(xdmf_data_format != "auto" && xdmf_data_format != "XML" &&
                    xdmf_data_format != "Binary" && xdmf_data_format != "HDF")
        << "Unknown \"xdmf_data_format\" setting \"" << xdmf_data_format
        << "\" (use \"auto\", \"XML\", \"Binary\" or \"HDF\")" << std::endl;

    const std::string entity_type = mParameters["entity_type"].GetString();
    KRATOS_ERROR_IF(entity_type != "automatic" && entity_type != "element" && entity_type != "condition")
        << "Unknown \"entity_type\" setting \"" << entity_type
        << "\" (use \"automatic\", \"element\" or \"condition\")" << std::endl;

    const std::string file_format = mParameters["file_format"].GetString();
    KRATOS_ERROR_IF(file_format != "default" && file_format != "ascii" && file_format != "binary")
        << "Unknown \"file_format\" setting \"" << file_format
        << "\" (use \"default\", \"ascii\" or \"binary\")" << std::endl;

    // Resolves the format name, throwing a descriptive error for unknown names
    // and formats compiled out of this build.
    FormatFromString(mParameters["format"].GetString());

    // Gauss point results extrapolated to the nodes are written as
    // non-historical nodal data (VtkOutput convention); the extrapolation
    // process itself is created lazily at the first write (it needs the
    // model part, which the IO does not hold at construction).
    for (const auto& r_name : mParameters["gauss_point_variables_extrapolated_to_nodes"].GetStringArray()) {
        mParameters["nodal_data_value_variables"].Append(r_name);
    }

    KRATOS_CATCH("")
}

MeshioPlusPlusIO::~MeshioPlusPlusIO()
{
    // The writers finalize themselves on destruction anyway; going through CloseOutput
    // keeps that in one place. A failure here can only be logged - a destructor must not
    // throw, which is exactly why CloseOutput is also public.
    try {
        CloseOutput();
    } catch (const std::exception& r_exception) {
        KRATOS_WARNING("MeshioPlusPlusIO")
            << "Could not finish the transient output of " << mFileName << ": "
            << r_exception.what() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::CloseOutput()
{
    KRATOS_TRY

    // Finalize before releasing: the destructor would do it too, but then a failure could
    // only be logged. Clearing the map is what actually hands the file back to the caller,
    // so a later delete is not undone by a writer still holding it.
    for (auto& r_entry : mXdmfWriters) {
        if (r_entry.second != nullptr) {
            r_entry.second->Finalize();
        }
    }
    mXdmfWriters.clear();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIO::GetDefaultParameters()
{
    return Parameters(R"({
        "format"                                      : "auto",
        "file_format"                                 : "default",
        "skin"                                        : true,
        "time_series"                                 : "automatic",
        "output_control_type"                         : "step",
        "output_precision"                            : 7,
        "label_precision"                             : 4,
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "entity_type"                                 : "automatic",
        "output_sub_model_parts"                      : false,
        "write_deformed_configuration"                : false,
        "write_ids"                                   : false,
        "xdmf_data_format"                            : "auto",
        "xdmf_auto_flush"                             : true,
        "xdmf_gzip_level"                             : -1,
        "nodal_solution_step_data_variables"          : [],
        "nodal_data_value_variables"                  : [],
        "nodal_flags"                                 : [],
        "element_data_value_variables"                : [],
        "element_flags"                               : [],
        "condition_data_value_variables"              : [],
        "condition_flags"                             : [],
        "gauss_point_variables_extrapolated_to_nodes" : [],
        "gauss_point_variables_in_elements"           : []
    })");
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> MeshioPlusPlusIO::GetSupportedFormats()
{
    std::set<std::string> names;
    for (const auto& r_entry : mio::registry_readers()) {
        names.insert(r_entry.first);
    }
    for (const auto& r_entry : mio::registry_writers()) {
        names.insert(r_entry.first);
    }
    return std::vector<std::string>(names.begin(), names.end());
}

std::vector<std::string> MeshioPlusPlusIO::GetSupportedReadFormats()
{
    std::vector<std::string> names;
    names.reserve(mio::registry_readers().size());
    for (const auto& r_entry : mio::registry_readers()) {
        names.push_back(r_entry.first);
    }
    return names;
}

std::vector<std::string> MeshioPlusPlusIO::GetSupportedWriteFormats()
{
    std::vector<std::string> names;
    names.reserve(mio::registry_writers().size());
    for (const auto& r_entry : mio::registry_writers()) {
        names.push_back(r_entry.first);
    }
    return names;
}

/***********************************************************************************/
/***********************************************************************************/

MeshioPlusPlusIO::Format MeshioPlusPlusIO::FormatFromString(const std::string& rFormatName)
{
    KRATOS_TRY

    std::string lower_name = rFormatName;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
                   [](unsigned char Character) { return std::tolower(Character); });

    if (lower_name.empty() || lower_name == "auto" || lower_name == "automatic") {
        return Format::AUTOMATIC;
    }

    const auto& r_format_map = GetFormatNameMap();
    const auto it = r_format_map.find(lower_name);
    if (it == r_format_map.end()) {
        std::ostringstream supported;
        for (const auto& r_name : GetSupportedFormats()) {
            supported << r_name << " ";
        }
        KRATOS_ERROR << "Unknown format \"" << rFormatName
                     << "\". Supported formats in this build: " << supported.str() << std::endl;
    }

    // Known format compiled out of this build: fail with the missing dependency.
    const char* p_missing_dependency = mio::registry_compiled_out(lower_name);
    KRATOS_ERROR_IF(p_missing_dependency != nullptr)
        << "Format \"" << lower_name << "\" is not available in this build (requires "
        << p_missing_dependency << ")" << std::endl;

    return it->second;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::string MeshioPlusPlusIO::FormatName(const Format Value)
{
    switch (Value) {
        case Format::AUTOMATIC: return "auto";
#define KRATOS_MESHIOPLUSPLUS_FORMAT_NAME(EnumName, Name) \
        case Format::EnumName: return Name;
        KRATOS_MESHIOPLUSPLUS_FORMATS(KRATOS_MESHIOPLUSPLUS_FORMAT_NAME)
#undef KRATOS_MESHIOPLUSPLUS_FORMAT_NAME
    }
    KRATOS_ERROR << "Unknown Format enum value" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

MeshioPlusPlusIO::Format MeshioPlusPlusIO::ResolveFormat(const std::filesystem::path& rPath)
{
    KRATOS_TRY

    std::string format_name;
    try {
        format_name = mio::resolve_format(rPath.string(), "");
    } catch (const std::exception& r_exception) {
        KRATOS_ERROR << "Cannot resolve a format from the extension of " << rPath
                     << ": " << r_exception.what() << std::endl;
    }

    const auto& r_format_map = GetFormatNameMap();
    const auto it = r_format_map.find(format_name);
    KRATOS_ERROR_IF(it == r_format_map.end())
        << "The extension of " << rPath << " resolves to unknown format \"" << format_name << "\"" << std::endl;
    return it->second;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool MeshioPlusPlusIO::IsFormatAvailable(const Format Value)
{
    if (Value == Format::AUTOMATIC) {
        return true;
    }
    const std::string name = FormatName(Value);
    if (mio::registry_compiled_out(name) != nullptr) {
        return false;
    }
    return mio::registry_readers().count(name) > 0 || mio::registry_writers().count(name) > 0;
}

/***********************************************************************************/
/***********************************************************************************/

std::string MeshioPlusPlusIO::ResolveEffectiveFormat(const bool CheckWritable) const
{
    KRATOS_TRY

    std::string format_name = mParameters["format"].GetString();
    std::transform(format_name.begin(), format_name.end(), format_name.begin(),
                   [](unsigned char Character) { return std::tolower(Character); });

    if (format_name.empty() || format_name == "auto" || format_name == "automatic") {
        try {
            format_name = mio::resolve_format(mFileName.string(), "");
        } catch (const std::exception& r_exception) {
            KRATOS_ERROR << "Cannot resolve a format from the extension of " << mFileName
                         << ": " << r_exception.what() << ". Set the \"format\" setting explicitly." << std::endl;
        }
    }

    const char* p_missing_dependency = mio::registry_compiled_out(format_name);
    KRATOS_ERROR_IF(p_missing_dependency != nullptr)
        << "Format \"" << format_name << "\" is not available in this build (requires "
        << p_missing_dependency << ")" << std::endl;

    if (CheckWritable) {
        KRATOS_ERROR_IF(mio::registry_writers().count(format_name) == 0)
            << "Format \"" << format_name << "\" is read-only or unknown - it cannot be written" << std::endl;
    } else {
        KRATOS_ERROR_IF(mio::registry_readers().count(format_name) == 0)
            << "Format \"" << format_name << "\" is write-only or unknown - it cannot be read" << std::endl;
    }

    return format_name;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::ReadModelPart(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rThisModelPart.IsDistributed())
        << "MeshioPlusPlusIO does not support distributed model parts" << std::endl;
    KRATOS_ERROR_IF(rThisModelPart.NumberOfNodes() > 0)
        << "The target model part \"" << rThisModelPart.FullName() << "\" must be empty" << std::endl;
    KRATOS_ERROR_IF_NOT(std::filesystem::exists(mFileName))
        << "The file " << mFileName << " does not exist" << std::endl;

    const std::string format_name = ResolveEffectiveFormat(false);

    // Read into the meshio++ Kratos backend: the materialized model part view
    // splits max-dimension cell blocks into elements, lower-dimension ones into
    // conditions, and turns integer tag arrays (gmsh physical groups, ...) into
    // sub model parts.
    mio::Mesh mesh = mio::registry_readers().at(format_name)(mFileName.string());

    // In .mdpa a "gmsh:physical" tag is how a Kratos properties id is stored, not a
    // physical group, so the automatic tag pass would synthesize a spurious
    // "gmsh_physical_<id>" sub model part beside the deck's real ones. In a genuine gmsh
    // file the very same key does mean a physical group, which must keep becoming a sub
    // model part - hence the exclusion is scoped to the format rather than applied
    // globally. Only the file-read path needs this: a mesh staged from a Kratos model
    // part never carries the key, because the backend reads "gmsh:physical" from staged
    // cell data and never synthesizes it from entity properties ids.
    if (format_name == "mdpa") {
        mesh.ExcludeTagSubModelPartKey("gmsh:physical");
    }

    mio::ModelPart& r_source = mesh.GetModelPart();

    // Bridge into the real Kratos model part (one bulk O(n) creation pass).
    // The four-argument overload also transfers the material data: meshio++ cannot fill a
    // Kratos::Properties itself (that needs KratosComponents, which it deliberately does
    // not link), so it hands each key/value pair over for the typed assignment.
    mio::to_model_part(
        r_source, rThisModelPart,
        [&rThisModelPart](std::size_t PropertiesId) {
            return rThisModelPart.HasProperties(PropertiesId)
                ? rThisModelPart.pGetProperties(PropertiesId)
                : rThisModelPart.CreateNewProperties(PropertiesId);
        },
        [](Properties::Pointer pProperties, const mio::PropertyValue& rValue) {
            Internals::ApplyMeshioProperty(pProperties, rValue);
        });

    KRATOS_INFO("MeshioPlusPlusIO") << "Read " << rThisModelPart.NumberOfNodes() << " nodes, "
        << rThisModelPart.NumberOfElements() << " elements and " << rThisModelPart.NumberOfConditions()
        << " conditions from " << mFileName << " (format \"" << format_name << "\")" << std::endl;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<double> MeshioPlusPlusIO::GetTimeValues() const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(std::filesystem::exists(mFileName))
        << "The file " << mFileName << " does not exist" << std::endl;

    const std::string format_name = ResolveEffectiveFormat(false);

    // Header-only summary (never the heavy arrays); falls back to a full read for formats
    // without a native metadata path, and reports no time values either way for a format
    // with no time-series concept - correct for every format, just not always as cheap.
    const mio::MeshMetadata metadata = mio::registry_read_metadata(mFileName.string(), format_name, mio::ReadOptions());
    return metadata.mTimeValues;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

int MeshioPlusPlusIO::GetNumberOfTimeSteps() const
{
    KRATOS_TRY

    return static_cast<int>(GetTimeValues().size());

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

int MeshioPlusPlusIO::GetTimeStepIndex(const double TimeValue) const
{
    KRATOS_TRY

    const std::vector<double> time_values = GetTimeValues();
    for (std::size_t i = 0; i < time_values.size(); ++i) {
        if (std::abs(time_values[i] - TimeValue) < 1e-9) {
            return static_cast<int>(i);
        }
    }
    return -1;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rThisModelPart.IsDistributed())
        << "MeshioPlusPlusIO does not support distributed model parts" << std::endl;

    // Extrapolate the requested gauss point results to nodal data first (the
    // process mutates nodal values, hence the contained const_cast - identical
    // effective behavior to VtkOutput, which holds a non-const reference).
    if (mParameters["gauss_point_variables_extrapolated_to_nodes"].size() > 0) {
        if (mpGaussToNodesProcess == nullptr) {
            Parameters gauss_parameters(R"({
                "echo_level"                 : 0,
                "area_average"               : true,
                "average_variable"           : "NODAL_AREA",
                "list_of_variables"          : [],
                "extrapolate_non_historical" : true
            })");
            gauss_parameters["list_of_variables"].SetStringArray(
                mParameters["gauss_point_variables_extrapolated_to_nodes"].GetStringArray());
            mpGaussToNodesProcess = std::make_unique<IntegrationValuesExtrapolationToNodesProcess>(
                const_cast<ModelPart&>(rThisModelPart), gauss_parameters);
        }
        mpGaussToNodesProcess->Execute();
    }

    WriteTarget(rThisModelPart, "");

    if (mParameters["output_sub_model_parts"].GetBool()) {
        for (const auto& r_sub_model_part : rThisModelPart.SubModelParts()) {
            WriteTarget(r_sub_model_part, "_" + rThisModelPart.Name() + "_" + r_sub_model_part.Name());
        }
    }

    ++mOutputStep;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::WriteTarget(
    const ModelPart& rModelPart,
    const std::string& rTargetSuffix
    )
{
    KRATOS_TRY

    const std::string format_name = ResolveEffectiveFormat(true);
    const std::string time_series = mParameters["time_series"].GetString();

    try {
        if (format_name == "xdmf" && time_series == "automatic") {
            // Transient XDMF: extend the target file's temporal collection in place.
            WriteXdmfStep(rModelPart, rTargetSuffix);
        } else if (time_series == "single_file") {
            WriteStatic(ComposeOutputPath(rTargetSuffix, ""), format_name, rModelPart);
        } else { // "file_series", or "automatic" for formats without in-file appending
            WriteStatic(ComposeOutputPath(rTargetSuffix, GetOutputLabel(rModelPart)), format_name, rModelPart);
        }
    } catch (const std::out_of_range& r_exception) {
        // meshio++ throws std::out_of_range when an entity references a node
        // that is not part of the written model part (possible for
        // inconsistent sub model parts when "output_sub_model_parts" is on).
        KRATOS_ERROR << "Could not write model part \"" << rModelPart.FullName()
                     << "\": an entity references a node the model part does not contain ("
                     << r_exception.what() << ")" << std::endl;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool MeshioPlusPlusIO::WritesElements() const
{
    return mParameters["entity_type"].GetString() != "condition";
}

bool MeshioPlusPlusIO::WritesConditions() const
{
    return mParameters["entity_type"].GetString() != "element";
}

/***********************************************************************************/
/***********************************************************************************/

std::filesystem::path MeshioPlusPlusIO::ComposeOutputPath(
    const std::string& rTargetSuffix,
    const std::string& rLabel
    ) const
{
    const std::filesystem::path directory = mFileName.parent_path();
    std::string name = mParameters["custom_name_prefix"].GetString()
        + mFileName.stem().string()
        + rTargetSuffix
        + mParameters["custom_name_postfix"].GetString()
        + (rLabel.empty() ? "" : "_" + rLabel)
        + mFileName.extension().string();
    return directory / name;
}

/***********************************************************************************/
/***********************************************************************************/

std::string MeshioPlusPlusIO::GetOutputLabel(const ModelPart& rModelPart) const
{
    const auto& r_process_info = rModelPart.GetProcessInfo();
    if (mParameters["output_control_type"].GetString() == "time") {
        const double time = r_process_info.Has(TIME) ? r_process_info[TIME] : static_cast<double>(mOutputStep);
        std::ostringstream label;
        label << std::fixed << std::setprecision(mParameters["label_precision"].GetInt()) << time;
        return label.str();
    }
    const int step = r_process_info.Has(STEP) ? r_process_info[STEP] : static_cast<int>(mOutputStep);
    return std::to_string(step);
}

/***********************************************************************************/
/***********************************************************************************/

double MeshioPlusPlusIO::GetOutputTimeValue(const ModelPart& rModelPart) const
{
    const auto& r_process_info = rModelPart.GetProcessInfo();
    if (mParameters["output_control_type"].GetString() == "time") {
        return r_process_info.Has(TIME) ? r_process_info[TIME] : static_cast<double>(mOutputStep);
    }
    return static_cast<double>(r_process_info.Has(STEP) ? r_process_info[STEP] : static_cast<int>(mOutputStep));
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<mio::XdmfTimeSeriesWriter::NamedArray> MeshioPlusPlusIO::CollectPointData(const ModelPart& rThisModelPart) const
{
    KRATOS_TRY

    std::vector<DataArray> point_data;
    const auto& r_nodes = rThisModelPart.Nodes();
    const std::size_t number_of_nodes = rThisModelPart.NumberOfNodes();

    // Historical nodal variables (validated against the solution step data)
    CollectVariableDataArrays(
        mParameters["nodal_solution_step_data_variables"].GetStringArray(), r_nodes, number_of_nodes,
        [](const auto& rNode, const auto& rVariable) -> decltype(auto) {
            return rNode.FastGetSolutionStepValue(rVariable);
        },
        [&rThisModelPart](const auto& rVariable) {
            KRATOS_ERROR_IF(!rThisModelPart.HasNodalSolutionStepVariable(rVariable))
                << "Variable " << rVariable.Name() << " is not a nodal solution step variable of model part \""
                << rThisModelPart.FullName() << "\"" << std::endl;
        },
        point_data);

    // Non-historical nodal variables (includes the extrapolated gauss point results)
    CollectVariableDataArrays(
        mParameters["nodal_data_value_variables"].GetStringArray(), r_nodes, number_of_nodes,
        [](const auto& rNode, const auto& rVariable) -> decltype(auto) { return rNode.GetValue(rVariable); },
        [](const auto&) {}, point_data);

    // Nodal flags (1/0, -1 when undefined - VtkOutput convention)
    CollectFlagDataArrays(mParameters["nodal_flags"].GetStringArray(), r_nodes, number_of_nodes, point_data);

    if (mParameters["write_ids"].GetBool()) {
        point_data.push_back(CollectIdsArray(r_nodes, number_of_nodes, "KRATOS_NODE_ID"));
    }

    return point_data;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<mio::XdmfTimeSeriesWriter::NamedArray> MeshioPlusPlusIO::CollectCellData(const ModelPart& rThisModelPart) const
{
    KRATOS_TRY

    const bool write_elements = WritesElements();
    const bool write_conditions = WritesConditions();
    const std::size_t number_of_elements = write_elements ? rThisModelPart.NumberOfElements() : 0;
    const std::size_t number_of_conditions = write_conditions ? rThisModelPart.NumberOfConditions() : 0;
    const bool write_ids = mParameters["write_ids"].GetBool();
    const auto gauss_point_variables = mParameters["gauss_point_variables_in_elements"].GetStringArray();

    auto non_historical_getter = [](const auto& rEntity, const auto& rVariable) -> decltype(auto) {
        return rEntity.GetValue(rVariable);
    };
    auto no_validation = [](const auto&) {};

    std::vector<DataArray> element_part;
    if (write_elements && number_of_elements > 0) {
        const auto& r_elements = rThisModelPart.Elements();
        CollectVariableDataArrays(mParameters["element_data_value_variables"].GetStringArray(), r_elements,
                                  number_of_elements, non_historical_getter, no_validation, element_part);
        CollectFlagDataArrays(mParameters["element_flags"].GetStringArray(), r_elements, number_of_elements,
                              element_part);
        if (write_ids) {
            element_part.push_back(CollectIdsArray(r_elements, number_of_elements, "KRATOS_ELEMENT_ID"));
            element_part.push_back(CollectPropertiesIdsArray(r_elements, number_of_elements));
        }
        if (!gauss_point_variables.empty()) {
            // CalculateOnIntegrationPoints is non-const, hence the contained
            // const_cast (identical effective behavior to VtkOutput)
            auto& r_mutable_elements = const_cast<ModelPart&>(rThisModelPart).Elements();
            CollectGaussPointDataArrays(gauss_point_variables, r_mutable_elements, number_of_elements,
                                        rThisModelPart.GetProcessInfo(), element_part);
        }
    }

    std::vector<DataArray> condition_part;
    if (write_conditions && number_of_conditions > 0) {
        const auto& r_conditions = rThisModelPart.Conditions();
        CollectVariableDataArrays(mParameters["condition_data_value_variables"].GetStringArray(), r_conditions,
                                  number_of_conditions, non_historical_getter, no_validation, condition_part);
        CollectFlagDataArrays(mParameters["condition_flags"].GetStringArray(), r_conditions, number_of_conditions,
                              condition_part);
        if (write_ids) {
            condition_part.push_back(CollectIdsArray(r_conditions, number_of_conditions, "KRATOS_CONDITION_ID"));
            condition_part.push_back(CollectPropertiesIdsArray(r_conditions, number_of_conditions));
        }
        if (!gauss_point_variables.empty()) {
            auto& r_mutable_conditions = const_cast<ModelPart&>(rThisModelPart).Conditions();
            CollectGaussPointDataArrays(gauss_point_variables, r_mutable_conditions, number_of_conditions,
                                        rThisModelPart.GetProcessInfo(), condition_part);
        }
    }

    return MergeCellDataParts(std::move(element_part), number_of_elements,
                              std::move(condition_part), number_of_conditions);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

mio::Mesh MeshioPlusPlusIO::BuildMeshWithData(const ModelPart& rThisModelPart) const
{
    KRATOS_TRY

    const bool write_elements = WritesElements();
    const bool write_conditions = WritesConditions();

    mio::Mesh mesh;
    Internals::FillMeshioModelPart(rThisModelPart, mesh.GetModelPart(), write_elements, write_conditions,
                        mParameters["write_deformed_configuration"].GetBool());

    // Nodal data set on the meshio++ model part becomes point data when the
    // staging is rebuilt from the model part view (same node container order).
    for (const auto& r_data : CollectPointData(rThisModelPart)) {
        mesh.GetModelPart().SetNodalData(r_data.mName, ToNDArray(r_data));
    }

    // Cell data: the combined arrays cover element rows first, then condition
    // rows; split them into the per-kind containers meshio++ restores per block.
    const std::size_t number_of_elements = write_elements ? rThisModelPart.NumberOfElements() : 0;
    const std::size_t number_of_conditions = write_conditions ? rThisModelPart.NumberOfConditions() : 0;
    for (const auto& r_data : CollectCellData(rThisModelPart)) {
        if (number_of_elements > 0) {
            DataArray element_slice;
            element_slice.mName = r_data.mName;
            element_slice.mNumComponents = r_data.mNumComponents;
            element_slice.mValues.assign(r_data.mValues.begin(),
                                        r_data.mValues.begin() + number_of_elements * r_data.mNumComponents);
            mesh.GetModelPart().SetElementalData(r_data.mName, ToNDArray(element_slice));
        }
        if (number_of_conditions > 0) {
            DataArray condition_slice;
            condition_slice.mName = r_data.mName;
            condition_slice.mNumComponents = r_data.mNumComponents;
            condition_slice.mValues.assign(r_data.mValues.begin() + number_of_elements * r_data.mNumComponents,
                                          r_data.mValues.end());
            mesh.GetModelPart().SetConditionalData(r_data.mName, ToNDArray(condition_slice));
        }
    }

    // The model part view was mutated directly: rebuild the point/cell staging
    // the format writers read from.
    mesh.InvalidateBlocks();
    return mesh;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::WriteStatic(
    const std::filesystem::path& rPath,
    const std::string& rFormatName,
    const ModelPart& rThisModelPart
    ) const
{
    KRATOS_TRY

    mio::Mesh mesh = BuildMeshWithData(rThisModelPart);


    // Honor an ascii/binary override where the format supports it
    const std::string file_format = mParameters["file_format"].GetString();
    const bool skin = mParameters["skin"].GetBool();
    if (file_format != "default") {
        if (WriteWithFileFormatOverride(rFormatName, file_format == "binary", skin, rPath, mesh)) {
            return;
        }
        KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "The \"file_format\" setting is not supported for format \""
            << rFormatName << "\" - using its default" << std::endl;
    }

    // The registry writers extract the boundary skin of volume meshes for the
    // surface-only formats; honor a "skin" opt-out by calling those writers
    // directly with their registry ascii/binary defaults.
    if (!skin) {
        if (rFormatName == "stl") {
            mio::write_stl(rPath.string(), mesh, /*binary=*/false, /*skin=*/false);
            return;
        }
        if (rFormatName == "ply") {
            mio::write_ply(rPath.string(), mesh, /*binary=*/true, /*skin=*/false);
            return;
        }
    }

    mio::registry_writers().at(rFormatName)(rPath.string(), mesh);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::WriteXdmfStep(
    const ModelPart& rThisModelPart,
    const std::string& rTargetSuffix
    )
{
    KRATOS_TRY

    auto& p_writer = mXdmfWriters[rTargetSuffix];
    if (p_writer == nullptr) {
        std::string data_format = mParameters["xdmf_data_format"].GetString();
        if (data_format == "auto") {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
            data_format = "HDF";
#else
            data_format = "Binary";
#endif
        }

        // "Append" continues a series left by a previous run instead of destroying it.
        // On a path with no file yet it is exactly "Truncate", so it is passed
        // unconditionally rather than probing the filesystem here.
        p_writer = std::make_unique<mio::XdmfTimeSeriesWriter>(
            ComposeOutputPath(rTargetSuffix, "").string(), data_format,
            mParameters["xdmf_gzip_level"].GetInt(), mio::XdmfSeriesMode::Append);

        // Rewrite the light data after every step, so a run that is killed or crashes
        // still leaves a readable .xdmf. Upstream re-serializes the whole document on each
        // flush, which is quadratic in the step count - hence the opt-out.
        p_writer->SetAutoFlush(mParameters["xdmf_auto_flush"].GetBool());

        // Resuming a series recovers the step count, the mesh grid and (since meshio++
        // 9.2.0) the point/cell counts the array overload validates against, so the static
        // grid must not be written again and every step can take the array fast path.
        if (p_writer->NumSteps() > 0) {
            KRATOS_INFO("MeshioPlusPlusIO")
                << "Continuing the existing XDMF time series \""
                << ComposeOutputPath(rTargetSuffix, "").string() << "\" at step "
                << p_writer->NumSteps() << std::endl;
        } else {
            mio::Mesh mesh;
            Internals::FillMeshioModelPart(rThisModelPart, mesh.GetModelPart(), WritesElements(), WritesConditions(),
                                mParameters["write_deformed_configuration"].GetBool());
            mesh.InvalidateBlocks();
            p_writer->WritePointsCells(mesh);
        }
    }

    // Only the values change from step to step, so the arrays are handed over directly
    // rather than re-staging the whole model part into a Mesh.
    p_writer->WriteData(GetOutputTimeValue(rThisModelPart), CollectPointData(rThisModelPart),
                        CollectCellData(rThisModelPart));

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::string MeshioPlusPlusIO::Info() const
{
    return "MeshioPlusPlusIO";
}

void MeshioPlusPlusIO::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "MeshioPlusPlusIO";
}

void MeshioPlusPlusIO::PrintData(std::ostream& rOStream) const
{
    rOStream << "File name: " << mFileName << "\n"
             << "Settings: " << mParameters.PrettyPrintJsonString() << std::endl;
}

} // namespace Kratos
