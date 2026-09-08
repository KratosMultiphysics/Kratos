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
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <utility>

// External includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/gid.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/operations/sniff.hpp"

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

using Internals::DataArray;

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
    X(GID, "gid")                        \
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

/***********************************************************************************/
/***********************************************************************************/

/// The compound extensions meshio++ registers, longest first. std::filesystem splits at the
/// *last* dot, which would tear "results.post.msh" into "results.post" + ".msh" - and a labelled
/// step built from those halves ("results.post_3.msh") no longer resolves to gid at all, since
/// resolve_format needs the whole ".post.msh" to beat ".msh"'s gmsh entry. Everything that
/// composes or matches an output name goes through SplitStemAndExtension instead.
const std::array<const char*, 4>& CompoundExtensions()
{
    static const std::array<const char*, 4> extensions = {".post.msh", ".post.res", ".post.bin", ".post.h5"};
    return extensions;
}

/// Splits a path into (stem, extension), honoring the compound extensions above and falling
/// back to std::filesystem's own last-dot split for every other format.
std::pair<std::string, std::string> SplitStemAndExtension(const std::filesystem::path& rPath)
{
    const std::string name = rPath.filename().string();
    for (const char* p_extension : CompoundExtensions()) {
        const std::string extension(p_extension);
        if (name.size() > extension.size() &&
            name.compare(name.size() - extension.size(), extension.size(), extension) == 0) {
            return {name.substr(0, name.size() - extension.size()), extension};
        }
    }
    return {rPath.stem().string(), rPath.extension().string()};
}

/// Resolves the GiD on-disk flavour. "auto" lets meshio++ infer it from the extension, except
/// that a "file_format" : "binary" - the ascii/binary axis every other format uses - is honored
/// as a request for the binary flavour rather than silently ignored.
mio::GidMode ResolveGidMode(
    const std::string& rGidMode,
    const std::string& rFileFormat
    )
{
    if (rGidMode != "auto") {
        return mio::gid_mode_from_name(rGidMode);
    }
    if (rFileFormat == "binary") {
        return mio::GidMode::Binary;
    }
    if (rFileFormat == "ascii") {
        return mio::GidMode::Ascii;
    }
    return mio::GidMode::Auto;
}

/// Reads the "provenance" setting. "default" leaves meshio++'s own default alone, which since
/// v10.17.0 is BestEffort - so the setting is an override, not the thing that turns the record
/// on. Anything else names a mode explicitly.
mio::detail::ProvenanceMode ResolveProvenanceMode(const std::string& rSetting)
{
    if (rSetting == "default") {
        return mio::detail::default_provenance_mode();
    }
    if (rSetting == "off") {
        return mio::detail::ProvenanceMode::Off;
    }
    if (rSetting == "best_effort") {
        return mio::detail::ProvenanceMode::BestEffort;
    }
    if (rSetting == "required") {
        return mio::detail::ProvenanceMode::Required;
    }
    KRATOS_ERROR << "Unknown \"provenance\" setting \"" << rSetting
                 << "\" (use \"default\", \"off\", \"best_effort\" or \"required\")" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

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

    // Throws by name on an unknown mode, so a typo fails at construction like the rest.
    ResolveProvenanceMode(mParameters["provenance"].GetString());

    // GiD's flavour axis is four-valued where every other format's is ascii/binary, so it gets
    // its own setting rather than overloading "file_format". mio::gid_mode_from_name throws on
    // an unknown name; catching it here keeps the failure at construction like the rest.
    const std::string gid_mode = mParameters["gid_mode"].GetString();
    if (gid_mode != "auto") {
        try {
            mio::gid_mode_from_name(gid_mode);
        } catch (const std::exception& r_exception) {
            KRATOS_ERROR << "Unknown \"gid_mode\" setting \"" << gid_mode
                         << "\" (use \"auto\", \"ascii\", \"binary\", \"hdf5\" or "
                         << "\"ascii_zipped\"): " << r_exception.what() << std::endl;
        }
    }

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

    // Both output kinds are drained here, and both must be released even if the other fails -
    // this is what a caller calls before deleting a run's output, so a writer still holding a
    // file after a partial failure is the one outcome to avoid. The first exception is
    // re-thrown once everything has been released.
    std::exception_ptr first_failure;

    // The GiD series exists only in memory until now, so this is where it becomes a file at all.
    try {
        FlushGidSeries();
    } catch (...) {
        mGidSeries.clear();
        first_failure = std::current_exception();
    }

    // Finalize before releasing: the destructor would do it too, but then a failure could
    // only be logged. Clearing the map is what actually hands the file back to the caller,
    // so a later delete is not undone by a writer still holding it.
    for (auto& r_entry : mXdmfWriters) {
        if (r_entry.second != nullptr) {
            try {
                r_entry.second->Finalize();
            } catch (...) {
                if (!first_failure) {
                    first_failure = std::current_exception();
                }
            }
        }
    }
    mXdmfWriters.clear();

    if (first_failure) {
        std::rethrow_exception(first_failure);
    }

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
        "write_mdpa_ids"                              : false,
        "provenance"                                  : "default",
        "gid_mode"                                    : "auto",
        "gid_analysis_name"                           : "Kratos",
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

std::string MeshioPlusPlusIO::SniffFormat(const std::string& rFileName)
{
    KRATOS_TRY

    return mio::sniff_format(rFileName);

    KRATOS_CATCH("")
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
    // meshio++ deliberately keeps "gid" out of registry_compiled_out() so a ".post.msh" path
    // always resolves to gid's own actionable error instead of falling through to gmsh - which
    // means the generic check below would call it available in a build with no gidpost. Ask the
    // format itself instead. Reported for the write side, as for every other format here; the
    // ascii reader is available even where the writer is not (mio::gid_readable).
    if (Value == Format::GID) {
        return mio::gid_available(mio::GidMode::Auto);
    }
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

    // A file series (every format but XDMF, by default): mFileName itself was never
    // written, so this is checked - and needs no format resolution - before the
    // existence guard and format-dependent lookup below.
    std::vector<double> series_time_values = DetectFileSeriesTimeValues();
    if (!series_time_values.empty()) {
        return series_time_values;
    }

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

Parameters MeshioPlusPlusIO::GetProvenance() const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(std::filesystem::exists(mFileName))
        << "The file " << mFileName << " does not exist" << std::endl;

    const std::string format_name = ResolveEffectiveFormat(false);
    const mio::MeshMetadata metadata =
        mio::registry_read_metadata(mFileName.string(), format_name, mio::ReadOptions());

    Parameters result(R"({})");
    // "recognised" is the honest distinction between "meshio++ wrote this" and "something left
    // a comment at the top of the file" - the lines are reported either way.
    result.AddBool("recognised", metadata.mProvenanceRecognised);
    result.AddStringArray("lines", metadata.mProvenance);
    return result;

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
        } else if (format_name == "gid" && time_series == "automatic") {
            // Transient GiD: buffer the step; CloseOutput writes the whole series at once.
            WriteGidSeriesStep(rModelPart, rTargetSuffix);
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
    const auto [stem, extension] = SplitStemAndExtension(mFileName);
    std::string name = mParameters["custom_name_prefix"].GetString()
        + stem
        + rTargetSuffix
        + mParameters["custom_name_postfix"].GetString()
        + (rLabel.empty() ? "" : "_" + rLabel)
        + extension;
    return directory / name;
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<double> MeshioPlusPlusIO::DetectFileSeriesTimeValues() const
{
    KRATOS_TRY

    const std::filesystem::path directory = mFileName.parent_path();
    const std::filesystem::path scan_directory = directory.empty() ? std::filesystem::path(".") : directory;

    std::error_code error_code;
    if (!std::filesystem::is_directory(scan_directory, error_code)) {
        return {};
    }

    // The ComposeOutputPath("", <label>) naming, minus the label itself: every root-target
    // step this IO (or another instance sharing the same prefix/postfix) wrote matches it.
    const auto [stem, extension] = SplitStemAndExtension(mFileName);
    const std::string required_prefix = mParameters["custom_name_prefix"].GetString()
        + stem
        + mParameters["custom_name_postfix"].GetString()
        + "_";
    const std::string required_suffix = extension;

    std::vector<double> time_values;
    for (const auto& r_entry : std::filesystem::directory_iterator(scan_directory, error_code)) {
        if (!r_entry.is_regular_file(error_code)) {
            continue;
        }
        const std::string name = r_entry.path().filename().string();
        if (name.size() <= required_prefix.size() + required_suffix.size() ||
            name.compare(0, required_prefix.size(), required_prefix) != 0 ||
            name.compare(name.size() - required_suffix.size(), required_suffix.size(), required_suffix) != 0) {
            continue;
        }

        const std::string label = name.substr(required_prefix.size(),
            name.size() - required_prefix.size() - required_suffix.size());
        try {
            std::size_t number_of_characters_parsed = 0;
            const double value = std::stod(label, &number_of_characters_parsed);
            if (number_of_characters_parsed == label.size()) {
                time_values.push_back(value);
            }
        } catch (const std::exception&) {
            // Not a numeric label (an unrelated file sharing the extension, say) - ignore it.
        }
    }

    std::sort(time_values.begin(), time_values.end());
    return time_values;

    KRATOS_CATCH("")
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

Internals::FieldDataSelection MeshioPlusPlusIO::GetFieldDataSelection() const
{
    Internals::FieldDataSelection selection;
    selection.NodalSolutionStepVariables = mParameters["nodal_solution_step_data_variables"].GetStringArray();
    selection.NodalDataValueVariables = mParameters["nodal_data_value_variables"].GetStringArray();
    selection.NodalFlags = mParameters["nodal_flags"].GetStringArray();
    selection.ElementDataValueVariables = mParameters["element_data_value_variables"].GetStringArray();
    selection.ElementFlags = mParameters["element_flags"].GetStringArray();
    selection.ConditionDataValueVariables = mParameters["condition_data_value_variables"].GetStringArray();
    selection.ConditionFlags = mParameters["condition_flags"].GetStringArray();
    selection.GaussPointVariables = mParameters["gauss_point_variables_in_elements"].GetStringArray();
    selection.WriteIds = mParameters["write_ids"].GetBool();
    selection.WriteMdpaIds = mParameters["write_mdpa_ids"].GetBool();
    return selection;
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<mio::XdmfTimeSeriesWriter::NamedArray> MeshioPlusPlusIO::CollectPointData(const ModelPart& rThisModelPart) const
{
    return Internals::CollectPointData(rThisModelPart, GetFieldDataSelection());
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<mio::XdmfTimeSeriesWriter::NamedArray> MeshioPlusPlusIO::CollectCellData(const ModelPart& rThisModelPart) const
{
    return Internals::CollectCellData(rThisModelPart, WritesElements(), WritesConditions(), GetFieldDataSelection());
}

/***********************************************************************************/
/***********************************************************************************/

mio::Mesh MeshioPlusPlusIO::BuildMeshWithData(const ModelPart& rThisModelPart) const
{
    return Internals::ModelPartToMeshWithData(
        rThisModelPart, WritesElements(), WritesConditions(),
        mParameters["write_deformed_configuration"].GetBool(), GetFieldDataSelection());
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

    // Names the model part the file came from, so the credit line meshio++ writes anyway also
    // says *what* was written. Non-copyable and non-movable, hence constructed in place.
    const mio::detail::ProvenanceScope provenance_scope(
        ResolveProvenanceMode(mParameters["provenance"].GetString()),
        mio::detail::current_provenance());
    mio::detail::provenance_set_source(rThisModelPart.FullName(), "kratos");

    mio::Mesh mesh = BuildMeshWithData(rThisModelPart);

    // GiD before the generic ascii/binary override: its flavour axis is four-valued and the
    // writer also takes an analysis name and a step value, none of which
    // WriteWithFileFormatOverride's (binary, skin) shape can carry.
    if (rFormatName == "gid") {
        mio::write_gid(rPath.string(), mesh,
                       ResolveGidMode(mParameters["gid_mode"].GetString(),
                                      mParameters["file_format"].GetString()),
                       mParameters["gid_analysis_name"].GetString(),
                       GetOutputTimeValue(rThisModelPart));
        return;
    }

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

    const mio::detail::ProvenanceScope provenance_scope(
        ResolveProvenanceMode(mParameters["provenance"].GetString()),
        mio::detail::current_provenance());
    mio::detail::provenance_set_source(rThisModelPart.FullName(), "kratos");

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

void MeshioPlusPlusIO::WriteGidSeriesStep(
    const ModelPart& rThisModelPart,
    const std::string& rTargetSuffix
    )
{
    KRATOS_TRY

    // meshio++'s write_gid_series pulls: it calls back for step 0, 1, ... until the callback
    // says stop, which it can only do once every step exists. So the step is staged and kept
    // rather than written, and CloseOutput does the writing. See mGidSeries on the cost.
    mGidSeries[rTargetSuffix].emplace_back(GetOutputTimeValue(rThisModelPart),
                                           BuildMeshWithData(rThisModelPart));

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::FlushGidSeries()
{
    KRATOS_TRY

    if (mGidSeries.empty()) {
        return;
    }

    const mio::detail::ProvenanceScope provenance_scope(
        ResolveProvenanceMode(mParameters["provenance"].GetString()),
        mio::detail::current_provenance());

    const mio::GidMode mode = ResolveGidMode(mParameters["gid_mode"].GetString(),
                                             mParameters["file_format"].GetString());
    const std::string analysis_name = mParameters["gid_analysis_name"].GetString();

    for (auto& r_entry : mGidSeries) {
        auto& r_steps = r_entry.second;
        if (r_steps.empty()) {
            continue;
        }
        const std::string path = ComposeOutputPath(r_entry.first, "").string();

        // The mesh is moved out rather than copied: the buffer is cleared right after, and a
        // staged mesh is the largest thing this class holds.
        mio::write_gid_series(
            path,
            [&r_steps](std::size_t Index, double& rTimeValue, mio::Mesh& rMesh) {
                if (Index >= r_steps.size()) {
                    return false;
                }
                rTimeValue = r_steps[Index].first;
                rMesh = std::move(r_steps[Index].second);
                return true;
            },
            mode, analysis_name);
    }
    mGidSeries.clear();

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
