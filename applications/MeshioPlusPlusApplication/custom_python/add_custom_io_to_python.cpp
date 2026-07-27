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
#include <filesystem>

// External includes
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

// Project includes
#include "custom_python/add_custom_io_to_python.h"
#include "custom_io/meshioplusplus_io.h"
#include "includes/model_part.h"

namespace Kratos::Python
{

void AddCustomIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto meshio_io = py::class_<MeshioPlusPlusIO, MeshioPlusPlusIO::Pointer, IO>(m, "MeshioPlusPlusIO",
        "Multi-format mesh IO based on the meshio++ library. Reads 41 and writes 43 mesh formats "
        "(vtu, vtk, gmsh, med, xdmf, abaqus, mdpa, ...) converting to/from Kratos model parts. "
        "Repeated WriteModelPart calls extend the output (XDMF temporal collection or file series).")
        .def(py::init<std::filesystem::path const&>(), py::arg("file_name"))
        .def(py::init<std::filesystem::path const&, Parameters>(), py::arg("file_name"), py::arg("parameters"))
        .def("ReadModelPart", &MeshioPlusPlusIO::ReadModelPart, py::arg("model_part"),
            "Reads the file into the given (empty) model part. Highest-dimension cells become "
            "elements, lower-dimension ones conditions, and named regions become sub model parts. "
            "Registered entity names are restored where the format carries them.")
        .def("WriteModelPart", py::overload_cast<const ModelPart&>(&MeshioPlusPlusIO::WriteModelPart), py::arg("model_part"),
            "Writes the model part. Repeated calls extend the current output: XDMF appends steps "
            "to the temporal collection of the file, other formats write a file series.")
        .def("GetNumberOfTimeSteps", &MeshioPlusPlusIO::GetNumberOfTimeSteps,
            "Number of transient steps this file's format reports, without a full read. "
            "0 for a format with no time-series concept.")
        .def("GetTimeValues", &MeshioPlusPlusIO::GetTimeValues,
            "The simulation time of each transient step; size matches GetNumberOfTimeSteps().")
        .def("GetTimeStepIndex", &MeshioPlusPlusIO::GetTimeStepIndex, py::arg("time_value"),
            "The 0-based index of the transient step whose time value matches time_value "
            "(within 1e-9), or -1 if none matches.")
        .def_static("GetDefaultParameters", &MeshioPlusPlusIO::GetDefaultParameters,
            "The default settings of the IO.")
        .def_static("GetSupportedFormats", &MeshioPlusPlusIO::GetSupportedFormats,
            "Names of every format this build can read or write.")
        .def_static("GetSupportedReadFormats", &MeshioPlusPlusIO::GetSupportedReadFormats,
            "Names of every format this build can read.")
        .def_static("GetSupportedWriteFormats", &MeshioPlusPlusIO::GetSupportedWriteFormats,
            "Names of every format this build can write.")
        .def_static("FormatFromString", &MeshioPlusPlusIO::FormatFromString, py::arg("format_name"),
            "Converts a format name to the Format enum ('gmsh' -> Format.GMSH). Throws on unknown names.")
        .def_static("FormatName", &MeshioPlusPlusIO::FormatName, py::arg("format"),
            "The canonical name of a Format value (Format.GMSH -> 'gmsh').")
        .def_static("ResolveFormat", &MeshioPlusPlusIO::ResolveFormat, py::arg("path"),
            "Resolves the format from a file path extension ('foo.vtu' -> Format.VTU).")
        .def_static("IsFormatAvailable", &MeshioPlusPlusIO::IsFormatAvailable, py::arg("format"),
            "False when the format was compiled out for lack of an optional dependency (HDF5/netCDF).")
        ;

    py::enum_<MeshioPlusPlusIO::Format>(meshio_io, "Format")
        .value("AUTOMATIC", MeshioPlusPlusIO::Format::AUTOMATIC)
        .value("ABAQUS", MeshioPlusPlusIO::Format::ABAQUS)
        .value("ANSYS", MeshioPlusPlusIO::Format::ANSYS)
        .value("ANSYSINP", MeshioPlusPlusIO::Format::ANSYSINP)
        .value("AVSUCD", MeshioPlusPlusIO::Format::AVSUCD)
        .value("CGNS", MeshioPlusPlusIO::Format::CGNS)
        .value("DEX", MeshioPlusPlusIO::Format::DEX)
        .value("DOLFIN", MeshioPlusPlusIO::Format::DOLFIN)
        .value("ENSIGHT", MeshioPlusPlusIO::Format::ENSIGHT)
        .value("EXODUS", MeshioPlusPlusIO::Format::EXODUS)
        .value("FLAC3D", MeshioPlusPlusIO::Format::FLAC3D)
        .value("FLUX", MeshioPlusPlusIO::Format::FLUX)
        .value("FREEFEM", MeshioPlusPlusIO::Format::FREEFEM)
        .value("GMSH", MeshioPlusPlusIO::Format::GMSH)
        .value("GMSH22", MeshioPlusPlusIO::Format::GMSH22)
        .value("H5M", MeshioPlusPlusIO::Format::H5M)
        .value("HMF", MeshioPlusPlusIO::Format::HMF)
        .value("IP", MeshioPlusPlusIO::Format::IP)
        .value("MDPA", MeshioPlusPlusIO::Format::MDPA)
        .value("MED", MeshioPlusPlusIO::Format::MED)
        .value("MEDIT", MeshioPlusPlusIO::Format::MEDIT)
        .value("MFF", MeshioPlusPlusIO::Format::MFF)
        .value("MFM", MeshioPlusPlusIO::Format::MFM)
        .value("MPHTXT", MeshioPlusPlusIO::Format::MPHTXT)
        .value("NASTRAN", MeshioPlusPlusIO::Format::NASTRAN)
        .value("NETGEN", MeshioPlusPlusIO::Format::NETGEN)
        .value("OBJ", MeshioPlusPlusIO::Format::OBJ)
        .value("OFF", MeshioPlusPlusIO::Format::OFF)
        .value("OPENFOAM", MeshioPlusPlusIO::Format::OPENFOAM)
        .value("PERMAS", MeshioPlusPlusIO::Format::PERMAS)
        .value("PLY", MeshioPlusPlusIO::Format::PLY)
        .value("STL", MeshioPlusPlusIO::Format::STL)
        .value("SU2", MeshioPlusPlusIO::Format::SU2)
        .value("SVG", MeshioPlusPlusIO::Format::SVG)
        .value("TECPLOT", MeshioPlusPlusIO::Format::TECPLOT)
        .value("TETGEN", MeshioPlusPlusIO::Format::TETGEN)
        .value("TIKZ", MeshioPlusPlusIO::Format::TIKZ)
        .value("TRIANGLE", MeshioPlusPlusIO::Format::TRIANGLE)
        .value("UGRID", MeshioPlusPlusIO::Format::UGRID)
        .value("UNV", MeshioPlusPlusIO::Format::UNV)
        .value("VTK", MeshioPlusPlusIO::Format::VTK)
        .value("VTP", MeshioPlusPlusIO::Format::VTP)
        .value("VTU", MeshioPlusPlusIO::Format::VTU)
        .value("WKT", MeshioPlusPlusIO::Format::WKT)
        .value("XDMF", MeshioPlusPlusIO::Format::XDMF)
        ;
}

} // namespace Kratos::Python
