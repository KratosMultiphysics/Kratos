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

// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/kahip_divide_heterogeneous_input_process.h"
#include "custom_processes/kahip_divide_heterogeneous_input_in_memory_process.h"

namespace Kratos::Python {

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // ── KaHIPDivideHeterogeneousInputProcess ─────────────────────────────────
    // Drop-in replacement for MetisDivideHeterogeneousInputProcess.
    // Positional constructor mirrors the MetisApplication API exactly so that
    // Python scripts can switch by changing only the import and class name.
    py::class_<
        KaHIPDivideHeterogeneousInputProcess,
        KaHIPDivideHeterogeneousInputProcess::Pointer,
        Process>(m, "KaHIPDivideHeterogeneousInputProcess")

        // ── Parameters-based constructors ────────────────────────────────────
        .def(py::init<IO&, unsigned int, Parameters>(),
             py::arg("io"), py::arg("number_of_partitions"), py::arg("settings"))
        .def(py::init<IO&, unsigned int, Parameters, bool>(),
             py::arg("io"), py::arg("number_of_partitions"), py::arg("settings"),
             py::arg("synchronize_conditions"))

        // ── Compatibility constructors (mirrors MetisApplication) ─────────────
        // KaHIPDivideHeterogeneousInputProcess(io, nparts)
        .def(py::init<IO&, unsigned int>())
        // KaHIPDivideHeterogeneousInputProcess(io, nparts, dimension)
        .def(py::init<IO&, unsigned int, int>())
        // KaHIPDivideHeterogeneousInputProcess(io, nparts, dimension, verbosity)
        .def(py::init<IO&, unsigned int, int, int>())
        // KaHIPDivideHeterogeneousInputProcess(io, nparts, dimension, verbosity, sync_cond)
        .def(py::init<IO&, unsigned int, int, int, bool>())
        ;

    // ── KaHIPDivideHeterogeneousInputInMemoryProcess ──────────────────────────
    // MPI in-memory variant: partitions on rank 0, scatters per-rank streams
    // to all ranks via Scatterv. Equivalent to
    // MetisDivideHeterogeneousInputInMemoryProcess.
    py::class_<
        KaHIPDivideHeterogeneousInputInMemoryProcess,
        KaHIPDivideHeterogeneousInputInMemoryProcess::Pointer,
        KaHIPDivideHeterogeneousInputProcess>(m, "KaHIPDivideHeterogeneousInputInMemoryProcess")

        // Parameters-based constructor
        .def(py::init<IO&, ModelPartIO&, const DataCommunicator&, Parameters>(),
             py::arg("io"), py::arg("serial_io"), py::arg("data_communicator"),
             py::arg("settings"))
        .def(py::init<IO&, ModelPartIO&, const DataCommunicator&, Parameters, bool>(),
             py::arg("io"), py::arg("serial_io"), py::arg("data_communicator"),
             py::arg("settings"), py::arg("synchronize_conditions"))

        // Compatibility constructors (mirrors MetisApplication)
        .def(py::init<IO&, ModelPartIO&, const DataCommunicator&>())
        .def(py::init<IO&, ModelPartIO&, const DataCommunicator&, int>())
        .def(py::init<IO&, ModelPartIO&, const DataCommunicator&, int, int>())
        .def(py::init<IO&, ModelPartIO&, const DataCommunicator&, int, int, bool>())
        ;
}

} // namespace Kratos::Python
