//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Jordi Cotela
//                   Carlos Roig
//


// System includes

// External includes

// Project includes
#include "custom_python/add_processes_to_python.h"

#include "custom_processes/metis_divide_heterogeneous_input_process.h"
#include "custom_processes/metis_divide_heterogeneous_input_in_memory_process.h"
#include "custom_processes/morton_divide_input_to_partitions_process.h"

#ifndef KRATOS_USE_METIS_5
#include "custom_processes/metis_divide_input_to_partitions_process.h"
#endif

namespace Kratos {
namespace Python {

void AddProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

#ifndef KRATOS_USE_METIS_5

    py::class_<MetisDivideInputToPartitionsProcess, MetisDivideInputToPartitionsProcess::Pointer, Process>(
        m,"MetisDivideInputToPartitionsProcess")
        .def(py::init<IO&, unsigned int, unsigned int>())
        .def(py::init<IO&, unsigned int>())
        ;
#endif
    py::class_<MetisDivideHeterogeneousInputProcess, MetisDivideHeterogeneousInputProcess::Pointer, Process>(
        m,"MetisDivideHeterogeneousInputProcess")
        .def(py::init<IO&, unsigned int>())
        .def(py::init<IO&, unsigned int, int>())
        .def(py::init<IO&, unsigned int, int, int>())
        .def(py::init<IO&, unsigned int, int, int, bool>())
        ;

    py::class_<MetisDivideHeterogeneousInputInMemoryProcess, MetisDivideHeterogeneousInputInMemoryProcess::Pointer, Process>(
        m,"MetisDivideHeterogeneousInputInMemoryProcess")
        .def(py::init<IO&, ModelPartIO&, unsigned int>())
        .def(py::init<IO&, ModelPartIO&, unsigned int, int>())
        .def(py::init<IO&, ModelPartIO&, unsigned int, int, int>())
        .def(py::init<IO&, ModelPartIO&, unsigned int, int, int, bool>())
        ;

    py::class_<MortonDivideInputToPartitionsProcess, MortonDivideInputToPartitionsProcess::Pointer, Process>(
        m,"MetisDivideNodalInputToPartitionsProcess")
        .def(py::init<IO&, std::size_t, int>())
        .def(py::init<IO&, std::size_t>())
        ;

}

} // namespace Python.
} // Namespace Kratos
