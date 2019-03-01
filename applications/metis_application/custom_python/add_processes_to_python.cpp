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
#include "includes/io.h"
#include "includes/model_part.h"
#include "includes/model_part_io.h"

#include "processes/process.h"

#include "custom_python/add_processes_to_python.h"

#include "custom_processes/metis_divide_heterogeneous_input_process.h"
#include "custom_processes/metis_divide_heterogeneous_input_in_memory_process.h"
#include "custom_processes/morton_divide_input_to_partitions_process.h"
#include "custom_processes/set_mpi_communicator_process.h"

#ifndef KRATOS_USE_METIS_5
#include "custom_processes/metis_partitioning_process.h"
#include "custom_processes/metis_divide_input_to_partitions_process.h"
#include "custom_processes/metis_contact_partitioning_process.h"
#include "custom_processes/metis_partitioning_process_quadratic.h"

#include "custom_processes/metis_scalar_reorder_process.h"
#endif

namespace Kratos {
namespace Python {

void AddProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

#ifndef KRATOS_USE_METIS_5
    py::class_<MetisScalarReorder, MetisScalarReorder::Pointer, Process>(m,"MetisScalarReorder")
        .def(py::init<ModelPart&>())
        ;

    py::class_<MetisPartitioningProcess, MetisPartitioningProcess::Pointer, Process>(m,"MetisPartitioningProcess")
        .def(py::init<ModelPart&, IO&, unsigned int, unsigned int>())
        .def(py::init<ModelPart&, IO&, unsigned int>())
        ;

    py::class_<MetisDivideInputToPartitionsProcess, MetisDivideInputToPartitionsProcess::Pointer, Process>(
        m,"MetisDivideInputToPartitionsProcess")
        .def(py::init<IO&, unsigned int, unsigned int>())
        .def(py::init<IO&, unsigned int>())
        ;

    py::class_<MetisContactPartitioningProcess, MetisContactPartitioningProcess::Pointer, MetisPartitioningProcess>(
        m, "MetisContactPartitioningProcess")
        .def(py::init<ModelPart&, IO&, unsigned int, std::vector<int>, unsigned int>())
        .def(py::init<ModelPart&, IO&, unsigned int, std::vector<int> >())
        ;

    py::class_<MetisPartitioningProcessQuadratic, MetisPartitioningProcessQuadratic::Pointer, MetisPartitioningProcess >(
        m, "MetisPartitioningProcessQuadratic")
        .def(py::init<ModelPart&, IO&, unsigned int, unsigned int>())
        .def(py::init<ModelPart&, IO&, unsigned int>())
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

    py::class_<SetMPICommunicatorProcess, SetMPICommunicatorProcess::Pointer, Process>(m,"SetMPICommunicatorProcess")
        .def(py::init<ModelPart&>())
        ;
}

} // namespace Python.
} // Namespace Kratos
