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
#include "includes/model_part.h"
#include "includes/io.h"
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

namespace Kratos
{

namespace Python
{



void AddProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;


#ifndef KRATOS_USE_METIS_5
    class_<MetisScalarReorder, MetisScalarReorder::Pointer, Process>(m,"MetisScalarReorder").def(init<ModelPart&>())
    ;

    class_<MetisPartitioningProcess, MetisPartitioningProcess::Pointer, Process>(m,"MetisPartitioningProcess")
    .def(init<ModelPart&, IO&, unsigned int, unsigned int>())
    .def(init<ModelPart&, IO&, unsigned int>())
    ;

    class_<MetisDivideInputToPartitionsProcess, MetisDivideInputToPartitionsProcess::Pointer, Process>(m,"MetisDivideInputToPartitionsProcess")
    .def(init<IO&, unsigned int, unsigned int>())
    .def(init<IO&, unsigned int>())
    ;

    class_<MetisContactPartitioningProcess, MetisContactPartitioningProcess::Pointer, MetisPartitioningProcess>(m, "MetisContactPartitioningProcess")
    .def(init<ModelPart&, IO&, unsigned int, std::vector<int>, unsigned int>())
    .def(init<ModelPart&, IO&, unsigned int, std::vector<int> >())
    ;

    class_<MetisPartitioningProcessQuadratic, MetisPartitioningProcessQuadratic::Pointer, MetisPartitioningProcess >(m, "MetisPartitioningProcessQuadratic")
    .def(init<ModelPart&, IO&, unsigned int, unsigned int>())
    .def(init<ModelPart&, IO&, unsigned int>())
    ;
#endif
    class_<MetisDivideHeterogeneousInputProcess, MetisDivideHeterogeneousInputProcess::Pointer, Process>(m,"MetisDivideHeterogeneousInputProcess")
    .def(init<IO&, unsigned int>())
            .def(init<IO&, unsigned int, int>())
            .def(init<IO&, unsigned int, int, int>())
            .def(init<IO&, unsigned int, int, int, bool>())
            ;

    class_<MetisDivideHeterogeneousInputInMemoryProcess, MetisDivideHeterogeneousInputInMemoryProcess::Pointer, Process>(m,"MetisDivideHeterogeneousInputInMemoryProcess")
    .def(init<IO&, unsigned int>())
            .def(init<IO&, unsigned int, int>())
            .def(init<IO&, unsigned int, int, int>())
            .def(init<IO&, unsigned int, int, int, bool>())
            ;

    class_<MortonDivideInputToPartitionsProcess, MortonDivideInputToPartitionsProcess::Pointer, Process>(m,"MetisDivideNodalInputToPartitionsProcess")
    .def(init<IO&, std::size_t, int>())
            .def(init<IO&, std::size_t>())
    ;

    class_<SetMPICommunicatorProcess, SetMPICommunicatorProcess::Pointer, Process>(m,"SetMPICommunicatorProcess")
    .def(init<ModelPart&>())
    ;

}

} // namespace Python.

} // Namespace Kratos
