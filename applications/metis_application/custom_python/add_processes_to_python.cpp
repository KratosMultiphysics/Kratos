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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/io.h"
#include "processes/process.h"
#include "custom_python/add_processes_to_python.h"
#include "python/vector_python_interface.h"

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



void AddProcessesToPython()
{
    using namespace boost::python;

    class_<std::vector<int> >("IndicesVector")
    .def(vector_indexing_suite<std::vector<int> >())
    ;
#ifndef KRATOS_USE_METIS_5
    class_<MetisScalarReorder, bases<Process> >("MetisScalarReorder",init<ModelPart&>())
    ;

    class_<MetisPartitioningProcess, bases<Process> >("MetisPartitioningProcess",
            init<ModelPart&, IO&, unsigned int, unsigned int>())
    .def(init<ModelPart&, IO&, unsigned int>())
    ;

    class_<MetisDivideInputToPartitionsProcess, bases<Process> >("MetisDivideInputToPartitionsProcess",
            init<IO&, unsigned int, unsigned int>())
    .def(init<IO&, unsigned int>())
    ;

    class_<MetisContactPartitioningProcess, bases<MetisPartitioningProcess> >("MetisContactPartitioningProcess",
            init<ModelPart&, IO&, unsigned int, std::vector<int>, unsigned int>())
    .def(init<ModelPart&, IO&, unsigned int, std::vector<int> >())
    ;

    class_<MetisPartitioningProcessQuadratic, bases<MetisPartitioningProcess> >("MetisPartitioningProcessQuadratic",
            init<ModelPart&, IO&, unsigned int, unsigned int>())
    .def(init<ModelPart&, IO&, unsigned int>())
    ;
#endif
    class_<MetisDivideHeterogeneousInputProcess, bases<Process> >("MetisDivideHeterogeneousInputProcess",
                                                                   init<IO&, unsigned int>())
            .def(init<IO&, unsigned int, int>())
            .def(init<IO&, unsigned int, int, int>())
            .def(init<IO&, unsigned int, int, int, bool>())
            ;

    class_<MetisDivideHeterogeneousInputInMemoryProcess, bases<Process> >("MetisDivideHeterogeneousInputInMemoryProcess",
                                                                   init<IO&, unsigned int>())
            .def(init<IO&, unsigned int, int>())
            .def(init<IO&, unsigned int, int, int>())
            .def(init<IO&, unsigned int, int, int, bool>())
            ;

    class_<MortonDivideInputToPartitionsProcess, bases<Process> >("MetisDivideNodalInputToPartitionsProcess",
                                                                   init<IO&, unsigned int, unsigned int>())
            .def(init<IO&, unsigned int>())
    ;

    class_<SetMPICommunicatorProcess, bases<Process> >("SetMPICommunicatorProcess",
            init<ModelPart&>())
    ;

}

} // namespace Python.

} // Namespace Kratos
