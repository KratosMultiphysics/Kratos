/*
==============================================================================
KratosPFEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-11-10 14:23:33 $
//   Revision:            $Revision: 1.2 $
//
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

