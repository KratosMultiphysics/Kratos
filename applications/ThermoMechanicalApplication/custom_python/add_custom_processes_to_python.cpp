/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-05-28 15:29:01 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes

// External includes
#include "pybind11/pybind11.h"


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "custom_processes/duplicate_interface_nodes_create_conditions_process.h"
#include "custom_processes/activation_deactivation_conditions_process.h"
#include "custom_processes/solidification_process.h"

#include "includes/node.h"

namespace Kratos
{

namespace Python
{
void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;


    py::class_<DuplicateInterfaceNodesCreateConditionsProcess, DuplicateInterfaceNodesCreateConditionsProcess::Pointer, Process>
    (m,"DuplicateInterfaceNodesCreateConditionsProcess")
    .def( py::init<ModelPart&, char* ,int, const Matrix >())
    .def("Execute", &DuplicateInterfaceNodesCreateConditionsProcess::Execute)
    .def("PairToId", &DuplicateInterfaceNodesCreateConditionsProcess::PairToId)
    .def("IdToPair", &DuplicateInterfaceNodesCreateConditionsProcess::IdToPair)
    ;
    py::class_<ActivationDeactivationConditionsProcess, ActivationDeactivationConditionsProcess::Pointer, Process>
    (m, "ActivationDeactivationConditionsProcess")
    .def( py::init<ModelPart& ,int, const Matrix >())
    .def("Execute", &ActivationDeactivationConditionsProcess::Execute)
    ;
    py::class_<SolidificationProcess, SolidificationProcess::Pointer, Process>
    (m, "SolidificationProcess")
    .def(py::init<ModelPart& ,const double  >())
    .def("Execute", &SolidificationProcess::Execute)
    ;

}

}  // namespace Python.

} // Namespace Kratos

