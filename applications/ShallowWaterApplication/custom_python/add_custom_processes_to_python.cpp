//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/elemental_refining_criteria_process.h"


namespace Kratos
{

namespace Python
{

    void  AddCustomProcessesToPython(pybind11::module& m)
    {
        namespace py = pybind11;

        py::class_<ElementalRefiningCriteriaProcess, ElementalRefiningCriteriaProcess::Pointer, Process>(m, "ElementalRefiningCriteriaProcess")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters>())
        .def(py::init<ModelPart&, Variable<double>, double, bool>())
        .def("Execute",&ElementalRefiningCriteriaProcess::Execute)
        ;
    }

}  // namespace Python.

} // Namespace Kratos
