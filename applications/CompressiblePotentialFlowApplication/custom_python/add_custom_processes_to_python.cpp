//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/kutta_condition_process.h"
#include "custom_processes/replace_elements_and_conditions_for_adjoint_problem_process.cpp"

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
	namespace py = pybind11;

    py::class_<KuttaConditionProcess, KuttaConditionProcess::Pointer, Process >
        (m, "KuttaConditionProcess")
        .def(py::init<ModelPart&>())
        ;

    py::class_<ReplaceElementsAndConditionsAdjointProcess, ReplaceElementsAndConditionsAdjointProcess::Pointer, Process >
        (m, "ReplaceElementsAndConditionsAdjointProcess")
        .def(py::init<ModelPart&>())
        .def("Execute",&ReplaceElementsAndConditionsAdjointProcess::Execute);
}

}  // namespace Python.

} // Namespace Kratos
