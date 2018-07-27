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
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/kutta_condition_process.h"



namespace Kratos
{

namespace Python
{
  void  AddCustomProcessesToPython(pybind11::module& m)
  {
	using namespace pybind11;

        class_<KuttaConditionProcess, KuttaConditionProcess::Pointer, Process >
        (m, "KuttaConditionProcess")
        .def(init<ModelPart&>())
        .def("Execute",&KuttaConditionProcess::Execute)
            ;
  }

}  // namespace Python.

} // Namespace Kratos
