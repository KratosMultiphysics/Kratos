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
#include "custom_processes/compute_lift_level_set_process.h"
#include "custom_processes/compute_gradient_adjoint_process.h"
#include "custom_processes/compute_gradient_numerical_process.h"

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
	using namespace pybind11;

        class_<KuttaConditionProcess, KuttaConditionProcess::Pointer, Process >
        (m, "KuttaConditionProcess")
        .def(init<ModelPart&>())
        .def("Execute",&KuttaConditionProcess::Execute);

        class_<ComputeLiftLevelSetProcess, ComputeLiftLevelSetProcess::Pointer, Process >
        (m, "ComputeLiftLevelSetProcess")
        .def(init<ModelPart&,Vector&>())
        .def("Execute",&ComputeLiftLevelSetProcess::Execute);

        class_<ComputeGradientAdjointProcess, ComputeGradientAdjointProcess::Pointer, Process >
        (m, "ComputeGradientAdjointProcess")
        .def(init<ModelPart&,Matrix&,Matrix&,Vector&>())
        .def("Execute",&ComputeGradientAdjointProcess::Execute);

        // class_<ComputeGradientNumericalProcess, ComputeGradientNumericalProcess::Pointer, Process >
        // (m, "ComputeGradientNumericalProcess")
        // .def(init<ModelPart&,Vector&>())
        // .def("Execute",&ComputeGradientNumericalProcess::Execute);
  }

}  // namespace Python.

} // Namespace Kratos
