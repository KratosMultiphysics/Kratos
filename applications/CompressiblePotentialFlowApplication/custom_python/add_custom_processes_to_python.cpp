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
#include "custom_processes/metrics_potential_hessian_process.h"

namespace Kratos {
namespace Python {

typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;

void  AddCustomProcessesToPython(pybind11::module& m)
{
	namespace py = pybind11;

    py::class_<KuttaConditionProcess, KuttaConditionProcess::Pointer, Process >
        (m, "KuttaConditionProcess")
        .def(py::init<ModelPart&>())
        ;
    // HESSIAN PROCESS
        py::class_<ComputePotentialHessianSolMetricProcess, ComputePotentialHessianSolMetricProcess::Pointer, Process>
        (m, "ComputePotentialHessianSolMetricProcess")
        .def(py::init<ModelPart&, Variable<double>&>())
        .def(py::init<ModelPart&, Variable<double>&, Parameters>())
        .def(py::init<ModelPart&, ComponentType&>())
        .def(py::init<ModelPart&, ComponentType&, Parameters>())
        ;

        m.attr("ComputeHessianSolMetricProcess2D") = m.attr("ComputeHessianSolMetricProcess");

}  // namespace Python.

} // Namespace Kratos
