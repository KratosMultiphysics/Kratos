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
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/elemental_refining_criteria_process.h"
#include "custom_processes/apply_perturbation_function_process.h"
#include "custom_processes/apply_sinusoidal_function_process.h"


namespace Kratos
{

namespace Python
{

    void  AddCustomProcessesToPython(pybind11::module& m)
    {
        namespace py = pybind11;

        py::class_<ElementalRefiningCriteriaProcess, ElementalRefiningCriteriaProcess::Pointer, Process>
        (m, "ElementalRefiningCriteriaProcess")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters>())
        .def(py::init<ModelPart&, const Variable<double>&, double, bool>())
        ;

        typedef ApplyPerturbationFunctionProcess<Variable<double>> ApplyPerturbationScalarFunctionProcess;
        py::class_<ApplyPerturbationScalarFunctionProcess, ApplyPerturbationScalarFunctionProcess::Pointer, Process>
        (m, "ApplyPerturbationFunctionToScalar")
        .def(py::init<ModelPart&, Node<3>::Pointer, Variable<double>&, Parameters&>())
        .def(py::init<ModelPart&, ModelPart::NodesContainerType&, Variable<double>&, Parameters&>())
        ;

        typedef ApplySinusoidalFunctionProcess<Variable<double>> ApplySinusoidalScalarFunctionProcess;
        py::class_<ApplySinusoidalScalarFunctionProcess, ApplySinusoidalScalarFunctionProcess::Pointer, Process>
        (m, "ApplySinusoidalFunctionToScalar")
        .def(py::init<ModelPart&, Variable<double>&, Parameters&>())
        ;

        typedef ApplySinusoidalFunctionProcess<Variable<array_1d<double,3>>> ApplySinusoidalVectorFunctionProcess;
        py::class_<ApplySinusoidalVectorFunctionProcess, ApplySinusoidalVectorFunctionProcess::Pointer, Process>
        (m, "ApplySinusoidalFunctionToVector")
        .def(py::init<ModelPart&, Variable<array_1d<double,3>>&, Parameters&>())
        ;

    }

}  // namespace Python.

} // Namespace Kratos
