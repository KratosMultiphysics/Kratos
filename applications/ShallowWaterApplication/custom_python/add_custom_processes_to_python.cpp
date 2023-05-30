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
#include "custom_processes/apply_perturbation_function_process.h"
#include "custom_processes/apply_sinusoidal_function_process.h"
#include "custom_processes/calculate_distance_to_boundary_process.h"
#include "custom_processes/depth_integration_process.h"
#include "custom_processes/write_from_sw_at_interface_process.h"


namespace Kratos
{

namespace Python
{

    void  AddCustomProcessesToPython(pybind11::module& m)
    {
        namespace py = pybind11;

        typedef ApplyPerturbationFunctionProcess<Variable<double>> ApplyPerturbationScalarFunctionProcess;
        py::class_<ApplyPerturbationScalarFunctionProcess, ApplyPerturbationScalarFunctionProcess::Pointer, Process>
        (m, "ApplyPerturbationFunctionToScalar")
        .def(py::init<ModelPart&, Node::Pointer, Variable<double>&, Parameters&>())
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

        py::class_<CalculateDistanceToBoundaryProcess, CalculateDistanceToBoundaryProcess::Pointer, Process>
        (m, "CalculateDistanceToBoundaryProcess")
        .def(py::init<Model&, Parameters>())
        .def(py::init<ModelPart&, ModelPart&, double>())
        ;

        py::class_<DepthIntegrationProcess<2>, DepthIntegrationProcess<2>::Pointer, Process>
        (m, "DepthIntegrationProcess2D")
        .def(py::init<Model&, Parameters>())
        ;

        py::class_<DepthIntegrationProcess<3>, DepthIntegrationProcess<3>::Pointer, Process>
        (m, "DepthIntegrationProcess3D")
        .def(py::init<Model&, Parameters>())
        ;

        py::class_<WriteFromSwAtInterfaceProcess<2>, WriteFromSwAtInterfaceProcess<2>::Pointer, Process>
        (m, "WriteFromSwAtInterfaceProcess2D")
        .def(py::init<Model&, Parameters>())
        ;

        py::class_<WriteFromSwAtInterfaceProcess<3>, WriteFromSwAtInterfaceProcess<3>::Pointer, Process>
        (m, "WriteFromSwAtInterfaceProcess3D")
        .def(py::init<Model&, Parameters>())
        ;

    }

}  // namespace Python.

} // Namespace Kratos
